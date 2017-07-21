#!/usr/bin/env python

# Process thread dumps into task/subtype/sort direction files and fit
# functions. Also work out the scaling of the costs that most closely
# fits the current functions to the data. XXX more details needed.

import matplotlib
matplotlib.use("Agg")

from scipy.optimize import curve_fit
import math
import numpy as np
import pylab as pl
import os
import sys

#  Our fitting functions.
def funcpair(x, a, b, c, d):
    return a + b*x[0]+ c*x[1] + d*x[0]*x[1]

def funcself(x, a, b, c):
    return a + b*x + c*x*x

def iszero(a):
    return 0

#  Utilities.
def gettask(line):
    # Lines look like ' task_type_none = 0,\n', ' task_type_sort,\n' or
    # 'task_type_count\n'. So find task_ and stop at space, comma or end of
    # line (although we don't expect to see "count").
    start = line.find("task_")
    for i in range(start,len(line)):
        if line[i] == " " or line[i] == "," or line[i] == "\n":
            return line[start:i]
    return None

#  Basic plot configuration.
PLOT_PARAMS = {"axes.labelsize": 10,
               "axes.titlesize": 10,
               "font.size": 12,
               "legend.fontsize": 12,
               "xtick.labelsize": 10,
               "ytick.labelsize": 10,
               "figure.figsize" : (12., 4.),
               "figure.subplot.left" : 0.1,
               "figure.subplot.right" : 0.9,
               "figure.subplot.bottom" : 0.1,
               "figure.subplot.top" : 0.85,
               "figure.subplot.wspace" : 0.,
               "figure.subplot.hspace" : 0.,
               "lines.markersize" : 6,
               "lines.linewidth" : 3.
               }
pl.rcParams.update(PLOT_PARAMS)

#  Tasks and subtypes extracted from ../src/task.h. XXX parameterise this.
TASKTYPES = []
SUBTYPES = []
with open(os.path.dirname(__file__) + "/../src/task.h") as task_file:
    for line in task_file:
        if "task_type_" in line:
            TASKTYPES.append(gettask(line))
        else:
            if "task_subtype_" in line:
                SUBTYPES.append(gettask(line))

#  Sorting directions in flags. Should give same results (0 to 12). See
#  sid_scale in scheduler.
SAMEFLAGS = [1, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2, 3]

#  Handle command-line.
if len(sys.argv) != 2:
    print "Usage: ", sys.argv[0], "input.dat"
    sys.exit(1)

#  Read input.
indata = pl.loadtxt(sys.argv[1])

#  Ignore global step lines (first for each rank).
indata = indata[indata[:,1] >= 0]

num_lines = pl.size(indata) / len(indata[0])
print "num_lines: ", num_lines

#  Create hash of unique keys "type subtype flags".
keys = {}
for line in range(num_lines):
    ttype = TASKTYPES[int(indata[line,3])]
    subtype = SUBTYPES[int(indata[line,4])]
    flags = int(indata[line,12])
    if flags > 12:
	flags = 12

    # Flags is a tag for these.
    if "recv" in ttype or "send" in ttype:
	flags = 0

    # All sorts are gathered together (multiple sorting directions).
    elif "sort" in ttype:
        flags = 0
    else:
        flags = SAMEFLAGS[flags]

    key = ttype + " " + subtype + " " + str(flags)
    if not key in keys:
        keys[key] = []
    keys[key].append(indata[line,:])

print "number of keys: ", len(keys)

# Create a file to contain all the solutions.
fsol = open( "swift-task-fitted-costs.txt", 'w')
fsol.write( "# task  subtask sortdir c1 c2 c3\n")

# Output the data to be fitted into files with the names "type.subtype.flags".dat.
# Also create fits and store these in "type.subtype.flags".fit and plot the fit...

for key in keys:
    print key
    fdat = open( key + ".dat", 'w' )
    ffit = open( key + ".fit", 'w' )
    fdat.write( "# rank otherrank rid itype isubtype self tic toc cicount cjcount cigcount cjgcount flags cost dt type subtype\n" )

    ydata = []
    sigma = []
    cidata = []
    cjdata = []
    cost = []

    ttype = TASKTYPES[int(keys[key][0][3])]
    subtype = SUBTYPES[int(keys[key][0][4])]

    num_lines = 0
    for line in keys[key]:
        res = ""
        for item in line:
            res = res + " " + str(int(item)) + " "
        dt = int(line[7]) - int(line[6])
        res = res + " " + str(dt) + " "
        res = res + " " + ttype + " "
        res = res + " " + subtype + "\n"
        fdat.write( res )

        cidata.append(line[8])
        cjdata.append(line[9])
        cost.append(line[14])
        ydata.append(dt)
        if ttype == "sort":
            #  Lots of pre-sorted data, try not to prefer that.
            sigma.append(1.0)
        else:
            sigma.append(math.sqrt(dt))

        num_lines = num_lines + 1

    ydata = np.array(ydata)
    sigma = np.array(sigma)
    cidata = np.array(cidata)
    cjdata = np.array(cjdata)
    cost = np.array(cost)

    # Self data has different fit. MPI is has cjdata, but really only like a self.
    if min(cjdata) == 0 or "_recv" in ttype or "_send" in ttype:
        isself = 1
    else:
        isself = 0

    if isself:
        ffit.write( "# cicount model dt dmodel fmodel scost\n" )
        minfunc = funcself
        xdata = cidata
        p0 = [0.,1.,1.]
    else:
        ffit.write( "# cicount cjcount model dt dmodel fmodel scost\n" )
        minfunc = funcpair
        xdata = [cidata,cjdata]
        p0 = [0.,1.,1.,1.]

    if len(cidata) < len(p0):
        print "Not enough data to fit '" + ttype + "." + subtype + "' skipping"
        print cidata
        continue
    popt, pcov = curve_fit(minfunc, xdata, ydata, p0=p0, sigma=sigma)
    print popt

    #  Format popt.
    fopt = []
    for i in range(len(popt)):
        fopt.append("{:f}".format(popt[i]))

    # Scale cost to solution for comparison.
    if isself:
        ymax = funcself(max(xdata), popt[0], popt[1], popt[2])
    else:
        ymax = funcpair([max(xdata[0]), max(xdata[1])], popt[0], popt[1], popt[2], popt[3])
    maxcost = float(max(cost))
    if maxcost > 0:
        scost = cost * ymax / maxcost
    else:
        scost = cost

    # Evaluate and output/plot.
    y = []
    fy = []
    x = []
    if isself:
        for i in range(num_lines):
            y.append(funcself(xdata[i], popt[0], popt[1], popt[2]))
            fy.append(funcself(xdata[i], 0, popt[1], popt[2]))
            ffit.write(str(xdata[i]) + " " + str(y[i]) + " " + str(ydata[i]) + " " + str(ydata[i] - y[i]) + " " + str(fy[i]) + " " + str(scost[i]) + "\n")
        x = xdata
        solution = fopt[0] + " + " + fopt[1] + " * cicount " + " + " + fopt[2] + " * cicount * cicount"
        fsol.write( key + " " + fopt[1] + " " + fopt[2] + "\n")
    else:
        for i in range(num_lines):
            y.append(funcpair([xdata[0][i], xdata[1][i]], popt[0], popt[1], popt[2], popt[3]))
            fy.append(funcpair([xdata[0][i], xdata[1][i]], 0, popt[1], popt[2], popt[3]))
            ffit.write(str(xdata[0][i]) + " " + str(xdata[1][i]) + " " + str(y[i]) + " " + str(ydata[i]) + " " + str(ydata[i] - y[i]) + " " + str(fy[i]) + " " + str(scost[i]) + "\n")
        x = xdata[0][:]
        solution = fopt[0] + " + " + fopt[1] + "* cicount" + " + " + fopt[2] + " * cjcount " + " + " + fopt[3] + " * cicount * cjcount"
        fsol.write( key + " " + fopt[1] + " " + fopt[2] + " " + fopt[3] + "\n")

    ffit.write("# Solution: " + solution + "\n")
    if maxcost > 0:
        ffit.write("# cost scale: " + str(ymax/maxcost) + "\n")
    else:
        ffit.write("# cost scale: 0\n")

    fig = pl.figure()
    ax = fig.add_subplot(1,1,1)
    pl.scatter(x, ydata, c="blue", label="data")
    pl.scatter(x, fy, c="green", label="model")
    pl.scatter(x, scost, c="cyan", label="costs")
    pl.scatter(x, y, c="red", label="fit")
    ax.set_title( "task.subtype.flags: " + key + "\n" + solution)
    ax.set_xlabel("ci count")
    ax.set_ylabel("ticks")

    ax.legend()

    pl.savefig(key + ".png")
    pl.show()
    pl.close("all")


    ffit.close()
    fdat.close()
fsol.close()
