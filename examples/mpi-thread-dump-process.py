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
import sys

#  Our fitting function.
def funcpair(x, a, b, c, d):
    return a + b*x[0]+ c*x[1] + d*x[0]*x[1]

def funcself(x, a, b, c):
    return a + b*x + c*x*x

def iszero(a):
    return 0

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

#  Tasks and subtypes. Indexed as in tasks.h.
TASKTYPES = ["none",
             "sort",
             "self",
             "pair",
             "sub_self",
             "sub_pair",
             "init",
             "ghost",
             "extra_ghost",
             "drift",
             "kick1",
             "kick2",
             "timestep",
             "send",
             "recv",
             "grav_gather_m",
             "grav_fft",
             "grav_mm",
             "grav_up",
             "grav_external",
             "cooling",
             "sourceterms",
             "count"]

SUBTYPES = ["none",
            "density",
            "gradient",
            "force",
            "grav",
            "external_grav",
            "tend",
            "xv",
            "rho",
            "gpart",
            "spart",
            "count"]

#  Sorting directions in flags. Should give same results (0 to 12). See
#  sid_scale in scheduler.
SAMEFLAGS = [1, 2, 1, 2, 3, 2, 1, 2, 1, 2, 3, 2, 3]

#  Handle command-line.
if len( sys.argv ) != 2:
    print "Usage: ", sys.argv[0], "input.dat"
    sys.exit(1)

#  Read input.
indata = pl.loadtxt( sys.argv[1] )

num_lines = pl.size(indata) / len(indata[0])
print "num_lines: ", num_lines

#  Create hash of unique keys "type.subtype.flags".
keys = {}
for line in range(num_lines):
    ttype = TASKTYPES[int(indata[line,3])]
    subtype = SUBTYPES[int(indata[line,4])]
    flags = int(indata[line,12])
    if flags > 12:
	flags = 12

    # Flags is a tag for these.
    if ttype == "recv" or ttype == "send":
	flags = 0

    # All sorts are gathered together (multiple sorting directions).
    if ttype == "sort":
        flags = 0

    flags = SAMEFLAGS[flags]
    key = ttype + "." + subtype + "." + str(flags)
    if not key in keys:
        keys[key] = []
    keys[key].append(indata[line,:])

print "number of keys: ", len(keys)

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
    if ttype == "none":
        continue

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
        cost.append(line[13])
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

    # Self data has different fit.
    if min(cjdata) == 0:
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
    solution = "Solution: "
    if isself:
        for i in range(num_lines):
            y.append(funcself(xdata[i], popt[0], popt[1], popt[2]))
            fy.append( funcself(xdata[i], 0, popt[1], popt[2]))
            ffit.write(str(xdata[i]) + " " + str(y[i]) + " " + str(ydata[i]) + " " + str(ydata[i] - y[i]) + " " + str(fy[i]) + " " + str(scost[i]) + "\n")
        x = xdata
        solution = solution + " " + fopt[0] + " + " + fopt[1] + " * cicount " + " + " + fopt[2] + " * cicount * cicount"

    else:
        for i in range(num_lines):
            y.append(funcpair([xdata[0][i], xdata[1][i]], popt[0], popt[1], popt[2], popt[3]))
            fy.append(funcpair([xdata[0][i], xdata[1][i]], 0, popt[1], popt[2], popt[3]))
            ffit.write(str(xdata[0][i]) + " " + str(xdata[1][i]) + " " + str(y[i]) + " " + str(ydata[i]) + " " + str(ydata[i] - y[i]) + " " + str(fy[i]) + " " + str(scost[i]) + "\n")
        x = xdata[0][:]
        solution = solution + " " + fopt[0] + " + " + fopt[1] + "* cicount" + " + " + fopt[2] + " * cjcount " + " + " + fopt[3] + " * cicount * cjcount"
    ffit.write("# " + solution + "\n")
    if maxcost > 0:
        ffit.write("# cost scale: " + str(ymax/maxcost) + "\n")
    else:
        ffit.write("# cost scale: 0\n")


    fig = pl.figure()
    ax = fig.add_subplot(1,1,1)
    pl.scatter(x, ydata, c="blue")
    pl.scatter(x, fy, c="green")
    pl.scatter(x, scost, c="cyan")
    pl.scatter(x, y, c="red")
    ax.set_title( "task.subtype.flags: " + key + "\n" + solution)
    ax.set_xlabel("ci count")
    ax.set_ylabel("ticks")

    pl.savefig(key + ".png")
    pl.show()


    ffit.close()
    fdat.close()
