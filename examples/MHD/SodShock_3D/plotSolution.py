###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
#               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

from scipy import stats
import h5py
import matplotlib
import makeIC as ic
import sys
import numpy as np
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Computes the analytical solution of the Sod shock and plots the SPH answer

# Generates the analytical  solution for the Sod shock test case
# The script works for a given left (x<0.5) and right (x>0.5) state and computes the solution at a later time t.
# This follows the solution given in (Brio & Wu, 1988)

limit_axis = True

# Plot parameters
params = {
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'font.size': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.figsize': (9.90, 6.45),
    'figure.subplot.left': 0.045,
    'figure.subplot.right': 0.99,
    'figure.subplot.bottom': 0.05,
    'figure.subplot.top': 0.99,
    'figure.subplot.wspace': 0.15,
    'figure.subplot.hspace': 0.12,
    'lines.markersize': 6,
    'lines.linewidth': 3.,
    'text.latex.unicode': True
}
matplotlib.rcParams.update(params)
matplotlib.rc(
    'font', **{
        'family': 'sans-serif',
        'sans-serif': ['Times']
    })

snap = int(sys.argv[1])


# Read the simulation data
sim = h5py.File("sodShock_%04d.hdf5" % snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

x = sim["/PartType0/Coordinates"][:, 0]
vx = sim["/PartType0/Velocities"][:, 0]
vy = sim["/PartType0/Velocities"][:, 1]
u = sim["/PartType0/InternalEnergies"][:]
S = sim["/PartType0/Entropies"][:]
P = sim["/PartType0/Pressures"][:]
rho = sim["/PartType0/Densities"][:]
By = sim["/PartType0/MagneticFields"][:, 1]


x += ic.x_min
N = 1000

# Bin te data
x_bin_edge = np.arange(-0.6, 0.6, 0.02)
x_bin = 0.5*(x_bin_edge[1:] + x_bin_edge[:-1])
rho_bin, _, _ = stats.binned_statistic(x, rho, statistic='mean',
                                       bins=x_bin_edge)
vx_bin, _, _ = stats.binned_statistic(x, vx, statistic='mean', bins=x_bin_edge)
vy_bin, _, _ = stats.binned_statistic(x, vy, statistic='mean', bins=x_bin_edge)
P_bin, _, _ = stats.binned_statistic(x, P, statistic='mean', bins=x_bin_edge)
S_bin, _, _ = stats.binned_statistic(x, S, statistic='mean', bins=x_bin_edge)
u_bin, _, _ = stats.binned_statistic(x, u, statistic='mean', bins=x_bin_edge)
By_bin, _, _ = stats.binned_statistic(x, By, statistic='mean', bins=x_bin_edge)
rho2_bin, _, _ = stats.binned_statistic(x, rho**2, statistic='mean',
                                        bins=x_bin_edge)
vx2_bin, _, _ = stats.binned_statistic(x, vx**2, statistic='mean',
                                       bins=x_bin_edge)
vy2_bin, _, _ = stats.binned_statistic(x, vy**2, statistic='mean',
                                       bins=x_bin_edge)
P2_bin, _, _ = stats.binned_statistic(x, P**2, statistic='mean',
                                      bins=x_bin_edge)
S2_bin, _, _ = stats.binned_statistic(x, S**2, statistic='mean',
                                      bins=x_bin_edge)
u2_bin, _, _ = stats.binned_statistic(x, u**2, statistic='mean',
                                      bins=x_bin_edge)
By2_bin, _, _ = stats.binned_statistic(x, By**2, statistic='mean',
                                       bins=x_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin**2)
vx_sigma_bin = np.sqrt(vx2_bin - vx_bin**2)
vy_sigma_bin = np.sqrt(vy2_bin - vy_bin**2)
P_sigma_bin = np.sqrt(P2_bin - P_bin**2)
S_sigma_bin = np.sqrt(S2_bin - S_bin**2)
u_sigma_bin = np.sqrt(u2_bin - u_bin**2)
By_sigma_bin = np.sqrt(By2_bin - By_bin**2)


# get the solution from ramses
filename = "data.csv"
ramses = np.genfromtxt(filename, names=True)
x_s = ramses["x"] + ic.x_min
rho_s = ramses["d"]
vx_s = ramses["u"]
vy_s = ramses["v"]
P_s = ramses["P"]
u_s = P_s / (rho_s * (ic.gamma - 1.))
By_s = ramses["B"]


# Plot the interesting quantities
plt.figure()

# Density profile --------------------------------
plt.subplot(231)
plt.plot(x, rho, '.', color='r', ms=0.5, alpha=0.2)
plt.plot(x_s, rho_s, '--', color='k', alpha=0.8, lw=1.2)
plt.errorbar(x_bin, rho_bin, yerr=rho_sigma_bin,
             fmt='.', ms=8.0, color='b', lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Density}}~\\rho$", labelpad=0)
if limit_axis:
    plt.xlim(-0.5, 0.5)
    plt.ylim(0.05, 1.1)

# Velocity profile vx --------------------------------
plt.subplot(232)
plt.plot(x, vx, '.', color='r', ms=0.5, alpha=0.2)
plt.plot(x_s, vx_s, '--', color='k', alpha=0.8, lw=1.2)
plt.errorbar(x_bin, vx_bin, yerr=vx_sigma_bin,
             fmt='.', ms=8.0, color='b', lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Velocity}}~v_x$", labelpad=0)
if limit_axis:
    plt.xlim(-0.5, 0.5)
    plt.ylim(-0.5, 0.95)

# Velocity profile vy --------------------------------
plt.subplot(233)
plt.plot(x, vy, '.', color='r', ms=0.5, alpha=0.2)
plt.plot(x_s, vy_s, '--', color='k', alpha=0.8, lw=1.2)
plt.errorbar(x_bin, vy_bin, yerr=vy_sigma_bin,
             fmt='.', ms=8.0, color='b', lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Velocity}}~v_x$", labelpad=0)
if limit_axis:
    plt.xlim(-0.5, 0.5)
    plt.ylim(-1.5, 0.1)

# Internal energy profile -------------------------
plt.subplot(234)
plt.plot(x, u, '.', color='r', ms=0.5, alpha=0.2)
plt.plot(x_s, u_s, '--', color='k', alpha=0.8, lw=1.2)
plt.errorbar(x_bin, u_bin, yerr=u_sigma_bin,
             fmt='.', ms=8.0, color='b', lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)
if limit_axis:
    plt.xlim(-0.5, 0.5)
    plt.ylim(0.5, 2.5)

# Pressure profile --------------------------------
plt.subplot(235)
plt.plot(x, P, '.', color='r', ms=0.5, alpha=0.2)
plt.plot(x_s, P_s, '--', color='k', alpha=0.8, lw=1.2)
plt.errorbar(x_bin, P_bin, yerr=P_sigma_bin,
             fmt='.', ms=8.0, color='b', lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Pressure}}~P$", labelpad=0)
if limit_axis:
    plt.xlim(-0.5, 0.5)
    plt.ylim(0.01, 1.1)

# Magnetic field profile --------------------------------
plt.subplot(236)
plt.plot(x, By, '.', color='r', ms=0.5, alpha=0.2)
plt.plot(x_s, By_s, '--', color='k', alpha=0.8, lw=1.2)
plt.errorbar(x_bin, By_bin, yerr=By_sigma_bin,
             fmt='.', ms=8.0, color='b', lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Magnetic Field}}~B_y$", labelpad=0)
if limit_axis:
    plt.xlim(-0.5, 0.5)
    plt.ylim(-1.1, 1.1)

# # Information -------------------------------------
# plt.subplot(236, frameon=False)

# plt.text(-0.49, 0.9, "Sod shock with  $\\gamma=%.3f$ in 3D at $t=%.2f$"
#          % (ic.gamma, time), fontsize=10)
# plt.text(-0.49, 0.8, "Left:~~ $(P_L, \\rho_L, v_L) = (%.3f, %.3f, %.3f)$"
#          % (ic.P_L, ic.rho_L, ic.v_L), fontsize=10)
# plt.text(-0.49, 0.7, "Right: $(P_R, \\rho_R, v_R) = (%.3f, %.3f, %.3f)$"
#          % (ic.P_R, ic.rho_R, ic.v_R), fontsize=10)
# plt.plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
# plt.text(-0.49, 0.5, "$\\textsc{Swift}$ %s"
#          % git, fontsize=10)
# plt.text(-0.49, 0.4, scheme, fontsize=10)
# plt.text(-0.49, 0.3, kernel, fontsize=10)
# plt.text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"
#          % (neighbours, eta), fontsize=10)
# plt.xlim(-0.5, 0.5)
# plt.ylim(0, 1)
# plt.xticks([])
# plt.yticks([])


plt.savefig("SodShock.png", dpi=200)
