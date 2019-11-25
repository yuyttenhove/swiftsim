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

limit_axis = False

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

# Analytic solution
c_L = np.sqrt(ic.gamma * ic.P_L / ic.rho_L)   # Speed of the rarefaction wave
c_R = np.sqrt(ic.gamma * ic.P_R / ic.rho_R)   # Speed of the shock front

# Helpful variable
Gama = (ic.gamma - 1.) / (ic.gamma + 1.)
beta = (ic.gamma - 1.) / (2. * ic.gamma)


# Characteristic function and its derivative, following Toro (2009)
def compute_f(P_3, P, c):
    u = P_3 / P
    if u > 1:
        term1 = ic.gamma*((ic.gamma+1.)*u + ic.gamma-1.)
        term2 = np.sqrt(2./term1)
        fp = (u - 1.)*c*term2
        dfdp = c*term2/P + (u - 1.) * c/term2 * (-1./term1**2) \
            * ic.gamma * (ic.gamma + 1.) / P
    else:
        fp = (u**beta - 1.)*(2.*c/(ic.gamma-1.))
        dfdp = 2.*c/(ic.gamma-1.)*beta*u**(beta-1.)/P
    return (fp, dfdp)


# Solution of the Riemann problem following Toro (2009)
def RiemannProblem(rho_L, P_L, v_L, rho_R, P_R, v_R):
    P_new = ((c_L + c_R + (v_L - v_R) * 0.5 * (ic.gamma-1.))
             / (c_L / P_L**beta + c_R / P_R**beta))**(1./beta)
    P_3 = 0.5*(P_R + P_L)
    f_L = 1.
    while np.fabs(P_3 - P_new) > 1e-6:
        P_3 = P_new
        (f_L, dfdp_L) = compute_f(P_3, P_L, c_L)
        (f_R, dfdp_R) = compute_f(P_3, P_R, c_R)
        f = f_L + f_R + (v_R - v_L)
        df = dfdp_L + dfdp_R
        dp = -f/df
        P_new = P_3 + dp
    v_3 = v_L - f_L
    return (P_new, v_3)


# Solve Riemann problem for post-shock region
(P_3, v_3) = RiemannProblem(ic.rho_L, ic.P_L, ic.v_L,
                            ic.rho_R, ic.P_R, ic.v_R)

# Check direction of shocks and wave
shock_R = (P_3 > ic.P_R)
shock_L = (P_3 > ic.P_L)

# Velocity of shock front and and rarefaction wave
if shock_R:
    v_right = ic.v_R + c_R**2*(P_3/ic.P_R - 1.)/(ic.gamma*(v_3-ic.v_R))
else:
    v_right = c_R + 0.5*(ic.gamma+1.)*v_3 - 0.5*(ic.gamma-1.)*ic.v_R

if shock_L:
    v_left = ic.v_L + c_L**2*(P_3/ic.p_L - 1.)/(ic.gamma*(v_3-ic.v_L))
else:
    v_left = c_L - 0.5*(ic.gamma+1.)*v_3 + 0.5*(ic.gamma-1.)*ic.v_L

# Compute position of the transitions
x_23 = -np.fabs(v_left) * time
if shock_L:
    x_12 = -np.fabs(v_left) * time
else:
    x_12 = -(c_L - ic.v_L) * time

x_34 = v_3 * time

x_45 = np.fabs(v_right) * time
if shock_R:
    x_56 = np.fabs(v_right) * time
else:
    x_56 = (c_R + ic.v_R) * time


# Prepare arrays
delta_x = (ic.x_max - ic.x_min) / N
x_s = np.arange(ic.x_min, ic.x_max, delta_x)
rho_s = np.zeros(N)
P_s = np.zeros(N)
vx_s = np.zeros(N)
vy_s = np.zeros(N)

# Compute solution in the different regions
for i in range(N):
    if x_s[i] <= x_12:
        rho_s[i] = ic.rho_L
        P_s[i] = ic.P_L
        vx_s[i] = ic.v_L
    if x_s[i] >= x_12 and x_s[i] < x_23:
        if shock_L:
            rho_s[i] = ic.rho_L*(Gama + P_3/ic.P_L)/(1. + Gama * P_3/ic.P_L)
            P_s[i] = P_3
            vx_s[i] = v_3
        else:
            rho_s[i] = ic.rho_L*(Gama * (0. - x_s[i])/(c_L * time) + Gama
                                 * ic.v_L/c_L + (1.-Gama))**(2./(ic.gamma-1.))
            P_s[i] = ic.P_L*(rho_s[i] / ic.rho_L)**ic.gamma
            vx_s[i] = (1.-Gama)*(c_L - (0. - x_s[i]) / time) + Gama*ic.v_L
    if x_s[i] >= x_23 and x_s[i] < x_34:
        if shock_L:
            rho_s[i] = ic.rho_L*(Gama + P_3/ic.P_L)/(1+Gama * P_3/ic.p_L)
        else:
            rho_s[i] = ic.rho_L*(P_3 / ic.P_L)**(1./ic.gamma)
        P_s[i] = P_3
        vx_s[i] = v_3
    if x_s[i] >= x_34 and x_s[i] < x_45:
        if shock_R:
            rho_s[i] = ic.rho_R*(Gama + P_3/ic.P_R)/(1. + Gama * P_3/ic.P_R)
        else:
            rho_s[i] = ic.rho_R*(P_3 / ic.P_R)**(1./ic.gamma)
        P_s[i] = P_3
        vx_s[i] = v_3
    if x_s[i] >= x_45 and x_s[i] < x_56:
        if shock_R:
            rho_s[i] = ic.rho_R
            P_s[i] = ic.P_R
            vx_s[i] = ic.v_R
        else:
            rho_s[i] = ic.rho_R*(Gama*(x_s[i])/(c_R*time) - Gama
                                 * ic.v_R/c_R + (1.-Gama))**(2./(ic.gamma-1.))
            P_s[i] = ic.p_R*(rho_s[i]/ic.rho_R)**ic.gamma
            vx_s[i] = (1.-Gama)*(-c_R - (-x_s[i])/time) + Gama * ic.v_R
    if x_s[i] >= x_56:
        rho_s[i] = ic.rho_R
        P_s[i] = ic.P_R
        vx_s[i] = ic.v_R


# Additional arrays
u_s = P_s / (rho_s * (ic.gamma - 1.))  # internal energy
s_s = P_s / rho_s**ic.gamma  # entropic function

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
# plt.plot(x_s, P_s, '--', color='k', alpha=0.8, lw=1.2)
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
