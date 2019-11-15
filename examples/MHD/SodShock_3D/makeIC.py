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

import h5py
import numpy as np

# Generates a swift IC file for the 3D Sod Shock in a periodic box

# Parameters
gamma = 2.          # Gas adiabatic index
x_min = -0.5
x_max = 0.5
rho_L = 1.             # Density left state
rho_R = 0.125          # Density right state
v_L = 0.               # Velocity left state
v_R = 0.               # Velocity right state
P_L = 1.               # Pressure left state
P_R = 0.1              # Pressure right state
By_L = 1.              # Magnetic field (y axis) left state
By_R = -1.             # Magnetic field (y axis) left state
Bz_L = 0.              # Magnetic field (z axis) left state
Bz_R = 0.              # Magnetic field (z axis) left state
Bx = 0.75
fileName = "sodShock.hdf5"


# ---------------------------------------------------
boxSize = (x_max - x_min)

glass_L = h5py.File("glassCube_64.hdf5", "r")
glass_R = h5py.File("glassCube_32.hdf5", "r")

pos_L = glass_L["/PartType0/Coordinates"][:, :]
pos_R = glass_R["/PartType0/Coordinates"][:, :]
h_L = glass_L["/PartType0/SmoothingLength"][:]
h_R = glass_R["/PartType0/SmoothingLength"][:]
ind_L = np.logical_and(pos_L[:, 1] < 0.5, pos_L[:, 2] < 0.5)
ind_L = np.logical_and(ind_L, pos_L[:, 0] < 0.5)
ind_R = np.logical_and(pos_R[:, 1] < 0.5, pos_R[:, 2] < 0.5)
ind_R = np.logical_and(ind_R, pos_R[:, 0] >= 0.5)
pos_L = pos_L[ind_L]
pos_R = pos_R[ind_R]
h_L = h_L[ind_L]
h_R = h_R[ind_R]

# Merge things
pos = np.append(pos_L, pos_R, axis=0)
h = np.append(h_L, h_R)

numPart_L = np.size(h_L)
numPart_R = np.size(h_R)
numPart = np.size(h)

vol_L = 0.5
vol_R = 0.5

# Generate extra arrays
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)
B = np.zeros((numPart, 3))
B[:, 0] = Bx

ind = pos[:, 0] < 0.5
u[ind] = P_L / (rho_L * (gamma - 1.))
m[ind] = rho_L * vol_L / numPart_L
v[ind, 0] = v_L
B[ind, 1] = By_L
B[ind, 2] = Bz_L

ind = pos[:, 0] >= 0.5
u[ind] = P_R / (rho_R * (gamma - 1.))
m[ind] = rho_R * vol_R / numPart_R
v[ind, 0] = v_R
B[ind, 1] = By_R
B[ind, 2] = Bz_R

# File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxSize, 0.5, 0.5]
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

# Particle group
grp = file.create_group("/PartType0")
grp.create_dataset('Coordinates', data=pos, dtype='d')
grp.create_dataset('Velocities', data=v, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('SmoothingLength', data=h, dtype='f')
grp.create_dataset('InternalEnergy', data=u, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')
grp.create_dataset("MagneticFields", data=B, dtype="f")


file.close()
