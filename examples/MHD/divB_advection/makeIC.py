#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

import numpy as np
from h5py import File
from copy import deepcopy

# Generates the initial condition for the divB cleaning test
# based on Tricco 2015 (Thesis)

# parameters

L = 2.  # Box size
N = 50  # number of particles per dimension
rho = 1.  # density
vx = 1.  # velocity along x axis
vy = 1.  # velocity along y axis
vz = 0.  # velocity along z axis
Bx = 1. / np.sqrt(4. * np.pi)  # magnetic field along x axis
By = 0.  # magnetic field along y axis
Bz = 1. / np.sqrt(4. * np.pi)  # magnetic field along z axis
p = 6.  # hydrodynamic pressure
gamma = 5. / 3.  # adiabatic index
r0 = 1. / np.sqrt(8.)  # radius of the initial perturbation


# Computation
def getB(x):
    r2 = np.sum(x**2, axis=1) / r0**2
    B = np.zeros(x.shape)
    ind = r2 < 1
    B[ind, 0] = Bx * (r2[ind]**4 - r2[ind]**2 + 1)
    B[~ind, 0] = 0.
    B[:, 1] = By
    B[:, 2] = Bz
    return B


def getInternalEnergy(B):
    B2 = np.sum(B**2, axis=1)
    v2 = vx**2 + vy**2 + vz**2
    return p / (gamma - 1) + 0.5 * rho * v2 + 0.5 * B2


def generate_square(num_on_side, side_length=1.0):
    """
    Generates a cube
    """

    values = np.linspace(0.0, side_length, num_on_side + 1)[:-1]

    positions = np.empty((num_on_side**2, 3), dtype=float)

    for x in range(num_on_side):
        for y in range(num_on_side):
            index = x + y * num_on_side

            positions[index, 0] = values[x]
            positions[index, 1] = values[y]
            positions[index, 2] = 0.5

    return positions


if __name__ == "__main__":
    x = generate_square(N, L)
    x -= 0.5

    # generate the arrays
    pos = deepcopy(x) + 0.5

    v = np.zeros(x.shape)
    v[:, 0] = vx
    v[:, 1] = vy
    v[:, 2] = vz

    m = np.ones(N**2) * rho * L**2 / N**2
    h = 2. * np.ones(N**2) * L / N
    B = getB(x)
    u = getInternalEnergy(B)
    ids = np.arange(N**2)

    f = File("divB.hdf5", "w")
    # Header
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [L, L, 0.]
    grp.attrs["NumPart_Total"] = [N**2, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [N**2, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFilesPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = 0
    grp.attrs["Dimension"] = 2

    # Units
    grp = f.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = 1.
    grp.attrs["Unit mass in cgs (U_M)"] = 1.
    grp.attrs["Unit time in cgs (U_t)"] = 1.
    grp.attrs["Unit current in cgs (U_I)"] = 1.
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.

    # Particle group
    grp = f.create_group("/PartType0")
    grp.create_dataset('Coordinates', data=pos, dtype='d')
    grp.create_dataset('Velocities', data=v, dtype='f')
    grp.create_dataset('Masses', data=m, dtype='f')
    grp.create_dataset('SmoothingLength', data=h, dtype='f')
    grp.create_dataset('InternalEnergy', data=u, dtype='f')
    grp.create_dataset('ParticleIDs', data=ids, dtype='L')
    grp.create_dataset("MagneticFields", data=B, dtype="f")


    f.close()
