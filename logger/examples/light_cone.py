#!/usr/bin/env python3
"""
Produce a light cone image from the logger file.
Example: python3 light_cone.py ../../examples/SedovBlast_3D/index_*dump
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import z_at_value
import unyt
import argparse
from swiftsimio.visualisation.projection import backends
from swiftsimio.visualisation.smoothing_length_generation import generate_smoothing_lengths
from scipy.interpolate import interp1d
sys.path.append("../.libs/")

import liblogger as logger

scatter_backend = backends["fast"]
kernel_gamma = 1.825742
res = 1080
save_file = "positions.npy"

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Produce a light cone image from the logger')

    default_files = "../../examples/SmallCosmoVolume/SmallCosmoVolume_DM/index_0000.dump"

    parser.add_argument("-n", '--number_steps', dest='n',
                        type=int, default=1000,
                        help='Number of slice used for the light cone')
    parser.add_argument("-l", '--length_unit', dest='l',
                        type=str, default="Mpc",
                        help='Length unit (unyt)')
    parser.add_argument("-v", '--velocity_unit', dest='v',
                        type=str, default="km / s",
                        help='Velocity unit (unyt)')
    parser.add_argument("-b", '--box_size', dest='b',
                        type=float, default=142.24,
                        help='Box size of the simulation')
    parser.add_argument("--evolution", dest='e',
                        action="store_true",
                        help='Disable projection along time and project along positions')
    parser.add_argument("-a", '--scale_factor', dest='a',
                        type=float, default=1.,
                        help='Maximal scale factor for the light cone')
    parser.add_argument('files', metavar='filenames', type=str, nargs="*",
                        help='The filenames of the logfiles')

    args = parser.parse_args()
    if len(args.files) == 0:
        args.files = default_files
    return args


# Read the arguments
args = parse_arguments()
print("basename: %s" % args.files)

# check the input
filename = args.files
if filename.endswith(".dump"):
    filename = filename[:-5]
else:
    raise Exception("It seems that you are not providing a logfile (.dump)")

# Compute the constants
c = unyt.c
u_l = unyt.unyt_quantity(1., args.l)
u_v = unyt.unyt_quantity(1., args.v)

def z(a):
    """
    Convert a scale factor into a redshift
    """
    return 1. / a - 1.


def get_min_scale_factor():
    """
    Compute the minimal scale factor from the speed of light,
    the box size and the maximal scale factor.
    """
    t = args.b * u_l / c
    t_now = cosmo.age(z(args.a))
    t = t_now - t.to_astropy()
    z_max = z_at_value(cosmo.age, t)
    return 1. / (z_max + 1.)

# Compute the maximal redshift
a_min = get_min_scale_factor()

def get_slice(a0, a1):
    """
    Compute the positions of a slice from the interval in scale factors.
    """
    z0 = z(a0)
    z1 = z(a1)
    z_max = z(a_min)

    t0 = cosmo.age(z0)
    t1 = cosmo.age(z1)
    t_min = cosmo.age(z_max)

    d0 = (t0 - t_min) * c.to_astropy()
    d1 = (t1 - t_min) * c.to_astropy()
    return d0.to(args.l).value, d1.to(args.l).value


# Create the list of output time
scale_factors = np.linspace(a_min, args.a, args.n + 1)

# Read the positions
if not os.path.isfile(save_file):
    # read the logger
    positions = np.empty((0, 3), dtype=np.float32)
    with logger.Reader(filename, verbose=0) as reader:
        a0, a1 = reader.get_time_limits()

        # Check the time interval
        if a0 > a_min or args.a > a1:
            raise Exception("Cannot process the required redshift.")

        # Generate the slices
        for i in range(args.n):
            # This could be improved by searching the exact time of the intersection
            # between the position and the light cone.
            a = 0.5 * (scale_factors[i] + scale_factors[i+1])
            pos, = reader.get_particle_data(
                ["Coordinates"], a)

            # Get the particles inside the slice
            d0, d1 = get_slice(scale_factors[i], scale_factors[i+1])
            ind = np.logical_and(pos[:, 2] > d0, pos[:, 2] <= d1)
            positions = np.append(positions, pos[ind, :], axis=0)

    # Save the file for later
    np.save(save_file, positions)
else:
    # Read the file
    print("Restoring positions from file")
    positions = np.load(save_file)

def periodic(pos):
    """
    Apply the periodic boundary conditions.
    """
    print("Applying periodic condition")
    for i in range(3):
        ind = pos[:, i] < 0
        pos[ind, i] += args.b
        ind = pos[:, i] > args.b
        pos[ind, i] -= args.b
    return pos

def change_ticks():
    """
    Change the x-ticks from position to redshift
    """
    ticks = plt.gca().get_xticks()
    x = np.linspace(0, args.b, args.n + 1)
    f = interp1d(x, scale_factors, fill_value="extrapolate")
    ticks = f(ticks)
    ticks = z(ticks)
    ticks = ["%.2f" % i for i in ticks]
    plt.gca().set_xticklabels(ticks)

print("Reading done")
# Compute secondary variables
positions = periodic(positions)
h = generate_smoothing_lengths(
    positions * u_l, args.b * u_l, kernel_gamma)
m = np.ones(positions.shape[0])

# Apply transformation
positions /= args.b
h /= args.b

# Make the image
if args.e:
    img = scatter_backend(positions[:, 1], positions[:, 2],
                          m, h, res)
else:
    img = scatter_backend(positions[:, 0], positions[:, 1],
                          m, h, res)

# Apply gamma
img **= 1. / 2.2

# Plot everything
plt.imshow(img, interpolation="hermite", extent=[0, args.b, 0, args.b],
           cmap="plasma")

if args.e:
    plt.xlabel("Redshift")
    plt.ylabel("Position [Mpc]")
    change_ticks()
else:
    plt.xlabel("Position [Mpc]")
    plt.ylabel("Position [Mpc]")


plt.show()
