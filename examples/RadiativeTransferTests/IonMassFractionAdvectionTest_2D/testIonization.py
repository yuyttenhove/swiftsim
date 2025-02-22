#!/usr/bin/env python3

# ----------------------------------------------------
# Check that the total amount of ionized species
# remains constant
# ----------------------------------------------------

import sys
import os
import swiftsimio
import numpy as np
import gc
import unyt
from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import SymLogNorm

snapshot_base = "output"


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshot(s) that are to be plotted 
    and return their names as list
    """

    snaplist = []

    dirlist = os.listdir()
    for f in dirlist:
        if f.startswith(snapshot_basename) and f.endswith("hdf5"):
            snaplist.append(f)

    snaplist = sorted(snaplist)

    return snaplist


def compare_data(snaplist):
    """
    Create and save the plot
    """

    HI = []
    HII = []
    HeI = []
    HeII = []
    HeIII = []

    for filename in snaplist:
        data = swiftsimio.load(filename)

        mXHI = data.gas.ion_mass_fractions.HI * data.gas.masses
        mXHII = data.gas.ion_mass_fractions.HII * data.gas.masses
        mXHeI = data.gas.ion_mass_fractions.HeI * data.gas.masses
        mXHeII = data.gas.ion_mass_fractions.HeII * data.gas.masses
        mXHeIII = data.gas.ion_mass_fractions.HeIII * data.gas.masses

        HI.append(mXHI.sum())
        HII.append(mXHII.sum())
        HeI.append(mXHeI.sum())
        HeII.append(mXHeII.sum())
        HeIII.append(mXHeIII.sum())

    plt.figure()
    plt.plot(range(len(snaplist)), HI, label="HI total mass")
    plt.plot(range(len(snaplist)), HII, label="HII total mass")
    plt.plot(range(len(snaplist)), HeI, label="HeI total mass")
    plt.plot(range(len(snaplist)), HeII, label="HeII total mass")
    plt.plot(range(len(snaplist)), HeIII, label="HeIII total mass")
    plt.legend()

    #  plt.show()
    plt.tight_layout()
    plt.savefig("total_abundancies.png", dpi=200)

    return


if __name__ == "__main__":

    snaplist = get_snapshot_list(snapshot_base)
    compare_data(snaplist)
