from __future__ import print_function, division, absolute_import

import numpy as np

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

except RuntimeError:
    pass

from .general import *

from .. import utils
from .. import config


def plot_sfh(sfh, show=True, save=False, from_bigbang=False):
    """ Make a quick plot of an individual sfh. """

    update_rcParams()

    fig = plt.figure(figsize=(12, 4))
    ax = plt.subplot()

    add_sfh(sfh, ax, from_bigbang=from_bigbang)

    if save:
        plt.savefig("model_sfh.pdf", bbox_inches="tight")
        plt.close(fig)

    if show:
        plt.show()
        plt.close(fig)

    return fig, ax


def add_sfh(sfh, ax, from_bigbang=False, zorder=4, color="black", z_axis=True, lw=2, alpha=1, ls="-", label=None):
    """ Creates a plot of sfr(t) for a given star-formation history. """

    # Set limits and x
    if from_bigbang:
        x = (sfh.age_of_universe - sfh.ages)*10**-9
        ax.set_xlim(sfh.age_of_universe*10**-9, 0)
    else:
        x = (sfh.ages)*10**-9
        ax.set_xlim(0,sfh.age_of_universe*10**-9)
    
    # Plot the sfh
    ax.plot(x, sfh.sfh, color=color, zorder=zorder, lw=lw, alpha=alpha, ls=ls, label=label)
    ax.set_ylim(bottom=0.)

    # Add redshift axis along the top
    if z_axis:
        zvals = [0,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30]
        z_axis = add_z_axis(ax, zvals=zvals, from_bigbang=from_bigbang)

    # Add labels
    if tex_on:
        ax.set_ylabel("$\\mathrm{SFR\\ /\\ M_\\odot\\ \\mathrm{yr}^{-1}}$")
        if from_bigbang:
            ax.set_xlabel("$\\mathrm{Age\\ of\\ Universe\\ /\\ Gyr}$")
        else:            
            ax.set_xlabel("$\\mathrm{Lookback \\ time\\ /\\ Gyr}$")

    else:
        ax.set_ylabel("SFR / M_sol yr^-1")
        if from_bigbang:
            ax.set_xlabel("Age of Universe / Gyr")
        else:
            ax.set_xlabel("Lookback time / Gyr")
            


    return z_axis
