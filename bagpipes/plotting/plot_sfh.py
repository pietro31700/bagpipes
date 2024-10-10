from __future__ import print_function, division, absolute_import

import numpy as np

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

except RuntimeError:
    pass

from .general import *


def plot_sfh(sfh, show=True, save=False, from_bigbang=False, log_x=False, log_y=False):
    """ Make a quick plot of an individual sfh. """

    update_rcParams()

    fig = plt.figure(figsize=(12, 4))
    ax = plt.subplot()

    add_sfh(sfh, ax, from_bigbang=from_bigbang, log_x=log_x, log_y=log_y)

    if save:
        plt.savefig("model_sfh.pdf", bbox_inches="tight")
        plt.close(fig)

    if show:
        plt.show()
        plt.close(fig)

    return fig, ax


def add_sfh(sfh, ax, from_bigbang=False, zorder=4, color="black", z_axis=True, lw=2, alpha=1, ls="-", label=None, log_x=False, log_y=False):
    """ Creates a plot of sfr(t) for a given star-formation history. """
    
    if log_x==True:
        ax.set_xscale("log")
        age_start = 1e-3
    else:
        age_start = 0
        
    # Set limits and x
    if from_bigbang:
        x = (sfh.age_of_universe - sfh.ages)*10**-9
        ax.set_xlim(sfh.age_of_universe*10**-9, age_start)
    else:
        x = (sfh.ages)*10**-9
        ax.set_xlim(age_start,sfh.age_of_universe*10**-9)
    
    # Plot the sfh
    ax.plot(x, sfh.sfh, color=color, zorder=zorder, lw=lw, alpha=alpha, ls=ls, label=label)
    
    if log_y==True:
        ax.set_yscale("log")
        ax.set_ylim(bottom=1.e-2)
    else:
        ax.set_ylim(bottom=0.)

    # Add redshift axis along the top
    if z_axis:
        zvals = [0,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30]
        z_axis = add_z_axis(ax, zvals=zvals, from_bigbang=from_bigbang,log_scale=log_x)
        
    if not label==None:
        ax.legend(frameon=False)
        
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
        
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
