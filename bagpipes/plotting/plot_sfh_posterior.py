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


def plot_sfh_posterior(fit, show=False, save=True, log_x=False, log_y=False, mean=False,from_bigbang=False, colorscheme="bw"):
    """ Make a plot of the SFH posterior. """

    update_rcParams()

    fig = plt.figure(figsize=(12, 4))
    ax = plt.subplot()

    add_sfh_posterior(fit, ax, log_x=log_x, log_y=log_y, mean=mean,colorscheme=colorscheme,from_bigbang=from_bigbang)

    if save:
        plotpath = "pipes/plots/" + fit.run + "/" + fit.galaxy.ID + "_sfh.pdf"            
        plt.savefig(plotpath, bbox_inches="tight")
        plt.close(fig)

    if show:
        plt.show()
        plt.close(fig)

    return fig, ax


def add_sfh_posterior(fit, ax, log_x=False, log_y=False, mean=False, from_bigbang=False, colorscheme="bw", z_axis=True, zorder=4,
                      label=None):

    color1 = "black"
    color2 = "gray"
    alpha = 0.6

    if colorscheme == "irnbru":
        color1 = "darkorange"
        color2 = "navajowhite"
        alpha = 0.6

    if colorscheme == "purple":
        color1 = "purple"
        color2 = "purple"
        alpha = 0.4

    if colorscheme == "blue":
        color1 = "dodgerblue"
        color2 = "dodgerblue"
        alpha = 0.7

    # Calculate median redshift and median age of Universe
    if "redshift" in fit.fitted_model.params:
        redshift = np.median(fit.posterior.samples["redshift"])

    else:
        redshift = fit.fitted_model.model_components["redshift"]

    age_of_universe = np.interp(redshift, utils.z_array, utils.age_at_z)

    # Calculate median and confidence interval for SFH posterior
    post = np.percentile(fit.posterior.samples["sfh"], (16, 50, 84), axis=0).T
    
    
    if log_x==True:
        ax.set_xscale("log")
        age_start = 1e-3
    else:
        age_start = 0
        
    # Plot the SFH
    if from_bigbang:
        x = age_of_universe - fit.posterior.sfh.ages*10**-9
        ax.set_xlim(age_of_universe, age_start)
    else:
        x = fit.posterior.sfh.ages*10**-9
        ax.set_xlim(age_start,age_of_universe)
        
        
    ax.plot(x, post[:, 1], color=color1, zorder=zorder+1,label="Median")
    ax.fill_between(x, post[:, 0], post[:, 2], color=color2,
                    alpha=alpha, zorder=zorder, lw=0, label=label)
    #if maxL ==True:
    #    max_likelihood_index = np.argmax(fit.results["lnlike"][fit.posterior.indices])
    #    maxl_sfh   = fit.posterior.samples["sfh"][max_likelihood_index]
    #    ax.plot(x, maxl_sfh, color="blue", alpha=0.5, zorder=zorder+1,label="MaxL")
        
    if mean ==True:
        mean_sfh = np.mean(fit.posterior.samples["sfh"],axis=0)
        
        #binning=5
        #len_binned = len(x)//binning
        #x_binned = np.zeros(len_binned)
        #y_binned = np.zeros(len_binned)
        #for ii in range(len_binned):
        #    x_binned[ii] = np.mean(x[binning*ii:binning*(ii+1)])      
        #    y_binned[ii] = np.mean(mean_sfh[binning*ii:binning*(ii+1)])      
        ax.plot(x,mean_sfh, color=color1, zorder=zorder+1,alpha = 0.5,label="Mean")
        
        
    if log_y==True:
        ax.set_yscale("log")
        ax.set_ylim(1e-2, np.max([ax.get_ylim()[1], 1.1*np.max(post[:, 2])]))
    else:
        ax.set_ylim(0., np.max([ax.get_ylim()[1], 1.1*np.max(post[:, 2])]))
        
    plt.legend(frameon=False)
    # Add redshift axis along the top
    if z_axis:
        zvals = [0,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30]
        ax2 = add_z_axis(ax, zvals=zvals, from_bigbang=from_bigbang, log_scale=log_x)

    # Set axis labels
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
            

    if z_axis:
        return ax2
