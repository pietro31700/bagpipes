from __future__ import print_function, division, absolute_import

import numpy as np

try:
    import corner
    import matplotlib.pyplot as plt

except RuntimeError:
    pass

from .general import *


def plot_corner(fit, show=False, save=True, bins=25, type="fit_params"):
    """ Make a corner plot of the fitted parameters. """

    update_rcParams()

    names = fit.fitted_model.params.copy()
    samples = np.copy(fit.posterior.samples2d)
    
    #print total mass formed for step function
    if "step" in fit.fit_instructions.keys():
        sfrs_indexes = [i for i in range(len(names)) if "sfr" in names[i]]
        bin_edges = fit.fit_instructions["step"]["bin_edges"]
        age_widths = np.diff(bin_edges)
        mass_formed = np.log10(np.sum(age_widths*1.e6*samples[:,sfrs_indexes],axis=1))
        
        names.append("step:massformed")
        samples = np.c_[samples,mass_formed]
        
    # Set up axis labels
    if tex_on:
        labels = fix_param_names(names)

    else:
        labels = names

    # Log any parameters with log_10 priors to make them easier to see
    for i in range(fit.fitted_model.ndim):
        if fit.fitted_model.pdfs[i] == "log_10":
            samples[:, i] = np.log10(samples[:, i])

            if tex_on:
                labels[i] = "$\\mathrm{log_{10}}(" + labels[i][1:-1] + ")$"

            else:
                labels[i] = "log_10(" + labels[i] + ")"

    # Replace any r parameters for Dirichlet distributions with t_x vals
    j = 0
    for i in range(fit.fitted_model.ndim):
        if "dirichlet" in fit.fitted_model.params[i]:
            comp = fit.fitted_model.params[i].split(":")[0]
            n_x = fit.fitted_model.model_components[comp]["bins"]
            t_percentile = int(np.round(100*(j+1)/n_x))

            samples[:, i] = fit.posterior.samples[comp + ":tx"][:, j]
            j += 1

            if tex_on:
                labels[i] = "$t_{" + str(t_percentile) + "}\\ /\\ \\mathrm{Gyr}$"

            else:
                labels[i] = "t" + str(t_percentile) + " / Gyr"

    # Make the corner plot (copied from corner package)
    factor = 2.0  # size of one side of one panel
    lbdim = 0.5 * factor  # size of left/bottom margin
    trdim = 0.2 * factor  # size of top/right margin
    whspace = 0.05  # w/hspace size
    plotdim = factor * len(labels) + factor * (len(labels) - 1.0) * whspace
    dim = lbdim + plotdim + trdim

    # Create a new figure if one wasn't provided.
    fig= plt.figure(figsize=(dim, dim), num=1, clear=True)
    fig = corner.corner(samples, labels=labels, quantiles=[0.16, 0.5, 0.84],
                        show_titles=True, title_kwargs={"fontsize": 13},
                        smooth=1., smooth1d=1., bins=bins, fig=fig)

    # Save the corner plot to file
    if save:
        plotpath = ("pipes/plots/" + fit.run + "/" + fit.galaxy.ID
                    + "_corner.pdf")

        plt.savefig(plotpath, bbox_inches="tight")
        fig.clear()
        plt.close(fig)

    # Alternatively show the corner plot
    if show:
        plt.show()
        fig.clear()
        plt.close(fig)

    return fig
