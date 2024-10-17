from __future__ import print_function, division, absolute_import

import numpy as np

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

except RuntimeError:
    pass

from .general import *
from .plot_galaxy import plot_galaxy


def plot_spectrum_posterior(fit, show=False, save=True):
    """ Plot the observational data and posterior from a fit object. """

    fit.posterior.get_advanced_quantities()

    update_rcParams()

    # First plot the observational data
    fig, ax, y_scale = plot_galaxy(fit.galaxy, show=False, return_y_scale=True)

    if fit.galaxy.spectrum_exists:
        add_spectrum_posterior(fit, ax[0], zorder=6, y_scale=y_scale[0])

    if fit.galaxy.photometry_exists:
        add_photometry_posterior(fit, ax[-1], zorder=2, y_scale=y_scale[-1])

    if save:
        plotpath = "pipes/plots/" + fit.run + "/" + fit.galaxy.ID + "_fit.pdf"
        plt.savefig(plotpath, bbox_inches="tight")
        plt.close(fig)

    if show:
        plt.show()
        plt.close(fig)

    return fig, ax


def add_photometry_posterior(fit, ax, zorder=4, y_scale=None, color1=None,
                             color2=None, skip_no_obs=False,
                             background_spectrum=True, label=None):

    if color1 == None:
        color1 = "navajowhite"

    if color2 == None:
        color2 = "darkorange"

    mask = (fit.galaxy.photometry[:, 1] > 0.)
    upper_lims = fit.galaxy.photometry[:, 1] + fit.galaxy.photometry[:, 2]
    ymax = 1.05*np.max(upper_lims[mask])

    if not y_scale:
        y_scale = float(int(np.log10(ymax))-1)

    # Calculate posterior median redshift.
    if "redshift" in fit.fitted_model.params:
        redshift = np.median(fit.posterior.samples["redshift"])

    else:
        redshift = fit.fitted_model.model_components["redshift"]

    # Plot the posterior photometry and full spectrum.
    log_wavs = np.log10(fit.posterior.model_galaxy.wavelengths*(1.+redshift))
    log_eff_wavs = np.log10(fit.galaxy.filter_set.eff_wavs)

    if background_spectrum:
        spec_post = np.percentile(fit.posterior.samples["spectrum_full"],
                                  (16, 84), axis=0).T*10**-y_scale

        spec_post = spec_post.astype(float)  # fixes weird isfinite error

        ax.plot(log_wavs, spec_post[:, 0], color=color1,
                zorder=zorder-1, label=label)

        ax.plot(log_wavs, spec_post[:, 1], color=color1,
                zorder=zorder-1)

        ax.fill_between(log_wavs, spec_post[:, 0], spec_post[:, 1],
                        zorder=zorder-1, color=color1, linewidth=0)

    phot_post = np.percentile(fit.posterior.samples["photometry"],
                              (16, 84), axis=0).T

    for j in range(fit.galaxy.photometry.shape[0]):

        if skip_no_obs and fit.galaxy.photometry[j, 1] == 0.:
            continue

        phot_band = fit.posterior.samples["photometry"][:, j]
        mask = (phot_band > phot_post[j, 0]) & (phot_band < phot_post[j, 1])
        phot_1sig = phot_band[mask]*10**-y_scale
        wav_array = np.zeros(phot_1sig.shape[0]) + log_eff_wavs[j]

        if phot_1sig.min() < ymax*10**-y_scale:
            ax.scatter(wav_array, phot_1sig, color=color2,
                       zorder=zorder, alpha=0.05, s=100, rasterized=True)

def add_spectrum_posterior(fit, ax, zorder=4, y_scale=None):

    ymax = 1.05*np.max(fit.galaxy.spectrum[:, 1])

    if not y_scale:
        y_scale = float(int(np.log10(ymax))-1)

    wavs = fit.galaxy.spectrum[:, 0]
    spec_post = np.copy(fit.posterior.samples["spectrum"])

    if "calib" in list(fit.posterior.samples):
        spec_post /= fit.posterior.samples["calib"]

    if "noise" in list(fit.posterior.samples):
        spec_post += fit.posterior.samples["noise"]

    post = np.percentile(spec_post, (16, 50, 84), axis=0).T*10**-y_scale

    ax.plot(wavs, post[:, 1], color="sandybrown", zorder=zorder, lw=1.5)
    ax.fill_between(wavs, post[:, 0], post[:, 2], color="sandybrown",
                    zorder=zorder, alpha=0.75, linewidth=0)


def plot_spectrum_beautiful(fit, line_close_up =False, show=False, save=True):
    fit.posterior.get_advanced_quantities()

    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({"font.size":17})

    wavs = fit.galaxy.spectrum[:, 0]
    spec_post = np.copy(fit.posterior.samples["spectrum"])

    if "calib" in list(fit.posterior.samples):
        spec_post /= fit.posterior.samples["calib"]

    if "noise" in list(fit.posterior.samples):
        spec_post += fit.posterior.samples["noise"]

    # Calculate median redshift
    if "redshift" in fit.fitted_model.params:
        redshift = np.median(fit.posterior.samples["redshift"])

    else:
        redshift = fit.fitted_model.model_components["redshift"]

    
    post = np.percentile(spec_post, (16, 50, 84), axis=0).T

    ymax = 1.2*np.max(post[:, 1])
    y_scale = float(int(np.log10(ymax))-1)

    colorDATA = "#000000"
    colorERROR = "#8ebad9"
    colorFIT = "#ed2024"

    fig, axes = plt.subplots(2,1,figsize=(9, 5.5),sharex=True,height_ratios=(3,1))
    fig.subplots_adjust(hspace=0)
    plt.subplot(211)
    for jj in range(2):
        axes[jj].yaxis.set_ticks_position('both')
        axes[jj].xaxis.set_ticks_position('both')
        axes[jj].tick_params(axis="both", direction='in', which="major",length=5,zorder=10)
        axes[jj].tick_params(axis="both", direction='in', which="minor",length=2,zorder=10)

    l_data1, = axes[0].step(fit.galaxy.spectrum[:,0],fit.galaxy.spectrum[:,1]/10**y_scale, label="Data",color=colorDATA,where="mid",alpha=1,lw=1)
    l_data2 = axes[0].fill_between(fit.galaxy.spectrum[:,0],(fit.galaxy.spectrum[:,1]-fit.galaxy.spectrum[:,2])/10**y_scale,(fit.galaxy.spectrum[:,1]+fit.galaxy.spectrum[:,2])/10**y_scale,color=colorERROR,alpha=1,step="mid",lw=0)
    
    l_fit1, = axes[0].step(fit.galaxy.spectrum[:,0], post[:, 1]/10**y_scale, color=colorFIT, label="SED fitting",alpha=0.9,lw=1.5,where="mid")
    
    legend  = axes[0].legend([(l_data2, l_data1), l_fit1], ["NIRSpec data", "SED best-fit"],frameon=True,loc="upper left",fontsize=14)
    legend.get_frame().set_linewidth(0.0)
    legend.get_frame().set_alpha(0.4)

    lines = {
        "Civ": {"wl0": 1549, "tex": r"$\textsc{C\,iv}$"},
        "Ciii": {"wl0": 1908, "tex": r"$\textsc{C\,iii]}$"},
        "Mgii": {"wl0": 2798, "tex": r"$\textsc{M{\rm g}\,ii}$","offset":0.1},
        "Oii": {"wl0": 3727, "tex": r"$\textsc{[O\,ii]}$"},
        "Hd": {"wl0": 4101.2, "tex": r"$\rm H\delta$", "offset": 0.2},
        "Hgamma": {"wl0": 4340, "tex": r"$\rm H\gamma$", "offset": 0.1},
        "Hbeta": {"wl0": 4861.3, "tex": r"$\rm H\beta$"},
        "OIIIb": {"wl0": 4958.9, "tex": r"$\textsc{[O\,iii]}$"},
        "OIIIa": {"wl0": 5006.8, "tex": r"$\textsc{[O\,iii]}$"},
        "Halpha": {"wl0": 6562.8, "tex": r"$\rm H\alpha$"},
        "SIIa": {"wl0": 6720, "tex": r"$\textsc{[S\,ii]}$"},
        "SIIb": {"wl0": 6733, "tex": r""},
        "SIIIa": {"wl0": 9075, "tex": r"$\textsc{[S\,iii]}$"},
        "SIIIb": {"wl0": 9535, "tex": r"$\textsc{[S\,iii]}$"},
        "HeI": {"wl0": 10836, "tex": r"$\textsc{H{\rm e}\,i}$"},
        "Pagamma": {"wl0": 10944, "tex": r"$\rm Pa\gamma$", "offset": 0.3},
        "Feiia": {"wl0": 12574, "tex": r"$\textsc{[F{\rm e}\,ii]}$"},
        "Pabeta": {"wl0": 12824, "tex": r"$\rm Pa\beta$"},
        "Feiib": {"wl0": 16444, "tex": r"$\textsc{[F{\rm e}\,ii]}$"},
        "Paalpha": {"wl0": 18760, "tex": r"$\rm Pa\alpha$"},
    }

    if line_close_up==False:
        axes[0].set_xlim(fit.galaxy.spectrum[0,0],fit.galaxy.spectrum[-1,0])
    else:

        wl_c = lines[line_close_up]["wl0"]* (1+redshift)
        lims_x = (max(wl_c*0.75,fit.galaxy.spectrum[0,0]),min(wl_c*1.25,fit.galaxy.spectrum[-1,0]))
        in_x_lim = (fit.galaxy.spectrum[:,0] > lims_x[0]) *  (fit.galaxy.spectrum[:,0] < lims_x[1])
        axes[0].set_xlim(lims_x)
        y_max_spectrum  = np.maximum(post[:, 1],fit.galaxy.spectrum[:,1])/10**y_scale
        axes[0].set_ylim( 0.9*y_max_spectrum[in_x_lim].min(), 1.2*y_max_spectrum[in_x_lim].max()) 


    lines = list(lines.values())
    means_p = [lines[ii]["wl0"] * (1+redshift) for ii in range(len(lines))]

    heigh_name_line = [0.1+((fit.galaxy.spectrum[:,1]/10**y_scale)[np.argmin(np.abs(means_p[jj]-fit.galaxy.spectrum[:,0]))]) for jj in range(len(means_p))]
    for ii in range(len(lines)):
        if means_p[ii]>1.05*fit.galaxy.spectrum[2,0] and  means_p[ii]<0.95*fit.galaxy.spectrum[-1,0]:
            axes[0].axvline(means_p[ii],ls="--",color="k",alpha=0.5,linewidth=1)
            if "offset" in lines[ii].keys():
                position = heigh_name_line[ii] + lines[ii]["offset"]*(max(heigh_name_line)-min(heigh_name_line))
            else:
                position = heigh_name_line[ii]
            t = axes[0].text(means_p[ii],position,lines[ii]["tex"],rotation='vertical',fontsize=13,ha='right',va="bottom",clip_on=True)#,clip_box = Bbox([[0,0], [axes[0].get_xlim()[1],1.2*axes[0].get_ylim()[1]]]))
            t.set_bbox(dict(facecolor='white', alpha=0.5,lw=0))




    ylim_now = axes[0].get_ylim()
    axes[0].set_ylim(top = ymax/10**y_scale)
    if ylim_now[0]<0:
        axes[0].set_ylim(bottom = 0)


    axes[1].axhline(0,ls="--",alpha=0.5,color="k")
    relative_errors = (fit.galaxy.spectrum[:,1]-post[:, 1])/fit.galaxy.spectrum[:,2]
    axes[1].set_ylabel(r"$\chi$",fontsize=15)
    axes[1].scatter(fit.galaxy.spectrum[:,0], relative_errors,color="k",marker=".",s=5,alpha=0.5)
    if not line_close_up:
        y_lim = np.percentile(relative_errors,[1,99])
        if y_lim[0]>-5 and y_lim[1]<5:
            axes[1].set_yticks([-3,0,3])
            axes[1].set_ylim((-5,5))
        else:
            axes[1].set_ylim(y_lim)
    else:
        axes[1].plot(fit.galaxy.spectrum[:,0], relative_errors,color=colorERROR,alpha=0.25,lw=1)
        axes[1].set_ylim(1.1*relative_errors[in_x_lim].min(),1.1*relative_errors[in_x_lim].max())

    axes[0].set_ylabel("$\\mathrm{f_{\\lambda}}\\ \\mathrm{/\\ 10^{"
                                + str(int(y_scale))
                                + "}\\ erg\\ s^{-1}\\ cm^{-2}\\ \\AA^{-1}}$",fontsize=15)
    axes[0].set_xlabel(r"$\lambda\,\rm[\AA]$")

    if save:
        if not line_close_up:
            plotpath = "pipes/plots/" + fit.run + "/" + fit.galaxy.ID + "_fit_b.pdf"
        else:
            plotpath = "pipes/plots/" + fit.run + "/" + fit.galaxy.ID + "_fit_b_" + line_close_up+ ".pdf"
        fig.savefig(plotpath, bbox_inches="tight")
        fig.clear()
        plt.close(fig)

    if show:
        plt.show()
        fig.clear()
        plt.close(fig)

    return fig, axes