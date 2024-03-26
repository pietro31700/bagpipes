from __future__ import print_function, division, absolute_import

import numpy as np

from .. import config
from .. import utils


class nebular(object):
    """Allows access to and maniuplation of nebular emission models.
    These must be pre-computed using Cloudy and the relevant set of
    stellar emission models. This has already been done for the default
    stellar models.

    Parameters
    ----------

    wavelengths : np.ndarray
        1D array of wavelength values desired for the stellar models.
    """

    def __init__(self, wavelengths, velshift):
        self.wavelengths = wavelengths
        self.velshift = velshift
        self.combined_grid, self.line_grid = self._setup_grids()

    def _setup_grids(self):
        """Loads Cloudy nebular continuum grid and resamples to the
        input wavelengths. Loads nebular line grids and adds line fluxes
        to the correct pixels in order to create a combined grid."""

        comb_grid = np.zeros(
            (
                self.wavelengths.shape[0],
                config.metallicities.shape[0],
                config.logU.shape[0],
                config.densities.shape[0],
                config.neb_ages.shape[0],
            )
        )

        line_grid = np.zeros(
            (
                config.line_wavs.shape[0],
                config.metallicities.shape[0],
                config.logU.shape[0],
                config.densities.shape[0],
                config.neb_ages.shape[0],
            )
        )

        for i in range(config.metallicities.shape[0]):
            for j in range(config.logU.shape[0]):
                for l in range(config.densities.shape[0]):

                    hdu_index = (
                        config.metallicities.shape[0] * config.logU.shape[0] * l
                        + config.metallicities.shape[0] * j
                        + i
                        + 1
                    )
                    raw_cont_grid = config.cont_grid[hdu_index].data
                    raw_line_grid = config.line_grid[hdu_index].data

                    line_grid[:, i, j, l, :] = raw_line_grid[1:, 1:].T

                    for k in range(config.neb_ages.shape[0]):
                        comb_grid[:, i, j, l, k] = np.interp(
                            self.wavelengths,
                            config.neb_wavs,
                            raw_cont_grid[k + 1, 1:],
                            left=0,
                            right=0,
                        )

        # Add the nebular lines to the resampled nebular continuum grid.
        for i in range(config.line_wavs.shape[0]):
            line_wav_shift = config.line_wavs[i] * (1 + (self.velshift / (3 * 10**5)))
            ind = np.abs(self.wavelengths - line_wav_shift).argmin()
            if ind != 0 and ind != self.wavelengths.shape[0] - 1:
                width = (self.wavelengths[ind + 1] - self.wavelengths[ind - 1]) / 2
                comb_grid[ind, :, :, :, :] += line_grid[i, :, :, :, :] / width

        return comb_grid, line_grid

    def spectrum(self, sfh_ceh, t_bc, logU, n_e):
        """Obtain a 1D spectrum for a given star-formation and
        chemical enrichment history, ionization parameter and t_bc.

        parameters
        ----------

        sfh_ceh : numpy.ndarray
            2D array containing the desired star-formation and
            chemical evolution history.

        logU : float
            Log10 of the ionization parameter.

        n_e  : float
            Electron density in HII regions

        t_bc : float
            The maximum age at which to include nebular emission.
        """

        return self._interpolate_grid(self.combined_grid, sfh_ceh, t_bc, logU, n_e)

    def line_fluxes(self, sfh_ceh, t_bc, logU, n_e):
        """Obtain line fluxes for a given star-formation and
        chemical enrichment history, ionization parameter and t_bc.

        parameters
        ----------

        sfh_ceh : numpy.ndarray
            2D array containing the desired star-formation and
            chemical evolution history.

        logU : float
            Log10 of the ionization parameter.

        n_e  : float
            Electron density in HII regions

        t_bc : float
            The maximum age at which to include nebular emission.
        """

        return self._interpolate_grid(self.line_grid, sfh_ceh, t_bc, logU, n_e)

    def _interpolate_grid(self, grid, sfh_ceh, t_bc, logU, n_e):
        """Interpolates a chosen grid in logU and collapses over star-
        formation and chemical enrichment history to get 1D models."""

        t_bc *= 10**9

        if logU == config.logU[0]:
            logU += 10**-10
        if n_e == config.densities[0]:
            n_e += 10**-5

        spectrum_low_logU_low_density = np.zeros_like(grid[:, 0, 0, 0, 0])
        spectrum_high_logU_low_density = np.zeros_like(grid[:, 0, 0, 0, 0])
        spectrum_low_logU_high_density = np.zeros_like(grid[:, 0, 0, 0, 0])
        spectrum_high_logU_high_density = np.zeros_like(grid[:, 0, 0, 0, 0])

        logU_ind = config.logU[config.logU < logU].shape[0]
        logU_weight = (config.logU[logU_ind] - logU) / (
            config.logU[logU_ind] - config.logU[logU_ind - 1]
        )

        density_ind = config.densities[config.densities < n_e].shape[0]
        density_weight = (config.densities[density_ind] - n_e) / (
            config.densities[density_ind] - config.densities[density_ind - 1]
        )

        index = config.age_bins[config.age_bins < t_bc].shape[0]
        weight = 1 - (config.age_bins[index] - t_bc) / config.age_widths[index - 1]

        for i in range(config.metallicities.shape[0]):
            if sfh_ceh[i, :index].sum() > 0.0:
                sfh_ceh[:, index - 1] *= weight

                spectrum_low_logU_low_density += np.sum(
                    grid[:, i, logU_ind - 1, density_ind - 1, :index]
                    * sfh_ceh[i, :index],
                    axis=1,
                )

                spectrum_high_logU_low_density += np.sum(
                    grid[:, i, logU_ind, density_ind - 1, :index] * sfh_ceh[i, :index],
                    axis=1,
                )

                spectrum_low_logU_high_density += np.sum(
                    grid[:, i, logU_ind - 1, density_ind, :index] * sfh_ceh[i, :index],
                    axis=1,
                )

                spectrum_high_logU_high_density += np.sum(
                    grid[:, i, logU_ind, density_ind, :index] * sfh_ceh[i, :index],
                    axis=1,
                )

                sfh_ceh[:, index - 1] /= weight
            print(density_weight,logU_weight)
        spectrum = (
            spectrum_low_logU_low_density *logU_weight * density_weight
            + spectrum_low_logU_high_density * logU_weight * (1-density_weight)
            + spectrum_high_logU_low_density * (1-logU_weight) * density_weight
            + spectrum_high_logU_high_density * (1-logU_weight) * (1-density_weight)
        )
        return spectrum
