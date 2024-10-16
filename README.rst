Bayesian Analysis of Galaxies for Physical Inference and Parameter EStimation
=============================================================================

Please refer to official documentation `bagpipes.readthedocs.io <http://bagpipes.readthedocs.io>`_. Last bagpipes version 1.0.4

Installation
------------

If bagpipes were previously installed please uninstall it:

.. code-block:: console

    pip uninstall bagpipes

Custom bagpipes must be installed in the following way (Linux & mac):

Go into the directory where you want to install bagpipes, and download this file:
`Grids file <https://mega.nz/file/di5FCTgI#M8-v6kFjj_aHanPqPLWC6XafvPtZyCWz1NLK9Deg5VI>`_
then run the following commands:

.. code-block:: console

    git clone https://github.com/pietro31700/bagpipes.git
    tar -xvf grids.tar.gz -C ./bagpipes/bagpipes/models/grids/
    pip install ./bagpipes/


To fit models to data with the code you will also need to install the `MultiNest <https://github.com/JohannesBuchner/MultiNest>`_ code. For more information please see the `official bagpipes documentation <http://bagpipes.readthedocs.io>`_.


The custom branch presents the following changes from the original package:
--------------------------------------------------------------------------

+ The nebular grids are built for logU up to +0.5, and max_redshift=15.

+ You can use BPASS by changing the variable BPASS to True in the config.py file before the installation (enabled by default). If you want to change the model after the installation, you need to first delete the file models/grids/d_igm_grid_inoue14.fits. In this way, the absorption grids will be built again. Even for Bpass, the nebular grids are built for logU up to +0.5

+ When fitting the SFH:

  + Two new SFH shapes are introduced:

    + **step**: This is simply a step function with each step independent from the others. The steps' edges are defined by you, than for each step you need to select the prior on the star formation rate. An example with uniform priors is the following

      .. code-block:: python

        step = {}
        step["metallicity"] = (0.05,0.3)
        step["metallicity_prior"] = "log_10"
        step["bin_edges"] = [0,10,50,200,600] #Myr from the galaxy time looking behind

        for i in range(1, len(step["bin_edges"])):
          # step["sfr1"] specifies the SFR of the bin closer to the galaxy time (new star formation)
          step["sfr" + str(i)] = (0, 50)
          step["sfr" + str(i) + "_prior"] = "uniform"

      **Note**: Step is the only SFH shape that does not use the key "massformed". If specified it will be ignored. However, when *step* is used, in the corner plot the total mass formed is also displayed as a derived quantity.

    + **relative_step**: This is similar to step, but each bin of star formation is relative to the SFR of the first bin (the newest star formation). In this case "massformed" is required as always

      .. code-block:: python

        rel_step = {}
        rel_step["metallicity"] = (0.05,0.3)
        rel_step["metallicity_prior"] = "log_10"
        rel_step["bin_edges"] = [0,10,50,200,600] #Myr from the galaxy time looking behind

        for i in range(1, len(step["bin_edges"])-1):
          # rel_step["rsfr1"] specifies the SFR of the second bin 
          # respect to the "0th" bin: the one closer to the galaxy time (newest star formation)
          rel_step["rsfr" + str(i)] = (0.5, 2)
          rel_step["rsfr" + str(i) + "_prior"] = "uniform"

  + Two new priors have been introduced:

    + "discrete" is a discrete prior with "n" equally weighted values between the limits. In the example below the first valid value  for the R_curve_multiplier variable is 1 and the last is 1.9

      .. code-block:: python

        fit_instructions["R_curve_multiplier"] = (0.,2.)
        fit_instructions["R_curve_multiplier_prior"] = "discrete"
        fit_instructions["R_curve_multiplier_prior_n"] = 20

      The discrete distribution can help to speed up the convergence of very informative likelihoods that do not need to be found so precisely.

    + "hyperbolic" can be used to overpopulate low values and keep a flat distribution elsewhere. "hyperbolic" has a parameter "eta" which selects the *knee* of the distribution. Above this value, the distribution is almost flat, below there are more occurrences.

      .. code-block:: python

        step = {}
        step["bin_edges"] = isolight_steps(n_bins=10,redshift=8,redshift_end=30)

        for i in range(1, len(step["bin_edges"])):
          step["sfr" + str(i)] = (0, 50)
          step["sfr" + str(i) + "_prior"] = "hyperbolic"
          step["sfr" + str(i) + "_prior_eta"] = 5


      **Note**: In the case "eta" becomes many times larger than the prior width the distribution becomes a square root.

+ When plotting the SFH:

  + You can select if also to plot the mean SFR value (instead of only the median SFR + 1σ CI) and if plot the SFH in log scale (both x and y axis)
    ``plot_sfh_posterior]`` has three new boolean parameters: ``mean``, ``log_x`` and ``log_y``. To enable the new options use the following arguments:

    .. code-block:: python

      plot_sfh_posterior(save=True,show=False,log_x=True,log_y=True,mean=True)
    
  + By default, the x-axis is written as lookback time from the time of the galaxy. To revert this option use:
      
    .. code-block:: python

      plot_sfh_posterior(save=True,show=False,from_bigbang=True)

    Moreover, more redshift values are printed on the second x-axis of the SFH plot.

+ I introduced a more beautiful spectrum plot which includes residual plot and line naming (WARNING: NIRSpec oriented plot). You can call this feature by invoking :

  .. code-block:: python

    fit.plot_spectrum_beautiful(save=True, show=True)
    fit.plot_spectrum_beautiful(line_close_up="OIIIa", save=True, show=False)

+ The h5 file contains more information in the attributes:

  + ``<h5 file>.attrs["parameter_names"]`` gives the ordered (as the samples2d in the same file) list of the names of the free parameters in the fit
  + ``<h5 file>.attrs["maxl_model"]`` gives the ready-to-use complete model of the galaxy as fitted. It is a dictionary. Import it with:

    .. code-block:: python

      maxl_params = eval(<h5 file>.attrs["maxl_model"].replace("array", "np.array").replace("float", "np.float"))
    
    The two ``replace`` must be used when dealing with R_curve. Also ``<h5 file>.attrs["fitted_model"]`` containing the fit parameters with the priors selected should be opened in the same way.

+ When modeling a galaxy:

  + A new key has been introduced to allow the addition of noise to the spectrum.

    .. code-block:: python

     model_components["flux_sensitivity"] = np.c_[wavelengths,sensitivity]

    where wavelengths must be in Angstrom and sensitivity in erg/(s*AA*cm^2). If "R_curve" is also provided to the model, the noise is added to the spectrum and then convolved with "R_curve" specifications.
  
  + The new key "R_curve_multiplier" has been added to increase the resolving power for targets that do not fill the slit.
    The resolving power curve provided is multiplied by this coefficient (probably it should greater than 1). It can be also be left as a free parameter.


Any previous Python file written for the standard bagpipes package works as usual.

Published papers and citing the code
------------------------------------

Bagpipes is described primarily in Section 3 of `Carnall et al. (2018) <https://arxiv.org/abs/1712.04452>`_, with further development specific to spectroscopic fitting described in Section 4 of `Carnall et al. (2019b) <https://arxiv.org/abs/1903.11082>`_. These papers are the best place to start if you want to understand how the code works.

If you make use of Bagpipes, please include a citation to `Carnall et al. (2018) <https://arxiv.org/abs/1712.04452>`_ in any publications. You may also consider citing `Carnall et al. (2019b) <https://arxiv.org/abs/1903.11082>`_, particularly if you are fitting spectroscopy.

Please note development of the code has been ongoing since these works were published, so certain parts of the code are no longer as described. Please inquire if in doubt.
