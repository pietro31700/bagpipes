**Bayesian Analysis of Galaxies for Physical Inference and Parameter EStimation**

Please refer to official documentation `bagpipes.readthedocs.io <http://bagpipes.readthedocs.io>`_. Last bagpipes version 1.0.4




**The custom branch present some changes from the original package:**

* The grids are built for logU up to +0.5
* In the h5 file are saved more information in the attributes:
    * :code:`<opened h5 file>.attrs["parameter_names"]` gives the ordered (as the samples2d in the same file) list of the names of the free parameters in the fit
    * :code:`<opened h5 file>.attrs["maxl_model"]` gives the ready-to-use complete model of the galaxy as fitted. It is a dictionary. Import it with 
      
      .. code::
         
         maxl_params = eval(<opened h5 file>.attrs["maxl_model"].replace("array", "np.array").replace("float", "np.float"))

      The two :code:`replace` must be used when dealing with R_curve. Also :code:`<opened h5 file>.attrs["fitted_model"]` that contains the parameters of the fit with the priors selected should be opened in the same way.

* When fitting the SFH
    * you can select if also to plot the mean SFR value (instead of only the median SFR + 1σ CI) and if plot the SFH in log scale
    :code:`plot_sfh_posterior` has two new boolean parameters: :code:`mean` and :code:`log_scale`. For enabling the new options use:
    .. code::

        plot_sfh_posterior(save=True,show=False,log_scale=True,mean=True)
    
    * By default the x-axis is written as time from the observed time of the galaxy. To revert this option use:
    .. code::

        plot_sfh_posterior(save=True,show=False,from_bigbang=True)
    Moreover, more redshift values are printed on the second x-axis o the SFH plot.

Any previous python file written for the standard bagpipes package works as usual.


ADD STEP, RELATIVE_STEP, IPERB
**Installation**

If bagpipes was previously installed please uninstall it:

.. code::

    pip uninstall bagpipes

Custom bagpipes must be installed as (linux & mac):

Go into the directory where you want to install bagpipes, download this file:
`Grids file <https://mega.nz/file/U65QWByS#WhU0ScTbRoO0wWeVt7ZAxJh9Iom_IOjGUV1RO2U6SCM>`_
than run the following commands:
.. code::

    git clone https://github.com/pietro31700/bagpipes.git 
    tar -xvf grids.tar.gz -C ./bagpipes/bagpipes/models/
    pip install ./bagpipes/


To fit models to data with the code you will also need to install the `MultiNest <https://github.com/JohannesBuchner/MultiNest>`_ code. For more information please see the `official bagpipes documentation <http://bagpipes.readthedocs.io>`_.

**Published papers and citing the code**

Bagpipes is described primarily in Section 3 of `Carnall et al. (2018) <https://arxiv.org/abs/1712.04452>`_, with further development specific to spectroscopic fitting described in Section 4 of `Carnall et al. (2019b) <https://arxiv.org/abs/1903.11082>`_. These papers are the best place to start if you want to understand how the code works.

If you make use of Bagpipes, please include a citation to `Carnall et al. (2018) <https://arxiv.org/abs/1712.04452>`_ in any publications. You may also consider citing `Carnall et al. (2019b) <https://arxiv.org/abs/1903.11082>`_, particularly if you are fitting spectroscopy.

Please note development of the code has been ongoing since these works were published, so certain parts of the code are no longer as described. Please inquire if in doubt.
