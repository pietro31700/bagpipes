**Bayesian Analysis of Galaxies for Physical Inference and Parameter EStimation**

Please refer to official documentation `bagpipes.readthedocs.io <http://bagpipes.readthedocs.io>`_. Last bagpipes version 1.0.4

The custom branch present some changes from the original package:

* In the h5 file are saved more information in the attributes:
    * :code:`<opened h5 file>.attrs["parameter_names"]` gives the ordered (as the samples2d in the same file) list of the names of the free parameters in the fit
    * :code:`<opened h5 file>.attrs["maxl_model"]` gives the ready-to-use complete model of the galaxy as fitted. It is a dictionary. Import it with 
      
      :code:`maxl_params = eval(<opened h5 file>.attrs['maxl_model'].replace("array", "np.array").replace("float", "np.float"))`

      The two :code:`replace` must be used when dealing with R_curve. Also :code: `<opened h5 file>.attrs["fitted_model"]`that contains the parameters of the fit with the priors selected should be opened in the same way

* When fitting the sfh you can select if also to plot the mean SFR value (instead of only the median SFR + 1σ CI) and if plot the SFH in log scale

**Installation**

Custom bagpipes must be installed as:

.. code::

    pip install bagpipes

To fit models to data with the code you will also need to install the `MultiNest <https://github.com/JohannesBuchner/MultiNest>`_ code. For more information please see the `official bagpipes documentation <http://bagpipes.readthedocs.io>`_.

**Published papers and citing the code**

Bagpipes is described primarily in Section 3 of `Carnall et al. (2018) <https://arxiv.org/abs/1712.04452>`_, with further development specific to spectroscopic fitting described in Section 4 of `Carnall et al. (2019b) <https://arxiv.org/abs/1903.11082>`_. These papers are the best place to start if you want to understand how the code works.

If you make use of Bagpipes, please include a citation to `Carnall et al. (2018) <https://arxiv.org/abs/1712.04452>`_ in any publications. You may also consider citing `Carnall et al. (2019b) <https://arxiv.org/abs/1903.11082>`_, particularly if you are fitting spectroscopy.

Please note development of the code has been ongoing since these works were published, so certain parts of the code are no longer as described. Please inquire if in doubt.