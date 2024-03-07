DIAMONDS - Bayesian Software
============================

.. image:: https://img.shields.io/badge/GitHub-DIAMONDS-yellow
    :target: https://github.com/EnricoCorsaro/DIAMONDS
.. image:: https://img.shields.io/badge/license-MIT-blue
    :target: https://github.com/EnricoCorsaro/DIAMONDS/blob/master/LICENSE.txt
.. image:: https://img.shields.io/badge/DOI-10.1051%2F0004--6361%2F201424181-blueviolet
    :target: https://www.aanda.org/articles/aa/abs/2014/11/aa24181-14/aa24181-14.html
.. image:: https://img.shields.io/badge/ASCL-1410.001-red
    :target: https://ascl.net/1410.001
.. image:: https://readthedocs.org/projects/diamonds/badge/?version=latest
    :target: https://diamonds.readthedocs.io/en/latest/?badge=latest
.. image:: https://img.shields.io/github/issues-closed/EnricoCorsaro/DIAMONDS
    :target: https://github.com/EnricoCorsaro/DIAMONDS/issues
.. image:: https://raw.githubusercontent.com/EnricoCorsaro/DIAMONDS/master/docs/figures/DIAMONDS_LOGO_WHITE.png
    :width: 500 px


Authors
^^^^^^^
- `Enrico Corsaro <mailto:enrico.corsaro@inaf.it>`_
- `Joris De Ridder <mailto:joris.deridder@kuleuven.be>`_

Description
^^^^^^^^^^^
The **DIAMONDS** (high-DImensional And multi-MOdal NesteD Sampling) code presented in this website is developed in ``C++11`` and structured in classes in order to be as much flexible and configurable as possible. The working scheme from the main function is as follow:

- Read an input dataset
- Set up model, likelihood, and prior distributions to be used in the Bayesian inference
- Set up a drawing algorithm for the nested sampling
- Configure and start the nested sampling
- Compute and print the final results

The code can be used for any application involving Bayesian parameter estimation and/or model selection problems. Users can supply new models, likelihood functions, and prior distributions whenever needed by taking the advantage of ``C++`` class polymorphism and inheritance. Any new model, likelihood, and prior distribution can be defined and implemented upon a basic template. In addition, it is possible to feed the basic nested sample with different drawing algorithms based on clustering, as well as different clustering algorithms can in principle be used.

The original paper publication that presents DIAMONDS and its application to asteroseismology is available at 
`E. Corsaro & J. De Ridder 2014 A&A, 571, 71 <https://www.aanda.org/articles/aa/abs/2014/11/aa24181-14/aa24181-14.html>`_

Navigation
^^^^^^^^^^
.. toctree::
   :maxdepth: 2
   
   package_content
   installation
   user_guide_manual
   publications
   events
   logo
