.. _user_guide_manual:

User Guide Manual
=================
The user guide manual is a comprehensive text to guide you to a proper usage of the DIAMONDS code. It is available for download in the docs/ folder of the public GitHub repository of DIAMONDS `here <https://github.com/EnricoCorsaro/DIAMONDS/blob/master/docs/pdf/DIAMONDS_UserGuide_2018.pdf>`_.

The manual is usually referred to the current version of DIAMONDS. It is updated from time to time as soon as significant changes to the code are made (your feedback is of course important!). The User Guide Manual is intended as a complement to what described in the original paper publication of the code.

The manual is structured in four different chapters and provides all the required information to set up and properly configure the code for your own application.
Some very useful code internal calibrations and diagrams obtained from existing published results and from new experimental testing are also included to help the user in configuring the code. The manual provides a comprehensive troubleshooting of the main computational errors and an essential guide on how to interpret the results and understand their reliability.

A table of content is given in the following.

1. Getting Started
    - Setting up the local working path
    - The input dataset
    - The model
    - The likelihood function
    - The prior probability distribution
    - Configuring the cluster algorithm
    - Configuring the nested sampling algorithm
    - What output is obtained and how to use it
         
2. The enlargement fraction *f* of the sampling ellipsoid
    - The number of clusters *N*:subscript:`clust`
    - The number of live points *N*:subscript:`live`
    - The shrinking rate :math:`\alpha`
    - The initial enlargement fraction *f*:subscript:`0`
    - Computational times *t*:subscript:`comp` and *t*:subscript:`norm`
    - How *N*:subscript:`clust` can change the enlargement fraction *f*
         
3. Tackling incomplete computations
    - Assertion failure
    - No better likelihood points found
    - Ellipsoid matrix decomposition failed
    - Computation unable to start with NaN values
    - Computation unable to start with zero evidence
         
4. Checking the results and understanding their reliability
    - How *f*:subscript:`0` and :math:`\alpha` can change the MPD
    - False multimodal MPD
    - False spike-like MPD
    - Truncated MPD and the role of uniform priors
    - How to extend the methodology to multiple analyses
