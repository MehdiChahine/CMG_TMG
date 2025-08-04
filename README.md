# CMG_TMG
This directory includes the Matlab functions implementing the Continuous Gaussian Mixture Approach to Sample Multivariate Gaussians constrained by linear inequalities.

This code is based on the worked introduced in the folowing paper: 
"A Continuous Gaussian Mixture Approach to Sample Multivariate Gaussians constrained by linear inequalities" by Mehdi AMROUCHE, Jérôme IDIER and Hervé CARFANTAN

DOI: hal-05195044
URL: https://hal.science/hal-05195044

Date of current version: July. 29,2025

Authors of the codes:
Mehdi AMROUCHE, Jérôme IDIER, and Hervé CARFANTAN 

(a) Description: see above.

(b) Platform: any platform on which Matlab is installed. No special toolbox is needed.

(c) Environment: Matlab (version 8 or more recent versions).

(d) Major component description:

    --- The directory provides three example scripts
        * SPX_example.m provides an example of sampling a TMG in the case of Single Type Constrains (STC) only. The illustrative example is the unit simplex.
        * BOX_example.m provides an example of sampling a TMG in the case of Box Type Constrains (BTC) only. The illustrative example is the unit ell_1 ball.
        * TMG_example.m provides an exmaple of sampling a TMG in the case of mixed STC and BTC constraints.

    --- And three functions to sample from TMGs
        * STC_CMG_TMG.m function that provides TMG samples in the case of Single Type Constrains (STC) only.
        * BTX_CMG_TMG.m function that provides TMG samples in the case of Box Type Constrains (BTC) only.
        * CMG_TMG.m     function that provides TMG samples in the case of mixed STC and BTC constraints.

