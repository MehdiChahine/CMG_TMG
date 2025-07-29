# CMG_TMG
This directory includes the Matlab functions implementing the Continuous Gaussian Mixture Approach to Sample Multivariate Gaussians constrained by linear inequalities.
The 
"A Continuous Gaussian Mixture Approach to Sample Multivariate Gaussians constrained by linear inequalities" by Mehdi AMROUCHE, Jérôme IDIER and Hervé CARFANTAN

DOI: XXX
Springer URL
HAL URL

Date of current version: July. 29,2025

Authors of the codes:
Mehdi AMROUCHE, Hervé CARFANTAN and Jérôme IDIER

(a) Description: see above.

(b) Platform: any platform on which Matlab is installed. No special toolbox is needed._

(c) Environment: Matlab (version 8 or more recent versions).

(d) Major component description:
    --- Major parts of the algorithm are divided into different Matlab scripts in the /Samplers directory:
        * init_PCGS.m initializes the variables of the PCGS sampler
        * sample_q_w.m joint sampling of (q,w) from the marginalized posterior following the RJMCMC framework
        * sample_x.m joint sampling of vector x (multidimensional Gaussian)
        * sample_h.m sampling of the parameter related to the impulse response (not used in the default setting)
        * sample_hyp.m sampling the hyper-parameters lambda (probability of having q_k =1) and sigma2k (the noise variance)

    --- Test program:
        * PCGS_BGM_S.m Sparse deconvolution problem using Bernoulli-Laplace prior, as described in section V-B of the corresponding paper.

(e) Detailed run instructions:
    Run PCGS_BGM_S.m calling

    > PCGS_BGM_S

    in Matlab.

(f) Output description:
    --- figure(1) & Data.pdf file show the observed data and the impulse response 
    --- figure(2) & qPM file show the posterior mean of Bernoulli variables q in comparison to the ground truth
    --- figure(3) & Estimation.pdf file show the estimation of the sparse signal x in comparison to the ground truth
    --- figure(4) & NoiseVarianceHistogram.pdf file show the histogram of the noise variance (hyper-parameter)
    --- figure(5) & LambdaHistogram.pdf file show the histogram of the hyper-parameter lambda, i.e., the probability of having q == 1
    --- figure(6) & ScaleHistogram.pdf file show the histogram of the scale hyper-parameter

(g) External materials:
    This MATLAB implementation uses the external function "rtnorm.m":
    --- "rtnorm.m" is a MATLAB function (under GNU/GPL license) developed by Vincent Mazet (ICube, CNRS/Université de Strasbourg)) which allows efficient sampling of random variables from a truncated Gaussian distribution
