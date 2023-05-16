# TT-Rare (Matlab)
Deep Inverse Rosenblatt Transports (DIRT) + MCMC sampling using Tensor Train (TT) approximation for rare event simulation. This folder implements examples from the paper [T. Cui, S. Dolgov, R. Scheichl, Deep importance sampling using tensor trains with application to a priori and a posteriori rare events](https://arxiv.org/abs/2209.01941).


## Rare event sampling in diffusion flow

This example benchmarks DIRT on sampling rare events of a tracer particle following a diffusion flow and escaping the domain faster than a chosen minimal time threshold `thres`. This model is described in Section 6 of the paper.

### Running the experiments

The files `test_diffusion_*.m` run benchmarks of the corresponding algorithms:
   - `test_diffusion_rare_dirt.m`  DIRT experiments with the diffusion particle tracer
   - `test_diffusion_rare_CE.m`  Cross Entropy experiments with the diffusion particle tracer

Each script will ask interactively (or take from *varargin*) the following parameters (case sensitive).
#### PDE/Prior parameters (all tests)
 * *sigma* Variance of the affine diffusion coefficient expansion (or of **log** of the coefficient for log-uniform and log-normal fields)
 * *corr_length* Correlation length of the expansion
 * *nu* Decay rate
 * *tol_kle* Relative tolerance for the KLE eigenvalue truncation

#### Inverse problem parameters (all tests)
 * *sigma_n* Variance of the observation noise
 * *m0* Number of observation points in each variable (total number of observations is m0*m0)
 * *y0* Synthetic truth value of the random variables. Special values are `'b'` indicating a high-contrast barrier inside the domain, and `'c'`, indicating aa high-contrast channel facilitating fast flow of the particle.
 * *Nsamples* Number of samples used for estimating the normalising constant and probability
 * *thres* Escape time threshold defining the event
 * *gamma* Smoothing parameter (width coefficient in the sigmoid approximation of the indicator function of the event)
 * *runs* Number of runs (replicas) of the experiment

#### DIRT approximation parameters (`test_diffusion_rare_dirt` only)

 * *npi* Number of discretization points in each random variable in the ratio functions on each level of DIRT
 * *rpi* TT rank of each DIRT ratio function (fixed-rank TT decompositions are used)
 * *beta* A vector of tempering powers. Should be in increasing order, and the last value should be 1.
 * *pp* Power of the prior tempering, so the actual tempering coefficient in front of the log-prior is beta^pp

#### Cross Entropy approximation parameters (`test_diffusion_rare_CE` only)

 * *Ng* Number of Gaussians constituting the mixture density fitted with cross entropy
 * *em_iter* Maximal number of iterations
 * *rho* 1 minus the resampling quantile


### Output variables

Each script creates certain variables that remain accessible in the main Matlab workspace.
Use the `who` command for the full list.
Some of the interesting variables are:

 * *Z_post* Normalising constant of the posterior
 * *ttimes_post* CPU time of posterior DIRT approximation and sampling
 * *evalcnt_post* Number of model evaluations during posterior approximation (for DIRT)
 * *n_iter_post* Number of iterations in Cross Entropy for posterior normalising constant approximation
 * *ess_post* N/ESS of the posterior
 * *tau_post* IACT of the Metropolised chain sampling from the posterior (DIRT only)
 * *hell_post* Hellinger-distance error of the posterior approximation (DIRT only)
 * *ttimes_event* CPU time of event DIRT approximation and sampling
 * *evalcnt_event* Number of model evaluations during event approximation (for DIRT)
 * *n_iter_event* Number of iterations in Cross Entropy for event probability approximation
 * *ess_event* N/ESS of the event biasing function
 * *tau_event* IACT of the Metropolised chain sampling the event (DIRT only)
 * *hell_event* Hellinger-distance error of the event indicator approximation (DIRT only)
 * *P_event* Probability of the event t<thres

Each of these variables is a vector of size *runs* x 1, so statistics over experiments can be computed. Each script will also print average values.

### Function files

 * `parse_diffusion_inputs.m`    A function for requesting/extracting model parameters
 * `build_kle_eig.m`             Discretization of the diffusion equation
 * `diffusion_QoI.m`             Solution of the diffusion equation
 * `loglike_time.m`              Log-biasing function of the event
 * `escape_time_vec.m`           Computes the particle escape time
 * `flux_interpolate_vec.m`      Computes the flux (interpolated using piecewise linear basis used for solving the PDE)

