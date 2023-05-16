# TT-Rare (Matlab)
Deep Inverse Rosenblatt Transports (DIRT) + MCMC sampling using Tensor Train (TT) approximation for rare event simulation. This repository implements examples from the paper [T. Cui, S. Dolgov, R. Scheichl, Deep importance sampling using tensor trains with application to a priori and a posteriori rare events](https://arxiv.org/abs/2209.01941).


## Rare event sampling in SIR model

This example benchmarks DIRT on sampling rare events of the maximal number of infected individuals in a compartmental SIR model exceeding a threshold `Imax`. This model is described in Section 5 of the paper.

### Running the experiments

The files `test_sir_*.m` run benchmarks of the corresponding algorithms and topologies of the compartments:
   - `test_sir_chain.m`        DIRT experiments with the chain topology of the compartments
   - `test_sir_austria.m`      DIRT experiments with the Austrian road topology of the compartments
   - `test_sir_CE.m`           Cross Entropy experiments with the chain topology of the compartments

Each script will ask interactively (or take from *varargin*) the following parameters (case sensitive).
#### ODE/Prior parameters (all tests)
 * *K* Number of compartments

#### Inverse problem parameters (all tests)
 * *sigma_n* Variance of the observation noise
 * *Imax* Threshold for the number of infected defining the event
 * *gammaev* Smoothing parameter (width coefficient in the sigmoid approximation of the indicator function of the event)
 * *Nsamples* Number of samples used for estimating the normalising constant and probability
 * *runs* Number of runs (replicas) of the experiment

#### DIRT approximation parameters (`test_sir_chain` and `test_sir_austria` only)

 * *npi* Number of discretization points in each random variable in the ratio functions on each level of DIRT
 * *rpi* TT rank of each DIRT ratio function (fixed-rank TT decompositions are used)
 * *beta* A vector of tempering powers. Should be in increasing order, and the last value should be 1.

#### Cross Entropy approximation parameters (`test_sir_CE` only)

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
 * *P_event* Probability of the event I>Imax

Each of these variables is a vector of size *runs* x 1, so statistics over experiments can be computed. Each script will also print average values.


### Function files

 * `parse_sir_inputs.m`    A function for requesting/extracting model parameters
 * `sir_rhs.m`             Right hand side of the SIR ODE
 * `sir_ll.m`              Log-likelihood of SIR (since the prior is uniform, this is also unnormalised log-posterior)
 * `ll_logsigmoid.m`       Log-biasing function of the event
 * `InfectedFun.m`         Computes the number of infected individuals for Cross Entropy

