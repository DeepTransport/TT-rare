# TT-Rare (Matlab)
Deep Inverse Rosenblatt Transports (DIRT) + MCMC sampling using Tensor Train (TT) approximation for rare event simulation. This folder implements examples from the paper [T. Cui, S. Dolgov, R. Scheichl, Deep importance sampling using tensor trains with application to a priori and a posteriori rare events](https://arxiv.org/abs/2209.01941).


## Cross Entropy method

This folder implements a cross entropy method using a Gaussian mixture. The main file to estimate the normalising constant of a given density (which can be used to estimate the unnormalised probability of the event, followed by dividing it by the normalising constant of the posterior, estimated separately) is

   - `cross_entropy.m`

which is accompanied by the following auxiliary function files:

   - `gm_density.m`     Evaluates each component of the mixture on given samples
   - `gm_init.m`        Initializes the mixture density
   - `gm_iter.m`        One iteration of the mixture update
   - `gm_samples.m`     iid samples from each component of the mixture density

