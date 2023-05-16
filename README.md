# TT-Rare (Matlab)
Deep Inverse Rosenblatt Transports (DIRT) + MCMC sampling using Tensor Train (TT) approximation for rare event simulation. This folder implements examples from the paper [T. Cui, S. Dolgov, R. Scheichl, Deep importance sampling using tensor trains with application to a priori and a posteriori rare events](https://arxiv.org/abs/2209.01941).

## Installation

The package is based on [TT-IRT](https://github.com/dolgov/TT-IRT), the main DIRT package using discrete TT decomposition, which in turn depends on [TT-Toolbox](https://github.com/oseledets/TT-Toolbox).

It should be sufficient to just run `install` script in Matlab. This will check the prerequisites and attempt to download them automatically, if Matlab has access to the Internet. For further information, as well as if the automatic download fails, please navigate to [TT-IRT](https://github.com/dolgov/TT-IRT) and follow the installation instructions for the **matlab/** section in the [`main README`](https://github.com/dolgov/TT-IRT/tree/master/README.md) file there.

## Examples

All files for running experiments start with a `test_` prefix. You may click on each item to open an individual directory or file.

 * [`sir/`](https://github.com/dolgov/TT-IRT/tree/master/matlab/examples/rare_events/sir)   Compartmental Susceptible-Infectious-Removed ODE model
   - `test_sir_chain.m`        DIRT experiments with the chain topology of the compartments
   - `test_sir_austria.m`      DIRT experiments with the Austrian road topology of the compartments
   - `test_sir_CE.m`           Cross Entropy experiments with the chain topology of the compartments
 * [`diffusion/`](https://github.com/dolgov/TT-IRT/tree/master/matlab/examples/rare_events/diffusion)  Diffusion model of the contaminant transport in groundwater
   - `test_diffusion_rare_dirt.m`  DIRT experiments with the diffusion particle tracer
   - `test_diffusion_rare_CE.m`  Cross Entropy experiments with the diffusion particle tracer

Each test can be run without any arguments. In this case, they will be interactively asked from the user. For batch runs, parameters can be passed in pairs of inputs ``'param_name1'``, ``param_value1``, ``'param_name2'``, ``param_value2``, and so on. For example,
```
test_sir_chain('K', 1, 'sigma_n', 1, 'Imax', 80, 'gammaev', 3e3/80, 'Nsamples', 2^14, 'runs', 1, 'npi', 17, 'rpi', 7, 'beta', 10.^(-4:1/3:0))
```
will run the SIR test with all default parameters. Only a subset of arguments can be given, e.g.
```
test_sir_chain('Nsamples', 1000, 'runs', 10)
```
will take the corresponding values from the inputs, and ask the user only for the remaining parameters.
Default parameters (corresponding to those in the paper) can be selected by entering empty input in the corresponding prompt.

Each test will print some statistical data (expected values, standard deviations, CPU times and so on).
Moreover, it will create all the variables in the main Matlab workspace, so they can be accessed afterwards.

### Further docs

Navigate to the [`main README`](https://github.com/dolgov/TT-IRT/tree/master/README.md) file of the TT-IRT repository and explore other folders.
In addition, each function file contains its own description in the first comment. See e.g.
```
help('tt_dirt_approx')
```
or open the file in the editor.

