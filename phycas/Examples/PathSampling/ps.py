# This example program performs a "polytomy" analysis similar to that used in
# the paper Lewis, P. O., M. T. Holder, and K. E. Holsinger. 2005.
# Polytomies and Bayesian phylogenetic inference. Systematic Biology 54:241-253.

from phycas import *
from math import exp

# Set up the substitution model
model.type      = 'hky'   # use the Hasegawa-Kishino-Yano (1985) model
model.num_rates = 4        # add discrete gamma rate heterogeneity with 4 rate categories

# Set prior on the gamma shape parameter to be an exponential distribution having mean 0.5
# Note that the value in parentheses is the inverse of the mean, not the mean.
model.gamma_shape_prior = ProbDist.Exponential(2.0)

# Set the prior for the base frequency parameters. The base frequency parameters in Phycas
# are only proportional to the base frequencies (i.e. they often add up to a value greater
# than 1.0, unlike the actual relative frequencies). These are normalized prior to calculation
# of the likelihood. Placing a Gamma(a,1) prior on these base frequency parameters is
# equivalent to placing a Dirichlet(a,a,a,a) prior on the actual base frequencies. If you
# change the prior, you should keep the second one (the scale) set to 1.0 if you want the
# joint base frequency prior to remain Dirichlet.
model.base_freq_param_prior = ProbDist.Gamma(1.0, 1.0)

# Set the prior for kappa, the ratio of the rate of transitions to the rate of transversions
model.kappa_prior = ProbDist.Exponential(1.0)

# Use a hyperparameter to govern the mean of the branch length prior
# Instead of specifying, say, Exponential(10.0) for branch lengths, this
# approach effectively sets the prior to Exponential(mu) and lets mu be
# a free parameter in the model (that will tune itself during the MCMC run
# to hover around the mean branch length). This means we are not setting a prior
# directly on branch lengths, but we do need to specify one for the "hyperparameter"
# mu. For this "hyperprior" (the prior for a hyperparameter) we used an Inverse
# Gamma distribution having mean 1.0 and variance 10.0, as first suggested by
# Suchard et al. in their 2001 MBE (18:1001-1013) paper on Bayesian model selection.
model.edgelen_hyperprior = ProbDist.InverseGamma(2.1, 1.0/1.1)

# Tell phycas that we want to allow polytomies
mcmc.allow_polytomies = True

# Tell phycas that we want to use the "polytomy" prior, which specifies
# how much higher the prior probability is for a given tree compared to
# a tree that has one more internal node (i.e. is slightly more resolved)
# Here, we are setting the value of C used to determine the topology
# prior ratio to exp(1), which equals the famous constant e (approx. 2.71828).
# This means tree topologies with k internal nodes will be e times
# more probable a priori than tree topologies having k+1 internal nodes.
# The value e has no particular significance, but it is aesthetically
# pleasing in that a slightly-more-resolved tree (with k+1 internal
# nodes) needs to be more than one log-likelihood unit better than a
# tree with just k internal nodes in order to overcome this prior.
mcmc.polytomy_prior   = True
mcmc.topo_prior_C     = exp(1.0)

# Start with a random tree
mcmc.starting_tree_source = 'random'

# Specify the data file (here we specified the location relative to this file)
# Note that you should use forward slashes ('/') even if running in Windows.
mcmc.data_file_name = 'ShoupLewis.nex'

# Let slice sampler have maximum freedom to extend slice
mcmc.slice_max_units = 0

# Specify the names of the files that will store the trees and parameter values
mcmc.out.trees.prefix = 'trees'
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = 'params'
mcmc.out.params.mode = REPLACE

# Specify the names of the file that will store the output
mcmc.out.log.prefix = 'output'
mcmc.out.log.mode = REPLACE

# Choose a pseudorandom number seed so that we can later repeat the analysis exactly if we so desire
mcmc.random_seed = 13957

# Tell Phycas that we want to run the MCMC analysis for 1000 cycles per beta value,
# and that we want to burn-in for 1000 cycles before sampling the first beta value.
mcmc.burnin = 1000
mcmc.ncycles = 1000
mcmc.sample_every = 10    # save tree and parameters every 10 cycles

# Setup path sampling
ps.nbetavals = 1
ps.maxbeta = 0.0
ps.minbeta = 0.0
ps()