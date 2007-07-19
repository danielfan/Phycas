# This example reanalyzes the land plant data used in Yang, Z., and B. Rannala. 2005.
# Systematic Biology 54(3): 455-470. Specifically, it compares the two edge length
# prior scenarios considered in Yang and Rannala's figure 6: external edges have an
# exponential prior with mean 0.1 and internal edges have a prior with mean either
# 0.1 or 0.0001. Only one combination of edge length priors is used in each run of 
# this example (search for mu_internal and mu_external to change these values).
#
# This example differs from the analysis in the Yang and Rannala paper in two ways:
# (1) it evaluates the Gelfand-Ghosh criterion for the purpose of asking whether
# there is reason to prefer one of these priors over the other; and (2) it does not
# fix any of the model parameters (Yang and Rannala fixed kappa=3.0, gamma shape=0.29,
# and base frequencies piA=0.272, piC=0.188, piG=0.233, and piT=0.307). One further
# difference applies if you set use_hyperpriors = True, in which case the means of the
# internal and external edge length priors are governed by separate hyperparameters,
# which absolves the researcher of the responsibility for selecting prior means.
#
# The results of running this example are described in the paper Lewis, P. O. A
# minimum posterior-predictive loss approach to model selection in Bayesian
# phylogenetics. Systematic Biology. In prep.

from phycas import *
from math import exp

# These three variables are provided to make it easy to find and change the internal
# and external edge length prior means, or to use the edge length hyperprior approach.
# Note that if use_hyperpriors is True, mu_internal and mu_external will be ignored.
mu_internal = 0.0001
mu_external = 0.1
use_hyperpriors = False

# Create a Phycas object
phycas = Phycas()

# Set various options before starting the run. To see all available options, you
# can go to the source: take a look at the beginning part of the Phycas.py file
#   phycas/Examples/LandPlants/LandPlants.py <- you are here
#   phycas/Phycas/Phycas.py                  <- here is Phycas.py

# Set up the edge length prior means. You can force Phycas to use just one prior for
# all edge lengths (the default) by setting one of these to None instead of to a
# probability distribution
phycas.external_edgelen_dist  = ProbDist.ExponentialDist(1.0/mu_external)
phycas.internal_edgelen_dist  = ProbDist.ExponentialDist(1.0/mu_internal)

# Specify that Gelfand-Ghosh measure should be calculated, the number of posterior-
# predictive simulations that should be performed per MCMC sample (in this case 2),
# and the name of the file where the Gelfand-Ghosh results will be saved (note how
# Python allows you to easily include the values of mu_internal and mu_external in
# the file name - the numbers replace the "%d" placeholders)
phycas.gg_do = True
phycas.gg_nreps = 2
phycas.gg_outfile = 'ggout.internal_%f.external_%f.txt' % (mu_internal, mu_external)

# Set up the substitution model
phycas.default_model = 'hky'   # use the Hasegawa-Kishino-Yano (1985) model
phycas.num_rates = 5           # add discrete gamma rate heterogeneity with 4 rate categories

# Set prior on the gamma shape parameter to have mean 0.5
phycas.gamma_shape_prior = ProbDist.ExponentialDist(2.0)

# Set the prior for the base frequency parameters. The base frequency parameters in Phycas
# are only proportional to the base frequencies (i.e. they often add up to a value greater
# than 1.0, unlike the actual relative frequencies). These are normalized prior to calculation
# of the likelihood. Placing a Gamma(a,1) prior on these base frequency parameters is
# equivalent to placing a Dirichlet(a,a,a,a) prior on the actual base frequencies. If you
# change the prior, you should keep the second parameter of the Gamma distribution (the scale
# parameter) set to 1.0 if you want the joint base frequency prior to remain Dirichlet.
phycas.base_freq_param_prior = ProbDist.GammaDist(1.0, 1.0)

# Set the prior for kappa, the ratio of the rate of transitions to the rate of transversions
phycas.relrate_prior = ProbDist.ExponentialDist(1.0)

# Use a hyperparameter to govern the mean of the branch length prior
# Instead of specifying, say, ExponentialDist(10.0) for branch lengths, this
# approach effectively sets the prior to ExponentialDist(mu) and lets mu be
# a free parameter in the model (that will tune itself during the MCMC run
# to hover around the mean branch length). This means we are not setting a prior
# directly on branch lengths, but we do need to specify one for the "hyperparameter"
# mu. For this "hyperprior" (the prior for a hyperparameter) we used an Inverse
# Gamma distribution having mean 1.0 and variance 10.0, as first suggested by
# Suchard et al. in their 2001 MBE (18:1001-1013) paper on Bayesian model selection.
phycas.using_hyperprior = use_hyperpriors 
phycas.edgelen_hyperprior = ProbDist.InverseGammaDist(2.1, 1.0/1.1)

# Tell phycas that we do not want to allow polytomies
phycas.allow_polytomies = False

# Start with a random tree
phycas.starting_tree_source = 'random'

# Specify the data file (here we specified the location relative to this file)
# If you need to specify a data file in a directory that is not the same as the
# one in which this script resides, you can do so, but always use forward slashes
# ('/') to separate file path elements even if running under Windows.
phycas.data_file_name = 'Yang_and_Rannala_Karol.nex'

# Let slice sampler have maximum freedom to extend slice
phycas.slice_max_units = 0

# Tell phycas that we want to run the MCMC analysis for 20000 cycles.
# Note that a cycle in Phycas differs from a generation in MrBayes.
# A cycle involves updating each non-branch-length parameter in the model
# as well as a certain number of Metropolis-Hastings updates of branch
# lengths and tree topology.
phycas.ncycles = 100000
phycas.sample_every = 20    # save tree and parameters every 10 cycles

# Finally, call mcmc(), which prepares phycas for the MCMC analysis, taking
# account of the changed settings above, then performs the actual analysis.
phycas.mcmc()
