# This example program performs a path sampling analysis similar to that used in
# the paper Xie, W., P. O. Lewis, Y. Fan, L. Kuo and M.-H. Chen, which will be
# submitted to Systematic Biology before the end of 2008. This analysis computes
# the log of the marginal (model) likelihood using two methods:
# (1) the thermodynamic integration (path sampling) method described by Lartillot
#     and Phillippe (2006. Systematic Biology 55(2): 195-207)
# (2) the steppingstone sampling method described in Xie et al. (to be submitted)
# Both methods perform much better than the harmonic mean method in common use,
# which usually greatly overestimates the marginal likelihood of a model.

import math, sys
from phycas import *

# Specify an output file name prefix that will begin the name of all output files
fnprefix                     = 'gtrg'

# Settings relating to substitution model (GTR+G)
model.type                  = 'gtr'    # use 'jc' for JC,'hky' for HKY model, 'gtr' for GTR model
model.num_rates             = 4        # use number > 1 for discrete gamma rate heterogeneity (the "+G" part of "GTR+G")
model.pinvar_model          = False    # change to True for proportion invariable sites (pinvar) model

# Starting values for model parameters
model.gamma_shape           = 0.5      # starting value for discrete gamma rate heterogeneity shape parameter
model.relrates              = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0]    # starting relative rates in this order: AC, AG, AT, CG, CT, GT
model.state_freqs           = [0.2, 0.3, 0.3, 0.2]  # starting state frequencies in this order: A,C,G,T

# Prior distributions for model parameters
edgelen_prior_mean          = 0.1
model.edgelen_hyperprior    = None      # use a non-hierarchical model for edge lengths
model.edgelen_prior         = ProbDist.Exponential(1.0/edgelen_prior_mean)  # prior for each edge length parameter
model.relrate_prior         = ProbDist.Exponential(0.001)   # vague prior for the GTR relative rate parameters
model.state_freq_prior      = ProbDist.Dirichlet((1.0, 1.0, 1.0, 1.0))  # flat joint prior for the state frequencies
model.gamma_shape_prior     = ProbDist.Exponential(0.001)   # vague prior for the discrete gamma rate heterogeneity shape parameter

# General settings 
mcmc.random_seed            = 91753  # pseudorandom number generator seed (can be any positive integer, specify zero to choose automatically)
mcmc.burnin                 = 1000   # number of MCMC cycles before any sampling is done
mcmc.ncycles                = 10000  # number of MCMC cycles per value of beta
mcmc.sample_every           = 10     # a tree and parameter sample will be taken whenever cycle modulo sample_every equals 0
mcmc.report_every           = 100    # a brief progress report will be made whenever cycle modulo report_every equals 0
mcmc.adapt_first            = 2      # slice samplers will be adapted the first time after adapt_first cycles, then after
                                     # twice adapt_first samples, then after four times adapt_first cycles, and so on.
                                     # This adaptation schedule is reset for each new beta value because the marginal
                                     # posterior densities of parameters can change dramatically with a new beta value.
mcmc.verbose                = True   # print adaptation summaries and progress reports to the console

# Settings relating to the tree used
mcmc.starting_tree_source = TreeCollection(newick='(1:0.258008,2:0.086726,(3:0.113174,(4:0.086740,((5:0.078662,8:0.183070):0.051180, \
(6:0.087103,((7:0.046703,10:0.098804):0.021001,9:0.072890):0.052088):0.034684):0.033875):0.044762):0.026713)')

# Settings related to the data file
mcmc.data_source            = 'green.nex'

# Settings related to output files
mcmc.out.log.mode           = REPLACE   # automatically replace the log file if it exists
mcmc.out.log.prefix         = fnprefix  # the log file will start with fnprefix followed by ".txt"
mcmc.out.params.mode        = REPLACE   # automatically replace the sampled parameter file if it exists
mcmc.out.params.prefix      = fnprefix  # the sampled parameter file will start with fnprefix followed by ".p"
mcmc.out.trees.mode         = REPLACE   # automatically replace the sampled tree file if it exists
mcmc.out.trees.prefix       = fnprefix  # the sampled tree file will start with fnprefix followed by ".t"

# Settings related to Metropolis-Hastings proposals
mcmc.tree_scaler_weight       = 10      # Propose rescaling the entire tree 10 times per cycle
model.update_freqs_separately = False   # Update the state frequencies jointly using a Metropolis-Hastings proposal (not separately using separate slice samplers)
mcmc.state_freq_weight        = 10      # Propose new set of state frequencies 10 times per cycle
mcmc.fix_topology             = True    # Fix the tree topology (but not edge lengths)
mcmc.edge_move_weight         = 100     # Do 100 EdgeMove updates per cycle (makes sense to set this to some multiple of the number of edges)

# Settings specific to path sampling (note: mcmc.ncycles and mcmc.sample_every are used to determine number of samples)
ps.nbetavals                  = 17      # Number of beta values to sample (number of intervals is one fewer than this value)
ps.minbeta                    = 0.0     # Smallest beta value to be sampled (should always be set to 0 - see note below)
ps.maxbeta                    = 1.0     # Largest beta value to be sampled (should always be set to 1 - see note below)
ps.shape1                     = 0.3     # First shape parameter of the beta sampling distribution (0.3 is a good choice)
ps.shape2                     = 1.0     # Second shape parameter of the beta sampling distribution (1.0 is a good choice)

# The ps.minbeta and ps.maxbeta settings should not be changed from their defaults if your intention is to 
# perform path sampling or steppingstone sampling. If, however, you wished to obtain a sample from the prior,
# this could be easily done by setting ps.minbeta = 0.0 and ps.maxbeta = 0.0. 

# Carry out the path sampling analysis; the log-marginal-likelihood computed by both
# the thermodynamic integration and steppingstone sampling methods will be reported 
# at the end.
ps()
