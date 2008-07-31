# This example is designed to demonstrate that the MCMC code is working correctly by
# exploring only the priors. The posterior resulting from such an analysis should
# approximate the joint prior distribution.

from phycas import *

#model.using_hyperprior         = True
model.num_rates                = 4
model.use_flex_model           = True
model.flex_ncat_move_weight    = 10        # number of times each cycle to attempt an ncat move
model.flex_num_spacers         = 1         # number of fake rates between each adjacent pair of real rates
model.flex_phi                 = 0.5       # proportion of ncat moves in which ncat is incremented (ncat is decremented with probability 1 - flex_phi)
model.flex_L                   = 5.0       # upper bound of interval used for unnormalized relative rate parameter values
model.flex_lambda              = 4.0       # parameter of Poisson prior on the number of extra categories
model.flex_prob_param_prior    = Exponential(1.0)
model.edgelen_dist             = Exponential(1.0)
model.use_inverse_shape        = False

mcmc.data_source              = None
mcmc.out.trees.prefix         = 'nodata.nex'
mcmc.out.trees.mode           = REPLACE
mcmc.out.params.prefix        = 'nodata.nex'
mcmc.out.params.mode          = REPLACE
mcmc.ntax                     = 10 
mcmc.starting_tree_source     = randomtree(n_taxa=10)
mcmc.ncycles                  = 1000
mcmc.nchains                  = 1
mcmc.sample_every             = 20
mcmc.report_every             = mcmc.ncycles//100
mcmc.adapt_first              = 100
mcmc.random_seed              = 463859
#mcmc.edgelen_prior_mean       = 1.0
mcmc.verbose                  = True
mcmc.ls_move_weight           = 100
mcmc.slice_weight             = 1
mcmc.slice_max_units          = 0

mcmc()
