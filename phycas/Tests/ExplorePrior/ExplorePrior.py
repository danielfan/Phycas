# This example is designed to demonstrate that the MCMC code is working correctly by
# exploring only the priors. The posterior resulting from such an analysis should
# approximate the joint prior distribution.

from phycas import *
phycas = Phycas()

phycas.data_source              = None
phycas.outfile_prefix           = 'nodata.nex'  # used for creating *.p and *.t file names
phycas.ntax                     = 10 
phycas.starting_tree_source     = 'random'
phycas.ncycles                  = 1000
phycas.nchains                  = 1
phycas.sample_every             = 20
phycas.report_every             = phycas.ncycles//100
phycas.adapt_first              = 100
phycas.random_seed              = '463859'
phycas.default_model            = 'hky'
phycas.using_hyperprior         = True
phycas.num_rates                = 4
phycas.use_flex_model           = True
phycas.flex_ncat_move_weight    = 10        # number of times each cycle to attempt an ncat move
phycas.flex_num_spacers         = 1         # number of fake rates between each adjacent pair of real rates
phycas.flex_phi                 = 0.5       # proportion of ncat moves in which ncat is incremented (ncat is decremented with probability 1 - flex_phi)
phycas.flex_L                   = 5.0       # upper bound of interval used for unnormalized relative rate parameter values
phycas.flex_lambda              = 4.0       # parameter of Poisson prior on the number of extra categories
phycas.flex_prob_param_prior    = ProbDist.ExponentialDist(1.0)
#phycas.edgelen_prior_mean       = 1.0
phycas.edgelen_dist             = ProbDist.ExponentialDist(1.0)
phycas.verbose                  = True
phycas.ls_move_weight           = 100
phycas.slice_weight             = 1
phycas.slice_max_units          = 0
phycas.use_inverse_shape        = False

phycas.mcmc()
