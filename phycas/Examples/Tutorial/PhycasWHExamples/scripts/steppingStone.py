# This script will regenerate the information used for the points corresponding to 
# K = 2 in Figure 10 of this paper:
#
# Xie, W., P. O. Lewis, Y. Fan, L. Kuo and M.-H. Chen. 2009.
# Improving Marginal Likelihood Estimation for Bayesian Phylogenetic 
# Model Selection. Syst. Biol., submitted Feb. 6, 2009.
#
# This figure plots estimates of the marginal likelihood (y-axis) for several
# different values of K, where K is the number of beta intervals used. See
# the Xie et al. paper for details.
#
# Note that this same script can be used to create both the HM estimate (call "mcmc()")
# as well as the estimates for SS and PS (call "ss()"). 

import math, sys
from phycas import *

#Toggling the following variables determines whether or not you use:
# GTR vs HKY, and
# Stepping stone vs harmonic mean estimator of the marginal likelihood
use_gtr = True
using_stepping_stone = True


# specify the data file name
mcmc.data_source = 'simhky.nex'


# settings relating to substitution model

if use_gtr:
    model.type = 'gtr'
    model.relrates = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0] # AC, AG, AT, CG, CT, GT
else:
    model.type = 'hky'
    model.kappa = 4.0


model.num_rates = 1 # use number > 1 for discrete gamma rate heterogeneity
model.pinvar_model = False # use True for proportion invariable sites model
if model.num_rates > 1:
    model.gamma_shape = 0.5

model.state_freqs = [0.25, 0.25, 0.25, 0.25] # order is A,C,G,T

# settings relating to prior distributions
model.edgelen_hyperprior = None
model.edgelen_prior = ProbDist.Exponential(1.0)
model.relrate_prior = ProbDist.Dirichlet((1.0,1.0,1.0,1.0,1.0,1.0)) 
model.kappa_prior = Exponential(1.0)
model.state_freq_prior = Dirichlet((1.00000, 1.00000, 1.00000, 1.00000))

# general settings 
rng = Lot()
rng.setSeed(13957)
mcmc.rng = rng
mcmc.burnin = 100 # number of cycles before any sampling is done
mcmc.ncycles = 1000 # number of MCMC cycles per value of beta. This is not enough, We are just using this for demo purposes
mcmc.sample_every = 10 # log-likelihood will be sampled whenever cycle modulo ps_sample_every equals 0
mcmc.report_every = 100
mcmc.adapt_first = 2 
mcmc.verbose = True


# specify the starting tree
mcmc.starting_tree_source = TreeCollection(newick='(t1:0.074783,(((((((t2:0.017390,t3:0.013871):0.014610,t4:0.038516):0.023079,t5:0.033673):0.011254,((t6:0.000718,t7:0):0.001124,t8:0.008383):0.020560):0.022048,(t14:0.127203,t15:0.017124):0.024730):0.002278,((t9:0.048058,t10:0.059900):0.024875,((t11:0.038944,t12:0.025554):0.014442,t13:0.034559):0.022475):0.002477):0.006325,t16:0.023526):0.006762,t17:0.025302)')

# probably best to *not* change any of these
mcmc.nchains = 1 # multiple chain analyses are not yet working
mcmc.slice_weight = 1 # means one slice sampling update of each model parameter per cycle

# fix the tree topology and use EdgeMove (which only change edge lengths) instead
# of LargetSimonMove (which changes both edge lengths and the tree topology)
mcmc.fix_topology           = True
mcmc.edge_move_weight       = 1      # do 100 EdgeMove updates per cycle (makes sense to set this to some multiple of the number of edges)


# specify how many times per cycle to rescale the entire tree
mcmc.tree_scaler_weight     = 1

# specify whether to update relative rates individually using slice sampling (True) 
# or jointly using a Metropolis-Hastings proposal (False)
model.update_relrates_separately = False

# specify how many times per cycle to propose new relative rates vector
mcmc.rel_rate_weight        = 1

# specify whether to update base frequencies individually using slice sampling (True) 
# or jointly using a Metropolis-Hastings proposal (False)
model.update_freqs_separately = False

# specify how many times per cycle to propose a new vector of nucleotide frequencies
mcmc.state_freq_weight      = 1

# probably best to *not* change any of these
mcmc.nchains                = 1         # multiple chain analyses are not yet working
mcmc.slice_weight           = 1         # means one slice sampling update of each model parameter per cycle


# settings specific to path/steppingstone sampling
ss.nbetavals = 5 # number of beta values 
ss.minbeta = 0.0 # smallest beta value to be sampled
ss.maxbeta = 1.0 # largest beta value to be sampled
ss.shape1 = 0.3 # first shape parameter of the beta sampling distribution
ss.shape2 = 1.0 # second shape parameter of the beta sampling distribution

# specify whether to use Metropolis-Hastings/slice sampling when exploring the prior (False)
# or to sample directly from the prior (True)
mcmc.draw_directly_from_prior = True

# specify tuning parameters for the Metropolist-Hastings proposals
# Note: there are two versions of each tuning parameter, one for the posterior
# and one for the prior (denoted 0).
# When exploring distribtions intermediate between prior and posterior, an 
# intermediate tuning parameter is used
# Using bolder moves for distributions that are more similar to the prior improves mixing
mcmc.tree_scaler_lambda = 0.5 # min_lambda
mcmc.tree_scaler_lambda0 = 1.0 # max_lambda
mcmc.rel_rate_psi = 200.0 # max_psi
mcmc.rel_rate_psi0 = 2.0 # min_psi
mcmc.state_freq_psi = 200.0 # max_psi
mcmc.state_freq_psi0 = 2.0 # min_psi
mcmc.edge_move_lambda     = 0.5     # min_lambda
mcmc.edge_move_lambda0    = 1.0     # max_lambda


# This will be the prefix for the name of output files generated
if using_stepping_stone:
    fnprefix = model.type + 'stepstone'
else:
    fnprefix = model.type + 'harmonic'

# settings related to the log file
mcmc.out.log.mode = REPLACE
mcmc.out.log.prefix = fnprefix

# settings related to the file that saves sampled parameter values
mcmc.out.params.mode = REPLACE
mcmc.out.params.prefix = fnprefix

# settings related to the file that saves sampled trees
mcmc.out.trees.mode = REPLACE
mcmc.out.trees.prefix = fnprefix

sump.out.log.prefix = 'sumt' + fnprefix

if using_stepping_stone:
    # start the MCMC analysis that will generate samples appropriate
    # for both path sampling and steppingstone sampling
    ss()
else:
    mcmc.ncycles *= ss.nbetavals
    mcmc()


sump(burnin=1,file="%s.p" % fnprefix)
