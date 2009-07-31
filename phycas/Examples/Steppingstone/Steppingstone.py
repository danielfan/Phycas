# This script will regenerate the information used for the points corresponding to 
# K = 2 in Figure 10 of this paper:
#
#   Xie, W., P. O. Lewis, Y. Fan, L. Kuo and M.-H. Chen. 2009.
#   Improving Marginal Likelihood Estimation for Bayesian Phylogenetic 
#   Model Selection. Syst. Biol., submitted Feb. 6, 2009.
#
# This figure plots estimates of the marginal likelihood (y-axis) for several
# different values of K, where K is the number of beta intervals used. See
# the Xie et al. paper for details.
#
# Note that this same script can be used to create both the HM estimate (call "mcmc()")
# as well as the estimates for SS and PS (call "ss()"). 

import math, sys
from phycas import *

# This will be the prefix for the name of output files generated
fnprefix                     = 'two'

# set prior means
prior_mean = 0.1

# settings relating to substitution model
model.type                  = 'gtr'    # use 'jc' for JC,'hky' for HKY model, 'gtr' for GTR model
model.num_rates             = 4        # use number > 1 for discrete gamma rate heterogeneity
model.pinvar_model          = False    # use True for proportion invariable sites model
#model.pinvar                = 0.2
model.gamma_shape           = 0.5
model.relrates              = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0]    # AC, AG, AT, CG, CT, GT
model.state_freqs           = [0.25, 0.25, 0.25, 0.25]  # order is A,C,G,T
#model.kappa                 = 4.0

# settings relating to prior distributions
model.edgelen_hyperprior    = None
model.edgelen_prior         = ProbDist.Exponential(1.0)
model.relrate_prior         = ProbDist.Dirichlet((1.0,1.0,1.0,1.0,1.0,1.0)) 
model.state_freq_prior      = ProbDist.Dirichlet((1.0, 1.0, 1.0, 1.0))
model.gamma_shape_prior     = ProbDist.Exponential(1.0)
#model.pinvar_prior          = ProbDist.Beta(1.0, 1.0)

# general settings 
mcmc.random_seed            = 17539
mcmc.burnin                 = 1000    # number of cycles before any sampling is done
mcmc.ncycles                = 200000  # number of MCMC cycles per value of beta
mcmc.sample_every           = 100     # log-likelihood will be sampled whenever cycle modulo ps_sample_every equals 0
mcmc.report_every           = 1000
mcmc.adapt_first            = 2    
mcmc.verbose                = True

# specify the data file name
mcmc.data_source            = 'green.nex'

# settings related to the log file
mcmc.out.log.mode           = REPLACE
mcmc.out.log.prefix         = fnprefix

# settings related to the file that saves sampled parameter values
mcmc.out.params.mode        = REPLACE
mcmc.out.params.prefix      = fnprefix

# settings related to the file that saves sampled trees
mcmc.out.trees.mode         = REPLACE
mcmc.out.trees.prefix       = fnprefix

# specify the starting tree
mcmc.starting_tree_source = TreeCollection(newick='(1:0.258008,2:0.086726,(3:0.113174,(4:0.086740,((5:0.078662,8:0.183070):0.051180, \
(6:0.087103,((7:0.046703,10:0.098804):0.021001,9:0.072890):0.052088):0.034684):0.033875):0.044762):0.026713)')

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
ss.nbetavals                = 3      # number of beta values 
ss.minbeta                  = 0.0    # smallest beta value to be sampled
ss.maxbeta                  = 1.0    # largest beta value to be sampled
ss.shape1                   = 0.3    # first shape parameter of the beta sampling distribution
ss.shape2                   = 1.0    # second shape parameter of the beta sampling distribution

# specify whether to use Metropolis-Hastings/slice sampling when exploring the prior (False)
# or to sample directly from the prior (True)
mcmc.draw_directly_from_prior = True

# specify tuning parameters for the Metropolist-Hastings proposals
# Note: there are two versions of each tuning parameter, one for the posterior
# and one for the prior. When exploring distribtions intermediate between 
# prior and posterior, an intermediate tuning parameter is used

mcmc.tree_scaler_lambda   = 0.5     # min_lambda
mcmc.tree_scaler_lambda0  = 1.0     # max_lambda

mcmc.edge_move_lambda     = 0.5     # min_lambda
mcmc.edge_move_lambda0    = 1.0     # max_lambda

mcmc.state_freq_psi       = 200.0       # max_psi
mcmc.state_freq_psi0      = 2.0         # min_psi

mcmc.rel_rate_psi         = 200.0       # max_psi
mcmc.rel_rate_psi0        = 2.0         # min_psi

# start the MCMC analysis that will generate samples appropriate
# for both path sampling and steppingstone sampling
ss()

# uncomment the two lines below (and comment out the line above) to
# estimate the marginal likelihood using the harmonic mean method
#mcmc.ncycles *= ss.nbetavals
#mcmc()

# sump reads in the samples from the param file and computes the
# path sampling and steppingstone sampling estimates
sump(burnin=1,file="%s.p" % fnprefix)

# uncomment the two lines below only if you are running this script on a server
# that can send email, otherwise the script will hang at this point. Also, you
# will need to replace MYUNIVERSITY and MYEMAIL with appropriate values

# from phycas.Utilities.enotify import sendemail
# sendemail(smtphost='smtp.MYUNIVERSITY.edu', toaddr='MYEMAIL@MYUNIVERSITY.edu', fromaddr='MYEMAIL@MYUNIVERSITY.edu', subject='your run named "%s" has finished' % (fnprefix,))
