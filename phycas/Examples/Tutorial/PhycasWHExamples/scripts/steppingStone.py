# This script demonstrates the use of the steppingstone sampler developed by:
#
# Xie, W., P. O. Lewis, Y. Fan, L. Kuo and M.-H. Chen. 2009.
# Improving Marginal Likelihood Estimation for Bayesian Phylogenetic 
# Model Selection. Syst. Biol., submitted Feb. 6, 2009.

from phycas import *

# The following variables are python variables that are not used in the rest
#   of Phycas.  Rather they are used within this script to make it easier to 
#   tweak the analysis and run it in a different configuration.
use_gtr = True              # Set this to True for GTR and False for HKY 

using_steppingstone = True # Set this to True for the steppingstone sample
                            #   or False to run a typical MCMC analysis and 
                            #   report the harmonic mean.

# specify the data file name
mcmc.data_source = 'simhky.nex'




# Settings relating to substitution model
# This demonstrates how we can use an if/else conditional statement in Python
#   note that the indentation of the code matters here.
if use_gtr:
    model.type = 'gtr'
    model.relrates = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0] # These will be the initial
                                                    # exchangeability rates, in the 
                                                    # order AC, AG, AT, CG, CT, GT
    # now we need to specify a prior:
    model.relrate_prior = ProbDist.Dirichlet((1.0,1.0,1.0,1.0,1.0,1.0)) 
else:
    model.type = 'hky'
    model.kappa = 4.0 # this is initial value of the transitition/transversion rate ratio.
    model.kappa_prior = Exponential(0.3)


model.num_rates = 1        # We won't worry with gamma-distributed rate het in this example
model.pinvar_model = False # ... or the invariant sites model.

# Set the starting base frequencies. The order in the list  is A,C,G,T,
model.state_freqs = [0.25, 0.25, 0.25, 0.25] 
model.state_freq_prior = Dirichlet((1.00000, 1.00000, 1.00000, 1.00000))

# settings relating to prior distributions
model.edgelen_prior = ProbDist.Exponential(1.0) # We'll use a fixed prior on the branch lengths
model.edgelen_hyperprior = None                 # ... not a hierarchical model.

# Create a random number generator and pass it to the mcmc command so that 
#   we can recreate this run if we want to later.
rng = Lot()
mcmc.rng = rng


# The following are general MCMC settings (not steppingstone sampling related).
mcmc.ncycles = 10000    # run the sample 10,000 generations at each value of beta 
mcmc.sample_every = 20  # sample the tree every 10 cycles
mcmc.report_every = 100 # report the status every 100 cycles
mcmc.adapt_first = 20   # tune the slice sampler after 10 cycles
mcmc.verbose = True
mcmc.nchains = 1        # MCMCMC multiple chain analyses are not yet working


# specify the starting tree. We'll use the model tree just so we can get 
#   decent results from a very short run in the lab.
mcmc.starting_tree_source = TreeCollection(newick='(t1:0.074783,(((((((t2:0.017390,t3:0.013871):0.014610,t4:0.038516):0.023079,t5:0.033673):0.011254,((t6:0.000718,t7:0):0.001124,t8:0.008383):0.020560):0.022048,(t14:0.127203,t15:0.017124):0.024730):0.002278,((t9:0.048058,t10:0.059900):0.024875,((t11:0.038944,t12:0.025554):0.014442,t13:0.034559):0.022475):0.002477):0.006325,t16:0.023526):0.006762,t17:0.025302)')


# fix the tree topology and use EdgeMove (which only change edge lengths) instead
# of LargetSimonMove (which changes both edge lengths and the tree topology)
mcmc.fix_topology = True


# The next lines control how often different moves are performed during MCMC. 
#   The integers are the number of times the move is used per-cycle
mcmc.slice_weight = 1 
mcmc.edge_move_weight = 1
mcmc.tree_scaler_weight = 1

# specify whether to update relative rates individually using slice sampling (True) 
# or jointly using a Metropolis-Hastings proposal (False)
model.update_relrates_separately = False
mcmc.rel_rate_weight        = 1

# specify whether to update base frequencies individually using slice sampling (True) 
# or jointly using a Metropolis-Hastings proposal (False)
model.update_freqs_separately = False
mcmc.state_freq_weight      = 1


# settings specific to path sampling and steppingstone sampling
# We'll be running MCMC with the likelihood raised to a power of beta for 
#   several values of beta.
ss.nbetavals = 20 # number of beta values
ss.minbeta = 0.0 # smallest beta value to be sampled. 0 means sampling from the prior.
ss.maxbeta = 1.0 # largest beta value to be sampled. 1 means sampling from the posterior.

# We pick the values of beta to use by calculating the quantiles of a 
#   distribution (confusingly, we use a Beta distribution).  The
#   next parameters specify a Beta(0.3, 1.0) distribution which
ss.shape1 = 0.3 # first shape parameter of the beta sampling distribution
ss.shape2 = 1.0 # second shape parameter of the beta sampling distribution




# specify whether to use Metropolis-Hastings/slice sampling when exploring the prior (False)
# or to sample directly from the prior (True)
mcmc.draw_directly_from_prior = True


# specify tuning parameters for the Metropolist-Hastings proposals
# Note: there are two versions of each tuning parameter, one for the posterior
#   and one for the prior (denoted 0).
# When exploring distribtions intermediate between prior and posterior, an 
#   intermediate tuning parameter is used
# Using bolder moves for distributions that are more similar to the prior
#   improves mixing, so that we fully explore these `heated` landscapes
mcmc.tree_scaler_lambda = 0.5 # min_lambda
mcmc.tree_scaler_lambda0 = 1.0 # max_lambda
mcmc.rel_rate_psi = 200.0 # max_psi
mcmc.rel_rate_psi0 = 2.0 # min_psi
mcmc.state_freq_psi = 200.0 # max_psi
mcmc.state_freq_psi0 = 2.0 # min_psi
mcmc.edge_move_lambda     = 0.5     # min_lambda
mcmc.edge_move_lambda0    = 1.0     # max_lambda


# This will be the prefix for the name of output files generated.
#   In the prefix we`ll indicate the model and whether or not we are 
#   running a steppingstone analysis or a harmonic mean run so that we can tell
#   what we ran later.
if using_steppingstone:
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

sump.out.log.prefix = 'sump' + fnprefix

if using_steppingstone:
    # start the MCMC analysis that will generate samples appropriate
    #   for both path sampling and steppingstone sampling
    ss()
else:
    # If we want to compare to the harmonic mean estimator we should at least
    #   run that sampler an approximately similar amount of time.  Since we
    #   will use the mcmc() function rather than ss(), the sampler will
    #   not be examining different beta values.  It will only sample the 
    #   posterior.  So we'll multiple the number of cycles by our 
    #   nbetavals.  If you have not changed the number above, this means that our
    #   MCMC will be 20 times as many cycles as the MCMC run on the posterior that
    #   we conducted in the steppingstone version.
    mcmc.ncycles *= ss.nbetavals
    mcmc()


# The sump command actually does the steppingstone sampling (it does
#   importance sampling of the likelihoods for the stepping stones) and the
#   path sampling.
# It will perform these calculations if the .p files have the format that
#   indicates that we have done several beta values.
sump_filename = fnprefix + ".p"
sump(burnin=0,file=sump_filename)
