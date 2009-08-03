# This script will draw trees and model parameters from the same 
import math, sys
from phycas import *

#Toggling the following variables determines whether or not you use:
# GTR vs HKY, and
use_gtr = True
print_paup_commands = True
# specify the data file name
mcmc.data_source = 'simhky.nex'


# settings relating to substitution model

if use_gtr:
    model.type = 'gtr'
    relrate_prior = ProbDist.Dirichlet((1.0,1.0,1.0,1.0,1.0,1.0)) 
else:
    model.type = 'hky'
    kappa_prior = Exponential(1.0)

model.num_rates = 1 # use number > 1 for discrete gamma rate heterogeneity
model.pinvar_model = False # use True for proportion invariable sites model
if model.num_rates > 1:
    model.gamma_shape = 0.5

# settings relating to prior distributions
model.edgelen_hyperprior = None




state_freq_prior = Dirichlet((1.00000, 1.00000, 1.00000, 1.00000))


file_contents = readFile("simhky.nex")
like.data_source = file_contents.characters

rand_tree_source = randomtree(distribution='Equiprobable', edgelen_dist=ProbDist.Exponential(1.0))

n_trees = 10

if print_paup_commands:
    sys.stderr.write("Set storebr;\n")

for ind, tree in enumerate(rand_tree_source):
    like.tree_source = TreeCollection(trees=[tree])
    model.state_freqs = list(state_freq_prior.sample())
    if use_gtr:
        model.relrates = list(relrate_prior.sample())
    else:
        model.kappa = kappa_prior.sample()

    if print_paup_commands:
        sys.stderr.write("begin trees; tree tree" + str(ind) + " = [&U] " + tree.makeNewick() + ";\nend;\n")
        if use_gtr:
            gt = model.relrates[5]
            rr = [str(i/gt) for i in model.relrates[:5]]
        else:
            rr = ['1', str(model.kappa), '1', '1', str(model.kappa)]
        freq_sum = sum(model.state_freqs)
        fr = [str(i/freq_sum) for i in model.state_freqs[:3]]
        sys.stderr.write("lscore / userbr nst=6 rmat = (%s) basefreq = (%s);\n" % (" ".join(rr), " ".join(fr)))
        
    print like()
    if ind >= 10:
        break


x = """
# probably best to *not* change any of these
mcmc.nchains = 1 # multiple chain analyses are not yet working
mcmc.slice_weight = 1 # means one slice sampling update of each model parameter per cycle

# settings specific to path/steppingstone sampling
ss.nbetavals = 10 # number of beta values 
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
mcmc.ls_move_lambda = 0.5 # min_lambda
mcmc.ls_move_lambda0 = 1.0 # max_lambda

# it is pretty easy to infer the tree, so we'll turn the larget-simon move weight down
# from 100 (the default to 10)
mcmc.ls_move_weight = 10

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

sump.out.log.prefix = fnprefix

if using_stepping_stone:
    # start the MCMC analysis that will generate samples appropriate
    # for both path sampling and steppingstone sampling
    ss()
else:
    mcmc.ncycles *= ss.nbetavals
    mcmc()


sump(burnin=1,file="%s.p" % fnprefix)



"""
