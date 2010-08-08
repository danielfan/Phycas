################################################################################
# The code hear is the steppingstone sampling portion of a full Phycas analysis.
#
# It assumes that the data have been read and the mcmc output files have been
#   specified 
#
# The code conducts the stepping stone analysis, and runs sump to calculate the
#   approximate value for the marginal likelihood of the data.
#
# Note that this code makes NO assumptions about what model has been set up. 
#   this means that we can reuse this code for sampling under any model we are
#   curious about.
############


################################################################################
# Currently, the steppingstone method works whith fixed trees only!
############
mcmc.fix_topology = True
 
################################################################################
# mcmc.ncycles is the number of iterations that the MCMC sampler will be run for 
#   for each level of the beta power
#
# 50 is NOT a long enough run, we have chosen 50 to make the lab run in a 
#   reasonable amount of time!
############
mcmc.ncycles = 50
mcmc.sample_every = 1

################################################################################
# Because we are updating:
#   - the relative rates jointly (model.update_relrates_separately = False)
#   - the base frequencies jointly (model.update_freqs_separately = False)
# we are NOT using slice sampling for these parameters.  So we can adjust the
# magnitude of steps that we take in Metropolis-Hastings.  Large psi values lead
# to smaller steps.  The default (300) usually leads to steps that are too big,
# so here we will make less-drastic proposals.
#########
mcmc.rel_rate_psi = 3000
mcmc.state_freq_psi = 3000



################################################################################
# Now we'll use the analysis tag that was created when we partitioned the data
#   to set appropriately name mcmc output files.
############
mcmc.out.log = analysis_tag + '.log'
mcmc.out.trees = analysis_tag + '.t'
mcmc.out.params = analysis_tag + '.p'




################################################################################
# We raise the target distribution density to the power beta in each portion of
#   the sampler.
# Using 6 beta values with shape1 = 1 (and shape2 left at its default of 1.0)
#   samples the beta values from the quantiles of a Beta(1, 1.0) distribution.
#   This is simply a uniform distribution.
# This results in runs with:
#    beta = 1  (in which the posterior density is sampled)
#    beta = 0.8
#    beta = 0.6
#    beta = 0.4
#    beta = 0.2
#    beta = 0 (in which the reference distribution is sampled).
############
ss.nbetavals = 6

################################################################################
# We only conduct "burn-in" for the initial run of the sample ss.xcycles specifies
#   the length of this burnin.  
############
ss.xcycles = 50
ss()

################################################################################
# mcmc.ncycles is the number of iterations that the MCMC sampler will be run for 
#   for each level of the beta power
#
# 50 is NOT a long enough run, we have chosen 50 to make the lab run in a 
#   reasonable amount of time!
############
sump.out.log =  'sump_' + mcmc.out.log.filename
sump.out.log.mode = REPLACE
sump.file = mcmc.out.params.filename
sump()
