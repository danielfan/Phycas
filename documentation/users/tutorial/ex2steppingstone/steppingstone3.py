execfile('steppingstoneReadData.py')

ERROR!:  You need to pick what model to use based on the Bayes Factor that you
    can calculate from the previous analyses. To set the preferred 
    model, you will uncomment one of the two execfile statements below! and 
    remove this message from the file.
    
#execfile('setUpByCodonPosPartition.py')
#execfile('setUpByGeneAndCodonPosPartition.py')



mcmc.fix_topology = False
 
################################################################################
# mcmc.ncycles is the number of iterations that the MCMC sampler will be run for 
#   for each level of the beta power
#
# 1000 is NOT a long enough run, we have chosen 50 to make the lab run in a 
#   reasonable amount of time!
############
mcmc.out.params = analysis_tag + '.p'
mcmc.ncycles = 750
mcmc.sample_every = 5
mcmc.starting_tree_source = randomtree()
mcmc.allow_polytomies = False

analysis_tag = 'tree_estimation'

mcmc.out.log = analysis_tag + '.log'
mcmc.out.trees = analysis_tag + '.t'
mcmc()

sump.out.log =  'sump_' + mcmc.out.log.filename
sump.out.log.mode = REPLACE
sump.file = mcmc.out.params.filename
sump()

sumt.out.log = 'sumt_' + mcmc.out.log.filename
sumt.out.splits = 'sumt_splits_' + mcmc.out.log.filename
sumt.out.trees = 'sumt_trees_' + mcmc.out.log.filename
sumt.trees = mcmc.out.trees.filename
sumt()


