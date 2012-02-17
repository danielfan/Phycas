# This file does a short MCMC analysis under the codon model

from phycas import *
filename = getPhycasTestData('green.nex')
blob = readFile(filename)

setMasterSeed(12245)

model.type                      = 'jc'
model.num_rates                 = 1
model.pinvar_model              = False
model.edgelen_hyperprior        = None
model.edgelen_prior             = Exponential(1.0)

mcmc.nchains                    = 1
mcmc.ncycles                    = 20
mcmc.sample_every               = 10
mcmc.report_every               = 1

mcmc.ls_move_weight             = 100
mcmc.tree_scaler_weight         = 1
mcmc.slice_weight               = 1

mcmc.starting_tree_source       = TreeCollection(newick='(8:0.56880388,(((3:0.40888265,(1:1.03799510,2:0.41917430):0.03417782):0.16416599,4:0.29333306):0.14865078,(6:0.28599164,((7:0.14870266,10:0.32973086):0.06151508,9:0.24129778):0.17828009):0.11396143):0.15762955,5:0.29601916);') # randomtree(n_taxa=len(blob.taxon_labels), rng=rng)
mcmc.data_source                = blob.characters

mcmc.out.log                    = 'beagletest.txt'
mcmc.out.log.mode               = REPLACE

mcmc.out.params.prefix          = 'params'
mcmc.out.params.mode            = REPLACE

mcmc.out.trees.prefix           = 'trees'
mcmc.out.trees.mode             = REPLACE

mcmc.use_beaglelib              = True
mcmc()
