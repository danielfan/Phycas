from phycas import *

if False:
    model.type = 'jc'
    jc = model()
    
    model.type = 'hky'
    hky = model()
    
    model.type = 'gtr'
    gtr = model()
    
    partition.addSubset(subset(1,1296,3), jc, 'first')
    partition.addSubset(subset(2,1296,3), hky, 'second')
    partition.addSubset(subset(3,1296,3), gtr, 'third')
    partition()

setMasterSeed(13579)
filename = getPhycasTestData('green.nex')
blob = readFile(filename)
mcmc.data_source = blob.characters
mcmc.out.log = 'output.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'gtr_test'
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = 'gtr_test'
mcmc.out.params.mode = REPLACE
mcmc.nchains = 1
mcmc.ncycles = 100
mcmc.sample_every = 5
mcmc.report_every = 20
mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels))
mcmc.data_source = blob.characters
mcmc.debugging = True
mcmc()
