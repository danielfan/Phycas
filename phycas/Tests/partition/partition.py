from phycas import *

if True:
    print '==> Setting up partition model...'
    
    model.type = 'gtr'
    model.relrates = [0.126648, 0.077024, 0.053009, 0.052451, 0.660005, 0.030864]
    model.state_freqs = [0.204130, 0.257487, 0.255923, 0.282460]
    model.update_freqs_separately = True
    model.update_relrates_separately = True
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None
    m1 = model()
    
    model.type = 'gtr'
    model.relrates = [0.013481, 0.038521, 0.046651, 0.105936, 0.758869, 0.036542]
    model.state_freqs = [0.179386, 0.345455, 0.127428, 0.347730]
    model.update_freqs_separately = True
    model.update_relrates_separately = True
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None
    m2 = model()
    
    model.type = 'gtr'
    model.relrates = [0.067018, 0.346256, 0.055179, 0.070112, 0.402758, 0.058677]
    model.state_freqs = [0.340719, 0.103901, 0.115752, 0.439628]
    model.update_freqs_separately = True
    model.update_relrates_separately = True
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None
    m3 = model()
    
    partition.addSubset(subset(1,1296,3), m1, 'first')
    partition.addSubset(subset(2,1296,3), m2, 'second')
    partition.addSubset(subset(3,1296,3), m3, 'third')
    partition()
else:
    print '==> Setting up model...'
    
    model.type = 'gtr'
    model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
    model.update_freqs_separately = False
    model.update_relrates_separately = False
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None

setMasterSeed(13579)
filename = getPhycasTestData('green.nex')
blob = readFile(filename)

mcmc.data_source = blob.characters
mcmc.out.log = 'parttest.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'parttest'
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = 'parttest'
mcmc.out.params.mode = REPLACE
mcmc.nchains = 1
mcmc.ncycles = 100
mcmc.adapt_first = 101
mcmc.sample_every = 1
mcmc.report_every = 10
mcmc.ls_move_weight = 1
mcmc.subset_relrates_weight = 1
mcmc.tree_scaler_weight = 1
mcmc.slice_weight = 1
mcmc.starting_tree_source = TreeCollection(newick='((((((9:0.1,(7:0.1,10:0.1):0.1):0.1,6:0.1):0.1,(5:0.1,8:0.1):0.1):0.1,4:0.1):0.1,3:0.1):0.1,2:0.1,1:0.1);') 
mcmc()
