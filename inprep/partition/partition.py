# Unpartitioned:
#	Phycas	-7606.49122446	TL = 1.70000000
#	PAUP*	-7606.491		TL = 1.70000
#	MrBayes -7606.491		TL = 1.700
# Partitioned:
#   Phycas  -6983.86059941
#	MrBayes -6983.863
#   PAUP*   -6983.861 = -1668.026 + -1422.016 + -3893.819
#            patterns      80          46         302    

from phycas import *

if partitioning:
	print '==> Setting up partition model...'
	
	if False:
		model.type = 'gtr'
		model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
		model.update_freqs_separately = False
		model.update_relrates_separately = False
		model.edgelen_prior = Exponential(10.0)
		model.edgelen_hyperprior = None

	else:
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

if False:
	like.data_source = blob.characters
	like.tree_source = TreeCollection(newick='((((((9:0.1,(7:0.1,10:0.1):0.1):0.1,6:0.1):0.1,(5:0.1,8:0.1):0.1):0.1,4:0.1):0.1,3:0.1):0.1,2:0.1,1:0.1);') 
	like.starting_edgelen_dist = None
	like.store_site_likes = True
	print '==> calling like() <=='
	lnL = like()
	print 'log-likelihood =',lnL
else:
	mcmc.data_source = blob.characters
	mcmc.out.log = 'parttest.txt'
	mcmc.out.log.mode = REPLACE
	mcmc.out.trees.prefix = 'parttest'
	mcmc.out.trees.mode = REPLACE
	mcmc.out.params.prefix = 'parttest'
	mcmc.out.params.mode = REPLACE
	mcmc.nchains = 1
	mcmc.ncycles = 20000
	mcmc.adapt_first = 2
	mcmc.sample_every = 10
	mcmc.report_every = 100
	mcmc.ls_move_weight = 100
	mcmc.subset_relrates_weight = 1
	mcmc.tree_scaler_weight = 1
	mcmc.slice_weight = 1
	#mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels))
	#mcmc.starting_tree_source = TreeCollection(newick='(1:0.258008,2:0.086726,(3:0.113174,(4:0.086740,((5:0.078662,8:0.183070):0.051180,(6:0.087103,((7:0.046703,10:0.098804):0.021001,9:0.072890):0.052088):0.034684):0.033875):0.044762):0.026713);') 
	mcmc.starting_tree_source = TreeCollection(newick='((((((9:0.1,(7:0.1,10:0.1):0.1):0.1,6:0.1):0.1,(5:0.1,8:0.1):0.1):0.1,4:0.1):0.1,3:0.1):0.1,2:0.1,1:0.1);') 
	mcmc.debugging = False
	mcmc()
