from phycas import *

#nchains = 8
prior_mean = 0.0001

phycas = Phycas()
phycas.random_seed = 13579
phycas.data_file_name = '../Data/nyldna4.nomissambig.nex'
phycas.outfile_prefix = 'path_sampling'
#phycas.nchains = nchains
phycas.is_standard_heating = False
phycas.ncycles = 10000
phycas.sample_every = 100
phycas.report_every = 500
phycas.adapt_first = 10
phycas.verbose = True
phycas.ls_move_weight = 100
phycas.tree_scaler_weight = 1
phycas.slice_weight = 1
phycas.using_hyperprior = False
phycas.edgelen_dist = ProbDist.ExponentialDist(1.0/prior_mean)
phycas.starting_tree_source = 'random'
phycas.default_model = 'jc'
phycas.num_rates = 1
phycas.estimate_pinvar = False
phycas.use_flex_model = False
#phycas.mcmc()
phycas.ps_toward_posterior = False
phycas.ps_burnin = 1000
phycas.ps_Q = 100
phycas.ps_nbetaincr = 100
phycas.pathsampling()
#print 'nchains =',nchains
print 'prior_mean =',prior_mean
