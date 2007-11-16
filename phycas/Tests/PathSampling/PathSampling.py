from phycas import *

nchains = 4
prior_mean = 1.0

phycas = Phycas()
phycas.outfile_prefix = 'path_sampling'
phycas.nchains = 4
phycas.is_standard_heating = False
phycas.ncycles = 10000
phycas.sample_every = 5
phycas.report_every = 200
phycas.adapt_first = 10
phycas.verbose = True
phycas.ls_move_weight = 100
phycas.tree_scaler_weight = 1
phycas.slice_weight = 1
phycas.using_hyperprior = False
phycas.edgelen_dist = ProbDist.ExponentialDist(1.0)
phycas.starting_tree_source = 'random'
phycas.default_model = 'jc'
phycas.num_rates = 1
phycas.estimate_pinvar = False
phycas.use_flex_model = False
phycas.random_seed = 13579
phycas.data_file_name = '../Data/nyldna4.nex'
phycas.mcmc()
print 'nchains =',nchains
print 'prior_mean =',prior_mean
