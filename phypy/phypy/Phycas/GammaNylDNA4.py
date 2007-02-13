from Phycas import *

phycas = Phycas()

phycas.data_file_name = '../../phypy/nyldna4.nex'
phycas.starting_tree_source = 'random'
phycas.ncycles = 1500
phycas.sample_every = 30
phycas.adapt_first = 100
phycas.random_seed = '13579'
phycas.default_model = 'hky'
phycas.num_rates = 4
phycas.using_hyperprior = True
phycas.edgelen_prior_mean = 0.1
phycas.verbose = True
phycas.ls_move_weight = 300
phycas.slice_weight = 1
phycas.gg_do = False

phycas.setup()
phycas.run()
