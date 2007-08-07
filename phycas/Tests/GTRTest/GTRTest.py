# This file does a short MCMC analysis under the GTR+I+G model

from phycas import *

phycas = Phycas()
phycas.outfile_prefix = 'gtr_test'
phycas.ncycles = 100
phycas.sample_every = 5 
phycas.report_every = 10 
phycas.adapt_first = 2
phycas.verbose = True
phycas.ls_move_weight = 100
phycas.tree_scaler_weight = 1
phycas.slice_weight = 1
phycas.slice_max_units = 0
phycas.use_inverse_shape = True
phycas.using_hyperprior = True
phycas.starting_tree_source = 'random'
phycas.default_model = 'gtr'
phycas.num_rates = 4
phycas.estimate_pinvar = True
phycas.use_flex_model = False
phycas.random_seed = 13579
phycas.data_file_name = '../Data/green.nex'
phycas.mcmc()

