# This file does a short MCMC analysis under the GTR+I+G model

from phycas import *
#Phycas.PhycassertRaisesException = True

mcmc.outfile_prefix = 'gtr_test'
mcmc.nchains = 1
mcmc.ncycles = 100
mcmc.sample_every = 5
mcmc.report_every = 20
mcmc.adapt_first = 10
mcmc.verbose = True
mcmc.ls_move_weight = 100
mcmc.tree_scaler_weight = 1
mcmc.slice_weight = 1
mcmc.slice_max_units = 0
mcmc.starting_shape = 2.0
mcmc.use_inverse_shape = True
mcmc.using_hyperprior = True
mcmc.starting_tree_source = 'random'
mcmc.default_model = 'gtr'
mcmc.num_rates = 4
mcmc.estimate_pinvar = True
mcmc.use_flex_model = False
mcmc.random_seed = 13579
mcmc.data_file_name = '../Data/green.nex'

import sys,os
if os.path.basename(sys.executable) == 'python_d.exe':
    raw_input('debug stop')

mcmc()

