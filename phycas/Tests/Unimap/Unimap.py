from phycas import *

model.type                  = 'hky'

mcmc.use_unimap             = True
mcmc.outfile_prefix         = 'unitest'
mcmc.random_seed            = 19375
mcmc.data_file_name         = 'green.nex'
mcmc.starting_tree_source   = 'random'
mcmc.ncycles                = 100
mcmc.sample_every           = 1
mcmc.report_every           = 10
mcmc.mapping_move_weight    = 1
mcmc.unimap_nni_move_weight = 1
mcmc.tree_scaler_weight     = 0
mcmc.slice_weight           = 1
mcmc.debugging              = True

if False:
	import sys,os
	if os.path.basename(sys.executable) == 'python_d.exe':
		raw_input('debug stop')

mcmc()

