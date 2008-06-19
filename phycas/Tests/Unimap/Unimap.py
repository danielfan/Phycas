from phycas import *

phycas                       = Phycas()

phycas.outfile_prefix        = 'unitest'
phycas.default_model         = 'hky'
phycas.random_seed           = 19375
phycas.data_file_name        = 'green.nex'
phycas.starting_tree_source  = 'random'
phycas.ncycles               = 100
phycas.sample_every          = 1
phycas.report_every          = 1
phycas.ls_move_weight        = 10
phycas.tree_scaler_weight    = 0
phycas.nielsen_move_weight   = 1
phycas.slice_weight          = 1
phycas.debugging             = False
phycas.use_unimap            = True

import sys,os
if os.path.basename(sys.executable) == 'python_d.exe':
    raw_input('debug stop')

phycas.mcmc()

#phycas.unimap()
#phycas.data_file_name       = 'green10sites.nex'
#phycas.starting_tree_source = 'usertree'
#phycas.tree_topology        = '(1:0.133913,2:0.053526,(3:0.084353,(4:0.071498,((5:0.064439,8:0.125345):0.032875,(6:0.071272,((7:0.045377,10:0.072536):0.018078,9:0.057064):0.037863):0.024930):0.025436):0.029039):0.022489)'
