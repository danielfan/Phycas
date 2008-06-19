from phycas import *

phycas                              = Phycas()

phycas.nchains                      = 1

phycas.random_seed                  = 13579
phycas.default_model                = 'hky'

phycas.starting_kappa               = 4.0
phycas.fix_kappa                    = True
phycas.kappa_prior                  = ProbDist.ExponentialDist(1.0)

phycas.starting_freqs               = [0.1, 0.2, 0.3, 0.4]
phycas.fix_freqs                    = True
phycas.base_freq_param_prior        = ProbDist.ExponentialDist(1.0)

phycas.num_rates                    = 4
phycas.starting_shape               = 0.14
phycas.fix_shape                    = True
phycas.use_inverse_shape            = False
phycas.gamma_shape_prior            = ProbDist.ExponentialDist(1.0)

phycas.estimate_pinvar              = True
phycas.starting_pinvar              = 0.27
phycas.fix_pinvar                   = True
phycas.pinvar_prior                 = ProbDist.BetaDist(1.0, 1.0)

phycas.fix_edgelens                 = True
phycas.edgelen_dist                 = ProbDist.ExponentialDist(10.0)

phycas.using_hyperprior             = True
phycas.starting_edgelen_hyperparam  = 0.05
phycas.fix_edgelen_hyperparam       = True
phycas.edgelen_hyperprior           = ProbDist.InverseGammaDist(2.1, 0.9090909)

phycas.data_source                  = 'file'
phycas.data_file_name               = '../Data/nyldna4.nex'

phycas.starting_tree_source         = 'random'

phycas.outfile_prefix               = 'fixed'

phycas.ncycles                      = 2500

import sys,os
if os.path.basename(sys.executable) == 'python_d.exe':
    raw_input('debug stop')

phycas.mcmc()
