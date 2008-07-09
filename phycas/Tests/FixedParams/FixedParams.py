from phycas import *
#Phycas.PhycassertRaisesException = True

mcmc.nchains                      = 1

mcmc.random_seed                  = 13579
mcmc.default_model                = 'hky'

mcmc.starting_kappa               = 4.0
mcmc.fix_kappa                    = True
mcmc.kappa_prior                  = ProbDist.ExponentialDist(1.0)

mcmc.starting_freqs               = [0.1, 0.2, 0.3, 0.4]
mcmc.fix_freqs                    = True
mcmc.base_freq_param_prior        = ProbDist.ExponentialDist(1.0)

mcmc.num_rates                    = 4
mcmc.starting_shape               = 0.14
mcmc.fix_shape                    = True
mcmc.use_inverse_shape            = False
mcmc.gamma_shape_prior            = ProbDist.ExponentialDist(1.0)

mcmc.estimate_pinvar              = True
mcmc.starting_pinvar              = 0.27
mcmc.fix_pinvar                   = True
mcmc.pinvar_prior                 = ProbDist.BetaDist(1.0, 1.0)

mcmc.fix_edgelens                 = True
mcmc.edgelen_dist                 = ProbDist.ExponentialDist(10.0)

mcmc.using_hyperprior             = True
mcmc.starting_edgelen_hyperparam  = 0.05
mcmc.fix_edgelen_hyperparam       = True
mcmc.edgelen_hyperprior           = ProbDist.InverseGammaDist(2.1, 0.9090909)

mcmc.data_source                  = 'file'
mcmc.data_file_name               = '../Data/nyldna4.nex'

mcmc.starting_tree_source         = 'random'

mcmc.outfile_prefix               = 'fixed'

mcmc.ncycles                      = 2500

import sys,os
if os.path.basename(sys.executable) == 'python_d.exe':
    raw_input('debug stop')

mcmc()
