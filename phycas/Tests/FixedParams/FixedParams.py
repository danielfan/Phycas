from phycas import *
#Phycas.PhycassertRaisesException = True

model.type                        = 'hky'

model.kappa                        = 4.0
model.fix_kappa                    = True
model.kappa_prior                  = ProbDist.ExponentialDist(1.0)

model.base_freqs                   = [0.1, 0.2, 0.3, 0.4]
model.fix_freqs                    = True
model.base_freq_param_prior        = ProbDist.ExponentialDist(1.0)

model.num_rates                    = 4
model.gamma_shape                  = 0.14
model.fix_shape                    = True
model.use_inverse_shape            = False
model.gamma_shape_prior            = ProbDist.ExponentialDist(1.0)

model.pinvar_model                 = True
model.pinvar                       = 0.27
model.fix_pinvar                   = True
model.pinvar_prior                 = ProbDist.BetaDist(1.0, 1.0)

mcmc.nchains                      = 1
mcmc.ncycles                      = 2500
mcmc.outfile_prefix               = 'fixed'
mcmc.random_seed                  = 13579
mcmc.fix_edgelens                 = True
mcmc.edgelen_dist                 = ProbDist.ExponentialDist(10.0)
mcmc.using_hyperprior             = True
mcmc.starting_edgelen_hyperparam  = 0.05
mcmc.fix_edgelen_hyperparam       = True
mcmc.edgelen_hyperprior           = ProbDist.InverseGammaDist(2.1, 0.9090909)
mcmc.data_source                  = 'file'
mcmc.data_file_name               = getPhycasTestData('nyldna4.nex')
mcmc.starting_tree_source         = 'random'

import sys,os
if os.path.basename(sys.executable) == 'python_d.exe':
    raw_input('debug stop')

mcmc()
