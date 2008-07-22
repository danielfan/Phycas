from phycas import *
#phycassertRaisesException = True

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

model.fix_edgelens                 = True
model.edgelen_dist                 = ProbDist.ExponentialDist(10.0)
model.starting_edgelen_hyperparam  = 0.05
model.fix_edgelen_hyperparam       = True
model.edgelen_hyperprior           = ProbDist.InverseGammaDist(2.1, 0.9090909)

mcmc.out.log                      = 'output.txt'
mcmc.out.log.mode                 = REPLACE
mcmc.out.trees                    = 'fixed.t'
mcmc.out.trees.mode               = REPLACE
mcmc.out.params                   = 'fixed.p'
mcmc.out.params.mode              = REPLACE
mcmc.nchains                      = 1
mcmc.ncycles                      = 2500
mcmc.random_seed                  = 13579
mcmc.data_source                  = 'file'
mcmc.data_file_name               = getPhycasTestData('nyldna4.nex')
mcmc.starting_tree_source         = 'random'

import sys,os
if os.path.basename(sys.executable) == 'python_d.exe':
    raw_input('debug stop')

mcmc()
