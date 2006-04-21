# This example is designed to demonstrate that the MCMC code is working correctly by
# exploring only the priors. The posterior resulting from such an analysis should
# approximate the joint prior distribution.

import Conversions
from Phycas import *

if __name__ == "__main__":
    phycas = Phycas()

    phycas.no_data                  = True
    phycas.data_file_name           = 'nodata.nex'  # used for creating *.p and *.t file names
    phycas.ntax                     = 10 
    phycas.starting_tree_source     = 'random'
    phycas.ncycles                  = 200000
    phycas.sample_every             = 100
    phycas.report_every             = phycas.ncycles//100
    phycas.adapt_first              = 100
    phycas.random_seed              = '463859'
    phycas.model_type               = 'hky'
    phycas.using_hyperprior         = True
    phycas.num_rates                = 4
    phycas.gamma_shape_mean         = 1.0
    phycas.edgelen_prior_mean       = 1.0
    phycas.verbose                  = True
    phycas.metropolis_weight        = 100
    phycas.slice_weight             = 1
    phycas.slice_max_units          = 0
    phycas.gg_do                    = False
    phycas.use_inverse_shape        = False

    phycas.setup()
    phycas.run()

