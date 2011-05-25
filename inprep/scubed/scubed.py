from phycas import *

rng = ProbDist.Lot()
rng.setSeed(98765)
n = Newick('(1,2,(3,(4,((5,8),(6,((7,10),9))))))')

model.type                          = 'gtr'

# base freqs
model.update_freqs_separately       = False
model.state_freqs                   = [0.25, 0.25, 0.25, 0.25]
model.state_freq_prior              = ProbDist.Dirichlet((1.0, 1.0, 1.0, 1.0))

# relative rates
model.update_relrates_separately    = False
model.relrates                      = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0]
model.relrate_prior                 = ProbDist.Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))

# discrete gamma rate heterogeneity
model.num_rates                     = 3
model.gamma_shape                   = 0.5
model.gamma_shape_prior             = ProbDist.Exponential(1.0)

# proportion invariable sites
model.pinvar_model                  = True
model.pinvar                        = 0.2
model.pinvar_prior                  = ProbDist.Beta(1.0, 1.0)

# edge lengths
model.edgelen_hyperprior            = None
model.edgelen_prior                 = ProbDist.Exponential(10.0)

mcmc.rng                  = rng
mcmc.burnin               = 1000
mcmc.ncycles              = 10000
mcmc.sample_every         = 10
mcmc.report_every         = 100
mcmc.data_source          = getPhycasTestData('green.nex')
mcmc.starting_tree_source = randomtree(newick=n, rng=rng)
mcmc.fix_topology         = True

mcmc.out.log              = 'sss.log'
mcmc.out.log.mode         = REPLACE

mcmc.out.trees            = 'sss.t'
mcmc.out.trees.mode       = REPLACE

mcmc.out.params           = 'sss.p'
mcmc.out.params.mode      = REPLACE

mcmc.draw_directly_from_prior = True
mcmc.ls_move_weight       = 100
mcmc.ls_move_lambda       = 0.2
mcmc.ls_move_lambda0      = 0.2

mcmc.edge_move_weight     = 10
mcmc.edge_move_lambda     = 0.2
mcmc.edge_move_lambda0    = 0.2

mcmc.state_freq_weight    = 1
mcmc.state_freq_psi       = 1000.0
mcmc.state_freq_psi0      = 1000.0

mcmc.rel_rate_weight      = 1
mcmc.rel_rate_psi         = 1000.0 
mcmc.rel_rate_psi0        = 1000.0

mcmc.tree_scaler_weight   = 1
mcmc.tree_scaler_lambda   = 0.2
mcmc.tree_scaler_lambda0  = 0.2

ss.nbetavals              = 3
ss.maxbeta                = 1.0
ss.minbeta                = 0.0
ss.shape1                 = 1
ss.shape2                 = 1
ss.ti                     = False
#ss()

sump.file                 = 'sss.p'
sump.out.log              = 'sss.sump.log'
sump.out.log.mode         = REPLACE
sump()

