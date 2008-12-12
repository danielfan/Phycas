from phycas import *

rng = ProbDist.Lot()
rng.setSeed(98765)
n = Newick('(1,2,(3,(4,((5,8),(6,((7,10),9))))))')

model.type                = 'jc'
model.edgelen_hyperprior  = None
model.edgelen_prior       = ProbDist.Exponential(10.0)

mcmc.rng                  = rng
mcmc.burnin               = 100
mcmc.ncycles              = 100
mcmc.sample_every         = 10
mcmc.report_every         = 10
mcmc.data_source          = getPhycasTestData('green.nex')
mcmc.starting_tree_source = randomtree(newick=n, rng=rng)
mcmc.fix_topology         = True
mcmc.edge_move_weight     = 1

mcmc.out.log              = 'output.log'
mcmc.out.log.mode         = REPLACE

mcmc.out.trees            = 'trees.t'
mcmc.out.trees.mode       = REPLACE

mcmc.out.params           = 'params.p'
mcmc.out.params.mode      = REPLACE

mcmc.draw_directly_from_prior = True
mcmc.ls_move_lambda       = 0.2
mcmc.ls_move_lambda0      = 0.5
mcmc.edge_move_lambda     = 0.2
mcmc.edge_move_lambda0    = 0.5
mcmc.state_freq_psi       = 300.0
mcmc.state_freq_psi0      = 1.0
mcmc.rel_rate_psi         = 300.0 
mcmc.rel_rate_psi0        = 1.0
mcmc.tree_scaler_lambda   = 0.2
mcmc.tree_scaler_lambda0  = 0.5

ps.nbetavals              = 5
ps.maxbeta                = 1.0
ps.minbeta                = 0.0
ps.shape1                 = 1
ps.shape2                 = 1
ps()
