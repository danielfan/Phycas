from phycas import *

rng = ProbDist.Lot()
rng.setSeed(98765)
n = Newick('(1,2,(3,(4,((5,8),(6,((7,10),9))))))')

model.type                = 'jc'
model.edgelen_hyperprior  = None
model.edgelen_dist        = ProbDist.Exponential(10.0)

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

ps.nbetavals              = 5
ps.maxbeta                = 1.0
ps.minbeta                = 0.0
ps()
