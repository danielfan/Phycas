from phycas import *

setMasterSeed(192847)

# specify GTR+G model
model.type = 'gtr'
model.num_rates = 4
model.pinvar_model = False

# prior specification
model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))
model.gamma_shape_prior = Lognormal(0.0, 1.0)
model.edgelen_prior = Exponential(1.0)
model.edgelen_hyperprior = None

# state frequency updater details
model.update_freqs_separately = False
mcmc.state_freq_psi = 3000.0
mcmc.state_freq_psi0 = 3000.0

# GTR relative rate updater details
model.update_relrates_separately = False
mcmc.rel_rate_psi = 3000.0
mcmc.rel_rate_psi0 = 3000.0

# scale tree during mcmc
mcmc.tree_scaler_weight  = 1.0
mcmc.tree_scaler_lambda  = 0.5
mcmc.tree_scaler_lambda0 = 0.5

# fix tree topology to ML tree
mcmc.fix_topology = True
mcmc.starting_tree_source = TreeCollection(filename='landplants.t')
mcmc.ls_move_weight     = 0.0

# edge length move details
mcmc.edge_move_weight   = 100.0
mcmc.edge_move_lambda   = 0.2
mcmc.edge_move_lambda   = 0.2

# specify data file
mcmc.data_source = 'Yang_and_Rannala_Karol_nomissambig.nex'

# specify log file
mcmc.out.log = 'sslandplants.log'
mcmc.out.log.mode = REPLACE

# specify file for parameter samples
mcmc.out.trees.prefix = 'sslandplants'
mcmc.out.trees.mode = REPLACE

# specify file for tree samples
mcmc.out.params.prefix = 'sslandplants'
mcmc.out.params.mode = REPLACE

# mcmc details
mcmc.ncycles = 1000
mcmc.sample_every = 10
mcmc.report_every = 100

# ss details
ss.nbetavals = 31
ss.ti = False
ss.xcycles = 1000
ss()

sump.file = 'sslandplants.p'
sump.burnin = 1
sump.out.log.prefix = 'sslandplants_sump'
sump.out.log.mode = REPLACE
sump()

