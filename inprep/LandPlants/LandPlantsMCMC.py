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

# updater details
model.update_freqs_separately = False
model.update_relrates_separately = False
mcmc.state_freq_psi = 3000.0
mcmc.rel_rate_psi = 3000.0

# scale tree during mcmc
mcmc.tree_scaler_weight = 1.0

# specify data file
mcmc.data_source = 'Yang_and_Rannala_Karol_nomissambig.nex'

# specify log file
mcmc.out.log = 'landplants.log'
mcmc.out.log.mode = REPLACE

# specify file for parameter samples
mcmc.out.trees.prefix = 'landplants'
mcmc.out.trees.mode = REPLACE

# specify file for tree samples
mcmc.out.params.prefix = 'landplants'
mcmc.out.params.mode = REPLACE

# fix tree topology to ML tree
mcmc.fix_topology = True
mcmc.starting_tree_source = TreeCollection(filename='gtrgml.tre')
mcmc.ls_move_weight     = 0.0
mcmc.edge_move_weight   = 100.0

# mcmc details
mcmc.ncycles = 20000
mcmc.sample_every = 20
mcmc()
