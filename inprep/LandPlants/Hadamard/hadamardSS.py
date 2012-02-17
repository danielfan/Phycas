from phycas import *

setMasterSeed(284197)

# specify GTR+G model
model.type = 'gtr'
model.num_rates = 4
model.pinvar_model = False

# prior specification
model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))
model.gamma_shape_prior = Exponential(1.0)
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

# specify data file
mcmc.data_source = 'hadamard2.nex'

# specify log file
mcmc.out.log = 'sshad2.log'
mcmc.out.log.mode = REPLACE

# specify file for parameter samples
mcmc.out.trees.prefix = 'sshad2'
mcmc.out.trees.mode = REPLACE

# specify file for tree samples
mcmc.out.params.prefix = 'sshad2'
mcmc.out.params.mode = REPLACE

# fix tree topology to ML tree
mcmc.fix_topology = True
mcmc.starting_tree_source = TreeCollection(newick='((((2:0.42549850,6:0.85640725):0.28678670,5:0.25718078):0.09449950,4:0.40359241):0.22071671,3:0.21503818,1:0.67608638)')
mcmc.ls_move_weight     = 0.0

# edge length move details
mcmc.edge_move_weight   = 100.0
mcmc.edge_move_lambda   = 0.2
mcmc.edge_move_lambda   = 0.2

# mcmc details
mcmc.ncycles = 5000            # 2000
mcmc.sample_every = 50
mcmc.report_every = 100

# ss details
ss.nbetavals = 21
ss.ti = False
ss.xcycles = 1000                # 1000
ss()

sump.file = 'sshad2.p'
sump.burnin = 1
sump.out.log.prefix = 'sshad2_sump'
sump.out.log.mode = REPLACE
sump()

