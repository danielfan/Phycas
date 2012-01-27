from phycas import *

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

idr.burnin = 100
idr.params = 'landplants.p'
idr.trees = TreeCollection(filename='landplants.t')
idr.data_source = 'Yang_and_Rannala_Karol_nomissambig.nex'
idr.out.log.prefix = 'idrlog'
idr.out.log.mode = REPLACE
idr()
