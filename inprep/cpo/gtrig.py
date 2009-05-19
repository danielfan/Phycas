from phycas import *

rng = ProbDist.Lot()
rng.setSeed(13579)

model.type = 'gtr'
model.num_rates = 4
model.pinvar_model = True
model.use_inverse_shape = False

model.relrates = [1.0, 4.0, 1.0, 1.0, 8.0, 1.0]
model.state_freqs = [0.2, 0.3, 0.3, 0.2]
model.gamma_shape = 0.5
model.pinvar = 0.2

randomtree.rng = rng
randomtree.taxon_labels = ['t%d' % (x+1) for x in range(50)] 
randomtree.edgelen_dist = Exponential(20.0)
randomtree.n_trees = 0

sim.file_name = 'GTRIGtrue.nex'
sim.nchar = 30000
sim.rng = rng
sim.tree_source = randomtree()
sim.taxon_labels = randomtree.taxon_labels

sim()