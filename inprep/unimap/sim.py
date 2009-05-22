from phycas import *

setMasterSeed(13579)

model.type = 'hky'
model.kappa = 4.0
model.state_freqs = [0.3, 0.2, 0.2, 0.3]
model.num_rates = 10
model.gamma_shape = 0.1
model.pinvar_model = False

sim.taxon_labels = ['A','B','C','D','E','F','G','H','I','J']
sim.tree_source = TreeCollection(newick='(A:0.05,(B:0.05,(C:0.05,(D:0.05,E:0.05):0.05):0.05):0.05,(F:0.05,(G:0.05,(H:0.05,(I:0.05,J:0.05):0.05):0.05):0.05):0.05)')
sim.nchar = 1000
sim.file_name = 'muchoratehet.nex'
sim()
