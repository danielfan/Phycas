import math
from phycas import *

model.type = 'gtr'
model.pinvar_model = True
model.edgelen_hyperprior = None
model.state_freqs = [0.339271, 0.154491, 0.134649, 0.371589]
model.relrates    = [1.144048, 5.419204, 0.454958, 1.766404, 5.546350, 1.0]
model.gamma_shape = 0.906291 
model.pinvar      = 0.442154
model.num_rates   = 4

filename = getPhycasTestData('rbcL50.nex')
blob = readFile(filename)
like.data_source = blob.characters
like.tree_source = TreeCollection(filename='gtrig.rbcL50.best.tre')
like.starting_edgelen_dist = None
like.uf_num_edges = 10
like.store_site_likes = True
lnL = like()

outf = open('output.txt', 'w')
outf.write('correct lnL = -18117.830737\n        lnL = %.6f\n\n' % lnL)
outf.close()
