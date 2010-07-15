from phycas import *

mcmc.data_source = 'green.nex'

# Set up K80+G model 
model.type="hky" 
model.state_freqs = [0.25, 0.25, 0.25, 0.25] 
model.fix_freqs = True
model.kappa = 2.0 
model.kappa_prior = BetaPrime(1.0, 1.0)
model.num_rates = 4
model.gamma_shape = 0.5
model.gamma_shape_prior = Exponential(1.0)

# Save the K80+G model 
m1 = model() 
m2 = model()

# Set up and save the HKY+G model 
model.fix_freqs = False 
m3 = model()

# Define partition subsets 
first = subset(1, 1296, 3) 
second = subset(2, 1296, 3) 
third = subset(3, 1296, 3)

mcmc.fix_topology = True
mcmc.starting_tree_source = TreeCollection(newick='''(1:0.31,2:0.10,(3:0.13,(4:0.09,
((5:0.10,6:0.20):0.06,(7:0.10,((8:0.05,10:0.11):0.02,9:0.08):0.06):0.04):0.04):0.05):0.03)''')
 
model.edgelen_hyperprior = None
model.edgelen_prior = Exponential(1.0)

ss.nbetavals = 11
mcmc.ncycles = 1000
ss.xcycles = 1000
ss()
sump.file = 'params.p'
sump()
