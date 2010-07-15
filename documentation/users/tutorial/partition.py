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

# Start the run 
mcmc()

# Summarize the posterior 
sumt.trees = 'trees.t' 
sumt()
