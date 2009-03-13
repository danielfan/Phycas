# This example program performs the simplest possible Bayesian MCMC analysis,
# where simplest means that default settings are used.

# Import the phycas package into the Python interpreter
from phycas import *

# Specify the data file
mcmc.data_source = 'green.nex'

# Perform an MCMC analysis using all default settings
mcmc()
