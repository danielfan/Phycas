# This example program performs the simplest possible Bayesian MCMC analysis,
# where simplest means that default settings are used.

# Import the phycas package into the Python interpreter
from phycas import *

# Create a Phycas object
phycas = Phycas()

# Specify the data file
phycas.data_file_name = 'green.nex'

# Perform an MCMC analysis using all default settings
phycas.mcmc()
