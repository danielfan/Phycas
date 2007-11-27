import sys,os
from phycas import *

p = Phycas()
raw_input('debug stop')
p.data_file_name = "greenrbcl.nex"
p.default_model = 'jc'
p.ncycles = 1000000
p.random_seed = 15397
p.outfile_prefix = 'samc_output'
p.samc_temperature = 0.6
p.samc()