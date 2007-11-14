import sys,os
from phycas import *

# raw_input('debug stop')

p = Phycas()
try:
    d = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p.data_file_name = os.path.join(d,"Data", "green.nex")
except NameError:
    p.data_file_name = os.path.join("phycas", "Tests","Data", "green.nex")
p.default_model = 'jc'
p.ncycles = 10000
p.random_seed = 15397
p.outfile_prefix = 'samc_output'
p.samc_move_debug = True
p.samc_temperature = 0.1
p.samc()