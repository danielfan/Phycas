import sys,os
from phycas import *

p = Phycas()
try:
    d = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p.data_file_name = os.path.join(d,"Data", "green.nex")
except NameError:
    p.data_file_name = os.path.join("phycas", "Tests","Data", "green.nex")
p.default_model = 'jc'
p.ncycles = 100000
p.random_seed = 13579
p.samc()

