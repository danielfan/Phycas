#!/usr/bin/env python
from phycas import *
import sys
fn = sys.argv[1]
coeff_var = float(sys.argv[2])
for line in open(fn, 'rU'):
    ls = line.split('=')
    if len(ls) == 1:
        print line.strip()
    else:
        key= ls[0].strip()
        value = eval('='.join(ls[1:]).strip())
        mean = value.getMean()
        sd = coeff_var*mean
        new_dist = Gamma(mean=mean, std_dev=sd)
        print key,'=', new_dist
        
