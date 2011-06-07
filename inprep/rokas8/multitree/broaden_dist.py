#!/usr/bin/env python
from phycas import *
import sys
coeff_var = float(sys.argv[1])
if len(sys.argv) > 2:
    fn = sys.argv[2]
    stream = open(fn, 'rU')
else:
    stream = sys.stdin
for line in stream:
    ls = line.split('=')
    if len(ls) == 1:
        print line.strip()
    else:
        key= ls[0].strip()
        old_dist = eval('='.join(ls[1:]).strip())
        mean = old_dist.getMean()
        old_sd = old_dist.getStdDev()
        sd = coeff_var*mean
        if sd > old_sd:
            new_dist = Gamma(mean=mean, std_dev=sd)
            comment = '# mean kept at ' + str(mean) + ' std_dev changed from '+ str(old_sd) + ' to ' + str(sd)
        else:
            sd = old_sd
            new_dist = old_dist
            comment = ''
        print key,'=', new_dist, comment
        
