#!/usr/bin/env python

from phycas import *

print help()
print help(readFile)
print help(mcmc)
print help(readFile)
mcmc.current()
mcmc.data_source = "green.nex"
mcmc.current()
mcmc.out.log = "basic.log"
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = "green"
mcmc.out.params.prefix = "green"
mcmc.ncycles = 2000 
mcmc.sample_every = 10
mcmc.current()
mcmc()
sumt.current()
sumt.trees = "green.t"
sumt.current()
sumt()
