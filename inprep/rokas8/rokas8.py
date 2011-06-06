import sys
from phycas import *

jobid = int(sys.argv[1])
setMasterSeed(13*jobid)

model.edgelen_hyperprior = None
model.edgelen_prior = Exponential(10.0)

#model.update_freqs_separately = False
#model.update_relrates_separately = False

mcmc.data_source = 'rokas8T1000C.nex'
mcmc.out.log = 'output.%d.txt' % (jobid,)
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = 'trees.%d' % (jobid,)
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = 'params.%d' % (jobid,)
mcmc.out.params.mode = REPLACE
mcmc.ncycles = 10000
mcmc.sample_every = 10
mcmc.report_every = 1000

model.type = 'gtr'
model.state_freqs = [0.25, 0.25, 0.25, 0.25]
model.state_freq_prior = Dirichlet((1.0,1.0,1.0,1.0))
model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
model.relrates = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0]
model.num_rates = 4
model.gamma_shape = 0.5
model.gamma_shape_prior = Exponential(1.0)
model.pinvar_model = False
model.pinvar = 0.5
model.pinvar_prior = Beta(1.0, 1.0)

mcmc.fix_topology = True
mcmc.allow_polytomies = False

tree_descriptions = [
    TREE_DESCRIPTION_HERE
]
newick = tree_descriptions[0]
import re
p = re.compile(r'\(1,\((.*)\)\)')
m = p.match(newick)
assert(m)
g = m.group(1)
newick='(1,' + g + ')'
print 'newick='+newick
b = []
for letter in newick[:-1]:
    b.append(letter)
    if letter.isdigit() or letter == ')':
        b.append(':0.05')
b.append(')')
newick = ''.join(b)

mcmc.starting_tree_source = TreeCollection(newick=newick)

ss.nbetavals = 25
ss()

sump.out.log = 'ss.sump.%d.txt' % (jobid,)
sump.out.log.mode = REPLACE
sump.file = 'params.%d.p' % (jobid,)
sump()
