from phycas import *

data_file_name = 'data186.nex'
master_seed = 97531
unpartitioned_analysis = False

setMasterSeed(master_seed)
blob = readFile(data_file_name)
nchar = blob.characters.getMatrix().n_char

tree_descriptions = []
for t in blob.trees:
    tree_descriptions.append(t.newick)
assert(len(tree_descriptions) > 0)

if unpartitioned_analysis:
    model.type = 'gtr'
    model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
    model.update_freqs_separately = False
    model.update_relrates_separately = False
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None

else:
    model.type = 'gtr'
    model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
    model.state_freqs = [0.25, 0.25, 0.25, 0.25]
    model.update_freqs_separately = False
    model.update_relrates_separately = False
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None
    m1 = model()
    
    model.type = 'gtr'
    model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
    model.state_freqs = [0.25, 0.25, 0.25, 0.25]
    model.update_freqs_separately = False
    model.update_relrates_separately = False
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None
    m2 = model()
    
    half_way = nchar/2
    partition.addSubset(subset(1,half_way),         m1, 'first_half')
    partition.addSubset(subset(half_way + 1,nchar), m2, 'last_half')
    partition()

mcmc.data_source = blob.characters
mcmc.out.log.prefix = data_file_name
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = data_file_name
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = data_file_name
mcmc.out.params.mode = REPLACE
mcmc.nchains = 1
mcmc.ncycles = 5000             # was 10000
mcmc.adapt_first = 2
mcmc.sample_every = 10
mcmc.report_every = 100
mcmc.state_freq_weight = 10     # was 1
mcmc.rel_rate_weight = 10       # was 1
mcmc.subset_relrates_weight = 1
mcmc.tree_scaler_weight = 1
mcmc.slice_weight = 1
mcmc.starting_tree_source = TreeCollection(newick=tree_descriptions[0]) 
mcmc.debugging = False

#mcmc.ls_move_weight = 100
mcmc.edge_move_weight = 50      # was 100
mcmc.fix_topology = True

ss.nbetavals = 11
ss()

sump.file = data_file_name + '.p'
sump.out.log = data_file_name + '.sump.txt'
sump.out.log.mode = REPLACE
sump()