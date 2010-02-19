from phycas import *

data_file_name = 'data162.nex'
master_seed = 1627
partitioned_analysis = False
plus_gamma = False
scubed = True
num_beta_values = 11
num_cycles = 2000
num_cycles_per_sample = 10
estimate_subset_relative_rates = False

setMasterSeed(master_seed)
blob = readFile(data_file_name)
nchar = blob.characters.getMatrix().n_char

tree_descriptions = []
for t in blob.trees:
    tree_descriptions.append(t.newick)
assert(len(tree_descriptions) > 0)

if partitioned_analysis:
    model.type = 'gtr'
    model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
    model.state_freqs = [0.25, 0.25, 0.25, 0.25]
    model.update_freqs_separately = False
    model.update_relrates_separately = False
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None
    model.pinvar_model = False
    if plus_gamma:
        model.num_rates = 4
        model.gamma_shape = 1.0
        model.gamma_shape_prior = Exponential(1.0/100.0)
    else:
        model.num_rates = 1
    m1 = model()
    
    model.type = 'gtr'
    model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
    model.state_freqs = [0.25, 0.25, 0.25, 0.25]
    model.update_freqs_separately = False
    model.update_relrates_separately = False
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None
    model.pinvar_model = False
    if plus_gamma:
        model.num_rates = 4
        model.gamma_shape = 1.0
        model.gamma_shape_prior = Exponential(1.0/100.0)
    else:
        model.num_rates = 1
    m2 = model()
    
    half_way = nchar/2
    partition.subset_relrates = [1.0,1.0]
    if estimate_subset_relative_rates:
        partition.fix_subset_relrates = False
    else:
        partition.fix_subset_relrates = True
    partition.addSubset(subset(1,half_way),         m1, 'first_half')
    partition.addSubset(subset(half_way + 1,nchar), m2, 'last_half')
    partition()
	
else:
    model.type = 'gtr'
    model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
    model.update_freqs_separately = False
    model.update_relrates_separately = False
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None
    model.pinvar_model = False
    if plus_gamma:
        model.num_rates = 4
        model.gamma_shape = 1.0
        model.gamma_shape_prior = Exponential(1.0/100.0)
    else:
        model.num_rates = 1

if num_beta_values == 0:
    pfx = '%s_%s_%s_%s_%d' % (data_file_name,plus_gamma and 'plusG' or 'minusG',partitioned_analysis and 'partitioned' or 'unpartitioned',estimate_subset_relative_rates and 'ssrrvar' or 'ssrrfix',master_seed)
else:
    pfx = '%s_%s_%s_%s_%s_%d_%d' % (data_file_name,plus_gamma and 'plusG' or 'minusG',partitioned_analysis and 'partitioned' or 'unpartitioned',scubed and 'sss' or 'ss',estimate_subset_relative_rates and 'ssrrvar' or 'ssrrfix',num_beta_values,master_seed)
mcmc.data_source = blob.characters
mcmc.out.log.prefix = pfx
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = pfx
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = pfx
mcmc.out.params.mode = REPLACE
mcmc.nchains = 1
mcmc.adapt_first = 2
mcmc.sample_every = num_cycles_per_sample
mcmc.report_every = 100

mcmc.state_freq_weight = 10
if scubed:
	mcmc.state_freq_psi = 500.0
	mcmc.state_freq_psi0 = 500.0
else:
	mcmc.state_freq_psi = 500.0
	mcmc.state_freq_psi0 = 1.0
	
mcmc.rel_rate_weight = 10
if scubed:
	mcmc.rel_rate_psi = 500.0
	mcmc.rel_rate_psi0 = 500.0
else:
	mcmc.rel_rate_psi = 500.0
	mcmc.rel_rate_psi0 = 1.0

mcmc.subset_relrates_weight = 0
if scubed:
	mcmc.subset_relrates_psi = 300.0
	mcmc.subset_relrates_psi0 = 300.0
else:
	mcmc.subset_relrates_psi = 300.0
	mcmc.subset_relrates_psi0 = 1.0

mcmc.tree_scaler_weight = 1
if scubed:
	mcmc.tree_scaler_lambda = 0.5
	mcmc.tree_scaler_lambda0 = 0.5
else:
	mcmc.tree_scaler_lambda = 0.5
	mcmc.tree_scaler_lambda0 = 1.0

mcmc.slice_weight = 1

mcmc.starting_tree_source = TreeCollection(newick=tree_descriptions[0]) 
mcmc.debugging = False

#mcmc.ls_move_weight = 100
mcmc.edge_move_weight = 50
if scubed:
	mcmc.edge_move_lambda = 1.0
	mcmc.edge_move_lambda0 = 1.0
else:
	mcmc.edge_move_lambda = 1.0
	mcmc.edge_move_lambda0 = 2.0

mcmc.fix_topology = True

if num_beta_values == 0:
    mcmc.ncycles = 10*num_cycles
    mcmc()
else:
    mcmc.ncycles = num_cycles
    ss.nbetavals = num_beta_values
    if scubed:
        ss.shape1 = 1.0
        ss.shape2 = 1.0
        ss.scubed = True
    else:
        ss.shape1 = 0.3
        ss.shape2 = 1.0
        ss.scubed = False
    ss()

sump.file = pfx + '.p'
sump.out.log = pfx + '.sump.txt'
sump.out.log.mode = REPLACE
sump()

print 'Run %s has finished.' % pfx
