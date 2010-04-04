from phycas import *

large = subset(1, 5000)
small = subset(5001,5100)

user_tree_def = ' (1:0.00156293,(((((((2:0.00325591,16:0.00402707):0.02205143,((20:0.00630496,29:0.00695889):0.01138659,28:0.01818355):0.00809410):0.01335766,(26:0.01188256,27:0.01250127):0.02960036):0.00967768,((((((3:0.00860238,15:0.01056681):0.00220042,4:0.01193905):0.00408679,5:0.01787676):0.00155750,((17:0.01386172,18:0.01057494):0.00315644,19:0.01190134):0.00143905):0.00324977,30:0.01220535):0.01044754,22:0.01935610):0.02073536):0.00355082,(((7:0.07257248,(25:0.07364494,(31:0.12962943,32:0.09693177):0.08359470):0.00544096):0.03392672,(((((8:0.01400035,(10:0.00393309,12:0.00684284):0.00452834):0.02117048,9:0.02142141):0.00187281,14:0.01577341):0.00403985,11:0.02125429):0.00374345,21:0.05267665):0.01590936):0.00295028,13:0.03355061):0.00389731):0.00448117,23:0.02567495):0.01624640,6:0.01444287):0.01588953,24:0.00385810)'
data_file_name = 'concatenated.nex'
fnprefix = 'subset_model'
sample_freq = 10
print_freq = 100
burn_in = 100
num_cycles = 5000

# model type (jc, hky or gtr)
model.type = 'jc'

# discrete gamma rate heterogeneity
model.num_rates = 1
model.gamma_shape = 0.5
model.gamma_shape_prior = Exponential(1.0)

# proportion of invariable sites
model.pinvar_model = False
model.pinvar = 0.2
model.pinvar_prior = Beta(1.0,1.0)

# edge lengths
model.edgelen_hyperprior = None    
model.edgelen_prior = Exponential(10.0)

# partition by gene
m1 = model()
m2 = model()
partition.addSubset(large, m1, 'large_slow')
partition.addSubset(small, m2, 'small_fast')
partition.fix_subset_relrates = False
partition.subset_relrates = [1.0, 1.0]
partition()

blob = readFile(data_file_name)
mcmc.data_source = blob.characters

# set output file names
mcmc.out.log.prefix = fnprefix
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = fnprefix
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = fnprefix
mcmc.out.params.mode = REPLACE

mcmc.nchains = 1

mcmc.edge_move_weight = 1
mcmc.ls_move_weight = 100

mcmc.subset_relrates_weight = 10
mcmc.subset_relrates_psi  = 1000.0
mcmc.subset_relrates_psi0 = 1000.0

mcmc.tree_scaler_weight = 0

mcmc.slice_weight = 1

#mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels))
mcmc.starting_tree_source = TreeCollection(newick=user_tree_def)
mcmc.debugging = False
mcmc.fix_topology = True
mcmc.adapt_first = 2
mcmc.sample_every = sample_freq
mcmc.report_every = print_freq

mcmc.burnin = burn_in

mcmc.ncycles = num_cycles
mcmc.uf_num_edges = 100
mcmc()

sump.file         = fnprefix + '.p'
sump.out.log      = fnprefix + '.sump.txt'
sump.out.log.mode = REPLACE
sump()
