from phycas import *

data_file_name = 'simple.nex'
master_seed = 18625
num_cycles = 22000
num_cycles_per_sample = 10

setMasterSeed(master_seed)
blob = readFile(data_file_name)
nchar = blob.characters.getMatrix().n_char

tree_descriptions = []
for t in blob.trees:
    tree_descriptions.append(t.newick)
assert(len(tree_descriptions) > 0)

model.type = 'gtr'
model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
model.state_freqs = [0.25, 0.25, 0.25, 0.25]
model.update_freqs_separately = False
model.update_relrates_separately = False
model.edgelen_prior = Exponential(10.0)
model.edgelen_hyperprior = None
model.pinvar_model = False
model.num_rates = 4
model.gamma_shape = 1.0
model.gamma_shape_prior = Exponential(1.0/100.0)
m1 = model()

model.type = 'gtr'
model.relrates = [1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0]
model.state_freqs = [0.25, 0.25, 0.25, 0.25]
model.update_freqs_separately = False
model.update_relrates_separately = False
model.edgelen_prior = Exponential(10.0)
model.edgelen_hyperprior = None
model.pinvar_model = False
model.num_rates = 1
model.gamma_shape = 1.0
model.gamma_shape_prior = Exponential(1.0/100.0)
m2 = model()

half_way = nchar/2
partition.subset_relrates = [1.0,1.0]
partition.fix_subset_relrates = False
partition.addSubset(subset(1,half_way),         m1, 'first_half')
partition.addSubset(subset(half_way + 1,nchar), m2, 'last_half')
partition()

#like.preorder_edgelens = [0.1, 0.1, 0.1, 0.1, 0.1]
like.data_source = blob.characters
like.tree_source = TreeCollection(newick=tree_descriptions[0])
lnl = like()
print 'lnL =',lnl
