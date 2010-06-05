from phycas import *

data_file_name   = 'rbcl738.nex'
rnseed = 215915
num_cycles = 100000
sample_freq = 10
print_freq = 10

def createModel():    
    # model type (jc, hky or gtr)
    model.type = 'jc'
    
    # GTR relative rates
    #model.update_relrates_separately = False
    #model.relrates = relrate_MLEs
    #model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
    
    # HKY model
    #model.kappa = 4.0
    #model.kappa_prior = ProbDist.BetaPrime(1.0, 1.0)
    
    # state frequencies
    #model.update_freqs_separately = False
    #model.state_freqs = (0.25, 0.25, 0.25, 0.25)
    #model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))
    
    # discrete gamma rate heterogeneity
    model.num_rates = 1
    #model.gamma_shape = 0.5
    #model.gamma_shape_prior = Exponential(1.0)
    
    # proportion of invariable sites
    model.pinvar_model = False
    #model.pinvar = pinvar_MLE
    #model.pinvar_prior = Beta(1.0, 1.0)
    
    # edge lengths
    model.edgelen_hyperprior = None    
    model.edgelen_prior = Exponential(10.0)

    # Unpartitioned: (All nucleotide positions: GTR+I+G)
    partition()
        
if __name__ == '__main__':
    # Set the seed for analyses
    rng = setMasterSeed(rnseed)
    
    fnprefix = '_uni' 
    
    # set up the partition model
    createModel()
    
    # read in the data
    blob = readFile(data_file_name)
    mcmc.data_source = blob.characters
    
    # set output file names
    mcmc.out.log.prefix = fnprefix
    mcmc.out.log.mode = REPLACE
    mcmc.out.trees.prefix = fnprefix
    mcmc.out.trees.mode = REPLACE
    mcmc.out.params.prefix = fnprefix
    mcmc.out.params.mode = REPLACE
    
    mcmc.allow_polytomies = False
    #mcmc.polytomy_prior = False
    #mcmc.topo_prior_C = 1.0
    
    mcmc.nchains = 1

    # unimap-specific settings
    mcmc.use_unimap = True
    mcmc.mapping_move_weight = 0
    mcmc.unimap_ls_move_weight = 1
    mcmc.unimap_nni_move_weight = 0
    mcmc.unimap_edge_move_weight = 0
    mcmc.unimap_thread_count = 2
    mcmc.unimap_sample_ambig_move_weight = 1
    mcmc.tree_scaler_weight     = 0
    mcmc.slice_weight           = 1
    mcmc.debugging              = False
    
    #mcmc.edge_move_weight = 1
    #mcmc.ls_move_weight = 100
    
    #mcmc.state_freq_weight = 10
    #mcmc.state_freq_psi = 500.0
    
    #mcmc.rel_rate_weight = 10
    #mcmc.rel_rate_psi = 500.0
    
    mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels))
    #mcmc.starting_tree_source = TreeCollection(newick=garli_stepwise_tree)
    mcmc.fix_topology = False
    mcmc.adapt_first = 2
    mcmc.ncycles = num_cycles
    mcmc.sample_every = sample_freq
    mcmc.report_every = print_freq
    mcmc.uf_num_edges = 100

    mcmc.reference_tree_source = TreeCollection(filename='garli.green.JC.best.all.tre')
    #mcmc.reference_tree_source = TreeCollection(filename='starting.tre')
    
    mcmc()
