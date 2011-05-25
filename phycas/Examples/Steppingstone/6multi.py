from phycas import *

method             = 'ss'     # 'ss' (steppingstone) or 'hm' (harmonic mean)
partition_scheme   = 'unpartitioned'  # 'unpartitioned', 'gene', 'codon', or 'genecodon'
num_beta_values    = 5
burn_in            = 0
cycles_per_beta    = 1000
xtra = 1 - cycles_per_beta
sample_freq        = 1
print_freq         = 10
rseed              = 13579

data_file_name     = 'sixtaxon.nex'

# HKY+G ml tree
user_tree_def = '((1:0.00156293,2:0.00325591):.1,(3:0.00156293,4:0.00325591):.1,(5:0.00156293,6:0.00325591):.1);'

def createModelForScheme():
    # Here are the charset definitions
    COI           = subset(   1,  774)
    COIfirst      = subset(   1,  774, 3)
    COIsecond     = subset(   2,  774, 3)
    COIthird      = subset(   3,  774, 3)
    COII          = subset( 775, 1476)
    COIIfirst     = subset( 775, 1476, 3)
    COIIsecond    = subset( 776, 1476, 3)
    COIIthird     = subset( 777, 1476, 3)
    ATPase8       = subset(1477, 1627)
    ATPase8first  = subset(1477, 1627, 3)
    ATPase8second = subset(1478, 1627, 3)
    ATPase8third  = subset(1479, 1627, 3)
    ATPase6       = subset(1628, 2090)
    ATPase6first  = subset(1629, 2090, 3)
    ATPase6second = subset(1630, 2090, 3)
    ATPase6third  = subset(1628, 2090, 3)
    firsts        = COIfirst  + COIIfirst  + ATPase8first  + ATPase6first
    seconds       = COIsecond + COIIsecond + ATPase8second + ATPase6second
    thirds        = COIthird  + COIIthird  + ATPase8third  + ATPase6third
    
    # For starting values, use MLEs based on unpartitioned GTR+G analysis of user_tree_def
    relrate_MLEs    = [4.28396, 34.10179,  1.39027,  2.53545, 40.60048,  1.00000]     # see note 1 below
    state_freq_MLEs = [0.376859, 0.096002, 0.094516, 0.432623]
    shape_MLE       = 0.146955
    pinvar_MLE      = 0.621227
        
    # model type (jc, hky or gtr)
    model.type = 'gtr'
    
    # GTR relative rates
    model.update_relrates_separately = False
    model.relrates = [x/sum(relrate_MLEs) for x in relrate_MLEs]     # see note 1 below
    model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
    
    # state frequencies
    model.update_freqs_separately = False
    model.state_freqs = state_freq_MLEs
    model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))
    
    # discrete gamma rate heterogeneity
    model.num_rates = 4
    model.gamma_shape = shape_MLE
    model.gamma_shape_prior = Exponential(1.0)
    
    # proportion of invariable sites
    model.pinvar_model = False
    
    # edge lengths
    model.edgelen_hyperprior = None    
    model.edgelen_prior = Exponential(40.0)
    #model.internal_edgelen_prior = Exponential(40.0)
    #model.external_edgelen_prior = Exponential(40.0)

    if partition_scheme == 'unpartitioned':
        # unpartitioned
        partition()
    elif partition_scheme == 'gene':
        # partition by gene
        m1 = model()
        m2 = model()
        m3 = model()
        m4 = model()
        partition.addSubset(COI,     m1, 'COI')
        partition.addSubset(COII,    m2, 'COII')
        partition.addSubset(ATPase8, m3, 'ATPase8')
        partition.addSubset(ATPase6, m4, 'ATPase6')
        partition()
    elif partition_scheme == 'codon':
        # partition by codon
        m1 = model()
        m2 = model()
        m3 = model()
        partition.addSubset(firsts,  m1, 'First codon positions')
        partition.addSubset(seconds, m2, 'Second codon positions')
        partition.addSubset(thirds,  m3, 'Third codon positions')
        partition()
    elif partition_scheme == 'genecodon':
        # partition by gene and codon position
        m1  = model()
        m2  = model()
        m3  = model()
        m4  = model()
        m5  = model()
        m6  = model()
        m7  = model()
        m8  = model()
        m9  = model()
        m10 = model()
        m11 = model()
        m12 = model()
        partition.addSubset(COIfirst,       m1, 'COIfirst')
        partition.addSubset(COIsecond,      m2, 'COIsecond')
        partition.addSubset(COIthird,       m3, 'COIthird')
        partition.addSubset(COIIfirst,      m4, 'COIIfirst')
        partition.addSubset(COIIsecond,     m5, 'COIIsecond')
        partition.addSubset(COIIthird,      m6, 'COIIthird')
        partition.addSubset(ATPase8first,   m7, 'ATPase8first')
        partition.addSubset(ATPase8second,  m8, 'ATPase8second')
        partition.addSubset(ATPase8third,   m9, 'ATPase8third')
        partition.addSubset(ATPase6first,  m10, 'ATPase6first')
        partition.addSubset(ATPase6second, m11, 'ATPase6second')
        partition.addSubset(ATPase6third,  m12, 'ATPase6third')
        partition()
    else:
        print 'Sorry, unrecognized scheme (%s)' % partition_scheme
        sys.exit()
        
def runMCMC(fnprefix):
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
    
    mcmc.nchains = 1
    
    mcmc.edge_move_weight = 1
    mcmc.ls_move_weight = 100
    
    mcmc.state_freq_weight = 10
    mcmc.state_freq_psi  = 1000.0
    mcmc.state_freq_psi0 = 1000.0
    
    mcmc.rel_rate_weight = 10
    mcmc.rel_rate_psi = 1000.0
    mcmc.rel_rate_psi = 1000.0
    
    mcmc.subset_relrates_weight = 10
    mcmc.subset_relrates_psi  = 1000.0
    mcmc.subset_relrates_psi0 = 1000.0
    
    mcmc.tree_scaler_weight = 1
    
    mcmc.slice_weight = 1
    
    #mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels))
    mcmc.starting_tree_source = TreeCollection(newick=user_tree_def)
    mcmc.debugging = False
    mcmc.fix_topology = False
    mcmc.adapt_first = 2
    mcmc.sample_every = sample_freq
    mcmc.report_every = print_freq
    
    mcmc.burnin = burn_in
    mcmc.uf_num_edges = 100
    
    
    if method == 'hm':
        # harmonic mean method
        mcmc.ncycles = cycles_per_beta*num_beta_values + xtra
        mcmc()
    elif method == 'ss':
        # steppingstone sampling method
        mcmc.ncycles = cycles_per_beta
        ss.override_fixed_topology_restriction = True
        ss.nbetavals = num_beta_values
        ss.shape1 = 1.0
        ss.shape2 = 1.0
        ss.ti = False
        ss.xcycles = xtra
        ss.refdist_definition_file = '6taxon_ref_dist.txt'

        #ss.minsample = min_wp_sample_size
        ss()
    else:
        print "method must be either 'hm' or 'ss'"
        sys.exit()
    
    sump.file         = fnprefix + '.p'
    sump.out.log      = fnprefix + '.sump.txt'
    sump.out.log.mode = REPLACE
    sump()
            
if __name__ == '__main__':
    # Set the seed for analyses
    rng = setMasterSeed(rseed)
    
    fnprefix = '%s_%s_%d' % (partition_scheme, method, rseed)
    
    # set up the partition model
    createModelForScheme()
    
    # set up and carry out the MCMC analysis
    runMCMC(fnprefix)
