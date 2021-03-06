import os
os.environ['DYLD_LIBRARY_PATH'] = '/Users/plewis/Documents/Projects/pdev/trunk/phycas/Conversions'

import sys
sys.path.append('/Users/plewis/Documents/Projects/pdev/trunk')

from phycas import *

g_scheme_index     = 3
g_random_seed      = 13579
g_fixed_topology   = True
g_burn_in          = 500
g_num_cycles       = 2000
g_sample_freq      = 1
g_print_freq       = 100
g_data_file_name   = 'marshall_sansrna.nex'

# HKY+G ml tree
# this tree has polytomies user_tree_def    = '(1:0.00155349,((((((2:0.00313198,16:0.00393977):0.02205944,((20:0.00609907,29:0.00668693):0.01295303,28:0.01768460):0.00791425):0.01246952,(26:0.01287999,27:0.01262314):0.02946479):0.01254132,(((((3:0.00811268,15:0.01095734):0.00206885,4:0.01206558):0.00467611,5:0.01820641,((17:0.01336253,18:0.01060919):0.00301371,19:0.01227330):0.00150356):0.00369882,30:0.01212312):0.01038863,22:0.01914910):0.02284554,(7:0.07564186,25:0.07701641,(31:0.13105682,32:0.09451660):0.09029667):0.03628257,(((((8:0.01356633,(10:0.00382838,12:0.00672068):0.00445306):0.02049844,9:0.02080556):0.00181409,14:0.01538416):0.00385776,11:0.02083608):0.00401140,21:0.05064474):0.02132007,13:0.03572793):0.00518702,23:0.02468321):0.01596404,6:0.01404921):0.01540777,24:0.00372350);'
g_user_tree_def    = '(1:0.00156293,(((((((2:0.00325591,16:0.00402707):0.02205143,((20:0.00630496,29:0.00695889):0.01138659,28:0.01818355):0.00809410):0.01335766,(26:0.01188256,27:0.01250127):0.02960036):0.00967768,((((((3:0.00860238,15:0.01056681):0.00220042,4:0.01193905):0.00408679,5:0.01787676):0.00155750,((17:0.01386172,18:0.01057494):0.00315644,19:0.01190134):0.00143905):0.00324977,30:0.01220535):0.01044754,22:0.01935610):0.02073536):0.00355082,(((7:0.07257248,(25:0.07364494,(31:0.12962943,32:0.09693177):0.08359470):0.00544096):0.03392672,(((((8:0.01400035,(10:0.00393309,12:0.00684284):0.00452834):0.02117048,9:0.02142141):0.00187281,14:0.01577341):0.00403985,11:0.02125429):0.00374345,21:0.05267665):0.01590936):0.00295028,13:0.03355061):0.00389731):0.00448117,23:0.02567495):0.01624640,6:0.01444287):0.01588953,24:0.00385810);'

def abort(message):
    print 'Aborting...'
    print 'Reason:',message
    import sys
    sys.exit(0)
    
def leftouts(v, start, stop):
    ref = range(start, stop + 1, 1)
    lo = []
    for i in v:
        if not i in ref:
            lo.append(i)
    return lo

def createModelForScheme(scheme_index):
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
    
    # check for sites left out of subset definitions
    left_out = leftouts(COI + COII + ATPase8 + ATPase6, 1, 2090)
    if len(left_out) > 0:
        abort('These sites are not assigned to any partition by the COI, COII, tRNA, ATPase8, and ATPase6 subsets: %s' % ','.join(['%d' % x for x in left_out]))
    left_out = leftouts(firsts + seconds + thirds, 1, 2090)
    if len(left_out) > 0:
        abort('These sites are not assigned to any partition by the firsts, seconds, thirds and tRNA subsets: %s' % ','.join(['%d' % x for x in left_out]))

    # MLEs based on unpartitioned GTR+G analysis of g_user_tree_def
    # -ln L      10191.56
    # Base frequencies:
    #   A        0.376859
    #   C        0.096002
    #   G        0.094516
    #   T        0.432623
    # Rate matrix R:
    #   AC        4.28396
    #   AG       34.10179
    #   AT        1.39027
    #   CG        2.53545
    #   CT       40.60048
    #   GT        1.00000
    # Shape      0.146955
    relrate_MLEs    = [4.28396, 34.10179,  1.39027,  2.53545, 40.60048,  1.00000]
    state_freq_MLEs = [0.376859, 0.096002, 0.094516, 0.432623]
    shape_MLE       = 0.146955
    pinvar_MLE      = 0.621227
        
    # model type (jc, hky or gtr)
    model.type = 'gtr'
    
    # GTR relative rates
    model.update_relrates_separately = False
    model.relrates = [x/sum(relrate_MLEs) for x in relrate_MLEs]
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

    if scheme_index == 0:
        # unpartitioned
        partition()
    elif scheme_index == 1:
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
    elif scheme_index == 2:
        # partition by codon
        m1 = model()
        m2 = model()
        m3 = model()
        partition.addSubset(firsts,  m1, 'First codon positions')
        partition.addSubset(seconds, m2, 'Second codon positions')
        partition.addSubset(thirds,  m3, 'Third codon positions')
        partition()
    elif scheme_index == 3:
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
        abort('Sorry, unrecognized scheme (%s)' % str(scheme))
        
def runMCMC(fnprefix):
    # read in the data
    blob = readFile(g_data_file_name)
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
    mcmc.starting_tree_source = TreeCollection(newick=g_user_tree_def)
    mcmc.debugging = False
    mcmc.fix_topology = g_fixed_topology
    mcmc.adapt_first = 2
    mcmc.sample_every = g_sample_freq
    mcmc.report_every = g_print_freq
    
    mcmc.burnin = g_burn_in
    mcmc.uf_num_edges = 100

    mcmc.ncycles = g_num_cycles
    cpo()
    
    sump.cpofile = 'sitelike.txt'
    sump()
                
def run():
    rng = setMasterSeed(g_random_seed)
    
    fnprefix = '%d_%d' % (g_scheme_index, g_random_seed)
    
    # set up the partition model
    createModelForScheme(g_scheme_index)
    
    # set up the MCMC analysis
    runMCMC(fnprefix)

if __name__ == '__main__':
    run()
