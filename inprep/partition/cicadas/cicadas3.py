from phycas import *

pid_offset = -1	# this amount added to pid 
nschemes = 4
nmethods = 2
nseeds   = 2

dry_run          = False
fixed_topology   = True
num_beta_values  = 3
burn_in          = 10000
cycles_per_beta  = 50000
sample_freq      = 100
print_freq       = 1000
data_file_name   = 'marshall.nex'

# HKY+G ml tree
user_tree_def    = '(1:0.00155349,((((((2:0.00313198,16:0.00393977):0.02205944,((20:0.00609907,29:0.00668693):0.01295303,28:0.01768460):0.00791425):0.01246952,(26:0.01287999,27:0.01262314):0.02946479):0.01254132,(((((3:0.00811268,15:0.01095734):0.00206885,4:0.01206558):0.00467611,5:0.01820641,((17:0.01336253,18:0.01060919):0.00301371,19:0.01227330):0.00150356):0.00369882,30:0.01212312):0.01038863,22:0.01914910):0.02284554,(7:0.07564186,25:0.07701641,(31:0.13105682,32:0.09451660):0.09029667):0.03628257,(((((8:0.01356633,(10:0.00382838,12:0.00672068):0.00445306):0.02049844,9:0.02080556):0.00181409,14:0.01538416):0.00385776,11:0.02083608):0.00401140,21:0.05064474):0.02132007,13:0.03572793):0.00518702,23:0.02468321):0.01596404,6:0.01404921):0.01540777,24:0.00372350);'

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
    tRNA          = subset(1477, 1538)
    ATPase8       = subset(1539, 1689)
    ATPase8first  = subset(1539, 1689, 3)
    ATPase8second = subset(1540, 1689, 3)
    ATPase8third  = subset(1541, 1689, 3)
    ATPase6       = subset(1690, 2152)
    ATPase6first  = subset(1691, 2152, 3)
    ATPase6second = subset(1692, 2152, 3)
    ATPase6third  = subset(1690, 2152, 3)
    firsts        = COIfirst  + COIIfirst  + ATPase8first  + ATPase6first
    seconds       = COIsecond + COIIsecond + ATPase8second + ATPase6second
    thirds        = COIthird  + COIIthird  + ATPase8third  + ATPase6third
    
    # check for sites left out of subset definitions
    left_out = leftouts(COI + COII + tRNA + ATPase8 + ATPase6, 1, 2152)
    if len(left_out) > 0:
        abort('These sites are not assigned to any partition by the COI, COII, tRNA, ATPase8, and ATPase6 subsets: %s' % ','.join(['%d' % x for x in left_out]))
    left_out = leftouts(firsts + seconds + thirds + tRNA, 1, 2152)
    if len(left_out) > 0:
        abort('These sites are not assigned to any partition by the firsts, seconds, thirds and tRNA subsets: %s' % ','.join(['%d' % x for x in left_out]))

    # MLEs based on unpartitioned GTR+G analysis of user_tree_def (max. lnL = -10407.54)
    relrate_MLEs    = [4.65042/88.22296, 35.65675/88.22296, 1.43228/88.22296, 2.61034/88.22296, 42.87317/88.22296, 1.00000/88.22296]
    state_freq_MLEs = [0.378622, 0.094271, 0.094319, 0.432788]
    shape_only_MLE  = 0.146081
    shape_pinv_MLE  = 1.862738
    pinvar_MLE      = 0.621227
        
    # model type (jc, hky or gtr)
    model.type = 'gtr'
    
    # GTR relative rates
    model.update_relrates_separately = False
    model.relrates = relrate_MLEs
    model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
    
    # state frequencies
    model.update_freqs_separately = False
    model.state_freqs = state_freq_MLEs
    model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))
    
    # discrete gamma rate heterogeneity
    model.num_rates = 4
    model.gamma_shape = shape_only_MLE
    model.gamma_shape_prior = Exponential(1.0)
    
    # proportion of invariable sites
    model.pinvar_model = False
    model.pinvar = pinvar_MLE
    model.pinvar_prior = Beta(1.0, 1.0)
    
    # edge lengths
    model.edgelen_hyperprior = None    
    model.edgelen_prior = Exponential(40.0)

    if scheme_index == 0:
        # unpartitioned
        partition()
    elif scheme_index == 1:
        # partition by codon (plus tRNA)
        m1 = model()
        m2 = model()
        m3 = model()
        m4 = model()
        partition.addSubset(firsts,  m1, 'First codon positions')
        partition.addSubset(seconds, m2, 'Second codon positions')
        partition.addSubset(thirds,  m3, 'Third codon positions')
        partition.addSubset(tRNA,    m4, 'tRNA')
        partition()
    elif scheme_index == 2:
        # partition by gene
        m1 = model()
        m2 = model()
        m3 = model()
        m4 = model()
        m5 = model()
        partition.addSubset(COI,     m1, 'COI')
        partition.addSubset(COII,    m2, 'COII')
        partition.addSubset(tRNA,    m3, 'rRNA')
        partition.addSubset(ATPase8, m4, 'ATPase8')
        partition.addSubset(ATPase6, m5, 'ATPase6')
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
        m13 = model()
        partition.addSubset(COIfirst,       m1, 'COIfirst')
        partition.addSubset(COIsecond,      m2, 'COIsecond')
        partition.addSubset(COIthird,       m3, 'COIthird')
        partition.addSubset(COIIfirst,      m4, 'COIIfirst')
        partition.addSubset(COIIsecond,     m5, 'COIIsecond')
        partition.addSubset(COIIthird,      m6, 'COIIthird')
        partition.addSubset(tRNA,           m7, 'rRNA')
        partition.addSubset(ATPase8first,   m8, 'ATPase8first')
        partition.addSubset(ATPase8second,  m9, 'ATPase8second')
        partition.addSubset(ATPase8third,  m10, 'ATPase8third')
        partition.addSubset(ATPase6first,  m11, 'ATPase6first')
        partition.addSubset(ATPase6second, m12, 'ATPase6second')
        partition.addSubset(ATPase6third,  m13, 'ATPase6third')
        partition()
    else:
        abort('Sorry, unrecognized scheme (%s)' % str(scheme))
        
def runMCMC(fnprefix, method_index):
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
    mcmc.fix_topology = fixed_topology
    mcmc.adapt_first = 2
    mcmc.sample_every = sample_freq
    mcmc.report_every = print_freq
    
    mcmc.burnin = burn_in

    if method_index == 0:
        # harmonic mean method
        mcmc.ncycles = cycles_per_beta*num_beta_values
        if not dry_run:
            mcmc()
    else:
        # steppingstone sampling method
        mcmc.ncycles = cycles_per_beta
        ss.nbetavals = num_beta_values
        ss.shape1 = 1.0
        ss.shape2 = 1.0
        ss.scubed = True
        if not dry_run:
            ss()
    
    if not dry_run:
        sump.file         = fnprefix + '.p'
        sump.out.log      = fnprefix + '.sump.txt'
        sump.out.log.mode = REPLACE
        sump()
        
    if dry_run:
        print 'dry run finished for', fnprefix
    
def run(cid, pid):
    # Want to spread runs over the following combinations:
    #   nschemes = 5 different partition schemes 
    #   nmethods = 2 different marginal likelihood estimation methods (hm vs. sss)
    #   nseeds   = 2 different runs per combination
    #   ms = nmethods*nseeds
    #
    #   pid    w = (pid % ms)      (pid // ms)     (w // nseeds)     (w % nseeds)  
    #   -----------------------------------------------------------------------------
    #     0          0                  0                0                0       
    #     1          1                  0                0                1       
    #   -----------------------------------------------------------------------------
    #     2          2                  0                1                0       
    #     3          3                  0                1                1       
    #   -----------------------------------------------------------------------------
    #     4          0                  1                0                0       
    #     5          1                  1                0                1       
    #   -----------------------------------------------------------------------------
    #     6          2                  1                1                0       
    #     7          3                  1                1                1       
    #   -----------------------------------------------------------------------------
    #     8          0                  2                0                0       
    #     9          1                  2                0                1       
    #   -----------------------------------------------------------------------------
    #    10          2                  2                1                0       
    #    11          3                  2                1                1       
    #   -----------------------------------------------------------------------------
    #    12          0                  3                0                0       
    #    13          1                  3                0                1       
    #   -----------------------------------------------------------------------------
    #    14          2                  3                1                0       
    #    15          3                  3                1                1       
    #   -----------------------------------------------------------------------------
    #    16          0                  4                0                0       
    #    17          1                  4                0                1       
    #   -----------------------------------------------------------------------------
    #    18          2                  4                1                0       
    #    19          3                  4                1                1       
    #   -----------------------------------------------------------------------------
    #                              scheme_index     method_index      seed_index
    
    pid += pid_offset
    
    # determine the file name prefix for output files
    ms = nmethods*nseeds
    w = pid % ms
    scheme_index = pid//ms
    method_index = w//nseeds
    seed_index   = w % nseeds

    # Set the seed for analyses
    rng = setMasterSeed(cid)
    for i in range(100*pid):
        rng.uniform()
    anal_seed = rng.getSeed()
    rng = setMasterSeed(anal_seed)
    
    fnprefix = '%d_%d_%s_%d' % (pid, scheme_index, (method_index == 0 and 'hm' or 'ss'), anal_seed)
    
    # set up the partition model
    createModelForScheme(scheme_index)
    
    # set up the MCMC analysis
    runMCMC(fnprefix, method_index)

if __name__ == '__main__':
    # A typical condor script might look like this:
    #
    # Executable  = /usr/bin/python
    # Requirements = ParallelSchedulingGroup == "stats group"
    # Universe   = vanilla
    # output     = out$(Cluster).$(Process).out
    # error      = out$(Cluster).$(Process).err
    # Log        = out$(Cluster).$(Process).log
    # Arguments  = "hummingbirds.py $(Cluster) $(Process)"
    # transfer_executable = false
    # transfer_input_files = hummingbirds.py,M3354.nex
    # should_transfer_files = YES
    # when_to_transfer_output = ON_EXIT
    # # notify_user = paul.lewis@uconn.edu
    # environment = "LD_LIBRARY_PATH=/scratch/phycas/Conversions PYTHONPATH=/scratch"
    # Queue 36
    #
    # To submit, save the above in a file named gohum.condor and type
    # condor_submit gohum.condor
    #
    # To kill cluster number 4971
    # condor_rm 4971
    #
    # To check status, use
    # condor_status -submitters
    #
    # condor_q -submitter plewis
    #
    # For help, go to http://www.stat.uconn.edu/www/resources/statcluster/overview.html

    # Get the cluster id and process id, both of which should have been provided on the command line
    assert len(sys.argv) > 2
    cluster_id = int(sys.argv[1])
    print 'cluster_id =',cluster_id
    process_id = int(sys.argv[2])
    print 'process_id =',process_id
    run(cluster_id, process_id)
