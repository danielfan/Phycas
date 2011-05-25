import sys
from phycas import *

# call as follows:
# python mth_test.py uni rbcl50.nex

setMasterSeed(9375)
model.type                  = 'jc'
model.edgelen_hyperprior = None
model.edgelen_prior = Exponential(10.0)

fixed_tree = 'FIXEDTREE' in sys.argv[1].upper()
fixed_topo = 'FIXEDTOPO' in sys.argv[1].upper()
do_mcmc = not sys.argv[1].upper().endswith('SUMT')
mcmc.use_unimap = sys.argv[1].upper().startswith('UNI')
fn = sys.argv[2]
file_pref = fn.split('.')[0]
mcmc.data_source = fn


if mcmc.use_unimap:
    prefix = 'uni'
else:
    prefix = 'fels'
if fixed_tree:
    suffix = 'fix_tree'
elif fixed_topo:
    suffix = 'fix_topo'
else:
    suffix = ''
mcmc.out.params = prefix + '_' + file_pref + '_' + suffix + '_' + 'params'
mcmc.out.trees = prefix + '_' + file_pref + '_' + suffix + '_' + 'trees'
mcmc.out.log = prefix + '_' + file_pref + '_' + suffix + '_' + 'log'

mcmc.out.log.mode = REPLACE
mcmc.out.params.mode = REPLACE
mcmc.out.trees.mode = REPLACE
mcmc.ncycles                = 10000
mcmc.sample_every           = 1
mcmc.report_every           = 10
mcmc.mapping_move_weight    = 1

if mcmc.use_unimap:
    if fixed_tree:
        mcmc.starting_tree_source = TreeCollection(newick='(1:0.05,2:0.05,(3:0.05,4:0.05):0.05);')
        mcmc.mapping_move_weight = 1
        mcmc.unimap_nni_move_weight = 0
        mcmc.unimap_sample_ambig_move_weight = 0
    elif fixed_topo:
        raise NotImplementedError()
    else:
        #mcmc.starting_tree_source = TreeCollection(newick='(1:0.12695,2:0.04969,(3:0.08372,(((6:0.07076,((9:0.05874,7:0.04764):0.01801,10:0.06761):0.03612):0.02586,(8:0.12347,5:0.06710):0.02890):0.02613,4:0.07260):0.02658):0.02408);')
        mcmc.mapping_move_weight = 1
        mcmc.unimap_ls_move_weight = 1
        mcmc.unimap_nni_move_weight = 0
        mcmc.unimap_edge_move_weight = 0
        mcmc.unimap_thread_count = 2
        mcmc.unimap_sample_ambig_move_weight = 1
        mcmc.ncycles = 50000
        mcmc.sample_every = 10
        mcmc.report_every = 10000
        mcmc.reference_tree_source = 'hky50taxaML.tre'        
else:
    if fixed_tree:
        raise NotImplementedError()
    elif fixed_topo:
        raise NotImplementedError()
    else:
        mcmc.ls_move_weight = 100
mcmc.tree_scaler_weight     = 0
mcmc.slice_weight           = 1
mcmc.debugging              = True

if False:
    import sys,os
    if os.path.basename(sys.executable) == 'python_d.exe':
        raw_input('debug stop')

if do_mcmc:
    print mcmc.curr
    mcmc()

sys.exit(0)

sumt.out.splits = prefix + '_' + file_pref + '_' + suffix + '_' + 'sumt_splits'
sumt.out.trees = prefix + '_' + file_pref + '_' + suffix + '_' + 'sumt_trees'
sumt.out.log = prefix + '_' + file_pref + '_' + suffix + '_' + 'sumt_log'

sumt.out.log.mode = REPLACE
sumt.out.splits.mode = REPLACE
sumt.out.trees.mode = REPLACE

sumt.trees = prefix + '_' + file_pref + '_' + suffix + '_' + 'trees'

sumt.tree_credible_prob = 1.0

#sumt()



