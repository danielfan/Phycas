import sys
from phycas import *

model = 'hky'
rnseed = 13975
debug = 0
nsites = 5000
ntax = 4

if __name__ == '__main__':
    if len(sys.argv) > 1:
        debug = 1
        
    phycas = Phycas()

    # Set data_source to None because we will not be reading data from a file
    phycas.data_source = None

    # Define the names of the taxa to use when the simulated data set is saved to a file
    if ntax == 4:
        phycas.sim_taxon_labels = ['A', 'B', 'C', 'D']
    else:
        phycas.sim_taxon_labels = [
        	'Chara',
        	'Conocephalum',
        	'Bazzania',
        	'Sphagnum',
        	'Osmunda',
        	'Picea',
        	'Iris',
        	'Asplenium',
        	'Nicotiana',
        	'Avena'
        ]

    # Create a model tree
    phycas.starting_tree_source = 'usertree'
    if ntax == 4:
        phycas.tree_topology = '(1:0.1,2:0.1,(3:0.1,4:0.1):0.1)'
    else:
        phycas.tree_topology = '(1:0.1,2:0.1,(3:0.1,(4:0.1,((5:0.1,8:0.1):0.1,(6:0.1,((7:0.1,10:0.1):0.1,9:0.1):0.1):0.1):0.1):0.1):0.1)'

    # Create a model
    assert model == 'jc' or model == 'hky', "model must be either 'jc' or 'hky'"
    if model == 'jc':
        phycas.default_model = 'jc'
    elif model == 'hky':
        phycas.default_model = 'hky'
        phycas.starting_kappa = 5.0
        phycas.starting_freqs = [0.25, 0.25, 0.25, 0.25]

    # Simulation settings
    phycas.random_seed   = rnseed
    #phycas.random_seed   = 0
    phycas.sim_file_name = 'simulated.nex'
    phycas.sim_nchar     = nsites

    if debug > 0:
        raw_input('debug stop')

    # Simulate data
    phycas.simulateDNA()

    # Now compute the likelihood of the model tree
    phycas.starting_tree_source    = 'random'
    phycas.data_source             = 'file'
    phycas.data_file_name          = 'simulated.nex'
    phycas.ls_move_weight          = 1 
    phycas.tree_scaler_weight      = 1
    phycas.nielsen_move_weight     = 1
    phycas.ncycles                 = 100
    phycas.sample_every            = 1
    phycas.report_every            = 10
    phycas.use_unimap              = True
    phycas.outfile_prefix          = 'doofus'
    phycas.mcmc()
    
    phycas.log_file_name           = 'logfile.txt'
    phycas.sumt_input_tree_file    = 'doofus.t'
    phycas.sumt_trees_prefix       = 'doofus_trees'
    phycas.sumt_splits_prefix      = 'doofus_splits'
    phycas.sumt_burnin             = 0
    phycas.sumt_outgroup_taxon     = phycas.sim_taxon_labels[0]
    phycas.sumt_output_replace     = True
    phycas.sumt_tree_credible_prob = 1.0
    phycas.sumt()

    
    # Add a PAUP block to the simulated.nex file to make it easy to check the results
    f = file('simulated.nex', 'a')
    f.write('\n')
    f.write('\nbegin paup;')
    f.write('\n  set criterion=likelihood storebrlen;')
    f.write('\nend;')
    f.write('\n')
    f.write('\nbegin trees;')
    f.write('\n  translate')
    for i,nm in enumerate(phycas.sim_taxon_labels):
        if nm.count(' ') > 0:
            f.write("\n    %d '%s'" % (i+1, nm))
        else:
            f.write("\n    %d %s" % (i+1, nm))
        if i < len(phycas.sim_taxon_labels) - 1:
            f.write(',')
        else:
            f.write(';')
    f.write('\n  utree simtree = %s;' % phycas.tree_topology)
    f.write('\nend;')
    f.write('\n')
    f.write('\nbegin paup;')
    f.write('\n  log start file=paup.log replace;')
    f.write('\n  set criterion=likelihood autoclose;')
    if model == 'hky':
        tratio = phycas.mcmc_manager.getColdChain().model.calcTRatio()
        f.write('\n  lset nst=2 variant=hky basefreq=(%.2f %.2f %.2f) tratio=%.1f rates=equal;' % (phycas.starting_freqs[0], phycas.starting_freqs[1], phycas.starting_freqs[2], tratio))
    else:
        f.write('\n  lset nst=1 basefreq=equal rates=equal;')
    f.write('\n  lscores 1 / userbrlen;')
    f.write('\n  lset nst=2 basefreq=estimate tratio=estimate rates=gamma shape=estimate;')
    f.write('\n  lscores 1 / nouserbrlen;')
    f.write('\n  describe 1 / brlens plot=phylogram;')
    f.write('\n  log stop;')
    f.write('\nend;')
    f.write('\n')
    f.close()
    