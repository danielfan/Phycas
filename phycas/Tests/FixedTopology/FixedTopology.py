import sys
from phycas import *

rnseed = 91375
nsites = 2000

if __name__ == '__main__':
    Phycas.PhycassertRaisesException = True

    # Set data_source to None because we will not be reading data from a file
    #simulate.data_source = None  # do not need this since it can be assumed for simulate

    # Define the names of the taxa to use when the simulated data set is saved to a file
    simulate.taxon_labels = [
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
    simulate.starting_tree_source = 'usertree'
    simulate.tree_topology = Newick('(1:0.1,2:0.1,(3:0.1,(4:0.1,((5:0.1,8:0.1):0.1,(6:0.1,((7:0.1,10:0.1):0.1,9:0.1):0.1):0.1):0.1):0.1):0.1)')

    # Create a model
    simulate.default_model = 'hky'
    simulate.kappa = 5.0
    simulate.base_freqs = [0.2, 0.3, 0.3, 0.2]

    # Simulation settings
    simulate.random_seed   = rnseed
    #simulate.random_seed   = 0
    simulate.file_name = 'simulated.nex'
    simulate.nchar     = nsites

    # Simulate data
    simulate()

    # Run MCMC analysis
    mcmc.starting_tree_source    = 'usertree'
    mcmc.data_source             = 'file'
    mcmc.data_file_name          = 'simulated.nex'
    mcmc.fix_topology            = 1 
    mcmc.edge_move_weight        = 17
    mcmc.tree_scaler_weight      = 1
    mcmc.ncycles                 = 100
    mcmc.sample_every            = 1
    mcmc.report_every            = 10
    mcmc.outfile_prefix          = 'fixdtree'
    
	if False:
		import sys,os
		if os.path.basename(sys.executable) == 'python_d.exe':
			raw_input('debug stop')
    
    mcmc()
    
    # Summarize the tree file
    if False:
        sumt.trees              = 'fixdtree.t'
        sumt.out.trees.prefix   = 'fixdtree_trees'
        sumt.out.splits.prefix  = 'fixdtree_splits'
        sumt.burnin             = 0
        sumt.outgroup_taxon     = phycas.sim_taxon_labels[0]
        sumt.tree_credible_prob = 1.0
        sumt()

    # Add a PAUP block to the simulated.nex file to make it easy to check the simulated data
    if False:
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
        tratio = phycas.mcmc_manager.getColdChain().model.calcTRatio()
        f.write('\n  lset nst=2 variant=hky basefreq=(%.2f %.2f %.2f) tratio=%.1f rates=equal;' % (phycas.base_freqs[0], phycas.base_freqs[1], phycas.base_freqs[2], tratio))
        f.write('\n  lscores 1 / userbrlen;')
        f.write('\n  lset nst=2 basefreq=estimate tratio=estimate rates=gamma shape=estimate;')
        f.write('\n  lscores 1 / nouserbrlen;')
        f.write('\n  describe 1 / brlens plot=phylogram;')
        f.write('\n  log stop;')
        f.write('\nend;')
        f.write('\n')
        f.close()
