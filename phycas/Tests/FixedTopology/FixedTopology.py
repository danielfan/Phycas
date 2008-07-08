import sys
from phycas import *

rnseed = 91375
nsites = 2000

if __name__ == '__main__':
    phycas = Phycas()

    # Set data_source to None because we will not be reading data from a file
    phycas.data_source = None

    # Define the names of the taxa to use when the simulated data set is saved to a file
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
    phycas.tree_topology = Newick('(1:0.1,2:0.1,(3:0.1,(4:0.1,((5:0.1,8:0.1):0.1,(6:0.1,((7:0.1,10:0.1):0.1,9:0.1):0.1):0.1):0.1):0.1):0.1)')

    # Create a model
    phycas.default_model = 'hky'
    phycas.starting_kappa = 5.0
    phycas.starting_freqs = [0.2, 0.3, 0.3, 0.2]

    # Simulation settings
    phycas.random_seed   = rnseed
    #phycas.random_seed   = 0
    phycas.sim_file_name = 'simulated.nex'
    phycas.sim_nchar     = nsites

    # Simulate data
    phycas.simulateDNA()

    # Run MCMC analysis
    phycas.starting_tree_source    = 'usertree'
    phycas.data_source             = 'file'
    phycas.data_file_name          = 'simulated.nex'
    phycas.fix_topology            = 1 
    phycas.edge_move_weight        = 17
    phycas.tree_scaler_weight      = 1
    phycas.ncycles                 = 100
    phycas.sample_every            = 1
    phycas.report_every            = 10
    phycas.outfile_prefix          = 'fixdtree'
    
    import sys,os
    if os.path.basename(sys.executable) == 'python_d.exe':
        raw_input('debug stop')
    
    phycas.mcmc()
    
    # Summarize the tree file
    if False:
        sumt.trees    = 'fixdtree.t'
        sumt.trees_prefix       = 'fixdtree_trees'
        sumt.splits_prefix      = 'fixdtree_splits'
        sumt.burnin             = 0
        sumt.outgroup_taxon     = phycas.sim_taxon_labels[0]
        sumt.output_replace     = True
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
        f.write('\n  lset nst=2 variant=hky basefreq=(%.2f %.2f %.2f) tratio=%.1f rates=equal;' % (phycas.starting_freqs[0], phycas.starting_freqs[1], phycas.starting_freqs[2], tratio))
        f.write('\n  lscores 1 / userbrlen;')
        f.write('\n  lset nst=2 basefreq=estimate tratio=estimate rates=gamma shape=estimate;')
        f.write('\n  lscores 1 / nouserbrlen;')
        f.write('\n  describe 1 / brlens plot=phylogram;')
        f.write('\n  log stop;')
        f.write('\nend;')
        f.write('\n')
        f.close()
