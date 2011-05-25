import sys
from phycas import *

model.type = 'hky'
rnseed = 13975
nsites = 2000
ntax = 10

if __name__ == '__main__':
    # Define the names of the taxa to use when the simulated data set is saved to a file
    if ntax == 4:
        sim.taxon_labels = ['A', 'B', 'C', 'D']
    else:
        sim.taxon_labels = [
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
    sim.tree_source = 'usertree'
    if ntax == 4:
        sim.tree_topology = Newick('(1:0.1,2:0.1,(3:0.1,4:0.1):0.1)')
    else:
        sim.tree_topology = Newick('(1:0.1,2:0.1,(3:0.1,(4:0.1,((5:0.1,8:0.1):0.1,(6:0.1,((7:0.1,10:0.1):0.1,9:0.1):0.1):0.1):0.1):0.1):0.1)')

    # Create a model
    assert model.type == 'jc' or model.type == 'hky', "model must be either 'jc' or 'hky'"
    if model.type == 'hky':
        model.kappa = 5.0
        model.state_freqs = [0.25, 0.25, 0.25, 0.25]

    # Simulation settings
    sim.random_seed = rnseed
    sim.file_name   = 'simulated.nex'
    sim.nchar       = nsites

    # Simulate data
    sim()

    # Run a unimap mcmc analysis
    mcmc.starting_tree_source    = 'random'
    mcmc.data_source             = 'simulated.nex'
    mcmc.tree_scaler_weight      = 1
    mcmc.fix_topology            = True
    mcmc.ncycles                 = 100
    mcmc.sample_every            = 1
    mcmc.report_every            = 10
    mcmc.use_unimap              = True
    mcmc.outfile_prefix          = 'doofus'
    mcmc.debugging               = True
    mcmc()
    
    # Summarize the tree file
    #sumt.trees              = 'doofus.t'
    #sumt.burnin             = 0
    #sumt.outgroup_taxon     = phycas.sim_taxon_labels[0]
    #sumt.tree_credible_prob = 1.0
    #sumt.out.tree.prefix    = "doofus_trees"
    #sumt.out.trees.mode     = REPLACE
    #sumt.out.splits.prefix  = "doofus_splits"
    #sumt.out.splits.mode    = REPLACE
    #sumt()

    # Add a PAUP block to the simulated.nex file to make it easy to check the simulated data
    f = file('simulated.nex', 'a')
    f.write('\n')
    f.write('\nbegin paup;')
    f.write('\n  set criterion=likelihood storebrlen;')
    f.write('\nend;')
    f.write('\n')
    f.write('\nbegin trees;')
    f.write('\n  translate')
    for i,nm in enumerate(sim.taxon_labels):
        if nm.count(' ') > 0:
            f.write("\n    %d '%s'" % (i+1, nm))
        else:
            f.write("\n    %d %s" % (i+1, nm))
        if i < len(sim.taxon_labels) - 1:
            f.write(',')
        else:
            f.write(';')
    f.write('\n  utree simtree = %s;' % sim.tree_topology)
    f.write('\nend;')
    f.write('\n')
    f.write('\nbegin paup;')
    f.write('\n  log start file=paup.log replace;')
    f.write('\n  set criterion=likelihood autoclose;')
    if model == 'hky':
        tratio = mcmc.mcmc_manager.getColdChain().model.calcTRatio()
        f.write('\n  lset nst=2 variant=hky basefreq=(%.2f %.2f %.2f) tratio=%.1f rates=equal;' % (mcmc.state_freqs[0], mcmc.state_freqs[1], mcmc.state_freqs[2], tratio))
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
    
