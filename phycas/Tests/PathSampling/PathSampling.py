import sys
from phycas import *

def analyze_results(phycas, rep, betavect, likevect):
    # open a file for writing (file name has replicate in it so each rep will not overwrite previous files)
    f = open('wangang-%d.txt' % (rep + 1), 'w')
    
    # line below produces 2-tuples where the first element is an index and the second is one of the values in betavect:
    # e.g. (0,0.9), (1,0.8), (2,0.7), (3,0.6), (4,0.5), (5,0.4), (6,0.3), (7,0.2), (8,0.1), (9,0.0)
    for i,b in enumerate(betavect):
        # line below sets like to each sampled likelihood in likevect[i]
        for like in likevect[i]:
            f.write('%d\t%.5f\t%.8f\n' % (i,b,like))
            
    # this line closes the file
    f.close()

prior_mean = 0.1

phycas = Phycas()

# general settings 
phycas.random_seed          = 13579
phycas.sim_taxon_labels     = ['A','B','C','D','E']

# settings relating to the tree used
phycas.starting_tree_source = 'usertree' # use 'random' for a random starting tree, 'usertree' for a specified topology
phycas.tree_topology        = Newick('(1:0.01,2:0.01,(3:0.01,(4:0.01,5:0.01):0.01):0.01)')

# settings relating to substitution model
phycas.default_model        = 'gtr'    # use 'jc' for JC,'hky' for HKY model, 'gtr' for GTR model
phycas.num_rates            = 4        # use number > 1 for discrete gamma rate heterogeneity
phycas.estimate_pinvar      = True     # use True for proportion invariable sites model
phycas.starting_pinvar      = 0.2
phycas.starting_shape       = 0.5
phycas.starting_relrates    = [1.0, 4.0, 1.0, 1.0, 4.0, 1.0]    # AC, AG, AT, CG, CT, GT
phycas.starting_freqs       = [0.2, 0.3, 0.3, 0.2]  # order is A,C,G,T

phycas.sim_file_name        = 'simulated.nex'
phycas.sim_nchar            = 1000

nreps = 1
for rep in range(nreps):
    # Simulate data
    phycas.data_source = None
    phycas.simulateDNA()

    phycas.data_source          = 'file'
    phycas.data_file_name       = 'simulated.nex'
    phycas.outfile_prefix       = 'junk'  # this will create junk.p and junk.t but these are not used, hence the prefix 'junk'
    phycas.is_standard_heating  = False   # True means heat posterior, False means heat only likelihood    

    # these settings specific to path sampling
    # number of cycles will be ps_burnin + (ps_Q*ps_nbetavals)
    phycas.ps_toward_posterior  = False # specify True for beta ps_minbeta -> ps_maxbeta, False means beta ps_maxbeta -> ps_minbeta
    phycas.ps_burnin            = 100   # number of cycles before any sampling is done
    phycas.ps_Q                 = 45    # number of MCMC cycles per value of beta
    phycas.ps_sample_every      = 5     # log-likelihood will be sampled whenever cycle modulo ps_sample_every equals 0
    phycas.ps_nbetavals         = 10    # number of beta values 
    phycas.ps_filename          = None  # set to file name to create file containing beta values and average likelihoods
    phycas.ps_minbeta           = 0.0   # smallest beta value to be sampled
    phycas.ps_maxbeta           = 0.9   # largest beta value to be sampled

    # mcmc settings
    # set sample_every so large that no samples will ever be saved in junk.t or junk.p
    phycas.sample_every         = phycas.ps_burnin + (phycas.ps_Q*phycas.ps_nbetavals) + 1
    phycas.report_every         = 100
    phycas.adapt_first          = 2     
    phycas.verbose              = True

    # settings relating to prior distributions
    phycas.edgelen_dist          = ProbDist.ExponentialDist(1.0/prior_mean)
    phycas.relrate_prior         = ProbDist.ExponentialDist(1.0) 
    phycas.gamma_shape_prior     = ProbDist.ExponentialDist(1.0)
    phycas.pinvar_prior          = ProbDist.BetaDist(1.0, 1.0)
    phycas.base_freq_param_prior = ProbDist.ExponentialDist(1.0)

    # try rescaling the entire tree at least once per cycle
    phycas.tree_scaler_weight    = 1
    
    # fix the tree topology, which means use EdgeMove (which only change edge lengths) instead
    # of LargetSimonMove (which changes both edge lengths and the tree topology)
    phycas.fix_topology          = True
    phycas.edge_move_weight      = 100      # do 100 EdgeMove updates per cycle (makes sense to set this to some multiple of the number of edges)

    # probably best to *not* change any of these
    phycas.nchains              = 1         # multiple chain analyses are not yet working
    phycas.slice_weight         = 1         # means one slice sampling update of each model parameter per cycle
    phycas.using_hyperprior     = False     # False means do not use hierarchical model for edge lengths
    phycas.use_flex_model       = False     # flex model is not yet published

    # this is the command that begins the analysis
    phycas.pathsampling()
    
    analyze_results(phycas, rep, phycas.wangang_sampled_betas, phycas.wangang_sampled_likes)
