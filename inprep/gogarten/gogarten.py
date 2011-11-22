from phycas import *

g_random_seed      = 13579
g_data_file_name   = 'gogarten01.nex'
g_tree_file_name   = 'gogarten.tre'
g_single_tree_file_name   = 'gogarten1000rooted.tre'

def runBinary(fnprefix):
    # set up model
    model.type = 'binary'
    model.fix_scaling_factor = False
    #model.scaling_factor_prior = ProbDist.Exponential(0.1)
    model.scaling_factor_prior = ProbDist.Lognormal(0.0,1.0)
    model.kappa_prior = ProbDist.Lognormal(0.0,1.0)
    #model.kappa_prior = ProbDist.BetaPrime(1.0,1.0)
    #model.kappa_prior = ProbDist.Exponential(0.1)
    
    jpg.out.details.prefix = 'binary.'+fnprefix
    jpg.out.details.mode = REPLACE
    jpg()
    
def runGain(fnprefix):
    # set up model
    model.type = 'gain'
    model.fix_scaling_factor = False
    #model.scaling_factor_prior = ProbDist.Exponential(0.1)
    model.scaling_factor_prior = ProbDist.Lognormal(0.0,1.0)

    jpg.out.details.prefix = 'gain.'+fnprefix
    jpg.out.details.mode = REPLACE
    jpg()

def runJPG(fnprefix):
    jpg.data_source = blob.characters    
    jpg.tree_source = TreeCollection(filename=g_tree_file_name)
    jpg.fromtree = 901
    jpg.totree = 1000
    jpg.nreps = 100
    
    runGain(fnprefix)
    runBinary(fnprefix)
    
def simJPG():
    sim.tree_source = TreeCollection(filename=g_single_tree_file_name)
    sim.taxon_labels = ['taxon-%d' % (i+1) for i in range(45)]
    sim.nchar = 100
    sim.edgelen_dist = Exponential(10.0)
    sim.file_name = 'doofus.nex'
    sim()
    
def runSS(fnprefix):
    # set up model
    model.type = 'binary'
    model.fix_scaling_factor = False
    #model.scaling_factor_prior = ProbDist.Exponential(0.1)
    model.scaling_factor_prior = ProbDist.Lognormal(0.0,1.0)
    model.kappa_prior = ProbDist.Lognormal(0.0,1.0)
    #model.kappa_prior = ProbDist.BetaPrime(1.0,1.0)
    #model.kappa_prior = ProbDist.Exponential(0.1)
    
    mcmc.data_source            = blob.characters    
    mcmc.starting_tree_source   = TreeCollection(filename=g_single_tree_file_name)
    mcmc.ncycles                = 1000
    mcmc.fix_topology           = True
    
    ss.ti           = False
    ss.maxbeta      = 1.0
    ss.minbeta      = 0.0
    ss.minsample    = 10
    ss.shape1       = 1.0
    ss.shape2       = 1.0
    ss.xcycles      = 0
    ss.nbetavals    = 21
    ss()

def runMCMC(fnprefix):
    global blob
    
    setMasterSeed(g_random_seed)
    
    # read in the data
    blob = readFile(g_data_file_name)
    
    #runJPG(fnprefix)
    simJPG()
    #runSS(fnprefix)
    
    
if __name__ == '__main__':
    rng = setMasterSeed(g_random_seed)
    fnprefix = 'out_gogarten_%d' % (g_random_seed,)    
    runMCMC(fnprefix)
