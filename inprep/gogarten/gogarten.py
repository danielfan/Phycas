from phycas import *

g_random_seed      = 13579
g_data_file_name   = 'gogarten01.nex'
g_tree_file_name   = 'gogarten.tre'

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

def runMCMC(fnprefix):
    setMasterSeed(g_random_seed)
    
    # read in the data
    blob = readFile(g_data_file_name)
    
    jpg.data_source = blob.characters    
    jpg.tree_source = TreeCollection(filename=g_tree_file_name)
    jpg.fromtree = 901
    jpg.totree = 1000
    jpg.nreps = 100
    
    runGain(fnprefix)
    runBinary(fnprefix)
    
if __name__ == '__main__':
    rng = setMasterSeed(g_random_seed)
    fnprefix = 'out_gogarten_%d' % (g_random_seed,)    
    runMCMC(fnprefix)
