from phycas import *

g_random_seed      = 13579
g_data_file_name   = 'gogarten01.nex'
g_tree_file_name   = 'gogarten.tre'

def runMCMC(fnprefix):
    setMasterSeed(g_random_seed)
    
    # read in the data
    blob = readFile(g_data_file_name)
    
    # set up model
    model.type = 'gain'
    model.fix_scaling_factor = False
    model.scaling_factor_prior = ProbDist.Exponential(0.1)
    
    jpg.data_source = blob.characters    
    jpg.tree_source = TreeCollection(filename=g_tree_file_name)
    jpg.out.log.prefix = fnprefix
    jpg.out.log.mode = REPLACE
    #jpg.scaling_factor_prior = ProbDist.Exponential(0.1)
    jpg.fromtree = 901
    jpg.totree = 1000
    jpg.nreps = 10
    jpg()
    
if __name__ == '__main__':
    rng = setMasterSeed(g_random_seed)
    fnprefix = 'out_gogarten_%d' % (g_random_seed,)    
    runMCMC(fnprefix)
