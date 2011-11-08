from phycas import *

g_random_seed      = 13579
g_data_file_name   = 'test01.nex'
g_tree_file_name   = 'test.tre'

# sf    avg. lnL
# 0.5     -3.07528004645
# 1.0     -2.46956902616
# 1.5     -2.15795241251
# 2.0     -1.96939280831
# 2.5     -1.84974940417
# 3.0     -1.7745875202
# 3.5     -1.73067197465
# 4.0     -1.7099700404
# 4.5  *  -1.70721002709
# 5.0     -1.71872906894
# 5.5     -1.74186917796
# 6.0     -1.77463469653
# 6.5     -1.81548565567
# 7.0     -1.86320677966
# 7.5     -1.91682098632
# 8.0     -1.97553028342
# 8.5     -2.0386741941
# 9.0     -2.10569977232
# 9.5     -2.17613950173
# 10.0    -2.24959469116

def abort(message):
    print 'Aborting...'
    print 'Reason:',message
    import sys
    sys.exit(0)
    
def runMCMC(fnprefix):
    setMasterSeed(g_random_seed)
    
    # read in the data
    blob = readFile(g_data_file_name)
    
    # set up model
    model.type = 'loss'
    model.scaling_factor = 10.0
    model.fix_scaling_factor = True
    
    jpg.data_source = blob.characters    
    jpg.tree_source = TreeCollection(filename=g_tree_file_name)
    jpg.fromtree = 1
    jpg.totree = 1
    jpg()
    
if __name__ == '__main__':
    rng = setMasterSeed(g_random_seed)
    fnprefix = 'out_gogarten_%d' % (g_random_seed,)    
    runMCMC(fnprefix)
