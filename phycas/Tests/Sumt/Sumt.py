from phycas import *

if __name__ == '__main__':
    phycas.log_file_name = 'logfile.txt'
    sumt.trees              = 'test.t'
    sumt.out.tree.prefix    = "trees"
    sumt.out.trees.mode     = REPLACE
    sumt.out.splits.prefix  = "splits"
    sumt.out.splits.mode    = REPLACE
    sumt.burnin             = 11
    sumt.outgroup_taxon     = '40 Cyanophora paradoxa'
    sumt.tree_credible_prob = 1.0
    sumt()
    
