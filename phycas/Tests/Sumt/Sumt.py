from phycas import *

if __name__ == '__main__':
    phycas.log_file_name = 'logfile.txt'
    sumt.trees = 'test.t'
    sumt.trees_prefix = 'trees'
    sumt.splits_prefix = 'splits'
    sumt.burnin = 11
    sumt.outgroup_taxon = '40 Cyanophora paradoxa'
    sumt.output_replace = True
    sumt.tree_credible_prob = 1.0
    sumt()
    
