from phycas import *

if __name__ == '__main__':
    phycas = Phycas()
    phycas.log_file_name = 'logfile.txt'
    phycas.sumt_input_tree_file = 'test.t'
    phycas.sumt_trees_prefix = 'trees'
    phycas.sumt_splits_prefix = 'splits'
    phycas.sumt_burnin = 11
    phycas.sumt_outgroup_taxon = '40 Cyanophora paradoxa'
    phycas.sumt_output_replace = True
    phycas.sumt()
    