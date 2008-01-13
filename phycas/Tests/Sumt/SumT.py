from phycas import *

if __name__ == '__main__':
    phycas = Phycas()
    phycas.log_file_name = 'logfile.txt'
    phycas.sumt_tfile_name = 'test.t'
    phycas.sumt_majrule_tree_file = 'contree.tre'
    phycas.sumt_burnin = 101
    #phycas.sumt_outgroup_taxon = '40 Cyanophora paradoxa'
    phycas.sumt()
    