from phycas import *

if __name__ == '__main__':
    # touch the output files so that they will always be replaced (otherwise
    #   the produced output on a new machine will not have the "replacing <filename>
    #   lines that the reference_output has and the test will fail.
    touch("splits.pdf")
    touch("trees.pdf")
    touch("trees.tre")
    
    phycas.log_file_name    = 'logfile.txt'
    sumt.trees              = 'test.t'
    sumt.out.trees.prefix   = "trees"
    sumt.out.trees.mode     = REPLACE
    sumt.out.splits.prefix  = "splits"
    sumt.out.splits.mode    = REPLACE
    sumt.burnin             = 11
    sumt.outgroup_taxon     = '40 Cyanophora paradoxa'
    sumt.tree_credible_prob = 1.0
    sumt.useGUI             = False
    sumt()
    
