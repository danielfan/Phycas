from phycas import *
sumt.trees = "HibbetGTRrun0.nex.t"
sumt.burnin = 1000
sumt.tree_credible_prob = 0.01
sumt.out.trees = "hibbet_sumt_trees"
sumt.out.splits.prefix = "hibbet_sumt_splits"
sumt()
