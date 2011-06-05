#!/usr/bin/env python
from phycas import Likelihood
from phycas import Phylogeny
from phycas.ProbDist import Gamma, Lot
from phycas.Utilities.io import TreeCollection
newick = "(1:1.0,2:1.0,(3:1.0,4:1.0):0.5);"

focal_tree = TreeCollection(newick=newick).trees[0]
ntips = focal_tree.getNObservables()
focal_tree.recalcAllSplits(ntips)
topo_ref_dist_calculator = Likelihood.FocalTreeTopoProbCalculatorBase(focal_tree)
split_dist_list = [('-***', Gamma(.1, .10)),
                   ('--*-', Gamma(.2, .10)),
                   ('-*--', Gamma(.3, .10)),
                   ('---*', Gamma(.4, .10)),
                   ('--**', Gamma(.5, .10)),]
rnseed     = 95629
rng = Lot()
rng.setSeed(rnseed)

for split_rep, dist in split_dist_list:
    s = Phylogeny.SplitBase()
    s.createFromPattern(split_rep)
    print split_rep, dist, s.createPatternRepresentation()
    dist.setLot(rng)
    topo_ref_dist_calculator.setEdgeLenDist(s, dist)
def_edge_len_dist = Gamma(.6, .10)
def_edge_len_dist.setLot(rng)
topo_ref_dist_calculator.setDefaultEdgeLenDist(def_edge_len_dist)

sampled_tree = TreeCollection(newick=newick).trees[0]



for i in range(100):
    topo_ref_dist_calculator.sampleTree(sampled_tree, rng)
    print 'tree', i, '= [&U]', sampled_tree.makeNewick(), ';'
    print 'lnprobs = ', topo_ref_dist_calculator.calcTopologyLnProb(sampled_tree, True)
    print
