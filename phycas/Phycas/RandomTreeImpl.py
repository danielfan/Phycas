import copy, os, sys, math, random, copy
from phycas.Utilities.PhycasCommand import _value_for_user
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas.Utilities.io import TreeCollection
from phycas import Phylogeny

class TreeSimulator(CommonFunctions, TreeCollection):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    A collection of simulated trees.  
    
    The n_trees attribute controls the size of the collection.  If it is set to
    0 then the collection will be bottomless and trees will not be stored (thus
    indexing the collection is not supported).
    
    If n_trees is > 0, then the trees will be simulated on demand, but stored 
    so that the same tree will be returned with subsequent indexing with the same
    number.
    """

    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes TreeSummarizer object by assigning supplied phycas object
        to a data member variable.
        
        """
        CommonFunctions.__init__(self, opts)
        self.r = self._getLot()
        self.opts.edgelen_dist.setLot(self.r)
        self._current = 0
        self.trees = []
        if not self.taxon_labels and self.n_taxa < 1:
            raise ValueError("TreeSimulator cannot be created with an empty list of taxon_labels and n_taxa set to 0")
        self.stdout.info("Generating %s" % str(self))

    def __str__(self):
        nt = self.n_trees
        n = nt > 0 and "%d " % nt or ""
        return "collection of %stree(s) simulated by the Yule process with edgelen_dist = %s" % (n, _value_for_user(self.opts.edgelen_dist))

    def writeTree(self, tree, name="", rooted=None):
        raise ValueError("TreeSimulator collection of trees cannot be used as an output destination")

    def finish(self):
        raise ValueError("TreeSimulator collection of trees cannot be used as an output destination")

    def __iter__(self):
        self._current = 0
        return self

    def _simulateTree(self):
        tl = self.taxon_labels
        if tl:
            t = Phylogeny.Tree(taxon_labels=self.taxon_labels)
            ntax = len(self.taxon_labels)
        else:
            t = Phylogeny.Tree()
            ntax = self.n_taxa
        tm = Phylogeny.TreeManip(t)
        tm.randomTree(ntax, self.r, self.opts.edgelen_dist, False)
        return t

    def next(self):
        n = self._current
        nt = self.n_trees
        if n >= nt and nt > 0:
            raise StopIteration()
        t = None
        if n < len(self.trees):
            t = self.trees[n]
        if t is None:
            t = self._simulateTree()
        self._current += 1
        if nt > 0:
            self._pad_to(self._current)
            self.trees[n] = t
        return t

    def __getitem__(self, k):
        nt = self.n_trees
        if nt < 1:
            raise TypeError("TreeSimulator without n_trees cannot return trees by index")
        if isinstance(k, slice):
            inds = k.indices(nt)
            l = []
            for i in range(inds[0], inds[1], inds[2]):
                l.append(self.__getitem__(i))
            return l
        if k < 0:
            k += nt
        if k > nt and nt > 0:
            raise IndexError("TreeSimulator index out of range")
        self._pad_to(k+1)
        t = self.trees[k]
        if t is None:
            t = self._simulateTree()
            self.trees[k] = t
        return t

    def _pad_to(self, k):
        diff_len = k - len(self.trees)
        if diff_len > 0:
            self.trees.extend([None]*diff_len)
        
            
    def __len__(self):
        if self.opts.n_trees < 1:
            raise TypeError("TreeSimulator without n_trees has no len()")
        return self.opts.n_trees

    def getNTrees(self):
        return self.opts.n_trees

    def setNTrees(self, x):
        self.opts.n_trees = x

    def getTaxa(self):
        return self.opts.taxon_labels

    def setTaxa(self, x):
        self.opts.taxon_labels = x

    def getNTaxa(self):
        return self.opts.n_taxa

    def setNTaxa(self, x):
        self.opts.n_taxa = x

    n_trees = property(getNTrees, setNTrees)
    n_taxa = property(getNTaxa, setNTaxa)
    taxon_labels = property(getTaxa, setTaxa)
