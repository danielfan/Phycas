from _Phylogeny import *

class TreeManip(TreeManipBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    A manipulator of Tree objects. Can be used to create a tree de novo,
    or rearrange the topology of an existing tree.

    """
    def __init__(self, t):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Constructs a TreeManip object that operates on the supplied tree t.
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> tm = Phylogeny.TreeManip(t)
        >>> d = ProbDist.ExponentialDist(10.0)
        >>> tm.starTree(5, d)
        >>> print t.walkPreorder()
        (0) -> [5] -> (4) -> (3) -> (2) -> (1)

        """
        TreeManipBase.__init__(self, t)

    # POLPY_NEWWAY        
    def buildTreeFromSplitVector(self, split_vect, edge_len_dist):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a tree having the splits defined in split_vect, which should
        be a list or tuple of string representations of splits. All splits in
        split_vect should be compatible. Later splits not compatible with
        earlier ones already in the tree will be ignored.
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> tm = Phylogeny.TreeManip(t)
        >>> v = ['-*****','--****', '--**--', '----**']
        >>> tm.buildTreeFromSplitVector(v, ProbDist.ExponentialDist(10.0))
        >>> print t.walkPreorder()
        (0) -> [6] -> (1) -> [4294967295] -> [4294967295] -> (3) -> (2) -> [4294967295] -> (5) -> (4)

        """
        TreeManipBase.buildTreeFromSplitVector(self, split_vect, edge_len_dist)

    def starTree(self, num_tips, edge_len_dist):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a star tree having num_tips tips and edge lengths drawn at
        random from the supplied edge_len_dist probability distribution.
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> tm = Phylogeny.TreeManip(t)
        >>> d = ProbDist.ExponentialDist(10.0)
        >>> tm.starTree(5, d)
        >>> print t.walkPreorder()
        (0) -> [5] -> (4) -> (3) -> (2) -> (1)

        """
        return TreeManipBase.starTree(self, num_tips, edge_len_dist)

    def yuleTree(self, num_tips, rng, lambd):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a Yule tree having num_tips tips, speciation rate lambd, and
        a topology and edge lengths determined by the supplied random number
        generator rng. num_tips should include the tip used to root the
        resulting tree.
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> rng = ProbDist.Lot(13579)
        >>> tm = Phylogeny.TreeManip(t)
        >>> tm.yuleTree(5, rng, 1.5)
        >>> print t.walkPreorder()
        (0) -> [5] -> [6] -> (1) -> (3) -> [7] -> (2) -> (4)
        >>> print t.makeNewick()
        (1:0.10962,(2:0.29077,4:0.29077):0.00340,(3:0.22779,5:0.22779):0.06638)

        """
        import phycas.ProbDist
        dist = phycas.ProbDist.ExponentialDist(lambd)
        dist.setLot(rng)
        TreeManipBase.randomTree(self, num_tips, rng, dist, True)

    def randTree(self, num_tips, rng, edge_len_dist):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a random tree topology having num_tips tips, with edge lengths
        determined by the supplied probability distribution edge_len_dist.
        num_tips should include the tip used to root the resulting tree.
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> rng = ProbDist.Lot(13579)
        >>> dist = ProbDist.ExponentialDist(0.5)
        >>> dist.setLot(rng)
        >>> tm = Phylogeny.TreeManip(t)
        >>> tm.randTree(5, rng, dist)
        >>> print t.walkPreorder()
        (0) -> [5] -> [6] -> [7] -> (1) -> (4) -> (3) -> (2)
        >>> print t.makeNewick()
        (1:0.90722,((2:0.22471,5:0.43341):0.36279,4:3.41679):0.56683,3:0.02039)

        """
        TreeManipBase.randomTree(self, num_tips, rng, edge_len_dist, False)

    def setRandomEdgeLengths(self, edge_len_dist):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets edge lengths using random draws from the probability distribution
        edge_len_dist.
        
        >>> from phycas import *
        >>> t = Phylogeny.Tree()
        >>> t.buildFromString('(1,(2,5),(3,4))')
        >>> rng = ProbDist.Lot(13579)
        >>> dist = ProbDist.ExponentialDist(0.5)
        >>> dist.setLot(rng)
        >>> tm = Phylogeny.TreeManip(t)
        >>> tm.setRandomEdgeLengths(dist)
        >>> print t.walkPreorder()
        1 -> [2] -> [0] -> 2 -> 5 -> [1] -> 3 -> 4
        >>> print t.makeNewick()
        (1:0.32887,(2:0.02039,5:1.18233):0.90722,(3:3.41679,4:0.72478):0.56683)

        """
        TreeManipBase.setRandomEdgeLens(self, edge_len_dist)

