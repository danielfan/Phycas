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
        
        >>> from ProbDist import *
        >>> from Phylogeny import *
        >>> t = Tree()
        >>> tm = TreeManip(t)
        >>> d = ExponentialDist(10.0)
        >>> tm.starTree(5, d)
        >>> print t.walkPreorder()
        (0) -> [5] -> (4) -> (3) -> (2) -> (1)

        """
        TreeManipBase.__init__(self, t)
        
    def starTree(self, num_tips, edge_len_dist):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a star tree having num_tips tips and edge lengths drawn at
        random from the supplied edge_len_dist probability distribution.
        
        >>> from ProbDist import *
        >>> from Phylogeny import *
        >>> t = Tree()
        >>> tm = TreeManip(t)
        >>> d = ExponentialDist(10.0)
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
        
        >>> from ProbDist import *
        >>> from Phylogeny import *
        >>> t = Tree()
        >>> rng = Lot(13579)
        >>> tm = TreeManip(t)
        >>> tm.yuleTree(5, rng, 1.5)
        >>> print t.walkPreorder()
        (0) -> [5] -> [6] -> (1) -> (3) -> [7] -> (2) -> (4)
        >>> print t.makeNewick()
        (1:0.10962,(2:0.29077,4:0.29077):0.00340,(3:0.22779,5:0.22779):0.06638)

        """
        import ProbDist
        dist = ProbDist.ExponentialDist(lambd)
        dist.setLot(rng)
        TreeManipBase.randomTree(self, num_tips, rng, dist, True)

    def randTree(self, num_tips, rng, edge_len_dist):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a random tree topology having num_tips tips, with edge lengths
        determined by the supplied probability distribution edge_len_dist.
        num_tips should include the tip used to root the resulting tree.
        
        >>> from ProbDist import *
        >>> from Phylogeny import *
        >>> t = Tree()
        >>> rng = Lot(13579)
        >>> dist = ExponentialDist(0.5)
        >>> dist.setLot(rng)
        >>> tm = TreeManip(t)
        >>> tm.randTree(5, rng, dist)
        >>> print t.walkPreorder()
        (0) -> [5] -> [6] -> [7] -> (1) -> (4) -> (3) -> (2)
        >>> print t.makeNewick()
        (1:0.90722,((2:0.22471,5:0.43341):0.36279,4:3.41679):0.56683,3:0.02039)

        """
        TreeManipBase.randomTree(self, num_tips, rng, edge_len_dist, False)

