from _Phylogeny import *

class Tree(TreeBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates the notion of a phylogenetic tree. This class has only
    methods necessary for representing the tree, copying trees from other
    trees, manipulating the tree, and representing the tree graphically.
    It does not perform any specialized activities, such as computing its
    own likelihood: these sorts of activities are left to member
    functions of other classes.

    """
    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls clear() to initialize data members.

        >>> from phycas.Phylogeny import *
        >>> t1 = Tree()
        >>> print t1.getNNodes()
        0

        """
        TreeBase.__init__(self)
        
    def __iter__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns nodes in preorder sequence.

        >>> from phycas.Phylogeny import *
        >>> tree = Tree()
        >>> tree.buildFromString('(a,b,(c,(d,e)x)y)z')
        >>> for nd in tree:
        ...     print nd.getNodeName(),
        ... 
        a z b y c x d e

        """
        TreeBase.refreshPreorder(self, None)
        nd = TreeBase.getFirstPreorder(self)
        while nd:
            yield nd
            nd = nd.getNextPreorder()
        
    def nodesWithEdges(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns nodes that have edges, which is all nodes in an unrooted tree
        except the root node. The nodes are visited in preorder sequence.

        >>> from phycas.Phylogeny import *
        >>> tree = Tree()
        >>> tree.buildFromString('(a,b,(c,(d,e)x)y)z')
        >>> for nd in tree.nodesWithEdges():
        ...     print nd.getNodeName(),
        ... 
        z b y c x d e

        """
        TreeBase.refreshPreorder(self, None)
        nd = TreeBase.getFirstPreorder(self)
        while nd:
            if not nd.isRoot():
                yield nd
            nd = nd.getNextPreorder()
        
    def tipNodes(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns nodes that correspond to tip nodes, where a tip node is any
        node (including the root node) that is of degree one (i.e. has only
        one edge attached to it). The nodes are visited in preorder sequence.

        >>> from phycas.Phylogeny import *
        >>> tree = Tree()
        >>> tree.buildFromString('(a,b,(c,(d,e)x)y)z')
        >>> for nd in tree.tipNodes():
        ...     print nd.getNodeName(),
        ... 
        a b c d e

        """
        TreeBase.refreshPreorder(self, None)
        nd = TreeBase.getFirstPreorder(self)
        while nd:
            if nd.isTip():
                yield nd
            nd = nd.getNextPreorder()
        
    def internalNodes(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns nodes that correspond to internal nodes, where an internal
        node is any node that has degree greater than two (i.e. has at least
        three edges attached to it). The nodes are visited in preorder
        sequence.

        >>> from phycas.Phylogeny import *
        >>> tree = Tree()
        >>> tree.buildFromString('(a,b,(c,(d,e)x)y)z')
        >>> for nd in tree.internalNodes():
        ...     print nd.getNodeName(),
        ... 
        z y x

        """
        TreeBase.refreshPreorder(self, None)
        nd = TreeBase.getFirstPreorder(self)
        while nd:
            if nd.isInternal():
                yield nd
            nd = nd.getNextPreorder()
        
    def buildFromString(self, newick):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Builds a tree from a newick tree description. If edge lengths are
        specified in the tree description, they must be specified for all
        edges. No checking of node names is performed as the tree is created.
        Before returning, roots tree at first tip node. Raises an XPhylogeny
        exception if a problem is encountered.

        >>> from phycas.Phylogeny import *
        >>> t2 = Tree()
        >>> t2.buildFromString('(c,d,(a, b))')
        >>> print t2.walkPreorder()
        c -> [5] -> d -> [4] -> a -> b

        """
        TreeBase.buildFromString(self ,newick)

    def clear(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns object to just-constructed state.

        >>> from phycas.Phylogeny import *
        >>> t3 = Tree()
        >>> t3.buildFromString('(c,d,(a, b))')
        >>> t3.getNNodes()
        6
        >>> t3.clear()
        >>> t3.getNNodes()
        0

        """
        TreeBase.clear(self)

    def getNInternals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of internal nodes currently composing the tree
        (i.e. getNInternals(). The number of internal nodes equals the number
        of nodes in the tree having degree 2 or higher. Its value is 1 for a
        (rooted or unrooted) star tree, n-2 for a fully-resolved unrooted
        tree, and n-1 for a fully-resolved rooted tree. Calls
        refreshNodeCounts() if node counts have been invalidated.

        >>> from phycas.Phylogeny import *
        >>> t4 = Tree()
        >>> t4.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t4.getNNodes()
        6
        >>> print t4.getNInternals()
        2

        """
        return TreeBase.getNInternals(self)

    def getNTips(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of leaf nodes currently composing the tree (i.e.
        getNLeaves(). The number of leaves equals the number of degree-1 nodes
        in the tree. Note that one of these degree-1 nodes is the root node.
        The number of leaves equals the number of taxa for unrooted trees.
        For rooted trees of n taxa, there are n+1 leaves because the root node
        is a leaf node but not a taxon. Calls refreshNodeCounts() if node
        counts have been invalidated.

        >>> from phycas.Phylogeny import *
        >>> t5 = Tree()
        >>> t5.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t5.getNNodes()
        6
        >>> print t5.getNTips()
        4

        """
        return TreeBase.getNTips(self)

    def getNNodes(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the total number of nodes currently composing the tree (i.e.
        getNLeaves() plus getNInternals(). For unrooted trees of n taxa, one
        of the tip nodes serves as the root node, so this number can range
        from n+1 (star tree) to 2n-2 (fully-resolved). For rooted trees, the
        range is n+2 (star tree) to 2n-1 (fully-resolved) because the root
        node exists but is not identical to one of the tips. Calls
        refreshNodeCounts() if node counts have been invalidated.

        >>> from phycas.Phylogeny import *
        >>> t6 = Tree()
        >>> t6.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t6.getNNodes()
        6

        """
        return TreeBase.getNNodes(self)

    def isRooted(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if tree is rooted, or False if the root is actually a
        tip.

        >>> from phycas.Phylogeny import *
        >>> t7 = Tree()
        >>> t7.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t7.isRooted()
        False

        """
        return TreeBase.isRooted(self)

    def hasEdgeLens(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if tree has edge lengths defined, False if edge lengths
        have not yet been set.

        >>> from phycas.Phylogeny import *
        >>> t8 = Tree()
        >>> t8.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t8.hasEdgeLens()
        False
        >>> t8.buildFromString('(fish:0.2,shark:0.25,(bird:0.1, mammal:0.15):0.09)')
        >>> print t8.hasEdgeLens()
        True
        
        """
        return TreeBase.hasEdgeLens(self)

    def walkPreorder(self, verbosity = 0):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Walks through the tree in preorder fashion (visits parents before
        descendants, and leftmost child before siblings), building up a string
        showing the path taken. Each unnamed internal node is represented by
        the node number (0..getNInternals() - 1) in square brackets, and each
        unnamed leaf node is represented by the node number (0..getNLeaves()
        - 1) in parentheses. If a name has been provided for a node, that
        name is used instead of the node number, both for internal and leaf
        nodes.

        >>> from phycas.Phylogeny import *
        >>> t9 = Tree()
        >>> t9.buildFromString('(fish,shark,(bird, mammal))')
        >>> print t9.walkPreorder()
        fish -> [5] -> shark -> [4] -> bird -> mammal

        """
        return TreeBase.debugWalkTree(self, True, verbosity)

    def walkPostorder(self, verbosity = 0):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Walks through the tree in postorder fashion (visiting children before
        parents and siblings on right before those on the left, building up a
        string showing the path taken. Each unnamed internal node is represented by
        the node number (0..getNInternals() - 1) in square brackets, and each
        unnamed leaf node is represented by the node number (0..getNLeaves()
        - 1) in parentheses. If a name has been provided for a node, that
        name is used instead of the node number, both for internal and leaf
        nodes.

        >>> from phycas.Phylogeny import *
        >>> tA = Tree()
        >>> tA.buildFromString('(fish,shark,(bird, mammal))')
        >>> print tA.walkPostorder()
        mammal -> bird -> [4] -> shark -> [5] -> fish

        """
        return TreeBase.debugWalkTree(self, False, verbosity)

    def rerootAtTip(self, num):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reroots the tree at the leaf node numbered num. An XPhylogeny
        exception is raised if a leaf node having number num cannot be found.

        >>> from phycas.Phylogeny import *
        >>> tC = Tree()
        >>> tC.buildFromString('((a,b),c,(d,e))')
        >>> print tC.walkPreorder()
        a -> [5] -> b -> [7] -> c -> [6] -> d -> e
        >>> tC.rerootAtTip(5)
        Traceback (most recent call last):
            ...
        Exception: there is no tip node having number 5
        >>> tC.rerootAtTip(4)
        >>> print tC.walkPreorder()
        e -> [6] -> d -> [7] -> c -> [5] -> b -> a

        """
        return TreeBase.rerootAtTip(self, num)

    def edgeLenSum(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sums all edge lengths in the tree. Raises an Exception if edge lengths
        were never provided for the tree.

        >>> from phycas.Phylogeny import *
        >>> tB = Tree()
        >>> tB.buildFromString('(fish,shark,(bird, mammal))')
        >>> print tB.edgeLenSum()
        Traceback (most recent call last):
            ...
        Exception: no edge lengths were specified for this tree
        >>> tB.buildFromString('(fish:0.2,shark:0.25,(bird:0.1, mammal:0.15):0.09)')
        >>> print round(tB.edgeLenSum(), 3)
        0.79

        """
        return TreeBase.edgeLenSum(self)

    def makeNewick(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string representing the tree in newick format (i.e. using
        nested parentheses). If a node is named, the name will be used in the
        tree description. Tip nodes are represented by node numbers if they
        do not have a name. If edge lengths are present, they will be shown.

        >>> from phycas.Phylogeny import *
        >>> t = Tree()
        >>> t.buildFromString('((a:0.11,b:0.12)x:0.10,c:0.21,(d:0.31,e:0.32)y:0.30)z')
        >>> print t.makeNewick()
        (a:0.11000,b:0.12000,(c:0.21000,(d:0.31000,e:0.32000)y:0.30000)z:0.10000)x

        """
        return TreeBase.makeNewick(self)

    def rectifyNumbers(self, taxon_names):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets tip node numbers in the tree to the index of the tip's name in
        the supplied object taxon_names, which should be a list or tuple of
        taxon names.

        >>> from phycas.Phylogeny import *
        >>> t = Tree()
        >>> t.buildFromString('((a,b),c,(d,e))')
        >>> print t.walkPreorder(verbosity=1)
        a (0) -> ? [5] -> b (1) -> ? [7] -> c (2) -> ? [6] -> d (3) -> e (4)
        >>> t.rectifyNumbers(['e','d','c','b','a'])
        >>> print t.walkPreorder(verbosity=1)
        a (4) -> ? [5] -> b (3) -> ? [7] -> c (2) -> ? [6] -> d (1) -> e (0)

        """
        return TreeBase.rectifyNumbers(self, taxon_names)

    def rectifyNames(self, taxon_names):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets tip node names in the tree to the names in the supplied object
        taxon_names, which should be a list or tuple of taxon names. The node
        numbers are treated as indices into the taxon_names list.

        >>> from phycas.Phylogeny import *
        >>> t = Tree()
        >>> t.buildFromString('((a,b),c,(d,e))')
        >>> print t.walkPreorder(verbosity=1)
        a (0) -> ? [5] -> b (1) -> ? [7] -> c (2) -> ? [6] -> d (3) -> e (4)
        >>> t.rectifyNames(['e','d','c','b','a'])
        >>> print t.walkPreorder(verbosity=1)
        e (0) -> ? [5] -> d (1) -> ? [7] -> c (2) -> ? [6] -> b (3) -> a (4)

        """
        return TreeBase.rectifyNames(self, taxon_names)

    def unselectAllNodes(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Turns the selected attribute off for all nodes in the tree.

        """
        return TreeBase.unselectAllNodes(self)

    def hasEdgeLens(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if edge lengths have been specified for the tree, False
        otherwise.

        """
        return TreeBase.hasEdgeLens(self)

    def getFirstPreorder(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the root node. Use the TreeNode getNextPreorder to walk
        through the other nodes in the tree in preorder fashion.

        """
        return TreeBase.getFirstPreorder(self)

    def getFirstPostorder(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the node at the upper right corner of the tree (the last
        node in the preorder sequence, which is also the first node in the
        postorder sequence). Use the TreeNode getNextPostorder to walk
        through the other nodes in the tree in postorder fashion.

        """
        return TreeBase.getFirstPostorder(self)