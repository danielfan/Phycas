from _Phylogeny import *

class TreeNode(TreeNodeBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates the notion of a node in a phylogenetic tree. 

    """
    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes data members.

        """
        TreeNodeBase.__init__(self)

    def getX(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current x-coordinate of this node. 

        """
        return TreeNodeBase.getX(self)
    
    def getY(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current y-coordinate of this node. 

        """
        return TreeNodeBase.getY(self)
    
    def setX(self, new_x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the x-coordinate of this node to the value new_x. 

        """
        TreeNodeBase.setX(self, new_x)
    
    def setY(self, new_y):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the y-coordinate of this node to the value new_y. 

        """
        TreeNodeBase.setY(self, new_y)
    
    def isSelected(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if the node is currently selected, False otherwise. 

        """
        return TreeNodeBase.isSelected(self)
    
    def selectNode(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Selects this node. 

        """
        TreeNodeBase.selectNode(self)
    
    def unselectNode(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Unselects this node. 

        """
        TreeNodeBase.unselectNode(self)
    
    def isRoot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns true if this node is currently serving as the root of the
        tree. 

        """
        return TreeNodeBase.isRoot(self)
    
    def isTip(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns true if this node is a tip node (i.e. a node having degree
        one).

        """
        return TreeNodeBase.isTip(self)
    
    def isInternal(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns true if this node is an internal node (i.e. a node having
        degree greater than one).

        """
        return TreeNodeBase.isInternal(self)
    
    def getLeftChild(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the left child node, or None if there is no left child.

        """
        return TreeNodeBase.getLeftChild(self)
    
    def getRightSib(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the right sibling node, or None if there is no right sib.

        """
        return TreeNodeBase.getRightSib(self)
    
    def getParent(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the parent node, or None if there is no parent.

        """
        return TreeNodeBase.getParent(self)
    
    def getNodeNumber(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the node number (as an int).

        """
        return TreeNodeBase.getNodeNumber(self)
    
    def getNodeName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the node name (as a string).

        """
        return TreeNodeBase.getNodeName(self)
    
    def getEdgeLen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the length of the edge subtending this node (as a float).

        """
        return TreeNodeBase.getEdgeLen(self)
    
    def setEdgeLen(self, new_edgelen):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the edge length to new_edgelen.

        """
        TreeNodeBase.setEdgeLen(self, new_edgelen)
    
    def getNextPreorder(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the next node in preorder sequence. If there is no next node,
        returns None.

        """
        return TreeNodeBase.getNextPreorder(self)

    def getNextPostorder(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the next node in postorder sequence. If there is no next node,
        returns None.

        """
        return TreeNodeBase.getNextPostorder(self)
