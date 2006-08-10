#if ! defined (PHO_TREENODE_H)
#define	PHO_TREENODE_H
#include <stack> //needed for some of the  templated functions

#include "phycas/trees/split.hpp"

//POL moved Score union and Length typedef to split.h

class TreeNodeAttr;
#if defined (SUSPECTING_FRIEND_NS_ISSUE)
namespace cipres
{
class NexmlNodeElementHandler;
}
#endif
class NodePlotData;
class DrawContext;
class TreeNodeAttribute;
class MultilineBound;
class PhoTaxaManager;
#include "phycas/trees/node_iterator_types.hpp"


/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Encapsulates a node (i.e. vertex) of a tree, with pointers to other nodes for navigation purposes, branch length
|	for purposes of drawing the tree and a pointer to a TreeNodeAttr that can be used to store arbitrary additional 
|	information (e.g. conditional likelihood arrays).
*/
class TreeNode
	{
		friend class HKYAdHocEvaluator; //POL@ used only for star tree paper
#	if defined(USE_PHOREST_PLOT_WINDOWS)
		friend class PhorestTreeCanvas;
#	endif

	public:
		~TreeNode();
		
		STATIC_CONST const double	edgelen_epsilon;	/* smallest allowable edge length */
		STATIC_CONST const double	edgelen_default;	/* default edge length */
		STATIC_DATA boost::shared_ptr<PhoTaxaManager> gTaxaMgr;

	protected:
		TreeNode		*lChild;			/* points to leftmost child */
		TreeNode		*par;				/* points to parent */
		TreeNode		*rSib;				/* points to sibling on right */
		TreeNode		*nextPreorder;		/* points to next node in pre-order sequence */
		TreeNode		*prevPreorder;		/* points to previous node in pre-order sequence */
		unsigned		nodeNum;			/* for tips, this is the taxon index, ranging from 0 to ntips-1 */
		Length			edgelen;			/* length of this node's edge */
		mutable Split	split;				/* used for identifying split defined by this node's edge */
		std::vector<TreeNodeAttribute *> attr;	/* attributes of this node that store the information needed to score a tree, acts as if the attributes were mutable*/
#		if !defined(USING_NODE_NUMBERS_ONLY)
			std::string		name;				/* name of node */
#		endif

	public:
		NodePlotData	*plotData;			/* only used when plotting tree */

		// Accessors
		//
		std::string				GetDebugDescription() const;
		Length					GetEdgeLen() const;
		double					GetFltEdgeLen() const;
		unsigned				GetIntEdgeLen() const;
		std::string				GetName() const;
		unsigned				GetNodeNumber() const;
		TreeNodeAttribute	   *GetTreeNodeAttribute(unsigned i) const;
		Split					GetSplit() const;
		
		// Queries
		//
		bool		IsLeafOrRoot() const;
		bool		IsShootTip() const;
		bool		IsRoot() const;
		bool		IsInternal() const {return !IsLeafOrRoot();}
		bool		EdgeIsInternal() const {return (IsInternal() && !(GetParent()->IsRoot()));}
		bool		IsPlotData() const;
		bool		IsFirstPostorderChild() const;
		bool 		IsSibling(const TreeNode* other) const;
		bool 		IsAdjacent(const TreeNode* other) const;
		bool		IsSelected() const;
		
		const 	TreeNode	*GetLeftChild() const;
		const 	TreeNode	*GetRightSib() const;
		const 	TreeNode	*GetParent() const;
		const 	TreeNode	*GetNextPreorder() const;
		const 	TreeNode	*GetNextPostorder() const;

				TreeNode	*GetLeftChild() {return lChild;}
				TreeNode	*GetRightSib()	{return rSib;}
				TreeNode	*GetParent()	{return par;}		
				TreeNode	*GetNextPreorder()	{return nextPreorder;}
				TreeNode	*GetNextPostorder() {return prevPreorder;}

		
		// Utilities
		//
		const TreeNode    * FindLeftSib() const					{   return FindLeftSibImpl_(); }
		const TreeNode    * FindRightChild() const				{   return FindRightChildImpl_(); }
		const TreeNode    * FindRightmostSibOrSelf() const		{   return FindRightmostSibOrSelfImpl_(); }
		const TreeNode	  * FindLastPreorderInClade() const		{   return FindLastPreorderInCladeImpl_(); }
		const TreeNode	  * FindFirstPreorderAfterClade() const {   return FindFirstPreorderAfterCladeImpl_(); }
		TreeNode		  * FindLeftSib()						{   return FindLeftSibImpl_(); }
		TreeNode		  * FindRightChild()					{   return FindRightChildImpl_(); }
		TreeNode		  * FindRightmostSib()					{   return FindRightmostSibOrSelfImpl_(); }
		TreeNode		  * FindLastPreorderInClade()			{   return FindLastPreorderInCladeImpl_(); }
		TreeNode		  * FindFirstPreorderAfterClade()		{   return FindFirstPreorderAfterCladeImpl_(); }
				unsigned	CountChildren() const;
				unsigned 	CountNumNodesBelow() const;
				unsigned 	CountNumLeafNodesBelow() const;
				void		ClassifyNodesBelow(unsigned * const nInternals, unsigned *const nTerminals) const;

		// Modifiers
		//
		void		RefreshSplit() const;
		void		SetNodeNumber(unsigned n, bool isTip); //if isTip, then the split will be reset.
		void		SetEdgeLen(Length x);
		void		SetFltEdgeLen(double x);
		void		SetIntEdgeLen(int x);
		void		SetInfo(const std::string &info);
		void		SetPct(double);
		void		SetSelectionStatus(bool s);

		// Drawing functions
		double		GetX() const;
		double		GetY() const;
		double		GetPct() const;

		void		AddToX(double);
		void		ScaleCoords(double xscaler, double yscaler);
		void		SetX(double newX);
		void		SetY(double newY);

		void		CopyPlotData(const TreeNode &r);
		void		RefreshPlotData();
		
		STATIC_DATA void DeleteSelfAndDescendants(TreeNode *);
		
		pre_iterator		clade_pre_begin();
		const_pre_iterator	clade_pre_begin() const;
		post_iterator		clade_post_begin();
		const_post_iterator	clade_post_begin() const;
		pre_iterator		clade_pre_end();
		const_pre_iterator	clade_pre_end() const;
		post_iterator		clade_post_end();
		const_post_iterator	clade_post_end() const;

	private:
		TreeNode(double edgeLen = DBL_MAX);
		TreeNode(unsigned ndNum, double edgeLen , bool intializeSplit);
		TreeNode(const TreeNode &); //only to be used by Tree does not copy the entire treestructure (only nodeNum and branch length)
		TreeNode &operator=(const TreeNode &); //not defined.  do not use.  declaration added to prevent compiler from supplying a definition

		TreeNode		  * FindLeftSibImpl_() const;
		TreeNode		  * FindRightChildImpl_() const;
		TreeNode		  * FindRightmostSibOrSelfImpl_() const;
		TreeNode		  * FindLastPreorderInCladeImpl_() const;
		TreeNode		  * FindFirstPreorderAfterCladeImpl_() const;

		void				AddLeftChild(TreeNode *lc);
		void				AddRightSib(TreeNode *lc);
		void				Clear();
		void				CreateSplit() const;
		void				ReplaceLeftChildHelper(TreeNode *nd); //called in tree rerooting 
		void				ReplaceRightSibHelper(TreeNode *nd); //called in tree rerooting 
		const TreeNode    * SearchForLastPreorderInClade() const;
		void				TipResetSplit() const;
		
#if !defined (SUSPECTING_FRIEND_NS_ISSUE)
	public:
#endif
		STATIC_DATA void DeleteSelfRightSibsAndDescendants(TreeNode *);
#		if ALT_TRAVERSE_TREE
			TreeNode	*SetPreorderPointers(TreeNode *theNewPreviousPreOrderForThisNode); //returns prevPreOrder for the next node
#		endif
		

		friend class Tree;	
		friend class ScoreableTree;	
		friend class TreeManip;
		friend class SplitManager;
		friend class EdgeLenAccumulator;
#if defined (SUSPECTING_FRIEND_NS_ISSUE)
		friend class cipres::NexmlNodeElementHandler;
#endif
	};
#if !defined(USING_NODE_NUMBERS_ONLY)											
/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Allows read access to protected data member `name'.
*/
inline std::string TreeNode::GetName() const
	{
	return name;
	}
#endif

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Sets `nodeNum' data member to `n'.
*/
inline void TreeNode::SetNodeNumber(unsigned n, bool setSplitByNum)
	{
	nodeNum = n;
	if (setSplitByNum)
		TipResetSplit();
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Allows read access to protected data member `nodeNum'.
*/
inline unsigned TreeNode::GetNodeNumber() const
	{
	return nodeNum;
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Allows read access to protected data member `nodeNum'.
*/
inline Split TreeNode::GetSplit() const
	{
	return split;
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Allows read access to protected data member `nextPreorder'.
*/
inline const TreeNode* TreeNode::GetNextPreorder() const
	{
	return nextPreorder;
	}
	
/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Returns `prevPreorder', which is the same thing as returning the next node in the post-order sequence.
*/
inline const TreeNode* TreeNode::GetNextPostorder() const
	{
	return prevPreorder;
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Allows read access to protected data member `par'.
*/
inline const TreeNode* TreeNode::GetParent() const
	{
	return par;
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Allows read access to protected data member `lchild'.
*/
inline const TreeNode* TreeNode::GetLeftChild() const
	{
	return lChild;
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	returns the first preorder node that is not a descendant of `this'
|	Returned node will be this's rSib, or par->rSib, or par->par->rSib, etc 
|	will return NULL for any "right" descendant of the root.
|	@impl (do not call directly)
*/
inline TreeNode *TreeNode::FindFirstPreorderAfterCladeImpl_() const
	{
	if (GetRightSib() != NULL)
		return const_cast<TreeNode *>(GetRightSib());
	return const_cast<TreeNode *>(SearchForLastPreorderInClade()->GetNextPreorder());
	}
	
/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	returns a pointer to this or thee rightmost descendant of this
|	@impl (do not call directly)
*/
inline TreeNode *TreeNode::FindLastPreorderInCladeImpl_() const
	{
	NXS_ASSERT(GetRightSib() == NULL || GetRightSib()->GetNextPostorder() != NULL);
	if (IsShootTip())
		return const_cast<TreeNode *>(this);
	if (GetRightSib() != NULL)
		return const_cast<TreeNode *>(GetRightSib()->GetNextPostorder());
	return const_cast<TreeNode *>(SearchForLastPreorderInClade());
	}
	
	
/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Allows read access to protected data member `rSib'.
*/
inline const TreeNode* TreeNode::GetRightSib() const
	{
	return rSib;
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Returns the EgdeLen union, interpreted as a Length
*/
inline Length TreeNode::GetEdgeLen() const
	{
	return edgelen;
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Returns the EgdeLen union, interpreted as a double
*/
inline double TreeNode::GetFltEdgeLen() const
	{
	return GetEdgeLen().f;
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Returns the EgdeLen union, interpreted as an unsigned
*/
inline unsigned TreeNode::GetIntEdgeLen() const
	{
	return GetEdgeLen().i;
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Allows write access to protected data member brlen.
*/
inline void TreeNode::SetEdgeLen(Length x)
	{
	edgelen = x;
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Allows write access to protected data member brlen. This function should always be used to set
|	the branch length, even when using a class that is a friend of TreeNode, because this provides a means of performing
|	validation of the new branch length.
*/
inline void TreeNode::SetFltEdgeLen(double x)
	{
	//@POL should we do this (make sure branch lengths are always at least TreeNode::edgelen_epsilon)?
	// Will doing this throw off the Hastings ratio?
	//
	Length newEdgeLen;
	newEdgeLen.f = (x < TreeNode::edgelen_epsilon ? TreeNode::edgelen_epsilon : x);
	SetEdgeLen(newEdgeLen);
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Allows write access to protected data member brlen. This function should always be used to set
|	the branch length, even when using a class that is a friend of TreeNode, because this provides a means of perform
|	validation of the new branch length.
*/
inline void TreeNode::SetIntEdgeLen(int x)
	{
	Length newEdgeLen;
	newEdgeLen.i = (x < 0 ? 0 : (unsigned) x);
	SetEdgeLen(newEdgeLen);
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Returns true if node is the root node, false otherwise.
*/
inline bool TreeNode::IsRoot() const
	{
	return (par == NULL);
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Returns true if node is a non-root leaf node, false otherwise.
*/
inline bool TreeNode::IsShootTip() const
	{
	return (GetLeftChild() == NULL);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Searches for rightmost child by starting with lChild and moving down rSib pointers until a node is found that has a
|	NULL rSib pointer. Returns NULL if this is a leaf node with no children.
|	@impl (do not call directly)
*/
inline TreeNode *TreeNode::FindRightChildImpl_() const
	{
	return const_cast<TreeNode *>(GetLeftChild() == NULL ? NULL : GetLeftChild()->FindRightmostSibOrSelf());
	}


/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Returns true if node is a tip node (the root or a leaf), false otherwise.
*/
inline bool TreeNode::IsLeafOrRoot() const
	{
	return (GetLeftChild() == NULL || GetParent() == NULL);
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Returns true if plotData has been allocated for this node, false if plotData equals NULL.
*/
inline bool TreeNode::IsPlotData() const
	{
	return (plotData != NULL);
	}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Returns true if the argument is a direct descendant or the anc of the node `other'.
*/
inline bool TreeNode::IsAdjacent(
  const TreeNode* other) const	
	{
	NXS_ASSERT(other);
	return (par == other || other->par == this);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if node is the first descendant of its immediate ancestor in the post-order traversal sequence to be 
|	visited. Example of where this is useful: when computing splits, it is the responsibility of the first descendant 
|	of an internal node to clear that internal node's split. Afterwards, each descendant adds its own bit(s) to the 
|	internal node's split.
*/
inline bool TreeNode::IsFirstPostorderChild() const
	{
	//return (!IsRoot() && !GetParent()->GetRightSib());
	return (!IsRoot() && !GetRightSib());	//@POL correct?
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns pointer to the requested TreeNodeAttribute structure (assumes that i is less than the number of stored 
|	attributes.
|	NOTE the attributes are "logically" mutable from the perspective of the node (so that scoring a tree is const) so
|	we return a non-const pointer from a const node
*/
inline  TreeNodeAttribute *TreeNode::GetTreeNodeAttribute(unsigned i) const
	{
	NXS_ASSERT(i < attr.size());
	return attr[i];
	}
	
// we might want to move the following class and TreeNode functions to a tree drawing header (e.g. tree/drawtree.h_
class NodePlotData
	{
	public:
						NodePlotData();

		bool			selected;		/* if true, this node and its branch have been selected by user */
		double			x;				/* x coordinate for use in drawing tree */
		double			y;				/* y coordinate for use in drawing tree */
		double			pct;			/* nodal support percentage */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes NodePlotData data members.
*/
inline NodePlotData::NodePlotData()
	{
	selected	= false;
	x			= 0.0;
	y			= 0.0;
	pct			= -1.0;	// negative value serves to indicate not yet set to anything meaningful
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets plotData->selected data member to the argument `s'. Assumes `plotData' is non-NULL.
*/
inline void TreeNode::SetSelectionStatus(
  bool s)
	{
	NXS_ASSERT(plotData != NULL);
	plotData->selected = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if plotData->selected is true, and returns false otherwise. Assumes `plotData' is non-NULL.
*/
inline bool TreeNode::IsSelected() const
	{
	NXS_ASSERT(plotData != NULL);
	return plotData->selected;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows read access to protected data member `x'.
*/
inline double TreeNode::GetX() const
	{
	NXS_ASSERT(plotData != NULL);
	return plotData->x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows read access to protected data member `y'.
*/
inline double TreeNode::GetY() const
	{
	NXS_ASSERT(plotData != NULL);
	return plotData->y;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if node is the root node, false otherwise.
*/
inline void TreeNode::SetPct(
  double p)
	{
	NXS_ASSERT(plotData != NULL);
	plotData->pct = p;
	}
		
/*----------------------------------------------------------------------------------------------------------------------
|	Allows read access to protected data member `pct'.
*/
inline double TreeNode::GetPct() const
	{
	NXS_ASSERT(plotData != NULL);
	return plotData->pct;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	sets the plotdata->x field
*/
inline void TreeNode::SetX(double newX)
	{
	NXS_ASSERT(plotData != NULL);
	plotData->x = newX;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	sets the plotdata->y field
*/
inline void TreeNode::SetY(double newY)
	{
	NXS_ASSERT(plotData != NULL);
	plotData->y = newY;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates a new NodePlotData structure if necessary
*/
inline void TreeNode::RefreshPlotData()
	{
	if (plotData == NULL)
		plotData = new NodePlotData;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `r->plotData' is not NULL, creates a new NodePlotData structure if necessary copies `r->plotData' to it.
*/
inline void TreeNode::CopyPlotData(const TreeNode &r)
	{
	if (r.IsPlotData())
		{
		if (plotData == NULL)
			{
			plotData = new NodePlotData;
			}
		plotData->x			= r.plotData->x;
		plotData->y			= r.plotData->y;
		plotData->pct		= r.plotData->pct;
		plotData->selected	= r.plotData->selected;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds lc as the lChild of `this' fixing this's lChild, nextPreorder pointers and lc's prevPreorder and par pointer
*/
inline void TreeNode::AddLeftChild(TreeNode *lc)
	{
	ReplaceLeftChildHelper(lc);
	lc->par = this;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Adds rs as the rSib of `this' fixing this's rSib pointer,
|	The last preorder descendant of `this's nextPreorder pointer, 
|	and rs's prevPreorder and par pointer
*/
inline void TreeNode::AddRightSib(TreeNode *rs)
	{
	ReplaceRightSibHelper(rs);
	rs->par = par;
	}
		
/*----------------------------------------------------------------------------------------------------------------------
|	Codes the current tree structure into the split datamember.  This function affects even the leaf nodes.
|	If the leaf nodes are still attached to the same taxa, but the tree structure has changed it is sufficient 
|	to call RefreshSplit();
|	Should be called in POSTORDER
*/
inline void TreeNode::CreateSplit() const
	{
	if (IsShootTip())
		TipResetSplit();
	else
		RefreshSplit();
	}

#endif
