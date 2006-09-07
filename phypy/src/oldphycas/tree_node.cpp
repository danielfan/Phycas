//#include "phycas/force_include.h"
#include "phypy/src/ncl/misc/string_extensions.hpp"
#include "phypy/src/oldphycas/taxa_manager.hpp"
#include "phypy/src/oldphycas/tree_node.hpp"
#include "phypy/src/ncl/taxa/nxs_taxa_manager.hpp"
#include "phypy/src/oldphycas/tree_node_attributes.hpp"
using std::vector;
using std::string;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::atof;
#endif

#if defined(ALLOW_NEGATIVE_EDGELENS) 
	//@ experimental may need to do more testing to get this option to work.  
	const double TreeNode::edgelen_epsilon	= -1.e300;
#else
	const double TreeNode::edgelen_epsilon	= 1.e-8;
#endif
const double TreeNode::edgelen_default	= 0.1;
boost::shared_ptr<PhoTaxaManager> TreeNode::gTaxaMgr;

std::string TreeNode::GetDebugDescription() const
	{
	std::string retStr;
	retStr << (IsRoot() ? "root node" : IsShootTip() ? "non-root leaf node" : "internal node") << ": ";
	if (IsLeafOrRoot())
		retStr << (GetName().length() > 0 ? GetName().c_str() : "<sin nombre>");
	retStr << " # = " << nodeNum << "\n  ";
	retStr << "edgelength = "<< GetFltEdgeLen() << "\n  ";
	retStr << "split = " << GetSplit().CreatePatternRepresentation();
	return retStr;
	}	
#if defined(USING_NODE_NUMBERS_ONLY)		
	/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
	|	Inlined function allowing read access to protected data member name.
	*/
	string TreeNode::GetName() const
		{
		PHYCAS_ASSERT(gTaxaMgr);
		return gTaxaMgr->GetLabel(GetNodeNumber());
		}

#endif

void TreeNode::Clear()
	{
	lChild = par = rSib = nextPreorder = prevPreorder = NULL;
	nodeNum = UINT_MAX;
	edgelen.f = edgelen_default;
#   if !defined(USING_NODE_NUMBERS_ONLY)
		name.clear();
#   endif
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Clears then sets the split field based on the node number
*/
void TreeNode::TipResetSplit() const
	{
	split.Clear();
	const unsigned n = GetNodeNumber();
	NXS_ASSERT(n != UINT_MAX);
	//if (n > 0)	//POL 21-Oct-2004  holdover from when we were assuming rooting was always at first taxon and that bit was not set?
		split.SetBit(n);	
	}


/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Zig-zags up the right side of the clade rooted at this node until it reaches a leaf. This leaf node is the last node in the clade according to the preorder 
|	sequence. Zig-zagging is necessary because there is no pointer to the right-most child of a node, so must use FindRightChild to move up the clade.
|	Cannot return NULL
*/
const TreeNode *TreeNode::SearchForLastPreorderInClade() const
	{
	const TreeNode *curr = this;
	const TreeNode *rChild = FindRightChild();
	while (rChild != NULL)
		{
		curr = rChild;
		rChild = curr->FindRightChild();
		}
	return curr;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	 Called during rerooting.  The function doesn't set all
| 	of the pointer used by tree nodes.  The caller needs to set the rest
|	Loses the previous rSib 
|	This function Sets: 	this->rSib = nd
|			LastPreorderInThisClade->nextPreorder = nd
|			nd->prevPreorder = LastPreorderInThisClade
|	Does NOT set 	nd->par = this->par
|					nd->LastPreorderInThisClade->nextPreorder = this->LastPreorderInThisClade->nextPreorder
|					(previous this->rSib node's pointer are not altered at all)
*/
void TreeNode::ReplaceRightSibHelper(TreeNode *nd)
	{
	TreeNode *lastPrInClade = FindLastPreorderInClade();
	rSib = nd;
	lastPrInClade->nextPreorder = nd;
	if (nd != NULL)
		nd->prevPreorder = lastPrInClade;
	}
		

/*----------------------------------------------------------------------------------------------------------------------
|	 Called during rerooting.  As the name implies, the function doesn't set all
| 	of the pointer used by tree nodes.  The caller needs to set the rest
|	Loses the previous lChild 
|	This function Sets: 	this->lChild = nd
|							this->nextPreorder = nd
|							nd->prevPreorder = this
|	Does NOT set 	nd->par = this
|					nd->rSib
|					(previous this->lChild node's pointer are not altered at all)
*/
void TreeNode::ReplaceLeftChildHelper(TreeNode *nd)
	{
	nextPreorder = lChild = nd;
	if (nd != NULL)
		nd->prevPreorder = this;
	}
		

/*----------------------------------------------------------------------------------------------------------------------
|	If called with an internal node this function recreates the split representation of the branch 
|	Should be called in POSTORDER.
|	Has no effect on leaves (if the leaves' splits aren't initialized call CreateSplit instead or TipResetSplit if you 
|	it is a tip)
*/
void TreeNode::RefreshSplit() const
	{
	const TreeNode *nd = GetLeftChild();
	if (nd != NULL)
		{
		split.Clear();
		for (; nd != NULL; nd = nd->GetRightSib())
			split.CombineWith(nd->split);
		}
	}

void TreeNode::DeleteSelfAndDescendants(TreeNode *toDie)
	{
	if (toDie != NULL)
		{
		DeleteSelfRightSibsAndDescendants(toDie->GetLeftChild());
		delete toDie;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recursively deletes all children and all right sibs (and their children)
*/
void TreeNode::DeleteSelfRightSibsAndDescendants(TreeNode *toDie)
	{
	if (toDie != NULL)
		{
		DeleteSelfRightSibsAndDescendants(toDie->GetLeftChild());
		DeleteSelfRightSibsAndDescendants(toDie->GetRightSib());
		delete toDie;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Default constructor for TreeNode (there are default values for both args). 
|	Initializes all pointers to NULL and brlen to edgeLen arg.
*/
TreeNode::TreeNode(unsigned ndNum, double edgeLen, bool setSplitBasedOnNum)
	: lChild(NULL),
	par(NULL),
	rSib(NULL),
	nextPreorder(NULL),
	prevPreorder(NULL),
	nodeNum(ndNum),
	plotData(NULL)
	{
	edgelen.f = edgeLen;	
	//if (setSplitBasedOnNum && ndNum != UINT_MAX && ndNum > 0)
	if (setSplitBasedOnNum && ndNum != UINT_MAX)	//POL 22-Oct-2004
		split.SetBit(ndNum);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Default constructor for TreeNode Initializes all pointers to NULL and brlen to edgeLen arg.
*/
TreeNode::TreeNode(double edgeLen)
	: lChild(NULL),
	par(NULL),
	rSib(NULL),
	nextPreorder(NULL),
	prevPreorder(NULL),
	nodeNum(UINT_MAX),
	plotData(NULL)
	{
	edgelen.f = edgeLen;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	copy constructor for TreeNode. Initializes all pointers to NULL
|	Only copies branchlength and nodenum
*/
TreeNode::TreeNode(const TreeNode &r)
	: lChild(NULL),
	par(NULL),
	rSib(NULL),
	nextPreorder(NULL),
	prevPreorder(NULL),
	nodeNum(r.nodeNum),
	edgelen(r.edgelen),
	plotData(NULL)
	{
	CopyPlotData(r);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor for TreeNode. The destructor just itself. Call DeleteSelfAndDescendants to delete a clade
*/
TreeNode::~TreeNode()
	{
	delete plotData;
	DeleteAndClearVecOfPtrs(attr);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Searches for sib on left by starting with lChild of parent and moving across rSib pointers until a node is found 
|	that has the current node as its rSib. Returns NULL if this node is the root of the tree, which has no parent and 
|	no siblings or if `this' is the left most child
|	@impl (do not call directly)
*/
TreeNode * TreeNode::FindLeftSibImpl_() const
	{
	if (par == NULL)
		return NULL;
	TreeNode * nd = par->GetLeftChild();
	if (nd == this)
		return NULL;
	PHYCAS_ASSERT(nd != NULL); // parent should have at least one child
	while (nd != NULL && nd->GetRightSib() != this)
		nd = nd->GetRightSib();
	return nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the rightmost node in this "level" of siblings.  Will return this (never returns NULL) if rSib if NULL.
|	or rSib
|	@impl (do not call directly)
*/
TreeNode *TreeNode::FindRightmostSibOrSelfImpl_() const 
	{
	TreeNode * nd = const_cast<TreeNode *>(this);
	while (nd->GetRightSib() != NULL)
		nd = nd->GetRightSib();
	return nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of immediate descendants (should be 0 or 2 for a tip node and 3 for the root).
*/
unsigned TreeNode::CountChildren() const
	{
	unsigned nDescendants = 0;
	for (const TreeNode *child = GetLeftChild(); child != NULL; child = child->GetRightSib())
		 ++nDescendants;
	return nDescendants;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of nodes (internal and terminal) below this node (returns 0 for leaf nodes)
*/
unsigned TreeNode::CountNumNodesBelow() const
	{
	unsigned nDescendants = 0;
	for (const TreeNode *child = GetLeftChild(); child != NULL; child = child->GetRightSib())
		 nDescendants += 1 + child->CountNumNodesBelow();
		 
	return nDescendants;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of nodes (internal and terminal) below this node (returns 0 for leaf nodes)
*/
unsigned TreeNode::CountNumLeafNodesBelow() const
	{
	unsigned nDescendants = 0;
	for (const TreeNode *child = GetLeftChild(); child != NULL; child = child->GetRightSib())
		nDescendants += (child->IsShootTip() ? 1 : child->CountNumLeafNodesBelow());
	return nDescendants;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	increments *nInternals and *nTerminals to determine the number of internal and leaf nodes below this node 
|	(calling this function for leaf nodes doesn't change nInternals or nTerminals)
*/
void TreeNode::ClassifyNodesBelow(unsigned * const nInternals, unsigned * const nTerminals) const
	{
	for (const TreeNode *child = GetLeftChild(); child != NULL; child = child->GetRightSib())
		{
		if (child->IsShootTip())
			*nTerminals += 1;
		else
			{
			*nInternals += 1;
			child->ClassifyNodesBelow(nInternals, nTerminals);
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Stores nodal information provided in string info. The info string normally represents a branch length read from a
|	newick-style tree definition, but could include some special comments. For now, this function treats info as simply
|	a branch length.
*/
void TreeNode::SetInfo(
  const string &info)	/* the branch length information; converted to float using atof and stored using SetFltEdgeLen member function */
	{
	PHYCAS_ASSERT(info.length() > 0);
	SetFltEdgeLen(atof(info.c_str()));
	}
	
#if 0
void TreeNode::Recurse(void &FuncObj(TreeNode*, int))
	{
	FuncObj(this, 1);
	if (!IsShootTip())
		{
		GetLeftChild()->Recurse(FuncObj);
		FuncObj(this, 2);
		}
	if (GetRightSib() != NULL)
		GetRightSib()->Recurse(FuncObj);
	FuncObj(this, 3);
	}
#endif

#if ALT_TRAVERSE_TREE
	TreeNode *TreeNode::SetPreorderPointers(TreeNode *p)
		{
		if (p != NULL)
			p->nextPreorder = this;
		prevPreorder = p;
		nextPreorder = NULL;
		TreeNode *const lastPreorderInThisClade = (lChild == NULL ? this : lChild->SetPreorderPointers(this));
		PHYCAS_ASSERT(p != NULL);
		return (rSib == NULL ? lastPreorderInThisClade : rSib->SetPreorderPointers(lastPreorderInThisClade));
		}
#endif
