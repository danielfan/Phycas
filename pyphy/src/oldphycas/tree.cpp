//#define POL_MESQUITE_COLORING
	
//#include "phycas/force_include.h"
#include <boost/function.hpp>
#include "pyphy/src/oldphycas/taxa_manager.hpp"
#include "phycas/trees/tree.hpp"
#include "phycas/trees/tree_node.hpp"
#include "phycas/trees/tree_manip.hpp"
#include "phycas/trees/tree_inl.hpp"
#include "ncl/nxs_token.hpp"
#include "ncl/taxa/nxs_taxa_manager.hpp"
#include "ncl/trees/full_tree_description.hpp"
#include "ncl/output/nxs_output.hpp"
#include "ncl/misc/string_extensions.hpp"
using std::set;
using std::ostream;
using std::string;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::printf;
	using std::sprintf;
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Unselects all nodes in the tree.
*/
void Tree::RefreshSelectionStatus()
	{
	for (TreeNode *nd = GetRoot(); nd != NULL; nd = nd->GetNextPreorder())
		{
		nd->RefreshPlotData();
		nd->SetSelectionStatus(false);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes ID for leaves as well as internals (RefreshID just does internals)
*/
void Tree::RefreshFullID()
	{
	GetRoot()->split.Clear();
	ForAllLeaves(boost::function1<void, TreeNode*>(&TreeNode::CreateSplit));
	RefreshID();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls FindTip to find the node corresponding to taxon `i', then calls Reroot to reroot the tree at that node. If 
|	node having number equal to `i' is not found, returns false; otherwise, returns true.
*/
bool Tree::RerootAtTip(
  unsigned i,			/*< is the index of the node to make the root */
  bool refresh_splits)	/*< if true, a traversal of the entire tree is done afterwards to ensure all splits are correct; specify false if the traversal is to be postponed (but in this case be sure to call Tree::RefreshID sometime soon!) */
	{
	TreeNode *nd = FindTip(i);
	if (nd == NULL)
		return false;
	Reroot(nd);
	if (refresh_splits)
		RefreshID();
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Finds the node corresponding to taxon `i', returning a pointer to the node if successful, or NULL if no tip node 
|	having that node number can be located.
*/
TreeNode *Tree::FindTip(unsigned tip_number)
	{
	TreeNode *nd = NULL;
	for (nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsLeafOrRoot() && nd->GetNodeNumber() == tip_number)
			break;
		}
	return nd;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Instantiates a new node (or returns a currently unused on).

TreeNode *Tree::CreateTreeNode()
	{
	TreeNode *newNd;
	if (nodeStorage.empty())
		newNd = new TreeNode();
	else
		{
		newNd = nodeStorage.top();
		nodeStorage.pop();
		}
	return newNd;
	}
*/
/*----------------------------------------------------------------------------------------------------------------------
|	Had been CreateTreeNode.  Name changed to imply that it is not a full copy of the node (the node copy constructor)
|	only copies the number and edgelen (and resets split field for tips)
*/
TreeNode* Tree::CreateSimilarTreeNode(const TreeNode &n)
	{
	return CreateTreeNode(n.GetNodeNumber(), n.GetFltEdgeLen(), n.IsLeafOrRoot());
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Instantiates a new node (or returns a currently unused on).
*/
TreeNode *Tree::CreateTreeNode(unsigned ndN, double edgeLen, bool setSplitByNum)
	{
	TreeNode *newNd;
	if (nodeStorage.empty())
		return new TreeNode(ndN, edgeLen, setSplitByNum);
	newNd = nodeStorage.top();
	newNd->Clear();
	nodeStorage.pop();
	newNd->SetNodeNumber(ndN, setSplitByNum);
	if (!setSplitByNum)
		newNd->split.Clear();
	newNd->SetFltEdgeLen(edgeLen);
	return newNd;
	}

Tree::Tree(const TreeID &inID)
 :firstPreorder(NULL), lastPreorder(NULL)
	{
	Clear();
	TreeManip treeMan(this);
	treeMan.SimpleBuildTreeFromID(Split::GetNTaxa(), inID);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Used in CopyTreeStructure to make a copy of a node that exists in another tree.
*/
void Tree::StoreTreeNode(TreeNode *u)
	{
	nodeStorage.push(u);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor to be used when there is no token available for specifying file position information.
*/
XBadTreeDef::XBadTreeDef(
  const char		*m) /* description of the failure */
  : msg(m) 
	{
	endpos = UINT_MAX;
	startpos = UINT_MAX;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Constructor requires a c_str message describing cause of exception and a reference to the NxsToken object that
|	holds the information about the point in the tree description string where the problem occurred.
*/
XBadTreeDef::XBadTreeDef(
  const char		*m, /* description of the failure */
  const NxsToken	&t)	/* reference to the Token that just caused the exception */
  : msg(m) 
	{
	endpos = (unsigned) t.GetFilePosition();
	startpos = endpos - t.GetTokenLength();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor to be used when there is no token available for specifying file position information.
*/
XBadTreeStructure::XBadTreeStructure(
  const char	*m,				/* description of the failure */
  unsigned		which_nodenum,	/* node being visited when problem was discovered */
  bool			is_leaf)		/* true if node being visited when problem was discovered was a leaf node, false otherwise */
	{
	msg		= m;
	nodeNum	= which_nodenum;
	isLeaf	= is_leaf;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor to be used when there is no token available for specifying file position information.
*/
XBadTreeStructure::XBadTreeStructure(
  const char	*m,				/* description of the failure */
  const TreeNode *n)		
	{
	msg		= m;
	nodeNum	= n->GetNodeNumber();
	isLeaf	= n->IsShootTip();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Initialize to initialize data members.
*/
Tree::Tree(): firstPreorder(NULL), lastPreorder(NULL)
	{
	Clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Standard CopyConstructor
*/
Tree::Tree(const Tree &r): firstPreorder(NULL), lastPreorder(NULL)
	{
	Clear();
	CopyTreeStructure(r);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Identical to building a tree with the default constructor and calling Tree::BuildTreeFromDescription
*/
Tree::Tree(const FullTreeDescription &d): firstPreorder(NULL), lastPreorder(NULL)
	{
	Clear();
	BuildTreeFromDescription(d);
	}
							
		
/*----------------------------------------------------------------------------------------------------------------------
|	Causes `this' tree to have the same topology as the argument
*/
void Tree::CopyTreeStructure(const Tree &r)
	{
	//@POL this does not yet copy the nodes in nodeStorage - should nodes in storage be copied?
	Clear();
	r.DebugCheckTreeStructure();
	nInternals = r.nInternals;
	nLeaves = r.nLeaves;
	preordersDirty = true;
	InvalidateID();
	hasEdgelens = r.hasEdgelens;
	max_x = r.max_x;
	max_y = r.max_y;
	isRooted = r.isRooted;
	const TreeNode *currNdToCopy = r.GetRoot();
	if (currNdToCopy != NULL)
		{
		TreeNode *nd = firstPreorder = CreateSimilarTreeNode(*currNdToCopy);
		bool movingDown = false;
		while (currNdToCopy != NULL)
			{
			if (!movingDown && currNdToCopy->GetLeftChild() != NULL)
				{
				currNdToCopy = currNdToCopy->GetLeftChild();
				nd->AddLeftChild( CreateSimilarTreeNode(*currNdToCopy));
				nd = nd->GetLeftChild();
				}
			else if (currNdToCopy->GetRightSib() != NULL)
				{
				movingDown = false;
				currNdToCopy = currNdToCopy->GetRightSib();
				TreeNode *temp = CreateSimilarTreeNode(*currNdToCopy);
				nd->AddRightSib(temp);
				nd = temp;
				}
			else 
				{
				movingDown = true;
				currNdToCopy = currNdToCopy->GetParent();
				nd = nd->GetParent();
				NXS_ASSERT((nd != NULL && currNdToCopy != NULL) ||(nd == NULL && currNdToCopy == NULL));
				}
			}
		lastPreorder = firstPreorder->FindLastPreorderInClade();
		}
	DebugCheckTreeStructure();
	}
		
/*----------------------------------------------------------------------------------------------------------------------
|	Calls Empty to delete all allocated memory.
*/
Tree::~Tree()
	{
	Clear();
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Deletes all of the nodes in the tree
|	calls recursive DeleteSelfAndDescendants
*/
void Tree::DeleteAllNodes()
	{
	TreeNode::DeleteSelfAndDescendants(firstPreorder);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes the node (added as a function so that the tree can recycle nodes)
*/
void Tree::DeleteSingleNode(TreeNode *nd)
	{
	delete nd; 
	}

	
/*----------------------------------------------------------------------------------------------------------------------
|	Called by the destructor. Deletes root (which deletes all other nodes recursively) and deletes workSpace c_string.
*/
void Tree::Clear()
	{
	DeleteAllNodes();
	firstPreorder		= NULL;
	lastPreorder		= NULL;
	nInternals			= 0;
	nLeaves				= 0;
	preordersDirty		= true;
	hasEdgelens			= false;
	max_x				= 0;
	max_y				= 0;
	isRooted			= false;

	InvalidateNodeCounts();
	InvalidateID();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	SHOULDN'T be called to create valid Trees!!! Use BuildTreeFromString for that!
|	Reads the tree, but does NOT reroot it at the lowest numbered node.
|	First calls FlushTree to eliminate any existing nodes, then builds the tree specified in s. If any problems 
|	arise in this process, Clear() is again called to clean things up and a XBadTreeDef exception is thrown. 
|	Does not translate taxon numbers to taxon names, and does not convert underscores found in taxon names to spaces. 
|	Before returning, calls TraverseTree to establish the nextPreorder links and RefreshID to build the TreeID data
|	member `id'.
|	Returns a pointer to the terminal with the lowest index (this node will root the tree if the tree is unrooted).
|	sets the following fields of the nodes (only):
|		par, rSib, lChild, nodeNum, edgeLen (name is set if USING_NODE_NUMBERS_ONLY is not defined)
|	
|	split, pre/post order pointers are NOT initialized.
*/
TreeNode *Tree::BuildUnfinishedTreeToCheckString(
  const char *s,					/* the string containing a valid Newick-style tree description (terminating semicolon Ok but not required) */
  bool translateNames,				/* if true, the names should be translated from taxon names */
  set<unsigned> *termNodeNumbers)	/* */	
	{
	assert(s != NULL);
	// Eliminate any existing nodes
	//
	Clear();

	unsigned lowestTipIndex = UINT_MAX;
	TreeNode *tipWithLowestIndex = NULL;
	try
		{
		unsigned nEdgeLengths = 0;
		bool expectingEdgeInfo = false;
		TreeNode *nd = CreateTreeNode(UINT_MAX, DBL_MAX, false);
		firstPreorder = nd;	//note we need to point the firstPreorder
		NxsToken token(s);
		enum {Prev_Tok_LParens = 0, Prev_Tok_RParens, Prev_Tok_Colon, Prev_Tok_Comma, Prev_Tok_Name, Prev_Tok_EdgeLen};
		unsigned previous = Prev_Tok_LParens;
		for (;;)
			{
			if (previous == Prev_Tok_Colon)
				{
				double tmp;
				token.ReadDoubleToken(&tmp);
				}
			else
				++token;
			if (token.AtEOF() && token.GetTokenReference().empty())
				break;
			char ch = token.GetTokenReference()[0];
			switch(ch)
				{
				case ';':
					assert(token.GetTokenLength() == 1);
					throw XBadTreeDef("premature semicolon", token);
					break;  

				case ')':
					assert(token.GetTokenLength() == 1);
					if (nd->IsRoot())
						throw XBadTreeDef("too many right parentheses", token);
					if (previous != Prev_Tok_RParens && previous != Prev_Tok_Name && previous != Prev_Tok_EdgeLen)
						throw XBadTreeDef("unexpected right parentheses", token);
					if (!nd->IsShootTip() && nd->GetNodeNumber() == UINT_MAX)
						nd->nodeNum = GetInternalNodeNum(string());
					nd = nd->GetParent(); // go down a level
					if (nd->GetLeftChild()->GetRightSib() == NULL)
						throw XBadTreeDef("internal node has only one child", token);
					previous = Prev_Tok_RParens;
					break;

				case ':':
					assert(token.GetTokenLength() == 1);
					if (previous != Prev_Tok_RParens && previous != Prev_Tok_Name && previous != Prev_Tok_EdgeLen)
						throw XBadTreeDef("unexpected colon", token);
					previous = Prev_Tok_Colon;
					assert(nd->GetFltEdgeLen() == DBL_MAX);
					//if(nd->GetFltEdgeLen() != DBL_MAX)
					// 	throw XBadTreeDef("multiple branch lengths specified for one node", token);
					break;

				case ',':
					assert(token.GetTokenLength() == 1);
					// Create sibling for nd
					//
					if (nd->IsRoot() || previous == Prev_Tok_LParens || previous == Prev_Tok_Comma || previous == Prev_Tok_Colon)
						throw XBadTreeDef("unexpected comma", token);

					if (!nd->IsShootTip() && nd->GetNodeNumber() == UINT_MAX)
						nd->nodeNum = GetInternalNodeNum(string());

					nd->rSib = CreateTreeNode(UINT_MAX, DBL_MAX, false);
					nd->rSib->par = nd->par;
					nd = nd->rSib;
					previous = Prev_Tok_Comma;

					break;

				case '(':
					assert(token.GetTokenLength() == 1);
					if (previous != Prev_Tok_LParens && previous != Prev_Tok_Comma)
						throw XBadTreeDef("too many left parentheses - not after a comma or left parentheses ", token);
					// Create new_node above and to the left of the current node
					//
					assert(nd->GetLeftChild() == NULL);
					nd->lChild = CreateTreeNode(UINT_MAX, DBL_MAX, false);
					nd->lChild->par = nd;
					nd = nd->lChild;
					previous = Prev_Tok_LParens;
					break;

				default:
					// Read taxon name or edge info
					//
					if (previous == Prev_Tok_Colon)
						{
						nd->SetInfo(token.GetTokenReference());
						++nEdgeLengths;
						previous = Prev_Tok_EdgeLen;
						}
					else
						{
						if (previous != Prev_Tok_RParens && previous != Prev_Tok_Comma && previous != Prev_Tok_LParens)
							{
							string s;
							s << "unexpected taxon name \"" << token.GetTokenReference() << '\"';
							throw XBadTreeDef(s.c_str(), token);
							}
						if (nd->IsShootTip())
							{
							nd->nodeNum = GetLeafNum(token.GetTokenReference(), translateNames);
							if (nd->nodeNum == UINT_MAX)
								{
								string s;
								s << token.GetTokenReference() << " is not a valid taxon name";
								throw XBadTreeDef(s.c_str(), token);	
								}
							if (termNodeNumbers->find(nd->nodeNum) != termNodeNumbers->end())
								{
								string s;
								s << "the taxon " << token.GetTokenReference() << " appears more than once in the tree description";
								throw XBadTreeDef(s.c_str(), token);
								}
							if (nd->nodeNum < lowestTipIndex)
								{
								tipWithLowestIndex = nd;
								lowestTipIndex = nd->nodeNum;
								}
							termNodeNumbers->insert(nd->nodeNum);
							}	
						else
							nd->nodeNum = GetInternalNodeNum(token.GetTokenReference());
						if (nd->nodeNum == UINT_MAX)
							throw XBadTreeDef("taxon unknown", token);
#						if !defined(USING_NODE_NUMBERS_ONLY)											
							nd->name = token.GetTokenReference();
#						endif
						previous = Prev_Tok_Name;
						}
					expectingEdgeInfo = false;
				}
			}

		if (firstPreorder->nodeNum == UINT_MAX)
			firstPreorder->nodeNum = GetInternalNodeNum(string());
		//	currently we're just throwing away any branch given to the root
		//
		//if (root->GetFltEdgeLen() != DBL_MAX)
		//	++nEdgeLengths;
		firstPreorder->SetFltEdgeLen(0.0);
		if (nEdgeLengths > 0)
			{
			if (nEdgeLengths < nLeaves + nInternals - 1)
				throw XBadTreeDef("some but not all edge lengths were specified");
			hasEdgelens = true;
			}
		else
			hasEdgelens = false;
		if (nd != GetRoot())
			throw XBadTreeDef("too many left parentheses");
		}
		
	catch(...)
		{
		Clear(); //@@@ we need to verify that this deletion of nodes on a bad tree isn't going to do bad things
		throw;
		}
	NXS_ASSERT(tipWithLowestIndex != NULL);
	return tipWithLowestIndex;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Makes nd's parent nd's rightmost child.  (nd->par pointer still points to this node, too).
|	Called by Reroot().
|	the Node's Parent's par pointer still points to the original "grandparent" node
|	Caller still has to fix nd->par, nd->prevPreorder
*/
void Tree::FlipNodeHelper(TreeNode *nd)
	{
	NXS_ASSERT(nd != NULL);
	NXS_ASSERT(nd->par != NULL);
	TreeNode *prevRightChild = nd->FindRightChild();
	if (prevRightChild == NULL)
		nd->ReplaceLeftChildHelper(nd->par);
	else
		prevRightChild->ReplaceRightSibHelper(nd->par);
	//	Now we erase ourselves from pointer chain of the previous parent's children
	//
	TreeNode *prevParChild = nd->par->lChild;
	if (prevParChild == nd)
		{
		NXS_ASSERT(nd->rSib != NULL || nd->par->par == NULL);	//we shouldn't be our old parent's left child unless we have a rSib or our parent is the root
		nd->par->ReplaceLeftChildHelper(nd->rSib);
		}
	else
		{
		while (prevParChild->rSib != nd)
			{
			NXS_ASSERT(prevParChild != NULL);	//if we trip this, nd, wasn't in the list of it's parent's children
			prevParChild = prevParChild->rSib;
			}	
		prevParChild->rSib = NULL;	//we set this to NULL, so that ReplaceRightSibHelper doesn't try to use nd's prevPreorder
		prevParChild->ReplaceRightSibHelper(nd->rSib);
		}
	nd->rSib = NULL;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Sets the tree's firstpreorder but not the nodes' 
|	splits are not initialized.
*/
void Tree::Reroot(
  TreeNode *nd) /*pointer to a tip that is attached to the tree */
	{
	if (nd->IsRoot())
		return;
	assert(nd->lChild == NULL && nd->par != NULL);	//assumes that we are rooting at a tip that had been attached 
	firstPreorder = nd;
	const bool treeHasEdgeLens = hasEdgelens;
	
	Length newEdgeLen;
	Length nextEdgeLen;
	newEdgeLen.f = DBL_MAX; //the new root get's this Length, //@@@ problem with union.  What do we do here? (we don't know whether or not the union is double )
	nextEdgeLen.f = DBL_MAX; //the new root get's this Length, 
	if (treeHasEdgeLens)
		nextEdgeLen = nd->GetEdgeLen();
	nd->prevPreorder = NULL;
	TreeNode *newParent = NULL;
	TreeNode *nextNodeToFlip = nd->par;
	while (nd->par->par != NULL)
		{
		FlipNodeHelper(nd);
		nd->par = newParent;
		if (treeHasEdgeLens)
			{
			nd->SetEdgeLen(newEdgeLen);
			newEdgeLen = nextEdgeLen;
			nextEdgeLen = nextNodeToFlip->GetEdgeLen();
			}
		newParent = nd;
		nd = nextNodeToFlip;
		nextNodeToFlip = nd->par;
		}
	//	we are near the root of tree.  If the previous root had been bifurcating (this
	//	could be true if we just read the tree in from a file) we need to delete
	//	the nextNodeToFlip, not flip it.
	//
	const unsigned nChildrenOrigRoot = nd->par->CountChildren();
	if (nChildrenOrigRoot != 2)
		{ 	
		assert(nChildrenOrigRoot != 0);
		FlipNodeHelper(nd);
		nd->par = newParent;
		if (treeHasEdgeLens)
			nd->SetEdgeLen(newEdgeLen);
		nextNodeToFlip->par = nd;	
		if (treeHasEdgeLens)
			nextNodeToFlip->SetEdgeLen(nextEdgeLen);
		}
	else
		{
		// had been rooted with a node of degree 2 - we need to delete the previous root
		// add our sibling as our rightmost (or only) child.
		TreeNode *formerSibling = (nd->rSib != NULL ? nd->rSib : nd->par->lChild);
		assert(formerSibling != NULL);
		if (formerSibling->rSib == nd)
			{
			NXS_ASSERT(nd->rSib == NULL);
			formerSibling->rSib = NULL;
			}
		else
			{
			NXS_ASSERT(nd->rSib == formerSibling && formerSibling->rSib == NULL);
			nd->rSib = NULL;
			}
		NXS_ASSERT(formerSibling != newParent);
		TreeNode *rightmostChild = nd->FindRightChild();
		if (rightmostChild != NULL)
			rightmostChild->ReplaceRightSibHelper(formerSibling);
		else
			nd->ReplaceLeftChildHelper(formerSibling);
		NXS_ASSERT(rightmostChild != NULL || newParent == NULL);
		formerSibling->par = nd;
		if (treeHasEdgeLens)
			{
			//the formerSibling gets all of the lenght of the basal branches
			nd->SetEdgeLen(newEdgeLen);
			newEdgeLen = formerSibling->GetEdgeLen();
			formerSibling->SetFltEdgeLen(newEdgeLen.f + nextEdgeLen.f);	//@@@ problem with using union.  What do we do here? (we don't know whether or not the union is double )
			}
		nd->par = newParent;
		--nInternals;
		StoreTreeNode(nextNodeToFlip);	
		}
	NXS_ASSERT(firstPreorder->GetLeftChild() != NULL && firstPreorder->GetLeftChild()->GetRightSib() == NULL);
	lastPreorder = (firstPreorder->GetLeftChild() != NULL ? firstPreorder->GetLeftChild()->FindLastPreorderInClade() : firstPreorder);
	lastPreorder->nextPreorder = NULL;
	}
	
#if defined(NCL_USE_NXS_CONSOLE_OUTPUT)
	void Tree::DebugRerootingTest(TreeNode *nd, NxsStdOutputStream &o)
		{
		Reroot(nd);
		DebugCheckTreeStructure();
		o << *this << ncl::endl;
		}
#endif

void Tree::BuildTreeFromDescription(const FullTreeDescription &d)
	{
	BuildTreeFromString(d.newick.c_str(), false, false );//@@@ should d.rooted as 3rd arg when rooted trees are working
	}

void Tree::BuildTreeFromString(const char *s, bool translateNames, bool readAsRooted)
	{
	set<unsigned> observerdTerminalIndices;
	TreeNode *rootToBe = BuildUnfinishedTreeToCheckString(s, translateNames, &observerdTerminalIndices);
	TraverseTree();
	isRooted = readAsRooted;
	if (!readAsRooted)
		Reroot(rootToBe);
	assert(!isRooted);
	RefreshFullID();
	DebugCheckTreeStructure();
	}

#if defined(USING_NODE_NUMBERS_ONLY)		
	/*----------------------------------------------------------------------------------------------------------------------
	|	Increments nLeaves and returns either: 1) the index of the taxon named in `nm' (if translateName true) or 2) the
	|	unsigned representation of `nm' (if translateName false). Used by BuildTreeFromString to number nodes. 
	|	Override in derived classes to set leaf numbers according to index into an array of valid leaf names.
	*/
	unsigned Tree::GetLeafNum(
		const std::string &nm,/* name of the leaf node for which a node number is needed */
	  bool translateName)	
		{
		unsigned node_number;
		++nLeaves;
		if (translateName && TreeNode::gTaxaMgr)
			{
			node_number = TreeNode::gTaxaMgr->FindIndex(nm);
			//@POL 19-Oct-2004 I suppose this dependency on taxaMgr is unavoidable, but perhaps 
			//one way around it is to pass in a pointer to a generic name-lookup object that
			//could be taxaMgr but the pointer could also be NULL, in which case this function
			//would act as if translateName == false
			//
			}
		else
			{
			bool success = IsAnUnsigned(nm, &node_number); 
			if ( !success ) node_number = UINT_MAX;
			}

		//POL 18-Oct-2004 
		if (node_number >= Split::splitNTax)
			{
			if (translateName)
				throw XBadTreeDef("taxon label in tree description not a valid taxon");
			else
				throw XBadTreeDef("taxon label in tree description too large (first taxon should be 0, not 1)");
			}

		return node_number;
		}
#else
	/*----------------------------------------------------------------------------------------------------------------------
	|	Returns nLeaves and then increments nLeaves. Used by BuildTreeFromString to number nodes. Virtual function that can
	|	be overridden in derived classes to set leaf numbers according to index into an array of valid leaf names.
	*/
	unsigned Tree::GetLeafNum(
		const std::string &nm,/* name of the leaf node for which a node number is needed */
	  bool )	
		{
	#	if defined(HAVE_PRAGMA_UNUSED)
	#		pragma unused(nm)
	#	endif
		return nLeaves++;
		}
#endif
/*----------------------------------------------------------------------------------------------------------------------
|	Returns nInternals and then increments nInternals. Used by BuildTreeFromString to number nodes. Virtual function
|	that can be overridden in derived classes number internal nodes according to some other scheme.
*/
unsigned Tree::GetInternalNodeNum(
  const string & nm)	/* name of the internal node for which a node number is needed */
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(nm)
#	endif

	return nInternals++;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides a crude representation of the tree structure for debugging purposes.
*/
void Tree::List(
  int nindent,	/* number of spaces to indent */
  TNodePtr p)	/* pointer to starting node (default is NULL, which means start at root node) */
	{
	TNodePtr q;

	if (p == NULL)
		p = GetRoot();

#define SHOW_XY
#	if defined(SHOW_XY)
	if (p->IsShootTip())
			printf("%*cNode %d \"%s\" (a=%d l=%d r=%d): x=%f y=%f\n", nindent, ' ', p->GetNodeNumber(), p->GetName().c_str(), Dex(p->GetParent()),
			  Dex(p->GetLeftChild()), Dex(p->GetRightSib()), p->GetX(), p->GetY());
		else
			printf("%*cNode %d (a=%d l=%d r=%d): x=%f y=%f\n", nindent, ' ', p->GetNodeNumber(), Dex(p->GetParent()), Dex(p->GetLeftChild()),
			  Dex(p->GetRightSib()), p->GetX(), p->GetY());
#	elif defined(VERBOSE)
	if (p->IsShootTip())
		printf("%*cNode %d \"%s\" (a=%d l=%d r=%d nxtPre=%d prvPre=%d)\n", nindent, ' ', Dex(p),
		  OrigTaxLabel(p->nodeNum), Dex(p->GetParent()), Dex(p->GetLeftChild()), Dex(p->GetRightSib()), Dex(p->GetNextPreorder()),
		  Dex(p->GetNextPostorder()));
	else
		printf("%*cNode %d (a=%d l=%d r=%d nxtPre=%d prvPre=%d)\n", nindent, ' ', Dex(p), Dex(p->GetParent()),
		  Dex(p->GetLeftChild()), Dex(p->GetRightSib()), Dex(p->GetNextPreorder()), Dex(p->GetNextPostorder()));
#	elif defined(TERSE)
	if (p->IsShootTip())
		printf("%*c%d\n", nindent, ' ', p->nodeNum);
	else
		printf("%*c.\n", nindent, ' ');
#	else
        if (p->IsShootTip())
			printf("%*cNode %d \"%s\" (a=%d l=%d r=%d)\n", nindent, ' ', p->nodeNum, p->name.c_str(), Dex(p->GetParent()),
			  Dex(p->GetLeftChild()), Dex(p->GetRightSib()));
		else
			printf("%*cNode %d (a=%d l=%d r=%d)\n", nindent, ' ', p->nodeNum, Dex(p->GetParent()), Dex(p->GetLeftChild()),
			  Dex(p->GetRightSib()));
#	endif
	//if (AbortRequested())
	//	return;
	q = p->GetLeftChild();
	while (q != NULL)
		{
		List(nindent + 1, q);
		//if (AbortRequested())
		//	return;
		q = q->GetRightSib();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If supplied node pointer is NULL, returns -1; otherwise returns node number.
*/
int Tree::Dex(
  TNodePtr p)	/* pointer to the node in question */
	{
	return (p == NULL) ? -1 : (int) p->GetNodeNumber();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recursively traverses the tree, starting at nd, checking all five navigational pointers (lChild, rSib, par, 
|	nextPreorder and prevPreorder) for each node. Any inconsistencies discovered result in an XBadTreeStructure 
|	exception being thrown. Should ordinarily be called without arguments so that the traversal begins at the root node.
*/
void Tree::DebugCheckTreeStructure(
  bool allow_polytomies,		/* if true, polytomies allowed; if false, every node should be of order 3 */
  const TreeNode *nd) const		/* the current node (NULL means start at root node) */
	{
	STATIC_DATA unsigned numinternals;
	STATIC_DATA unsigned numlvs;
	STATIC_DATA const TreeNode *last_node_visited;
	if (nd == NULL)
		{
		nd = GetRoot();
		numinternals = 0;
		numlvs = 0;
		last_node_visited = NULL;

		// Root node should be a tip node, have no par, no rSib, and only one child, which is an internal node
		//
		if (nd->GetParent() != NULL)
			throw XBadTreeStructure("root node has a parent", nd);
		if (nd->GetRightSib() != NULL)
			throw XBadTreeStructure("root node has a sibling", nd);
		if (nd->GetLeftChild() == NULL)
			throw XBadTreeStructure("root node has no children", nd);
		if (nd->GetLeftChild()->GetRightSib() != NULL)
			throw XBadTreeStructure("root node has more than one child", nd);
		if (nd->GetLeftChild()->IsShootTip())
			throw XBadTreeStructure("single child of root node should be an internal node", nd);

		// Trees's lastPreorder pointer should point to the last node in root's clade (in the preorder sequence)
		//
		const TreeNode *last = nd->FindLastPreorderInClade();
		if (lastPreorder != last)
			throw XBadTreeStructure("lastPreorder does not point to last preorder node in root's clade", nd);
		}

	// Update numlvs and numinternals
	//
	if (nd->IsRoot() || nd->IsShootTip())
		++numlvs;
	if (!nd->IsShootTip() && !nd->IsRoot())
		++numinternals;

	double edgeLength = (nd->IsRoot() ? 0.0 : nd->GetFltEdgeLen());
#	if !defined(ALLOW_NEGATIVE_EDGELENS)
		// Floating-point edge lengths should never be negative
		//
		if (edgeLength < 0.0)
			{
			const string msg = MakeStrPrintF("Edge length for node number %d is negative (%.6f)", nd->GetNodeNumber(), nd->GetFltEdgeLen());
			throw XBadTreeStructure(msg.c_str(), nd);
			}
#	endif

	if (edgeLength == DBL_MAX && HasEdgeLens())
		throw XBadTreeStructure("At least one floating point edge length has not been initialized", nd);

	
	// The prevPreorder pointer should match last_node_visited because we are traversing the tree in preorder fashion
	//
	if (nd->GetNextPostorder() != last_node_visited)
		throw XBadTreeStructure("prevPreorder pointer is incorrect", nd);

	// If not at root, last_node_visited should be non-NULL
	//
	if (!nd->IsRoot() && last_node_visited == NULL)
		throw XBadTreeStructure("last_node_visited is NULL but not at root", nd);

	// If not at root, the nextPreorder pointer of last_node_visited should be equal to nd
	//
	if (last_node_visited != NULL && last_node_visited->GetNextPreorder() != nd)
		throw XBadTreeStructure("nextPreorder pointer is incorrect", last_node_visited);

	// If nd is NOT a leaf, then lChild should be equal to nextPreorder for nd
	//
	const TreeNode *child = nd->GetLeftChild();
	if (!nd->IsShootTip() && child != nd->GetNextPreorder())
		throw XBadTreeStructure("lChild not same as nextPreorder", nd);

	// If nd IS a leaf, and if rSib is not NULL, then rSib should be equal to nextPreorder for nd
	//
	if (nd->IsShootTip() && nd->GetRightSib() != NULL && nd->GetRightSib() != nd->GetNextPreorder())
		throw XBadTreeStructure("rSib exists but is not same as nextPreorder for leaf node", nd);

	// If nd is not a leaf, all of its children should name nd as their parent
	//
	for (child = nd->GetLeftChild(); child != NULL; child = child->GetRightSib())
		{
		if (child->GetParent() != nd)
			throw XBadTreeStructure("not all children name this node as their parent", nd);
		}

	// If left sib exists, prevPreorder should point to the last node in the left sib's clade; 
	// otherwise, prevPreorder should be same as par
	//
	const TreeNode *lSib = nd->FindLeftSib();
	if (lSib != NULL && nd->GetNextPostorder() != lSib->FindLastPreorderInClade())
		throw XBadTreeStructure("prevPreorder does not point to last node in left sibling's clade", nd);
	if (lSib == NULL && nd->GetNextPostorder() != nd->GetParent())
		throw XBadTreeStructure("prevPreorder should equal par if there is no sibling to the left", nd);

	// nd->par->rSib should not equal nd
	//
	if (nd->GetParent() != NULL && nd->GetParent()->GetRightSib() == nd)
		throw XBadTreeStructure("nd->par->rSib should not equal nd", nd);

	// nd->par should not equal nd->nextPreorder
	//
	if (nd->GetParent() == nd->nextPreorder)
		throw XBadTreeStructure("nd->par should not equal nd->nextPreorder", nd);

	// nd->nextPreorder should not equal nd (and the same goes for nd->prevPreorder)
	//
	if (nd->GetNextPreorder() == nd)
		throw XBadTreeStructure("nd->nextPreorder should not equal nd", nd);
	if (nd->GetNextPostorder() == nd)
		throw XBadTreeStructure("nd->prevPreorder should not equal nd", nd);

	// check for polytomies if allow_polytomies is false
	//
	bool unallowed_polytomy = false;
	if (!nd->IsRoot() && !allow_polytomies)
		{
		const TreeNode *lc = nd->GetLeftChild();
		if (lc != NULL)
			{
			const TreeNode *rc = lc->GetRightSib();
			if (rc->GetRightSib() != NULL)
				unallowed_polytomy = true;
			}
		}

	if (unallowed_polytomy)
		{
		throw XBadTreeStructure("non-root node has more than 2 children (and polytomies are not allowed)", nd);
		}

	last_node_visited = nd;

	if (nd->GetLeftChild() != NULL)
		DebugCheckTreeStructure(allow_polytomies, nd->GetLeftChild());
	if (nd->GetRightSib() != NULL)
		DebugCheckTreeStructure(allow_polytomies, nd->GetRightSib());

	if (nd->IsRoot())
		{
		if (numlvs != GetNLeaves())
			{
			string msg = MakeStrPrintF("nLeaves is set to %d but there are really %d leaves in the tree", nLeaves, numlvs);
			throw XBadTreeStructure(msg.c_str(), nd);
			}

		if (numinternals != GetNInternals())
			{
			string msg = MakeStrPrintF("nInternals is set to %d but there are really %d internal nodes in the tree", nInternals, numinternals);
			throw XBadTreeStructure(msg.c_str(), nd);
			}
		}
	}

#if !ALT_TRAVERSE_TREE
	/*----------------------------------------------------------------------------------------------------------------------
	|	Recursively traverses the tree, starting at nd, setting the nextPreorder and prevPreorder traversal pointers. Also
	|	sets data members `nLeaves' and `nInternals'. Should ordinarily be called without arguments so that the traversal 
	|	begins at the root node.
	*/
	void Tree::TraverseTree(
		TreeNode *nd)				/* the current node (NULL means start at root node) */
		{
		STATIC_DATA TreeNode *preheld; // used to keep track of node whose nextPreorder pointer has not yet been assigned

		//
		// X   Y   V  U  Preorder:  W, Z, X, Y, V, U
		// |   |   |  |  Postorder: U, V, Y, X, Z, W (exactly the reverse of preorder sequence)
		// +-Z-+   |  |  
		//   |     |  |
		//   +-----W--+  
		//               

		if (nd == NULL)
			{
			nd						= GetRoot();
			lastPreorder			= NULL;
			preheld					= NULL;
			nLeaves					= 0;
			nInternals				= 0;
			}

		if (preheld != NULL && nd->GetParent() != NULL && nd->GetParent()->GetNextPreorder() != nd)
			{
			// preheld will be set to a node if that node is at the upper-right
			// corner of its lineage, because at that point it is not obvious
			// who its next preorder neighbor is going to be. We need to back
			// out of some recursions and go right once in order to find this 
			// next preorder node. Thus, preheld is used to hold our place 
			// until that time. If preheld is non-NULL, and if the current 
			// node is not the nextPreorder node for its ancestor, then preheld's 
			// nextPreorder pointer should be assigned to this node. 
			// After the assignment is made, preheld can be 
			// set back to NULL to indicate that no nodes are currently 
			// in limbo with respect to their nextPreorder pointers.
			//
			preheld->nextPreorder = nd;
			nd->prevPreorder = preheld;
			preheld = NULL;
			lastPreorder = NULL; // preheld was not the last node in the pre-order sequence after all
			}

		// Update nLeaves and nInternals
		//
		if (nd->IsRoot() && (nd->CountChildren() == 1))
			++nLeaves;
		else if (nd->IsShootTip())
			++nLeaves;
		else
			++nInternals;

		if (nd->GetLeftChild() != NULL)
			{
			nd->nextPreorder = nd->GetLeftChild();
			nd->GetLeftChild()->prevPreorder = nd;
			}
		else if (nd->GetRightSib() != NULL)
			{
			nd->nextPreorder = nd->GetRightSib();
			nd->GetRightSib()->prevPreorder = nd;
			}
		else
			{
			// We've reached the upper right corner of this clade
			//
			assert(preheld == NULL);
			preheld = nd;
			nd->nextPreorder = NULL;
			lastPreorder = preheld; // preheld might be last node in pre-order sequence
			}

		if (nd->GetLeftChild() != NULL)
			TraverseTree(nd->GetLeftChild());

		if (nd->GetRightSib() != NULL)
			TraverseTree(nd->GetRightSib());

		if (nd == GetRoot())
			{
			preordersDirty = false;
			node_counts_valid = true;
			}
		}
#else
	//@POL are we using this?
	void Tree::TraverseTree(
		TreeNode *nd)				/* the current node (NULL means start at root node) */
		{
		if (nd == NULL)
			{
			nd = GetRoot();
			if (nd == NULL)
				return;
			nd->SetPreorderPointers(NULL);
			preordersDirty = false;
			}
		else
			nd->SetPreorderPointers(nd->GetNextPostorder());
		}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Loops through tree in preorder fashion, showing for each node the node's name, whether it is a leaf or an internal
|	node, the node's number, and the node's floating-point edge length.
*/
void Tree::DebugShowPreorder(ostream & out)
	{
	if (preordersDirty)
		TraverseTree();
	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		out << nd->GetDebugDescription() << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Loops through tree in postorder fashion, showing for each node the node's name, whether it is a leaf or an internal
|	node, the node's number, and the node's floating-point edge length.
*/
void Tree::DebugShowPostorder(ostream & out)
	{
	if (preordersDirty)
		TraverseTree();
	for (TreeNode *nd = GetFirstPostorder(); nd != NULL; nd = nd->GetNextPostorder())
		out << nd->GetDebugDescription() << std::endl;
	}

std::string & AppendColonAndEdgeLen(std::string & s, const DblFormatter & formatter, bool fltlens, const TreeNode &nd);

inline std::string & AppendColonAndEdgeLen(std::string & s,  const DblFormatter & formatter, bool fltlens, const TreeNode &nd)
	{
#if defined(POL_MESQUITE_COLORING)
	s << ':';
	if (fltlens)
		formatter.FormatDouble(s, nd.GetFltEdgeLen());
	else
		s << nd.GetIntEdgeLen();

	// Mesquite colors: 
	//	5	= red
	//	11	= green
	//	14	= blue
	//	7	= yellow
	if (nd.IsSelected())
		s << "<color=5>";
#else
	s << ':';
	if (fltlens)
		return formatter.FormatDouble(s, nd.GetFltEdgeLen());
	s << nd.GetIntEdgeLen();
#endif
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Stores tree structure and edge lengths in the form of a Newick-style tree description and saves in the supplied
|	string `s'.
*/
std::string & Tree::AppendNewickRepresentation(
  std::string & s,			/**< is the NxsString object in which the representation will be saved */
  bool useTaxonNames,		/**< if true, inserts taxon labels stored in TreeNode::taxaMgr; if false, inserts node numbers (indices of taxa in taxaMgr) */
  bool edgelens,			/**< determines whether edge lengths will be saved (true) or omitted (false) from tree description */
  bool fltlens,				/**< determines whether floating point edge lengths will be saved (true) or the integer edge lengths will be saved (false), contingent upon `edgelens' being true */
  unsigned fltprec) const	/**< determines number of decimal places used to represent floating point edge lengths, contingent upon both `edgelens' and 'fltlens' being true */
	{
	//@POL at the moment, designed only to work only with unrooted trees rooted at one of the leaves
	DblFormatter df(UINT_MAX, fltprec);
	s << '(';
	unsigned openParens = 1;
	assert(TreeNode::gTaxaMgr);
	PhoTaxaManager & taxMgr = *TreeNode::gTaxaMgr;
	const TreeNode * nd = GetFirstPreorder();
	
	assert(nd->IsRoot());
	// MTH 13-Dec-2004 moving root handling out of the preorder loop, fixes bug that assigned the branch length to the root incorrectly.
	//POL 19-Oct-2004
	AppendTaxonLabel(s, nd->GetNodeNumber(), taxMgr, useTaxonNames);
	nd = nd->GetNextPreorder();
	if (nd == NULL) //one taxon tree, presumably this won't happen, but it is easy to deal with
		return s << ')'; 
	if (edgelens) // we give the root's child's edge length to the root because we print as a basal polytomy (not rooted at a leaf)
		AppendColonAndEdgeLen(s, df, fltlens, *nd);
	s << ',';
	if (nd->GetNextPreorder() == NULL)//two taxon tree, presumably this won't happen, but it is easy to deal with
		{ 
		AppendTaxonLabel(s, nd->GetNodeNumber(), taxMgr, useTaxonNames);
		if (edgelens)  // 0 branch length because we gave the root's child's edge length to the root because we print as a basal polytomy (not rooted at a leaf)
			return s << ":0)";
		else
			return s << ")";
		}
	nd = nd->GetNextPreorder(); // note: we don't add an open parentheses here, this is how we get the polytomy at the root
	assert(nd->GetNextPreorder()); // if we trip this we have a two taxon tree with 3 nodes (degree 2 node in the middle)
	for (; nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsShootTip())
			{
			AppendTaxonLabel(s, nd->GetNodeNumber(), taxMgr, useTaxonNames);
			if (edgelens)
				AppendColonAndEdgeLen(s, df, fltlens, *nd);
			if (nd->GetRightSib() != NULL)
				s << ',';
			else if (nd->nextPreorder != NULL)
				{
				const TreeNode *tempAncNode = nd;
				do
					{
					s << ')';
					--openParens;
					assert(openParens > 0);
					tempAncNode = tempAncNode->par;
					assert(tempAncNode != NULL);
					if (edgelens)
						AppendColonAndEdgeLen(s, df, fltlens, *tempAncNode);
					}
				while (tempAncNode->rSib != nd->nextPreorder);
				s << ',';
				}
			}
		else
			{
			// MTH 13-Dec-2004 if (!nd->GetParent()->IsRoot() test no longer needed because root handling is now done outside the loop
			//POL 18-Oct-2004
			//if (!nd->GetParent()->IsRoot())
			//	{
				// Only generate parentheses if this is not the special interior node just above the root node 
				// (i.e. the one with no siblings, only children). The reason for not generating parentheses in
				// this one case is to create a more normal looking unrooted tree description. For example, 
				// rather than (0,(1,2,3,4,5)) for a star tree with 6 taxa, spit out (0,1,2,3,4,5) instead.
				//
			s << '(';
			++openParens;
			//	}
			}
		}
	
	if (openParens > 0)
		{
		const TreeNode *tmpNode = GetLastPreorder();
		while (openParens > 0)
			{
			s << ')';
			--openParens;
			assert(tmpNode != NULL);
			tmpNode = tmpNode->par;
			if (!tmpNode->GetParent()->IsRoot() && edgelens)
				AppendColonAndEdgeLen(s, df, fltlens, *tmpNode);
			}
		}
	return s;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes nLeaves and nInternals to their correct values. Uses preorder pointers to walk through tree, so results
|	will be unpredictable if Tree::firstPreorder or any of the TreeNode::nextPreorder pointers along the way are not
|	set correctly.
*/
void Tree::RefreshNodeCounts() const
	{
	nLeaves = 0;
	nInternals = 0;
	for (TreeNode *nd = firstPreorder; nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsShootTip() || nd->IsRoot())
			++nLeaves;
		else
			++nInternals;
		}
	node_counts_valid = true;
	}

//#define DEBUG_REFRESHID
//#if defined(DEBUG_REFRESHID)
/*----------------------------------------------------------------------------------------------------------------------
|	Brings the `id' data member up to date. It is the responsibility of classes that change the structure of the tree 
|	to call this member function so that the id data member is always valid.
*/
void Tree::RefreshID() const
	{
	id.Clear();

#	if defined(DEBUG_REFRESHID)
		ofstream tmpf("createid.txt");
		tmpf << "split_set size = " << id.GetNumSplitsInID() << endl;
		tmpf << "ntips          = " << GetNLeaves() << endl;
		tmpf << "nnodes         = " << GetNNodes() << endl;
		tmpf << "Split::ntax          = " << Split::splitNTax << endl;
		tmpf << "Split::bits_per_unit = " << Split::bits_per_unit << endl;
		tmpf << "Split::nunits        = " << Split::nunits << endl;
		tmpf << endl;
#	endif

	string s;
	const TreeNode *const rootPtr = GetRoot(); //root is the last postorder
	for (const TreeNode *nd = GetFirstPostorder(); nd != rootPtr; nd = nd->GetNextPostorder())
		{
		nd->CreateSplit();
#		if defined(DEBUG_REFRESHID)
			tmpf << ((nd->IsRoot() || nd->IsShootTip()) ? "tip" : "internal");
			tmpf << " node " << nd->nodeNum << ":\n  split  : ";
			s.clear();
			nd->split.CreateAndAppendPatternRepresentation(&s);
			tmpf << s << "  " << nd->split.CreateIdRepresentation() << "\n\n";
#		endif

		if (!nd->IsShootTip())
			{
			// Add non-trivial splits to split manager the tree id
			//
			id.AddSplit(nd->split);
			
#			if defined(DEBUG_REFRESHID)
				tmpf << "id now has " << id.GetNumSplitsInID() << " splits\n\n";
#			endif
			}
		}

#	if defined(DEBUG_REFRESHID)
		tmpf.close();
#	endif

	id_valid = true;

	//debug_splits_broken = false;
	}
//#else
#if 0
	/*----------------------------------------------------------------------------------------------------------------------
	|	Brings the `id' data member up to date. It is the responsibility of classes that change the structure of the tree 
	|	to call this member function so that the id data member is always valid.
	*/
	void Tree::RefreshID()
		{
		id.Clear();
		const TreeNode *const rootPtr = GetRoot(); //root is the last postorder
		for (TreeNode *nd = GetFirstPostorder();; nd = nd->GetNextPostorder())
			{
			nd->RefreshSplit();
			// Add non-trivial splits to split manager the tree id
			//
			if (!nd->IsLeafOrRoot())
				id.AddSplit(nd->split);
			if (nd == rootPtr)
				break;
			}
		id_valid = true;

		//debug_splits_broken = false;
		}
#endif
