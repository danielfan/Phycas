//#include "phycas/force_include.h"
#include <set>
#include <list>
#include <boost/function.hpp>
#include <iostream>
#include <iomanip>
#include "pyphy/src/probability_distribution.hpp"
#include "pyphy/src/oldphycas/tree_manip.hpp"
#include "pyphy/src/oldphycas/split.hpp"
#include "pyphy/src/oldphycas/tree_id.hpp"
#include "pyphy/src/oldphycas/tree.hpp"
#include "pyphy/src/oldphycas/tree_node.hpp"
#include "pyphy/src/ncl/taxa/nxs_taxa_manager.hpp"
#include "pyphy/src/oldphycas/tree_inl.hpp"
using std::list;
using std::vector;
using std::string;


/*#include "phycas/tree/drawcontext.h"

#include "phycas/rand/lot.h"
#include <boost/bind.hpp>
*/
//#define DEBUG_BUILD_TREE_FROM_ID
#if defined(DEBUG_BUILD_TREE_FROM_ID)
#	include <iostream>
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Expects splits in tree_id to be ordered from smallest (least inclusive) to largest. All nodes in tree have branch 
|	lengths equal to 1.0 after this function. Use SplitManager::SetBrlensFromSplits function to assign branch lengths
|	based on splits in SplitManager. No TreeID is made here, so if SplitManager::SetBrlensFromSplits is not called 
|	(it creates a TreeID), you should call CreateID explicitly after building the tree. Assumes tree_id contains at 
|	least some splits.
|	Note that there is a nearly identical method called TreeManip::SimpleBuildTreeFromID that assumes that all taxa in 
|	the taxa manager are in the tree being built.
*/
void TreeManip::BuildTreeFromID(
  NxsTaxaManager *taxaMgr,	/**< is the taxa manager (knows which taxa are currently active) */
  const TreeID& tree_id,	/**< s the TreeID object specifying the splits to be used in constructing the tree */
  unsigned root_at)			/**< is the index of the tip node serving as the root */	//POL added 22-Oct-2004
	{
	assert(taxaMgr != NULL);
	unsigned nlvs = taxaMgr->GetNumActive();
	unsigned maxlvs = taxaMgr->GetNumTaxa();
	unsigned nextInternalNodeNum = maxlvs;
	const double commonDefEdgeLen = 1.0;

	SplitSet const &ss = tree_id.GetSplitSet();
	assert(!ss.empty());

#if defined(DEBUG_BUILD_TREE_FROM_ID)
	std::string tmps;
	std::ofstream tmpf("idbuild.txt");
	tmpf << "tree_id comprises these " << (unsigned)ss.size() << " splits:" << std::endl;
	for (SplitSet::const_iterator ssi = ss.begin(); ssi != ss.end(); ++ssi)
		{
		const Split &s = (*ssi);
		tmps.clear();
		s.CreateAndAppendPatternRepresentation(&tmps);
		tmpf << tmps << "\n";
		}
	tmpf << std::endl;
#endif

	typedef list<TreeNode *> NodeList;
	NodeList nodes;
	const NxsIndexSet &active_taxa = taxaMgr->GetActiveSet();
	NxsIndexSet::const_iterator iter = active_taxa.begin();
	for (; iter != active_taxa.end(); ++iter)
		{
		unsigned leaf = *iter;
		if (leaf != root_at)
			{
			TreeNode *leafNode = tree->CreateTreeNode(leaf, commonDefEdgeLen, true);
#if defined(DEBUG_BUILD_TREE_FROM_ID)
			tmpf << "Adding node number " << leaf << " to node list" << std::endl;
#endif
			nodes.push_back(leafNode);
			}
		}

	NodeList::iterator tmp;
	for (SplitSet::const_iterator ssi = ss.begin(); ssi != ss.end(); ++ssi)
		{
		const Split &s = (*ssi);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
		tmpf << "\n" << std::setw(6) << "s";
		tmpf << " --> ";
		tmps.clear();
		s.CreateAndAppendPatternRepresentation(&tmps);
		tmpf << tmps << std::endl;
#endif

		TreeNode *newnd = tree->CreateTreeNode(nextInternalNodeNum++, commonDefEdgeLen, false);
		NodeList::iterator it = nodes.begin();
		for (; it != nodes.end();)
			{
			TreeNode *child = (*it);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
			tmpf << std::setw(6) << child->GetNodeNumber() << " --> ";
			tmps.clear();
			child->split.CreateAndAppendPatternRepresentation(&tmps);
			tmpf << tmps;
#endif

			bool subsumed = child->split.SubsumedIn(s);
			if (subsumed)
				{
				InsertSubtree(child, newnd, TreeManip::kOnRight);
				newnd->split.CombineWith(child->split);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
				tmpf << " ==> subsumed in s, deleting node";
#endif

				tmp = it++;
				nodes.erase(tmp);
				}
			else
				++it;

#if defined(DEBUG_BUILD_TREE_FROM_ID)
			tmpf << std::endl;
#endif
			}

#if defined(DEBUG_BUILD_TREE_FROM_ID)
		tmpf << "\nAdding node number " << newnd->GetNodeNumber() << " to node list" << std::endl;
#endif

		nodes.push_back(newnd);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
		tmpf << "Node list now: ";
		for (tmp = nodes.begin(); tmp != nodes.end(); ++tmp)
			{
			tmpf << (*tmp)->GetNodeNumber() << " ";
			}
		tmpf << "\n" << std::endl;
#endif
		}

#if defined(DEBUG_BUILD_TREE_FROM_ID)
	tmpf << "\n\nAll remaining nodes ";
	for (NodeList::iterator it = nodes.begin(); it != nodes.end(); ++it)
		{
		tmpf << "[" << (*it)->GetNodeNumber() << "]";
		}
	tmpf << std::endl;
	//tmpf << " will be added to root (taxon 0):\n  ";
	tmpf << " will be added to root (taxon " << root_at << "):\n  "; //POL 22-Oct-2004

	tmpf.close();
#endif

	tmp = nodes.begin();
	TreeNode *lastNode = *tmp;
	//TreeNode *rootNode = tree->CreateTreeNode(0,0.0, false);
	TreeNode *rootNode = tree->CreateTreeNode(root_at, 0.0, false);	//POL 22-Oct-2004
	InsertSubtree(lastNode, rootNode, TreeManip::kOnRight);

	tree->firstPreorder = rootNode;
	tree->TraverseTree();

	std::cerr << "GetNLeaves() = " << tree->GetNLeaves() << std::endl;
	std::cerr << "nlvs = " << nlvs << std::endl;

	assert(tree->GetNLeaves() == nlvs);
	assert(rootNode == tree->GetRoot());

	tree->RefreshFullID();
	tree->DebugCheckTreeStructure();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Expects splits in tree_id to be ordered from smallest (least inclusive) to largest. All nodes in tree have branch 
|	lengths equal to 1.0 after this function. Use SplitManager::SetBrlensFromSplits function to assign branch lengths
|	based on splits in SplitManager. No TreeID is made here, so if SplitManager::SetBrlensFromSplits is not called 
|	(it creates a TreeID), you should call CreateID explicitly after building the tree.
|	Note that there is a nearly identical method called TreeManip::BuildTreeFromID that takes a reference to a
|	NxsTaxaManager as its first argument. This allows it to handle cases in which fewer taxa are in the tree 
|	specified by the tree_id than are in the data matrix. Also, handles case in which user has excluded some taxa.
|	Because it is never safe to assume that no taxa have been excluded, this version should be deprecated soon.
*/
void TreeManip::SimpleBuildTreeFromID(
  unsigned nlvs,			/* the number of leaves in the tree */
  const TreeID& tree_id,	/* the TreeID object specifying the splits to be used in constructing the tree */
  unsigned root_at)			/* the index of the tip node serving as the root */	//POL added 22-Oct-2004
	{
	const double commonDefEdgeLen = 1.0;
	unsigned nextInternalNodeNum = nlvs;

	SplitSet const &ss = tree_id.GetSplitSet();

#if defined(DEBUG_BUILD_TREE_FROM_ID)
	std::ofstream tmpf("idbuild.txt");
#endif

	if (ss.empty())
		{
#if defined(DEBUG_BUILD_TREE_FROM_ID)
		tmpf << "split set empty, building star tree" << std::endl;
#endif

		SimpleBuildStarTree(nlvs, root_at);	//POL 22-Oct-2004 added root_at
		return;
		}

	typedef list<TreeNode *> NodeList;
	NodeList nodes;
	for (unsigned leaf = 0; leaf < nlvs; ++leaf)
		{
		if (leaf != root_at)
			{
			TreeNode *leafNode = tree->CreateTreeNode(leaf, commonDefEdgeLen, true);
#if defined(DEBUG_BUILD_TREE_FROM_ID)
			tmpf << "Adding node number " << leaf << " to node list" << std::endl;
#endif
			nodes.push_back(leafNode);
			}
		}

	std::string tmps;
	TreeNode *lastNode = NULL;
	NodeList::iterator tmp;
	for (SplitSet::const_iterator ssi = ss.begin(); ssi != ss.end(); ++ssi)
		{
		const Split &s = (*ssi);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
		tmpf << "\n" << std::setw(6) << "s";
		tmpf << " --> ";
		tmps.clear();
		s.CreateAndAppendPatternRepresentation(&tmps);
		tmpf << tmps << std::endl;
#endif

		TreeNode *newnd = tree->CreateTreeNode(nextInternalNodeNum++, commonDefEdgeLen, false);
		NodeList::iterator it = nodes.begin();
		for (; it != nodes.end();)
			{
			TreeNode *child = (*it);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
			tmpf << std::setw(6) << child->GetNodeNumber() << " --> ";
			tmps.clear();
			child->split.CreateAndAppendPatternRepresentation(&tmps);
			tmpf << tmps;
#endif

			bool subsumed = child->split.SubsumedIn(s);
			if (subsumed)
				{
				InsertSubtree(child, newnd, TreeManip::kOnRight);
				newnd->split.CombineWith(child->split);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
				tmpf << " ==> subsumed in s, deleting node";
#endif

				tmp = it++;
				nodes.erase(tmp);
				}
			else
				++it;

#if defined(DEBUG_BUILD_TREE_FROM_ID)
			tmpf << std::endl;
#endif
			}

#if defined(DEBUG_BUILD_TREE_FROM_ID)
		tmpf << "\nAdding node number " << newnd->GetNodeNumber() << " to node list" << std::endl;
#endif

		lastNode = newnd;
		nodes.push_back(newnd);

#if defined(DEBUG_BUILD_TREE_FROM_ID)
		tmpf << "Node list now: ";
		for (tmp = nodes.begin(); tmp != nodes.end(); ++tmp)
			{
			tmpf << (*tmp)->GetNodeNumber() << " ";
			}
		tmpf << "\n" << std::endl;
#endif
		}

#if defined(DEBUG_BUILD_TREE_FROM_ID)
	tmpf << "\n\nAll remaining nodes ";
	for (NodeList::iterator it = nodes.begin(); it != nodes.end(); ++it)
		{
		tmpf << "[" << (*it)->GetNodeNumber() << "]";
		}
	tmpf << std::endl;
	//tmpf << " will be added to root (taxon 0):\n  ";
	tmpf << " will be added to root (taxon " << root_at << "):\n  "; //POL 22-Oct-2004

	tmpf.close();
#endif

#if 1
	// Create the root node
	//
	TreeNode *rootNode = tree->CreateTreeNode(root_at, 0.0, false);	//POL 22-Oct-2004

	// Want to add lastNode to rootNode.
	//
	InsertSubtree(lastNode, rootNode, TreeManip::kOnRight);
	rootNode->split.CombineWith(lastNode->split);
#else
	tmp = nodes.begin();
	TreeNode *lastNode = *tmp;
	//TreeNode *rootNode = tree->CreateTreeNode(0,0.0, false);
	TreeNode *rootNode = tree->CreateTreeNode(root_at, 0.0, false);	//POL 22-Oct-2004
	InsertSubtree(lastNode, rootNode, TreeManip::kOnRight);
#endif

	tree->firstPreorder = rootNode;
	tree->TraverseTree();

	std::cerr << "GetNLeaves() = " << tree->GetNLeaves() << std::endl;
	std::cerr << "nlvs = " << nlvs << std::endl;

	//assert(tree->GetNLeaves() == nlvs);
	assert(rootNode == tree->GetRoot());

	tree->RefreshFullID();
	tree->DebugCheckTreeStructure();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a star tree having `nlvs' leaves. Each edge has a length drawn from an exponential distribution with mean
|	`lambda' (hazard parameter 1/`lambda').
*/
void TreeManip::SimpleBuildStarTree(
  unsigned	nlvs,		/**< number of tip nodes (root plus leaves) */
  unsigned	root_at,		/**< node number of root node (should correspond to index of tip in taxa manager */	//POL added 22-Oct-2004
#if defined(POL_PHYCAS)
  phycas::Lot *		rnd,		/**< random number generator to be used */
#else
  Lot *		rnd,		/**< random number generator to be used */
#endif
  double	lambda)		/**< mean edge length */
	{
	assert(tree != NULL);
	tree->Clear();
	tree->hasEdgelens = true;
	const double commonEdgeLen = 1.0;	// used only if rnd == NULL

	ExponentialDistribution expDist;
	expDist.SetLot(rnd); //POL 20Jun2005
	expDist.SetMeanAndVariance(lambda, lambda*lambda); //@ is this right? check mean and variance of exponential distribution

	// Create the root node
	//
	//TreeNode *rootNd = tree->CreateTreeNode(0, 0.0, false);	
	TreeNode *rootNd = tree->CreateTreeNode(root_at, 0.0, false);	//POL 22-Oct-2004
	tree->firstPreorder = rootNd;

	// Create the hub node
	//
	const double hubLen = (rnd == NULL ? commonEdgeLen : expDist.Sample()); //POL 20Jun2005 was Sample(*rnd)
	TreeNode *hub = tree->CreateTreeNode(nlvs, hubLen, false);
	InsertSubtree(hub, rootNd, TreeManip::kOnLeft);

#if 1	//POL 22-Oct-2004
	for (unsigned i = 0; i < nlvs; ++i)
		{
		if (i != root_at)
			{
			const double ndEdgeLen = (rnd == NULL ? commonEdgeLen : expDist.Sample()); //POL 20Jun2005 was Sample(*rnd)
			TreeNode *nd = tree->CreateTreeNode(i, ndEdgeLen, true);
			InsertSubtree(nd, hub, TreeManip::kOnLeft);
			}
		}
#else
	for (unsigned i = 1; i < nlvs; ++i)
		{
		const double ndEdgeLen = (rnd == NULL ? commonEdgeLen : expDist.Sample()); //POL 20Jun2005 was Sample(*rnd)
		TreeNode *nd = tree->CreateTreeNode(i, ndEdgeLen, true);
		InsertSubtree(nd, hub, TreeManip::kOnLeft);
		}
#endif

	tree->TraverseTree();
	assert(tree->GetNLeaves() == nlvs);
	assert(tree->GetNInternals() == 1);
	assert(rootNd == tree->GetRoot());
	rootNd->split.Clear();
	tree->ForAllLeaves(boost::function1<void, TreeNode*>(&TreeNode::CreateSplit));
	tree->RefreshID();
	tree->DebugCheckTreeStructure();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Let Sk be the time between the kth and (k+1)st birth. Let Wk be the time of the kth birth. Let X(t) be the number of
|	lineages existing at time t. Each lineage is independent and splits with probability lambda dt in each infinitesimal
|	time interval.
|>
|	------------------------------------------ W4 = time of 4th birth = S0 + S1 + S2 + S3
|	 \             \        /     |
|	  \             \      /      |
|	   \             \    /       S3    X(t) = 3
|	    \             \  /        |
|	     \             \/         |
|	------------------------------------------ W3 = time of 3rd birth = S0 + S1 + S2
|	      \            /          |
|	       \          /           |
|	        \        /            |
|	         \      /             S2    X(t) = 2
|	          \    /              |
|	           \  /               |
|	            \/                |
|	------------------------------------------ W2 = time of 2nd birth = S0 + S1
|	            /                 |
|	           /                  |
|	          /                   S1    X(t) = 1
|	         /                    |
|	        /                     |
|	------------------------------------------ W1 = time of 1st birth = S0 = 0.0
|>
|	The cumulative distribution function of the sojourn time S1 given speciation rate lambda is
|>
|	Pr(S1 < t) = 1 - Pr{no speciation in lineage | lambda, t} = 1 - exp{-lambda t}
|>
|	Thus, S1 is exponentially distributed with hazard parameter lambda. The cdf of the sojourn time S2 is 
|>
|	Pr(S2 < t) = 1 - Pr{no speciation in either lineage | lambda, t}
|	           = 1 - Pr{no speciation in lineage 1} Pr{no speciation in lineage 2}
|	           = 1 - exp{-lambda t} exp{-lambda t}
|	           = 1 - exp{-2 lambda t}
|>
|	Thus, S2 is exponentially distributed with hazard parameter 2*lambda. In general, for j > 0, sojourn time Sj is
|	exponentially distributed with hazard parameter j*lambda, where j is the number of extant lineages.
*/
void TreeManip::SimpleBuildYuleTree(
#if defined(POL_PHYCAS)
  phycas::Lot	&rnd,		/* random number generator to be used */
#else
  Lot		&rnd,		/* random number generator to be used */
#endif
  double	lambda,		/* instantaneous speciation rate */
  unsigned	nlvs)		/* number of tip nodes (leaves) in final rooted tree */
	{
	//ASCIIDrawContext ascii(os);

	//@POL this function essentially builds an unrooted tree because node 0 is the root
	// need to decide how to deal with rooted trees - is the root node number 0 for rooted trees?
	assert(tree != NULL);
	tree->Clear();
	tree->hasEdgelens = true;

	unsigned i, j, k;
	double mean, variance, new_edgelen;
	TreeNode **ndarr;
	unsigned nextLeafNodeNum = 2; // 0 is reserved for root node, and 1 is given to the first tip created
	unsigned nextInternalNodeNum = nlvs;

	ExponentialDistribution expDist;
	expDist.SetLot(&rnd); //POL 20Jun2005

	// This vector keeps all nodes representing growing tips handy
	//
	vector<TreeNode *> shootTips;

	// Create the root node
	//
	TreeNode *rootNd = tree->CreateTreeNode(0, 0.0, false);
	tree->firstPreorder = rootNd;

	// Create the first leaf node
	//
	TreeNode *nd = tree->CreateTreeNode(1, 0.0, true);
	shootTips.push_back(nd);

	// Connect the root node to the leaf node
	//
	rootNd->prevPreorder = NULL;
	rootNd->nextPreorder = nd;
	rootNd->lChild = nd;
	rootNd->rSib = NULL;
	rootNd->par = NULL;

	nd->prevPreorder = rootNd;
	nd->nextPreorder = NULL;
	nd->lChild = NULL;
	nd->rSib = NULL;
	nd->par = rootNd;

	for (i = 2; i < nlvs; ++i)
		{
		// Choose one of the current tips to speciate
		//
		j = rnd.SampleUInt(i-1);
		nd = shootTips[j];

		// Create an internal node
		//
		TreeNode *parent = tree->CreateTreeNode(nextInternalNodeNum++, DBL_MAX, false);
		// parent becomes nd, nd can then remain a shoot tip
		//
		parent->prevPreorder = nd->prevPreorder;
		parent->nextPreorder = nd->nextPreorder; 
		parent->lChild = nd;
		parent->rSib = nd->rSib;
		parent->par = nd->par;

		// Check for prevPreorder pointer directed at nd, redirect to parent
		//
		if (nd->nextPreorder != NULL)
			nd->nextPreorder->prevPreorder = parent;

		// Check for nextPreorder pointer directed at nd, redirect to parent
		//
		if (nd->prevPreorder != NULL)
			nd->prevPreorder->nextPreorder = parent;

		// Check for lChild pointer directed at nd, redirect to parent
		//
		assert(nd->par != NULL);
		if (nd->par->lChild == nd)
			nd->par->lChild = parent;

		// Check for rSib pointer directed at nd, redirect to parent
		//
		TreeNode *lSib = nd->FindLeftSib();
		if (lSib != NULL)
			lSib->rSib = parent;

		nd->prevPreorder = parent;
		nd->nextPreorder = NULL;
		nd->lChild = NULL;
		nd->rSib = NULL;
		nd->par = parent;

		parent->SetFltEdgeLen(nd->GetFltEdgeLen()); //@POL will nd->GetFltEdgeLen() always return 0.0?
		nd->SetFltEdgeLen(0.0);

		// Create a new tip representing the sister taxon
		//
		TreeNode *sister = tree->CreateTreeNode(nextLeafNodeNum++, 0.0, true);
		nd->rSib = sister;
		nd->nextPreorder = sister;

		sister->prevPreorder = nd;
		sister->nextPreorder = parent->nextPreorder;
		sister->lChild = NULL;
		sister->rSib = NULL;
		sister->par = parent;

		parent->nextPreorder = nd;

		// Choose a new brlen
		//
		double mean = 1.0/(lambda*(i - 1));
		double variance = mean*mean;
		expDist.SetMeanAndVariance(mean, variance);
		new_edgelen = expDist.Sample(); // POL 20Jun2005 was Sample(rnd)

		// Add the new edge length to all existing tips
		//
		ndarr = &shootTips[0];
		for (k = 0; k < i - 1; ++k)
			{
			if (k == j)
				parent->SetFltEdgeLen(new_edgelen + parent->GetFltEdgeLen());
			else
				ndarr[k]->SetFltEdgeLen(new_edgelen + ndarr[k]->GetFltEdgeLen());
			}

		// Add sister to the vector of shoot tips
		//
		shootTips.push_back(sister);

		//tree->TraverseTree();
		//tree->Draw(ascii);
		}

	// Choose a final brlen
	//
	assert(i == nlvs);
	mean = 1.0/(lambda*nlvs);
	variance = mean*mean;
	expDist.SetMeanAndVariance(mean, variance);
	new_edgelen = expDist.Sample(); //POL 20Jun2005 was Sample(rnd)

	// Add the new edge length to all existing tips
	//
	ndarr = &shootTips[0];
	for (k = 0; k < nlvs - 1; ++k)
		{
		ndarr[k]->SetFltEdgeLen(new_edgelen + ndarr[k]->GetFltEdgeLen());
		}

	tree->TraverseTree();
	assert(tree->GetNLeaves() == nlvs);
	tree->RefreshFullID();
	tree->DebugCheckTreeStructure();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Detaches node 's' from tree, keeping intact subtree rooted at `s'. Assumes `s' is non-NULL and is neither the root
|	node nor the only child of the root node. Sets all pointers of `s' to NULL except lChild and nextPreorder.
*/
void TreeManip::DetachSubtree(
  TreeNode *s)	/* the root of the subtree to be detached */
	{
	assert(s != NULL);
	assert(!s->IsRoot());
	assert(!s->GetParent()->IsRoot());

	// Save pointers to relevant nodes
	//
	TreeNode *s_par					= s->par;
	TreeNode *s_lSib				= s->FindLeftSib();
	TreeNode *s_rSib				= s->rSib;
	TreeNode *s_prevPreorder		= s->prevPreorder;
	TreeNode *slast					= s->FindLastPreorderInClade();
	TreeNode *slast_nextPreorder	= slast->nextPreorder;

	// Completely detach s and seal up the wound
	//
	s->par = NULL;

	s->rSib = NULL;
	if (s_lSib == NULL)
		s_par->lChild = s_rSib;
	else
		s_lSib->rSib = s_rSib;

	s->prevPreorder = NULL;
	if (s_prevPreorder == NULL)
		tree->firstPreorder = slast_nextPreorder;
	else
		s_prevPreorder->nextPreorder = slast_nextPreorder;

	slast->nextPreorder = NULL;
	if (slast_nextPreorder == NULL)
		tree->lastPreorder = s_prevPreorder;
	else
		slast_nextPreorder->prevPreorder = s_prevPreorder;

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts node 's' as child of 'u' keeping intact subtree rooted at `s'. Assumes both `s' and `u' are non-NULL. If 
|	`targetSib' is specified, `s' will become an immediate sibling of `targetSib' (right or left sibling determined by 
|	`m'). If `targetSib' is not specified, `s' will be added as the rightmost or leftmost child of `u' (again,
|	depending on value of `m').
*/
void TreeManip::InsertSubtree(
  TreeNode *s,			/* node to be added as a child of u */
  TreeNode *u,			/* the new parent of s */
  InsertMode m,			/* if kOnRight, s will be added to the right of targetSib or, if targetSib is not specified, s will become rightmost child of u (vice versa if kOnLeft specified) */
  TreeNode *targetSib)	/* if non-NULL, s will become an immediate sibling of this node (right or left depending on value of `m') */
	{
	assert(u != NULL);
	assert(s != NULL);
	assert(targetSib == NULL || (targetSib->par == u));

	TreeNode *slast					= s->FindLastPreorderInClade();
	TreeNode *u_lChild				= u->lChild;
	TreeNode *u_rChild				= u->FindRightChild();
	TreeNode *ulast					= u->FindLastPreorderInClade();
	TreeNode *ulast_nextPreorder	= ulast->nextPreorder;
	TreeNode *u_nextPreorder		= u->nextPreorder;

	// Identify possible reductions in complexity
	//
	if (targetSib != NULL)
		{
		if ((m == TreeManip::kOnLeft) && (targetSib == u_lChild))
			targetSib = NULL;
		if ((m == TreeManip::kOnRight) && (targetSib == u_rChild))
			targetSib = NULL;
		if ((m == TreeManip::kOnLeft) && (targetSib != NULL))
			{
			targetSib = targetSib->FindLeftSib();
			m = TreeManip::kOnRight;
			}
		}

	// Make s a child of u
	//
	if (targetSib == NULL)
		{
		if (m == TreeManip::kOnRight)
			{
			s->rSib				= NULL;
			s->par				= u;
			s->prevPreorder 	= ulast;
			slast->nextPreorder = ulast_nextPreorder;
			ulast->nextPreorder	= s;

			if (ulast_nextPreorder == NULL)
				tree->lastPreorder = slast;
			else
				ulast_nextPreorder->prevPreorder = slast;

			if (u_lChild == NULL)
				u->lChild		= s;
			else
				u_rChild->rSib	= s;
			}
		else
			{
			s->rSib					= u_lChild;
			s->par					= u;
			s->prevPreorder 		= u;
			slast->nextPreorder		= u_nextPreorder;

			if (u_lChild != NULL)
				u_lChild->prevPreorder	= slast;
			u->lChild				= s;
			u->nextPreorder			= s;
			}
		}
	else
		{
		// Save pointers to relevant nodes
		//
		TreeNode *targetSib_rSib				= targetSib->rSib;
		TreeNode *targetSiblast					= targetSib->FindLastPreorderInClade();
		TreeNode *targetSiblast_nextPreorder	= targetSiblast->nextPreorder;

		// Can assume that targetSib is not rightmost child of u and also that m == TreeManip::kOnRight
		//
		s->rSib						= targetSib_rSib;
		s->par						= u;
		s->prevPreorder 			= targetSiblast;
		slast->nextPreorder 		= targetSiblast_nextPreorder;
		targetSiblast->nextPreorder	= s;
		targetSib->rSib				= s;
		}

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Removes `s' (a sibling of `u') and makes it a child of `u'. If `m' equals `TreeManip::kOnRight', `s' will be added
|	as the immediate right sibling of `targetSib' or, if `targetSib' is not specified, `s' will become the rightmost 
|	child of `u'. If `m' is `TreeManip::kOnLeft', `s' will be added as the immediate left sibling of `targetSib' or, 
|	if `targetSib' is not specified, `s' will become the leftmost child of `u'. Assumes that `u', `u->par',`s' and 
|	`s->par' are all non-NULL and that `u->par' equals `s->par'. Also assumes that targetSib is a child of u if 
|	targetSib is specified.
|>
|	Before:                After:
|	
|	a       b                s   a   b
|	 \     /		          \  |  /
|	  \   /		               \ | /
|	   \ /		                \|/
|	    u	c   s                u       c
|	     \  |  /		          \     /
|	      \ | /			           \   /
|	       \|/			            \ /
|	        v	                     v
|	       /				        /
|>
*/
void TreeManip::SibToChild(
  TreeNode *u,			/* this node's sibling will be removed and attached as the leftmost child of this node */
  TreeNode *s,			/* the sibling of u to be removed */
  InsertMode ,			/* if kOnRight, s will be added to the right of targetSib or, if targetSib is not specified, s will become rightmost child of u (vice versa if kOnLeft specified) */
  TreeNode *targetSib)	/* if specified, s will become and immediate sibling of this node (which should be a child of u) */
	{
	assert(u				!= NULL);
	assert(u->par			!= NULL);
	assert(s->par			!= NULL);
	assert(u->par == s->par);
	assert(targetSib == NULL || (targetSib->par == u));

	DetachSubtree(s);
	InsertSubtree(s, u, TreeManip::kOnRight, targetSib);

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Removes leaf node `u' and its associated edge.
*/
void TreeManip::DeleteLeaf(
  TreeNode *u)	/* this node's edge will be removed */
	{
	// Note which nodes have pointers potentially aimed at u (will need to be redirected)
	//
	TreeNode *u_par		= u->par;
	TreeNode *u_rSib	= u->rSib;
	TreeNode *u_lSib	= u->FindLeftSib();
	TreeNode *u_prevPre = u->prevPreorder;
	TreeNode *u_nextPre = u->nextPreorder;

	// Redirect pointers aimed at u
	//
	if (u_nextPre == NULL)
		tree->lastPreorder = u_prevPre;
	else
		u_nextPre->prevPreorder = u_prevPre;

	if (u_lSib == NULL)
		u_par->lChild = u_rSib;
	else
		u_lSib->rSib = u_rSib;

	assert(u_prevPre != NULL);
	u_prevPre->nextPreorder = u_nextPre;

	// Detach and then delete u
	//
	assert(u->lChild == NULL);
	u->rSib			= NULL;
	u->par			= NULL;
	u->prevPreorder = NULL;
	u->nextPreorder = NULL;
	tree->StoreTreeNode(u);

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Removes leftmost child of node `u' and makes it the immediate left sibling of w. Used in performing Larget-Simon
|	moves when polytomies are present. Does not necessarily leave tree in a consistent state (i.e. tree may contain
|	nodes of order 2 afterwards).
|>
|	Before:                After:
|	
|	      c       w              c   a   w
|	       \     /		          \  |  /
|	        \   /		           \ | /
|	         \ /		            \|/
|	  a   b   v			     b       v
|	   \  |  /			      \     /
|	    \ | /			       \   /
|	     \|/			        \ /
|	      u				         u
|	     /				        /
|>
*/
void TreeManip::LChildToLSib(
  TreeNode *u,	/* this node's leftmost child will be removed */
  TreeNode *w)	/* the removed node will become this node's left sibling */
	{
	assert(u				!= NULL);
	assert(u->lChild		!= NULL);	// must have a node to mode
	assert(w				!= NULL);
	assert(w->par			!= NULL);
	assert(w->par->lChild	!= NULL);

	// Make copies of the important pointers
	//
	TreeNode *v						= w->par;
	TreeNode *a						= u->lChild;
	TreeNode *a_last				= a->FindLastPreorderInClade();
	TreeNode *a_rSib				= a->rSib;
	TreeNode *a_prevPre				= a->prevPreorder;
	TreeNode *w_prevPre				= w->prevPreorder;
	TreeNode *w_lSib				= w->FindLeftSib();
	TreeNode *w_lSibLast			= (w_lSib == NULL ? NULL : w_prevPre);
	//TreeNode *v_lChild				= v->lChild;

	// Make the switch
	//
	a->rSib							= w;
	a->par							= v;
	a->prevPreorder					= w_prevPre;
	a_last->nextPreorder			= w;
	a_rSib->prevPreorder			= a_prevPre;
	u->lChild						= a_rSib;
	u->nextPreorder					= a_rSib;
	w->prevPreorder					= a_last;
	if (w_lSib == NULL)
		{
		v->lChild					= a;
		v->nextPreorder				= a;
		if (u == w)
			u->lChild				= a_rSib;
		}
	else
		{
		w_lSib->rSib				= a;
		w_lSibLast->nextPreorder	= a;
		}

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Removes rightmost child of node `u' and makes it the immediate right sibling of w. Used in performing Larget-Simon
|	moves when polytomies are present. Does not necessarily leave tree in a consistent state (i.e. tree may contain
|	nodes of order 2 afterwards). Assumes w either equals u or is above u in the tree as it is currently rooted.
|>
|	Before:                After:
|	
|	w       c                w   a   c
|	 \     /		          \  |  /
|	  \   /		               \ | /
|	   \ /		                \|/
|	    v	b  a                 v       b
|	     \  |  /		          \     /
|	      \ | /			           \   /
|	       \|/			            \ /
|	        u	                     u
|	       /				        /
|>
*/
void TreeManip::RChildToRSib(
  TreeNode *u,	/* this node's rightmost child will be removed */
  TreeNode *w)	/* the removed node will become this node's right sibling */
	{
	assert(u				!= NULL);
	assert(u->lChild		!= NULL);
	assert(w				!= NULL);
	assert(w->par			!= NULL);
	assert(w->par->lChild	!= NULL);

	// Make copies of the important pointers
	//
	TreeNode *w_par						= w->par;
	TreeNode *a							= u->FindRightChild();
	//	TreeNode *a_par						= a->GetParent();
	TreeNode *w_last					= w->FindLastPreorderInClade();
	TreeNode *w_last_nextPre			= w_last->nextPreorder;
	TreeNode *a_last					= a->FindLastPreorderInClade();
	TreeNode *a_last_nextPre			= a_last->nextPreorder;
	TreeNode *w_rSib					= w->rSib;
	TreeNode *a_lSib					= a->FindLeftSib();
	TreeNode *a_lSib_last				= a_lSib->FindLastPreorderInClade();

	// Make the switch
	//
	w->rSib								= a;
	a->par								= w_par;
	a->rSib								= w_rSib;
	a_lSib->rSib						= NULL;

	if (u == w)
		{
		a_lSib_last->nextPreorder			= a;
		a->prevPreorder						= a_lSib_last;
		}
	else
		{
		w_last->nextPreorder				= a;
		a->prevPreorder						= w_last;

		if (w_last_nextPre != a)
			{
			w_last_nextPre->prevPreorder	= a_last;
			a_last->nextPreorder			= w_last_nextPre;
			a_lSib_last->nextPreorder		= a_last_nextPre;
			}

		if (a_last_nextPre == NULL && a_lSib_last != w_last)
			tree->lastPreorder				= a_lSib_last;
		else if (w_last_nextPre != a)
			a_last_nextPre->prevPreorder	= a_lSib_last;

		if (w_rSib != NULL)
			w_rSib->prevPreorder			= a_last;
		}

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a specialized NNI swap for the following situation (where swap2 is the deepest node in the affected
|	subtree):
|>
|	      swap1   B
|	         \   /
|	          \ /
|	     A     V
|	      \   /
|	       \ /
|	        U
|	       /
|	      /
|	   swap2
|>
|	This move effectively moves V down (carrying B along with it) until V is between U and swap2 (thinking of it as
|	moving U - with A attached - up until it is between V and swap1 yields the same rearrangement). If polytomies
|	are allowed, V carries not only B but any other siblings of B (except swap1) as it slides along the path from
|	swap1 to swap2. Note: because swap2 must be the parent of U, there is no need to actually supply it to this
|	function.
*/
void TreeManip::NNISwapSpecial(
  TreeNode *swap1)
	{
	assert(swap1 != NULL);

	TreeNode *V = swap1->GetParent();
	assert(V != NULL);

	TreeNode *U = V->GetParent();
	assert(U != NULL);

	TreeNode *R = U->FindRightChild();
	if (R == V)
		R = swap1;

	TreeNode *L = U->GetLeftChild();
	if (L == V)
		L = swap1;

	while (U->GetLeftChild() != V)
		LChildToLSib(U, swap1);
	while (V->GetRightSib() != NULL)
		RChildToRSib(U, swap1);
	while (V->GetLeftChild() != L)
		LChildToLSib(V, V);
	while (R->GetRightSib() != NULL)
		RChildToRSib(V, V);

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Swaps the two nodes specified by swap1 and swap2. Assumes tree, tree->root, swap1 and swap2 are all non-NULL.
|	Polytomies are ok at u or v or both, but assumes that the tree is NOT rooted at either swap1 or swap2.
*/
void TreeManip::NNISwap(
  TreeNode *swap1, 
  TreeNode *swap2)
	{
	assert(tree != NULL);
	assert(tree->GetRoot() != NULL);
	assert(tree->GetRoot() != swap1);
	assert(tree->GetRoot() != swap2);
	assert(swap1 != NULL);
	assert(swap2 != NULL);
	assert(swap1->GetParent() != NULL);
	assert(swap2->GetParent() != NULL);

	// Let node v be the higher of the two parent nodes, with u being the lower of the two parent nodes, 
	// and let x be v's child and y be u's child:
	//
	//           \ /   /
	//            x   /
	//             \ /
	//        \ /   v
	//         y   /
	//          \ /
	//           u
	//          /
	//
	// Assume swap1 is the higher of the two swap nodes
	//
	TreeNode* x = swap1;
	TreeNode* v = swap1->GetParent();
	TreeNode* y = swap2;
	TreeNode* u = swap2->GetParent();

	if (u->GetParent() == v)
		{
		// swap2 was actually the higher of the two swap nodes
		//
		y = swap1;
		u = swap1->GetParent();
		x = swap2;
		v = swap2->GetParent();
		}

	assert(v->GetParent() == u);

	// +------------------------------------------------------------------+
	// | First order of business is to save the new navigational pointer  |
	// | targets before we go changing things                             |
	// +------------------------------------------------------------------+

	TreeNode *xlast 					= x->FindLastPreorderInClade();
	TreeNode *ylast 					= y->FindLastPreorderInClade();

	TreeNode *xlastnext 				= xlast->nextPreorder;
	TreeNode *ylastnext 				= ylast->nextPreorder;

	TreeNode *u_lChild					= u->lChild;
	TreeNode *v_lChild					= v->lChild;

	TreeNode *x_rSib					= x->rSib;
	TreeNode *y_rSib					= y->rSib;

	TreeNode *x_par						= x->par;
	TreeNode *y_par						= y->par;

	TreeNode *x_prevPreorder			= x->prevPreorder;
	TreeNode *y_prevPreorder			= y->prevPreorder;

	TreeNode *xlast_nextPreorder		= xlast->nextPreorder;
	TreeNode *ylast_nextPreorder		= ylast->nextPreorder;

	TreeNode *xlsib 					= x->FindLeftSib();
	TreeNode *ylsib 					= y->FindLeftSib();

	TreeNode *xlsiblast					= (xlsib == NULL ? NULL : xlsib->FindLastPreorderInClade());
	TreeNode *ylsiblast					= (ylsib == NULL ? NULL : ylsib->FindLastPreorderInClade());

	// +--------------------------------------------------------------------------+
	// | Now that we are done making pointer reassigments, we can make the switch |
	// +--------------------------------------------------------------------------+
	u->lChild				= (u_lChild == y ? x : u_lChild);
	v->lChild				= (v_lChild == x ? y : v_lChild);

	x->rSib					= y_rSib;
	y->rSib					= x_rSib;

	x->par					= y_par;
	y->par					= x_par;

	x->prevPreorder			= (y_prevPreorder == xlast ? ylast : y_prevPreorder);
	y->prevPreorder			= x_prevPreorder;

	xlast->nextPreorder		= ylast_nextPreorder;
	ylast->nextPreorder		= (xlast_nextPreorder == y ? x : xlast_nextPreorder);

	u->nextPreorder			= (u_lChild == y ? x : u_lChild);
	v->nextPreorder			= (v_lChild == x ? y : v_lChild);

	if (xlsib != NULL)
		{
		xlsib->rSib = y;
		xlsiblast->nextPreorder = y;
		}
	if (ylsib != NULL)
		{
		ylsib->rSib = x;
		if (ylsiblast != xlast)
			ylsiblast->nextPreorder = x;
		}

	if (xlastnext == NULL)
		tree->lastPreorder = ylast;
	else if (xlastnext != y)
		xlastnext->prevPreorder	= ylast;

	if (ylastnext == NULL)
		tree->lastPreorder = xlast;
	else if (ylastnext != x)
		ylastnext->prevPreorder	= xlast;

	// This rearrangement invalidates the TreeID and node counts
	//
	tree->InvalidateID();
	tree->InvalidateNodeCounts();
	}
