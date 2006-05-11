#include "phycas/force_include.h"
#include "pyphy/phylogeny/tree_manip.hpp"
#include "pyphy/prob_dist/basic_lot.hpp"
#include "phycas/rand/probability_distribution.hpp"
#include "pyphy/phylogeny/basic_tree.hpp"

using namespace phycas;

#if POLPY_NEWWAY
#include "pyphy/prob_dist/basic_lot.hpp"
#endif

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Assigns all edge lengths in the tree using independent draws from the ProbabilityDistribution object pointed to by
|	`d'.
*/
void TreeManip::setRandomEdgeLens(ProbDistShPtr d)
	{
	//@POL would like to do something like that shown below, but doesn't compile with d wrapped in var, and not using
	// var wrapper results in Sample() only getting called once (all edge lengths set to the same value)
	//std::for_each(tree->begin(), tree->end(), boost::lambda::bind(&TreeNode::SetEdgeLen, boost::lambda::_1, boost::lambda::var(d)->Sample()));
	for (preorder_iterator nd = tree->begin(); nd != tree->end(); ++nd)
		{
		nd->SetEdgeLen(d->Sample());
		}
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Begins with left child of parent of `start' and calls GetRightSib() until the left sibling of `start' is located.
|	Assumes `start' is non-NULL.
*/
TreeNode * TreeManip::FindLeftSib(
  TreeNode * start)	/**< is the node whose left sibling is sought */
	{
	assert(start != NULL);
	TreeNode * nd = start->GetParent();

	// If start has no parent, then there's no way it can have a left sibling
	if (nd == NULL)
		return NULL;	

	nd = nd->GetLeftChild();

	// Parent of start should have at least one child
	assert(nd != NULL); 

	// If left child of start's parent equals start, then start is an only child and has no siblings
	if (nd == start)
		return NULL;	

	TreeNode * leftsib = NULL;
	while (nd != start)
		{
		if (nd == NULL)
			{
			throw XPhylogeny("pointer inconsistency in FindLeftSib");
			}
		leftsib = nd;
		nd = nd->GetRightSib();
		}
	return leftsib;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Begins with left child of `start'. If left child is NULL, returns NULL, otherwise, returns the rightmost sibling of
|	the left child of `start'. Assumes `start' is non-NULL.
*/
TreeNode * TreeManip::FindRightmostChild(
  TreeNode * start)	/**< is the parent node whose children will be searched */
	{
	assert(start != NULL);
	TreeNode * curr = start->GetLeftChild();
	TreeNode * rightmost = NULL;
	while (curr != NULL)
		{
		rightmost = curr;
		curr = curr->GetRightSib();
		}
	return rightmost;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Zig-zags up the right side of the clade starting with the node `start' until reaching a tip. This tip node is the 
|	last node in the clade according to the preorder sequence. Zig-zagging is necessary because there is no pointer to
|	the right-most child of a node, so we must use FindRightmostChild to move up the clade. Assumes `start' is non-NULL.
*/
TreeNode * TreeManip::FindLastPreorderInClade(
  TreeNode * start)	/**< is the deepest node belonging to the clade in question */
	{
	assert(start != NULL);
	TreeNode * curr = start;
	TreeNode * rChild = FindRightmostChild(curr);
	while (rChild != NULL)
		{
		curr = rChild;
		rChild = FindRightmostChild(curr);
		}
	return curr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Inserts node 's' as child of 'u' keeping the subtree rooted at `s' intact. Assumes both `s' and `u' are non-NULL. If 
|	`targetSib' is specified, `s' will become an immediate sibling of `targetSib' (right or left sibling determined by 
|	`m'). If `targetSib' is not specified, `s' will be added as the rightmost or leftmost child of `u' (again, depending
|	on the value of `m').
*/
void TreeManip::InsertSubtree(
  TreeNode * s,			/**< is the node to be added as a child of u */
  TreeNode * u,			/**< is the new parent of s */
  InsertMode m,			/**< if kOnRight, s will be added to the right of targetSib or, if targetSib is not specified, s will become rightmost child of u (vice versa if kOnLeft specified) */
  TreeNode * targetSib)	/**< if non-NULL, s will become an immediate sibling of this node (right or left depending on value of `m') */
	{
	assert(u != NULL);
	assert(s != NULL);
	assert(targetSib == NULL || (targetSib->par == u));

	TreeNode * slast				= FindLastPreorderInClade(s);
	TreeNode * u_lChild				= u->lChild;
	TreeNode * u_rChild				= FindRightmostChild(u);
	TreeNode * ulast				= FindLastPreorderInClade(u);
	TreeNode * ulast_nextPreorder	= ulast->nextPreorder;
	TreeNode * u_nextPreorder		= u->nextPreorder;

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
			targetSib = FindLeftSib(targetSib);
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
		TreeNode * targetSib_rSib				= targetSib->rSib;
		TreeNode * targetSiblast				= FindLastPreorderInClade(targetSib);
		TreeNode * targetSiblast_nextPreorder	= targetSiblast->nextPreorder;

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
	tree->InvalidateTreeID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a star tree having `ntips' tips. Each edge has a length drawn from the supplied ProbabilityDistribution. It
|	is assumed that the random number generator has already been set for `edge_len_dist'.
*/
void TreeManip::starTree(
  unsigned		ntips,			/**< Number of tip nodes (including the tip at which the tree is rooted) */
  ProbDistShPtr	edge_len_dist)	/**< Probability distribution from which to draw the edge lengths */
	{
	assert(ntips > 0);
	assert(edge_len_dist);
	assert(tree);

	tree->Clear();
	tree->hasEdgeLens = true;

	// Create the root node
	//
	TreeNode * rootNd = tree->GetNewNode();
	rootNd->SetEdgeLen(0.0);
	rootNd->SetNodeNum(0);
	tree->firstPreorder = rootNd;

	// Create the hub node
	//
	const double hubLen = edge_len_dist->Sample();
	TreeNode * hub = tree->GetNewNode();
	hub->SetEdgeLen(hubLen);
	hub->SetNodeNum(ntips);
	InsertSubtree(hub, rootNd, TreeManip::kOnLeft);

	for (unsigned i = 1; i < ntips; ++i)
		{
		const double ndEdgeLen = edge_len_dist->Sample();
		TreeNode * nd = tree->GetNewNode();
		nd->SetEdgeLen(ndEdgeLen);
		nd->SetNodeNum(i);
		InsertSubtree(nd, hub, TreeManip::kOnLeft);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If the `yule' parameter is true, and the supplied probability distribution is ExponentialDist(lambda), this function
|	builds a Yule tree with speciation rate `lambda'. If the `yule' parameter is false, then the supplied probability
|	distribution is used to generated independent edge lengths and the resulting tree is not expected to be ultrametric.
|	If you want the same random number generator to be used for both edge lengths and branching order, call the setLot()
|	function of the probability distribution object, passing the same Lot object used in the call to 
|	TreeManip::randomTree. The resulting tree is unrooted, so the supplied `ntips' includes the tip serving as the root
|	node. A Yule tree is more appropriately interpreted as a rooted tree, but thus far Phycas only deals with unrooted 
|	trees so we are using this function as simply a way to generate a random topology. 
|	
|	Yule trees: Let Sk be the sojourn time between the kth and (k+1)st birth. Let Wk be	the time of the kth birth. Let 
|	X(t) be the number of lineages existing at time t. Each lineage is independent and splits with probability `lambda' 
|	dt in each infinitesimal time interval.
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
|	The cumulative distribution function of the sojourn time S1 given speciation rate `lambda' is
|>
|	Pr(S1 < t) = 1 - Pr{no speciation in lineage | lambda, t} = 1 - exp{-lambda t}
|>
|	Thus, S1 is exponentially distributed with hazard parameter `lambda'. The cumulative distribution function of the 
|	sojourn time S2 is 
|>
|	Pr(S2 < t) = 1 - Pr{no speciation in either lineage | lambda, t}
|	           = 1 - Pr{no speciation in lineage 1} Pr{no speciation in lineage 2}
|	           = 1 - exp{-lambda t} exp{-lambda t}
|	           = 1 - exp{-2 lambda t}
|>
|	Thus, S2 is exponentially distributed with hazard parameter 2*`lambda'. In general, for j > 0, sojourn time Sj is
|	exponentially distributed with hazard parameter j*`lambda', where j is the number of extant lineages.
|	
|	While finding the MLE of `lambda' is not relevant here, it might be useful to calling functions to obtain the 
|	MLE of `lambda' in order to supply that value to this function (to create a Yule tree that is similar to a certain
|	known tree, for example). If one is given T, the total tree length, and the number of non-root tip nodes n, the MLE
|	of `lambda' is simply n/T. Here is the derivation for a tree with 4 non-root tips as an example:
|>
|	Pr(S1) = lambda*exp(-lambda*S1)
|	Pr(S2) = 2*lambda*exp(-2*lambda*S2)
|	Pr(S3) = 3*lambda*exp(-3*lambda*S3)
|	Pr(S4) = 4*lambda*exp(-4*lambda*S4)
|>
|	The likelihood of `lambda' is thus proportional to the product of these four expressions:
|>
|	L = (4!)*(lambda^4)*exp(-lambda*T)
|>
|	The log-likelihood is thus
|>
|	lnL = n ln(lambda) - lambda*T
|>
|	The derivative of lnL with respect to lambda is n/lambda - T, the maximum of which is at n/T.
*/
void TreeManip::randomTree(
  unsigned ntips,				/**< is the number of tip nodes in the final tree (includes the tip node serving as the root) */
  LotShPtr	rng,				/**< is the random number generator used to determine the branching order (but not edge lengths) */
  ProbDistShPtr	edge_len_dist,	/**< is the probability distribution used to generate new edge lengths (should be exponential(lambda) if `yule' is true, where lambda is the instantaneous speciation rate) */
  bool yule)					/**< if true, a Yule tree will be generated; otherwise, edge lengths will be independent draws from `edge_len_dist' */
	{
	//@POL should assert that edge_len_dist is of type ExponentialDist if yule is true?
	assert(edge_len_dist);
	assert(tree);
	tree->Clear();
	tree->hasEdgeLens = true;

	unsigned i, j, k;
	double new_edgelen;
	unsigned nextTipNodeNum = 2; // 0 is reserved for root node, and 1 is given to the first tip created
	unsigned nextInternalNodeNum = ntips;

	// This vector keeps all nodes representing growing tips handy
	//
	std::vector<TreeNode *> shootTips;

	// Create the root node
	//
	TreeNode * rootNd = tree->GetNewNode();
	rootNd->SetEdgeLen(0.0);
	rootNd->SetNodeNum(0);
	tree->firstPreorder = rootNd;

	// Create the first leaf node (note: root is a tip but not a leaf, shoot tip is synonymous with leaf)
	//
	TreeNode * nd = tree->GetNewNode();
	nd->SetEdgeLen(yule ? 0.0 : edge_len_dist->Sample());
	nd->SetNodeNum(1);
	shootTips.push_back(nd);

	// Connect the root node to the leaf node
	//
	rootNd->lChild = nd;
	rootNd->rSib = NULL;
	rootNd->par = NULL;

	nd->lChild = NULL;
	nd->rSib = NULL;
	nd->par = rootNd;

	for (i = 2; i < ntips; ++i)
		{
		unsigned num_leaves = i - 1;

		// Choose one of the current leaves to speciate
		//
		j = rng->SampleUInt(num_leaves);
		nd = shootTips[j];

		// Create an internal node
		//
		TreeNode * parent = tree->GetNewNode();
		parent->SetEdgeLen(yule ? 0.0 : edge_len_dist->Sample()); // was DBL_MAX, can that be right?
		parent->SetNodeNum(nextInternalNodeNum++);

		// parent becomes nd, nd can then remain a leaf
		//
		parent->lChild = nd;
		parent->rSib = nd->rSib;
		parent->par = nd->par;

		// Check for rSib pointer directed at nd, redirect to parent
		//
		TreeNode * lSib = FindLeftSib(nd);
		if (lSib != NULL)
			lSib->rSib = parent;

		// Check for lChild pointer directed at nd, redirect to parent
		//
		assert(nd->par != NULL);
		if (nd->par->lChild == nd)
			nd->par->lChild = parent;

		nd->lChild = NULL;
		nd->rSib = NULL;
		nd->par = parent;

		if (yule)
			{
			parent->SetEdgeLen(nd->GetEdgeLen());
			nd->SetEdgeLen(0.0);
			}

		// Create a new tip representing the sister taxon
		//
		TreeNode * sister = tree->GetNewNode();
		sister->SetEdgeLen(yule ? 0.0 : edge_len_dist->Sample());
		sister->SetNodeNum(nextTipNodeNum++);
		nd->rSib = sister;

		sister->lChild = NULL;
		sister->rSib = NULL;
		sister->par = parent;

		// Choose a new edge length if building a Yule tree
		//
		if (yule)
			{
			new_edgelen = edge_len_dist->Sample()/num_leaves;

			// Add the new edge length to all existing tips
			//
			for (k = 0; k < num_leaves; ++k)
				{
				if (k == j)
					parent->SetEdgeLen(new_edgelen + parent->GetEdgeLen());
				else
					shootTips[k]->SetEdgeLen(new_edgelen + shootTips[k]->GetEdgeLen());
				}
			}

		// Add sister to the vector of shoot tips
		//
		shootTips.push_back(sister);
		}

	// Choose a final edge length if building a Yule tree
	//
	if (yule)
		{
		assert(i == ntips);
		new_edgelen = edge_len_dist->Sample()/ntips;

		// Add the new edge length to all existing tips
		//
		for (k = 0; k < ntips - 1; ++k)
			{
			shootTips[k]->SetEdgeLen(new_edgelen + shootTips[k]->GetEdgeLen());
			}
		}

	tree->RefreshPreorder();
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
|	This move effectively moves V down (carrying B along with it) until V is between U and swap2 (you can also think of 
|	it as moving U up - with A attached - until it is between V and swap1). If polytomies are allowed, V carries not 
|	only B but any other siblings of B (except swap1) as it slides along the path from swap1 to swap2. Note: because 
|	swap2 must be the parent of U, there is no need to actually supply it to this function.
*/
void TreeManip::NNISwapSpecial(
  TreeNode * swap1)	/**< is the node whose parent is slid down one node, carrying all of its children except swap1 */
	{
	assert(swap1 != NULL);

	TreeNode * V = swap1->GetParent();
	assert(V != NULL);

	TreeNode * U = V->GetParent();
	assert(U != NULL);

	TreeNode * R = FindRightmostChild(U);
	if (R == V)
		R = swap1;

	TreeNode * L = U->GetLeftChild();
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
	tree->InvalidateTreeID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Swaps the two nodes specified by swap1 and swap2. Assumes tree, tree->root, swap1 and swap2 are all non-NULL.
|	Polytomies are ok at u or v or both, but assumes that the tree is NOT rooted at either swap1 or swap2.
*/
void TreeManip::NNISwap(
  TreeNode * swap1,		/**< */
  TreeNode * swap2)		/**< */
	{
	assert(tree != NULL);
	assert(tree->GetFirstPreorder() != NULL);
	assert(tree->GetFirstPreorder() != swap1);
	assert(tree->GetFirstPreorder() != swap2);
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
	TreeNode * x = swap1;
	TreeNode * v = swap1->GetParent();
	TreeNode * y = swap2;
	TreeNode * u = swap2->GetParent();

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

	TreeNode * xlast 					= FindLastPreorderInClade(x);
	TreeNode * ylast 					= FindLastPreorderInClade(y);

	TreeNode * xlastnext 				= xlast->nextPreorder;
	TreeNode * ylastnext 				= ylast->nextPreorder;

	TreeNode * u_lChild					= u->lChild;
	TreeNode * v_lChild					= v->lChild;

	TreeNode * x_rSib					= x->rSib;
	TreeNode * y_rSib					= y->rSib;

	TreeNode * x_par					= x->par;
	TreeNode * y_par					= y->par;

	TreeNode * x_prevPreorder			= x->prevPreorder;
	TreeNode * y_prevPreorder			= y->prevPreorder;

	TreeNode * xlast_nextPreorder		= xlast->nextPreorder;
	TreeNode * ylast_nextPreorder		= ylast->nextPreorder;

	TreeNode * xlsib 					= FindLeftSib(x);
	TreeNode * ylsib 					= FindLeftSib(y);

	TreeNode * xlsiblast				= (xlsib == NULL ? NULL : FindLastPreorderInClade(xlsib));
	TreeNode * ylsiblast				= (ylsib == NULL ? NULL : FindLastPreorderInClade(ylsib));

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
	tree->InvalidateTreeID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Removes leftmost child of node `u' and makes it the immediate left sibling of `w'. Used in performing Larget-Simon
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
  TreeNode * u,	/* this node's leftmost child will be removed */
  TreeNode * w)	/* the removed node will become this node's left sibling */
	{
	assert(u				!= NULL);
	assert(u->lChild		!= NULL);	// must have a node to mode
	assert(w				!= NULL);
	assert(w->par			!= NULL);
	assert(w->par->lChild	!= NULL);

	// Make copies of the important pointers
	//
	TreeNode * v					= w->par;
	TreeNode * a					= u->lChild;
	TreeNode * a_last				= FindLastPreorderInClade(a);
	TreeNode * a_rSib				= a->rSib;
	TreeNode * a_prevPre			= a->prevPreorder;
	TreeNode * w_prevPre			= w->prevPreorder;
	TreeNode * w_lSib				= FindLeftSib(w);
	TreeNode * w_lSibLast			= (w_lSib == NULL ? NULL : w_prevPre);
	//TreeNode * v_lChild			= v->lChild;

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
	tree->InvalidateTreeID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Removes rightmost child of node `u' and makes it the immediate right sibling of w. Used in performing Larget-Simon
|	moves when polytomies are present. Does not necessarily leave tree in a consistent state (i.e. tree may contain
|	nodes of order 2 afterwards). Assumes `w' either equals `u' or is above `u' in the tree as it is currently rooted.
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
  TreeNode * u,	/**< this node's rightmost child will be removed */
  TreeNode * w)	/**< the removed node will become this node's right sibling */
	{
	assert(u				!= NULL);
	assert(u->lChild		!= NULL);
	assert(w				!= NULL);
	assert(w->par			!= NULL);
	assert(w->par->lChild	!= NULL);

	// Make copies of the important pointers
	//
	TreeNode * w_par					= w->par;
	TreeNode * a						= FindRightmostChild(u);
	//	TreeNode * a_par				= a->GetParent();
	TreeNode * w_last					= FindLastPreorderInClade(w);
	TreeNode * w_last_nextPre			= w_last->nextPreorder;
	TreeNode * a_last					= FindLastPreorderInClade(a);
	TreeNode * a_last_nextPre			= a_last->nextPreorder;
	TreeNode * w_rSib					= w->rSib;
	TreeNode * a_lSib					= FindLeftSib(a);
	TreeNode * a_lSib_last				= FindLastPreorderInClade(a_lSib);

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
	tree->InvalidateTreeID();
	tree->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Detaches node 's' from tree, keeping intact subtree rooted at `s'. Assumes `s' is non-NULL and is neither the root
|	node nor the subroot node (i.e. only child of the root node). Sets all pointers of `s' to NULL except `lChild' and 
|	`nextPreorder'.
*/
void TreeManip::DetachSubtree(
  TreeNode * s)	/**< is the root of the subtree to be detached */
	{
	assert(s != NULL);
	assert(!s->IsRoot());
	assert(!s->GetParent()->IsRoot());

	// Save pointers to relevant nodes
	//
	TreeNode * s_par				= s->par;
	TreeNode * s_lSib				= FindLeftSib(s);
	TreeNode * s_rSib				= s->rSib;
	TreeNode * s_prevPreorder		= s->prevPreorder;
	TreeNode * slast				= FindLastPreorderInClade(s);
	TreeNode * slast_nextPreorder	= slast->nextPreorder;

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
	TreeNode *u_lSib	= FindLeftSib(u);
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

// *********************************************************************************************************************
// *********************************************************************************************************************
// ************************************ Phycas not yet ready for functions below here **********************************
// *********************************************************************************************************************
// *********************************************************************************************************************

#if 0

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
	assert(rootNode == tree->GetFirstPreorder());

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
	assert(rootNode == tree->GetFirstPreorder());

	tree->RefreshFullID();
	tree->DebugCheckTreeStructure();
	}

#endif
