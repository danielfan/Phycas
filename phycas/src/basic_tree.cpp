/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string>
#include <cstdio>
#include <cctype>
#include "phycas/src/phycas_string.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/xphylogeny.hpp"
//#define TESTING_TOWARD_NODE_ITERATOR
#if defined(TESTING_TOWARD_NODE_ITERATOR)
#	include "phycas/src/edge_iterators.hpp"
#endif
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Push node onto end of `nodeStorage' vector (nodes are stored rather than deleted to save having to reallocate them
|	later.
*/
void Tree::StoreTreeNode(TreeNode * u)
	{
	nodeStorage.push_back(u);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Should only be called from within RerootAtThisTip() or RerootAtThisInternal. Makes node `m' (the mover) the 
|	rightmost child of node `t' (the target). Assumes both `m' and `t' are non-NULL. Because the purpose of this 
|	function is to move nodes below the new root node to the appropriate point above the new root node, this method 
|	also assumes that `t' is always a descendant (perhaps remote) of `m'.
*/
void Tree::RerootHelper(TreeNode * m, TreeNode * t)
	{
	PHYCAS_ASSERT(m != NULL);
	PHYCAS_ASSERT(t != NULL);

	// Save nodes to which m attaches
	TreeNode * m_lChild	= m->lChild;
	TreeNode * m_rSib	= m->rSib;
	TreeNode * m_par		= m->par;

	// Starting with t, walk down tree to identify x, the child of m that is on the path from m to t
	TreeNode * x = t;
	while (x->par != m)
		{
		x = x->par;
		PHYCAS_ASSERT(x != NULL);
		}
	TreeNode * x_rSib = x->rSib;

	// identify x_lSib, the immediate left sibling of x (will be NULL if x is lChild of m)
	TreeNode * x_lSib = NULL;
	if (x != m_lChild)
		{
		x_lSib = m_lChild;
		while (x_lSib->rSib != x)
			{
			x_lSib = x_lSib->rSib;
			PHYCAS_ASSERT(x_lSib != NULL);
			}
		}

	// identify m_lSib, the immediate left sibling of m (will be NULL if m is root node or is lChild of its parent)
	TreeNode *m_lSib = NULL;
	if (m_par != NULL && m != m_par->lChild)
		{
		m_lSib = m_par->lChild;
		while (m_lSib->rSib != m)
			{
			m_lSib = m_lSib->rSib;
			PHYCAS_ASSERT(m_lSib != NULL);
			}
		}

	// Put x where m is now
	preorderDirty = true;
	if (m_par == NULL)
		{
		// m is the root node
		PHYCAS_ASSERT(m_rSib == NULL);
		PHYCAS_ASSERT(m_lSib == NULL);
		x->rSib = NULL;
		x->par = NULL;
		if (x == m_lChild)
			m->lChild = x_rSib;
		else
			x_lSib->rSib = x_rSib;
		}
	else if (m == m_par->lChild)
		{
		// m is leftmost child of its parent
		x->rSib = m_rSib;
		x->par = m_par;
		m->rSib = NULL;
		m->par = NULL;
		m_par->lChild = x;
		if (x == m_lChild)
			m->lChild = x_rSib;
		else
			x_lSib->rSib = x_rSib;
		}
	else
		{
		// m is not leftmost child of its parent
		m_lSib->rSib = x;
		x->rSib = m_rSib;
		x->par = m_par;
		m->rSib = NULL;
		m->par = NULL;
		if (x == m_lChild)
			m->lChild = x_rSib;
		else
			x_lSib->rSib = x_rSib;
		}
	
	// Make m the new rightmost child of t
	m->par = t;
	if (t->lChild == NULL)
		t->lChild = m;
	else
		{
		// Find rightmost child of t
		m_lSib = t->lChild;
		while (m_lSib->rSib != NULL)
			m_lSib = m_lSib->rSib;

		// Make rightmost child of t the left sib of m
		m_lSib->rSib = m;
		}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Reroots the tree at the node specified by the TreeNode pointer `nd', which is assumed to be non-NULL and an internal
|	node. 
*/
void Tree::RerootAtThisInternal(TreeNode * nd)
	{
	//POL: could be merged with RerootAtThisTip, but note that the two use different root test functions: IsTipRoot vs. IsInternalRoot
	PHYCAS_ASSERT(!nd->IsTip());
	PHYCAS_ASSERT(nd != NULL);
	TreeNode * t = nd;
	TreeNode * m = nd->par;
	while (!nd->IsInternalRoot())
		{
		//std::cerr << "In RerootAt: t = " << t->GetNodeName() << ", m = " << m->GetNodeName() << std::endl; //POL-debug
		PHYCAS_ASSERT(nd->HasParent());

		// Begin by swapping the mover's edge length with nd's edge length
		if (HasEdgeLens())
			{
			double tmp_edgelen = m->edgeLen;
			m->edgeLen = nd->edgeLen;
			nd->edgeLen = tmp_edgelen;
			}

		// Make the "mover" node m (along with all of its children except nd)
		// the rightmost child of the "target" node t
		RerootHelper(m, t);

		// The next target node is always the previous mover, and the next mover node is always nd's parent
		t = m;
		m = nd->par;
		}

	firstPreorder = nd;
	nd->prevPreorder = NULL;

	//std::cerr << "\n\nWalking tree just before leaving RerootAt...\n"; //POL-debug
	//std::cerr << DebugWalkTree(true, 2) << std::endl; //POL-debug
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reroots the tree at the node specified by the TreeNode pointer `nd', which is assumed to be non-NULL and a tip node.
*/
void Tree::RerootAtThisTip(TreeNode * nd)
	{
	//POL: could be merged with RerootAtThisInternal, but note that the two use different root test functions: IsTipRoot vs. IsInternalRoot
	PHYCAS_ASSERT(nd->IsTip());
	PHYCAS_ASSERT(nd != NULL);
	TreeNode * t = nd;
	TreeNode * m = nd->par;
	while (!nd->IsTipRoot())
		{
		//std::cerr << "In RerootAt: t = " << t->GetNodeName() << ", m = " << m->GetNodeName() << std::endl; //POL-debug
		PHYCAS_ASSERT(nd->HasParent());

		// Begin by swapping the mover's edge length with nd's edge length
		if (HasEdgeLens())
			{
			double tmp_edgelen = m->edgeLen;
			m->edgeLen = nd->edgeLen;
			nd->edgeLen = tmp_edgelen;
			}

		// Make the "mover" node m (along with all of its children except nd)
		// the rightmost child of the "target" node t
		RerootHelper(m, t);

		// The next target node is always the previous mover, and the next mover node is always nd's parent
		t = m;
		m = nd->par;
		}

	firstPreorder = nd;
	nd->prevPreorder = NULL;

	//std::cerr << "\n\nWalking tree just before leaving RerootAt...\n"; //POL-debug
	//std::cerr << DebugWalkTree(true, 2) << std::endl; //POL-debug
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Searches for tip (degree = 1) node whose `nodeNum' data member equals `num', then reroots the tree at that node. If
|	a node by that number cannot be found, throws an XPhylogeny exception.
*/
void Tree::RerootAtTip(unsigned num)
	{
	//std::cerr << "*** in RerootAtTip ***" << std::endl; 

	TreeNode *nd = FindTipNode(num);
	if (nd == NULL)
		{
		workspace.clear();
		workspace << "there is no tip node having number ";
		throw XPhylogeny(append_unsigned(workspace, num));
		}
	else
		RerootAtThisTip(nd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes preorder pointers (`Tree::firstPreorder', `Tree::lastPreorder', `TreeNode::nextPreorder', 
|	`TreeNode::prevPreorder') by starting with the root node (`firstPreorder') and walking through the tree by visiting
|	`lChild', `rSib' and `par' pointers. Sets `preorderDirty' to false when done.
*/
void Tree::RefreshPreorder(TreeNode *nd)
	{
	if (!preorderDirty)
		return;

	if (nd == NULL)
		nd = firstPreorder;

	PHYCAS_ASSERT(nd != NULL);

	if (nd->lChild == NULL && nd->rSib == NULL)
		{
		// nd has no children and no siblings, so nextPreorder is the right sibling of 
		// the first ancestral node that has a right sibling.
		TreeNode *anc = nd->par;
		while ((anc != NULL) && (anc->rSib == NULL))
			anc = anc->par;
		if (anc == NULL)
			{
			// We descended all the way to the root without finding an ancestor with
			// a right sibling, so nd must be the upper-right-most node in the tree
			lastPreorder = nd;
			nd->nextPreorder = NULL;
			}
		else
			{
			// We found an ancestor with a right sibling without having to go all
			// the way to the root of the tree
			nd->nextPreorder = anc->rSib;
			anc->rSib->prevPreorder = nd;
			}
		}
	else if (nd->lChild == NULL && nd->rSib != NULL)
		{
		// nd has no children (it is a tip), but does have a sibling on its right
		nd->nextPreorder = nd->rSib;
		nd->rSib->prevPreorder = nd;
		}
	else if (nd->lChild != NULL && nd->rSib == NULL)
		{
		// nd has children (it is an internal node) but no siblings on its right
		nd->nextPreorder = nd->lChild;
		nd->lChild->prevPreorder = nd;
		}
	else
		{
		// nd has both children and siblings on its right
		nd->nextPreorder = nd->lChild;
		nd->lChild->prevPreorder = nd;
		}

	if (nd->lChild != NULL)
		{
		RefreshPreorder(nd->lChild);
		}

	if (nd->rSib != NULL)
		{
		RefreshPreorder(nd->rSib);
		}

	if (nd == firstPreorder)
		preorderDirty = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The supplied `name_vector' should be a vector of observable tip names. Assumes that observable nodes in the tree 
|	already have names, and the goal is to number them appropriately (i.e. set their `nodeNum' to the position of their 
|	name in the supplied `name_vector'. If an observable node is encountered that has no name, or if the number of 
|	observable nodes is not equal to the length of `name_vector', an XPhylogeny exception is thrown. Depends on the
|	data member `observable' having been set correctly for all nodes in the tree.
*/
void Tree::RectifyNumbers(std::vector<std::string> name_vector)
	{
	unsigned name_vector_len = (unsigned)name_vector.size();
	bool ok = (GetNObservables() == name_vector_len);
	if (!ok)
		throw XPhylogeny("number of names must equal number of observable tips");

	if (preorderDirty)
		RefreshPreorder();

	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsObservable())
			{
			std::vector<std::string>::iterator it = std::find(name_vector.begin(), name_vector.end(), nd->GetNodeName());
			if (it == name_vector.end())
				{
				std::string s = "node name (";
				s << nd->GetNodeName();
				s << ") not found in supplied list of names";
				throw XPhylogeny(s);
				}
			else
				{
				unsigned n = (unsigned)(it - name_vector.begin());
				nd->SetNodeNum(n);
				}
			} 
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The supplied `name_vector' should be a vector of tip names. Assumes that the tip nodes have numbers but that their 
|	names should be set according to the values supplied in `name_vector'; that is, a node whose number is k will have
|	its name set to `name_vector'[k]. Any tip having the default number TreeNode::nodeNumInitValue will cause an 
|	XPhylogeny exception to be thrown. An XPhylogeny exception will also be thrown if the number of tips in the tree 
|	does not match the length of `name_vector'.
*/
void Tree::RectifyNames(std::vector<std::string> name_vector)
	{
	//std::cerr << "*** in RectifyNames ***" << std::endl; 

	unsigned name_vector_len = (unsigned)name_vector.size();
	bool ok = (GetNObservables() == name_vector_len);
	if (!ok)
		throw XPhylogeny(str(boost::format("number of names (%d) not equal to the number of tips (%d)") % name_vector_len % GetNObservables()));

	if (preorderDirty)
		RefreshPreorder();

	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsObservable())
			{
			unsigned n = nd->GetNodeNumber();
			if (n >= 0 && n < name_vector_len)
				{
				nd->SetNodeName(name_vector[n]);
				}
			else
				{
				std::string s = "node number (";
				s << n;
				s << ") out of range";
				throw XPhylogeny(s);
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `selected' data member to false for all nodes in the tree.
*/
void Tree::UnselectAllNodes()
	{
	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		nd->UnselectNode();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `selected' data member to true for all nodes in the tree.
*/
void Tree::SelectAllNodes()
	{
	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		nd->SelectNode();
		}
	}

#if 0 // including this function results in this VC error: fatal error C1055: compiler limit : out of keys
/*----------------------------------------------------------------------------------------------------------------------
|	Returns pointer to TreeNode representing the deepest point of intersection between the lineage leading to `tip1' and
|	the lineage leading to `tip2'. Both `tip1' and `tip2' are assumed to be tip nodes, but one or the other could
|	currently be serving as the root. If either `tip1' or `tip2' is serving as the root, then the node returned is the 
|	immediate descendant of this node. Otherwise, the chain of parents of both tips is walked until the other tip shows
|	up as a member of the split of the other.
*/
TreeNode * Tree::FindMRCA(unsigned tip1, unsigned tip2)
	{
	assert(0);	//POL not yet written
	return NULL;
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Zig-zags up the right side of the clade starting with the node `start' until reaching a tip. This tip node is the 
|	last node in the clade according to the preorder sequence. Zig-zagging is necessary because there is no pointer to
|	the right-most child of a node, so we must use FindRightmostChild to move up the clade. Assumes `start' is non-NULL.
*/
TreeNode * Tree::FindLastPreorderInClade(
  TreeNode * start)	/**< is the deepest node belonging to the clade in question */
	{
	//@POL This function was stolen from the TreeManip class - makes more sense here
	PHYCAS_ASSERT(start != NULL);
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
|	Begins with left child of parent of `start' and calls GetRightSib() until the left sibling of `start' is located.
|	Assumes `start' is non-NULL.
*/
TreeNode * Tree::FindLeftSib(
  TreeNode * start)	/**< is the node whose left sibling is sought */
	{
	//@POL This function was stolen from the TreeManip class - makes more sense here
	PHYCAS_ASSERT(start != NULL);
	TreeNode * nd = start->GetParent();

	// If start has no parent, then there's no way it can have a left sibling
	if (nd == NULL)
		return NULL;	

	nd = nd->GetLeftChild();

	// Parent of start should have at least one child
	PHYCAS_ASSERT(nd != NULL); 

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
TreeNode * Tree::FindRightmostChild(
  TreeNode * start)	/**< is the parent node whose children will be searched */
	{
	//@POL This function was stolen from the TreeManip class - makes more sense here
	PHYCAS_ASSERT(start != NULL);
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
|	Inserts node 's' as child of 'u' keeping the subtree rooted at `s' intact. Assumes both `s' and `u' are non-NULL. If 
|	`targetSib' is specified, `s' will become an immediate sibling of `targetSib' (right or left sibling determined by 
|	`on_right'). If `targetSib' is not specified, `s' will be added as the rightmost or leftmost child of `u' (again,
|	depending on the value of `on_right').
*/
void Tree::InsertSubtree(
  TreeNode * s,			/**< is the node to be added as a child of u */
  TreeNode * u,			/**< is the new parent of s */
  bool on_right,		/**< if true, `s' will be added to the right of `targetSib' or, if` targetSib' is not specified, `s' will become rightmost child of `u' (vice versa if false specified) */
  TreeNode * targetSib)	/**< if non-NULL, `s' will become an immediate sibling of this node (right or left depending on value of `on_right') */
	{
	//@POL This function was stolen from the TreeManip class - makes more sense here
	PHYCAS_ASSERT(u != NULL);
	PHYCAS_ASSERT(s != NULL);
	PHYCAS_ASSERT(targetSib == NULL || (targetSib->par == u));

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
		if ((!on_right) && (targetSib == u_lChild))
			targetSib = NULL;
		if ((on_right) && (targetSib == u_rChild))
			targetSib = NULL;
		if ((!on_right) && (targetSib != NULL))
			{
			targetSib = FindLeftSib(targetSib);
			on_right = true;
			}
		}

	// Make s a child of u
	//
	if (targetSib == NULL)
		{
		if (on_right)
			{
			s->rSib				= NULL;
			s->par				= u;
			s->prevPreorder 	= ulast;
			slast->nextPreorder = ulast_nextPreorder;
			ulast->nextPreorder	= s;

			if (ulast_nextPreorder == NULL)
				this->lastPreorder = slast;
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

		// Can assume that targetSib is not rightmost child of u and also that on_right == true
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
	this->InvalidateTreeID();
	this->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Detaches node 's' from tree, keeping intact subtree rooted at `s'. Assumes `s' is non-NULL and is neither the root
|	node nor the subroot node (i.e. only child of the root node). Sets all pointers of `s' to NULL except `lChild' and 
|	`nextPreorder'.
*/
void Tree::DetachSubtree(
  TreeNode * s)	/**< is the root of the subtree to be detached */
	{
	//@POL This function was stolen from the TreeManip class - makes more sense here
	PHYCAS_ASSERT(s != NULL);
	PHYCAS_ASSERT(!s->IsTipRoot());
	PHYCAS_ASSERT(!s->GetParent()->IsTipRoot());
	if (s->GetParent()->IsTipRoot())
		{
		std::cerr << "oops" << std::endl;
		}

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
		this->firstPreorder = slast_nextPreorder;
	else
		s_prevPreorder->nextPreorder = slast_nextPreorder;

	slast->nextPreorder = NULL;
	if (slast_nextPreorder == NULL)
		this->lastPreorder = s_prevPreorder;
	else
		slast_nextPreorder->prevPreorder = s_prevPreorder;

	// This rearrangement invalidates the TreeID and node counts
	//
	this->InvalidateID();
	this->InvalidateNodeCounts();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs a postorder traversal and, for each internal node visited, sorts clades from smallest to largest. If the
|	`right' argument is true, then the largest clades will be on the right after ladderization. If the `right' argument
|	is false, then the leftmost child will be the largest clade.
*/
void Tree::Ladderize(
  bool right)	/**> if true, largest clades will be on right; if false, largest clades will be on left */
	{
	for (TreeNode * nd = GetLastPreorder(); nd != NULL; nd = nd->GetNextPostorder())
		{
		bool istip = nd->IsTip();
		TreeNode * par = nd->GetParent();
		TreeNode * rsib = nd->GetRightSib();

		if (par)
			{
			if (istip)
				{
				if (rsib)
					par->tmp += 1.0;
				else
					par->tmp = 1.0;
				}
			else
				{
				if (rsib)
					par->tmp += nd->tmp;
				else
					par->tmp = nd->tmp;
				}
			}
		
		// If node is an internal node, then reorder its children from smallest (left) to
		// largest (right) based on the tmp data member of each child
		if (!istip)
			{
			typedef std::multimap< unsigned, TreeNode * > ChildMap;
			typedef std::pair< unsigned, TreeNode * > ChildMapKeyValuePair;
			ChildMap childmap;
			TreeNode * child = nd->GetLeftChild();
			while (child != NULL)
				{
				unsigned clade_size = (unsigned)child->tmp;
				childmap.insert(ChildMapKeyValuePair(clade_size, child));
				TreeNode * next_child = child->GetRightSib();
				DetachSubtree(child);
				child = next_child;
				}

			// The childmap is ordered from smallest to largest already, so we need only 
			// to march down the map adding nodes as we go to the right of the last one
			if (right)
				{
				for (ChildMap::iterator childiter = childmap.begin(); childiter != childmap.end(); ++childiter)
					{
					InsertSubtree(childiter->second, nd, true);
					}
				}
			else
				{
				for (ChildMap::reverse_iterator childiter = childmap.rbegin(); childiter != childmap.rend(); ++childiter)
					{
					InsertSubtree(childiter->second, nd, true);
					}
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Follows preorder pointers until a tip (degree = 1) node is encountered with `nodeNum' equal to `num'. Returns a 
|	pointer to this node, or 0 if a tip node numbered `num' cannot be not found.
*/
TreeNode *Tree::FindTipNode(unsigned num)
	{
	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip() && num == nd->GetNodeNumber())
			return nd;
		}

	return 0;
	}

#if defined(TESTING_TOWARD_NODE_ITERATOR)
bool isValid(const TreeNode * refNd, const TreeNode * neighborCloserToEffectiveRoot)
	{
	return refNd->IsTip();
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Computes and returns the longest path from the root to any leaf node in the tree.
*/
double Tree::calcTotalHeight()
	{
	//@POL: this function could be made more general (i.e. right now it depends on the root being a tip)
	//@POL: this function could be made more efficient
	PHYCAS_ASSERT(GetFirstPreorder() != NULL);
	PHYCAS_ASSERT(GetFirstPreorder()->GetNextPreorder() != NULL);

	// Create a map (node number -> height) that will keep track of the maximum height above each internal node
	typedef std::map<unsigned, double> HeightMap;
	HeightMap hmap;

	// Walk through nodes in postorder fashion, for each internal node adding an entry to hmap
	for (TreeNode * nd = GetLastPreorder(); nd != NULL; nd = nd->GetNextPostorder())
		{
		if (nd->IsInternal())
			{
			// Figure out which path above this node is greatest and insert that value into hmap
			double max_child_height = 0.0;
			for (TreeNode * child = nd->GetLeftChild(); child != NULL; child = child->GetRightSib())
				{
				// height for child is really child's edge length plus the maximum height above child,
				// which is stored in hmap if child is an internal node
				double h = (hasEdgeLens ? child->GetEdgeLen() : 1.0);
				unsigned child_nodenum = child->GetNodeNumber();
				if (child->IsInternal())
					h += hmap[child_nodenum];
				if (h > max_child_height)
					max_child_height = h;
				}
			unsigned nd_nodenum = nd->GetNodeNumber();
			hmap.insert(HeightMap::value_type(nd_nodenum,max_child_height));
			}
		}

	// Maximum tree height is now equal to the maximum height above the subroot plus the subroot's edge length
	TreeNode * subroot = GetFirstPreorder()->GetNextPreorder();
	unsigned subroot_nodenum = subroot->GetNodeNumber();
	double subroot_height = (hasEdgeLens ? subroot->GetEdgeLen() : 1.0);
	double height_above_subroot = hmap[subroot_nodenum];
	return subroot_height + height_above_subroot;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Walks tree using preorder pointers, returning a string documenting the journey. For this tree description:
|	"(c:0.4,d:0.5,(a:0.1, b:0.2):0.3);", the returned string would be: "[1] -> c (0) -> d (1) -> [0] -> a (2) -> b (3)"
|	assuming that the `preorder' parameter was true. Note that node names are output if they were supplied, the internal
|	node numbers are inside square brackets, and the tip node numbers are inside parentheses. The tip nodes are, by 
|	default, numbered according to the preorder sequence (and the internal nodes are numbered in postorder sequence).
|	If a vector of tip names was supplied using RectifyNames(), then the tip node numbers will be equal to the 0-offset
|	index of the tip into that vector. This function is intended to be used for debugging purposes and is thus not 
|	designed to be particularly fast. Verbosity levels: 0 means minimal (names for named nodes, numbers otherwise);
|	1 means add node numbers and indicate using brackets or parentheses whether node is an [internal], (tip) or {root}
|	node; 2 means show nearly everything about each node (each node requires a paragraph to explain).
*/
std::string Tree::DebugWalkTree(bool preorder, unsigned verbosity)
	{
	std::string s;
#	if defined(TESTING_TOWARD_NODE_ITERATOR)
	NodeValidityChecker validFunctor = isValid;
	for (TreeNode * z = GetFirstPreorder(); z != NULL; z = z->GetNextPreorder())
		{
		if (z->IsTip())
			continue;
		s << "\nStarting with " << z->nodeName << "...\n";
		effective_postorder_edge_iterator i(z, validFunctor);
		effective_postorder_edge_iterator e;
		for (; i != e; ++i)
#	else
	TreeNode *nd = (preorder ? GetFirstPreorder() : GetLastPreorder());
	for (; nd != NULL; nd = (preorder ? nd->GetNextPreorder() : nd->GetNextPostorder()))
#	endif	
		{
#		if defined(TESTING_TOWARD_NODE_ITERATOR)
			TreeNode * nd = i->first;
#		endif
		if (verbosity == 2)
			{
			s << "\n";
			nd->AppendNodeInfo(s);
			}
		else
			{
			// Gather information about node name and number
			const std::string &nm = nd->GetNodeName();
			std::string tmpstr = str(boost::format("%d") % nd->GetNodeNumber());

			// Output the node name
			if (verbosity == 0)
				{
				if (nm.length() > 0)
					{
					// output node name
					s << nm;
					}
				else
					{
					// output node number
					if (nd->IsTip())
						{
						s << "(";
						s << tmpstr;
						s << ")";
						}
					else
						{
						s << "[";
						s << tmpstr;
						s << "]";
						}
					}
				}
			else // verbosity == 1
				{
				if (nm.length() > 0)
					{
					// output node name
					s << nm;
					}
				else
					{
					s << "?";
					}

				// Output the node number
				if (nd->IsTip())
					{
					s << " (";
					s << tmpstr;
					s << ")";
					}
				else
					{
					s << " [";
					s << tmpstr;
					s << "]";
					}
				}

			// Output an arrow if not the last node in the series
			bool show_arrow = false;
			if (preorder)
				show_arrow = (nd != lastPreorder);
			else
				show_arrow = (nd != firstPreorder);
			if (show_arrow)
				s << " -> ";
			}
		}
#	if defined(TESTING_TOWARD_NODE_ITERATOR)
	}
#	endif

	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the sum of all edge lengths in the tree.  Uses preorder pointers to walk through tree, calling 
|	RefreshPreorder() if these pointers have never been set or have been invalidated.
*/
double Tree::EdgeLenSum()
	{
	if (!HasEdgeLens())
		throw XPhylogeny("no edge lengths were specified for this tree");

	//std::cerr << "\nEntering EdgeLenSum()" << std::endl; //POL-debug

	// todo: use algorithm accumulate for this, which entails writing node iterator
	//       this would make this function inlineable, I think
	if (preorderDirty)
		RefreshPreorder();

	TreeNode *nd = GetFirstPreorder();
	double sum = 0.0;

	for (nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		//std::cerr << "nd->GetNodeNumber() = " << nd->GetNodeNumber() << ", sum = " << sum << std::endl; //POL-debug
		sum += nd->GetEdgeLen();
		}
	//std::cerr << "sum = " << sum << std::endl; //POL-debug

	//std::cerr << "\nLeaving EdgeLenSum()" << std::endl; //POL-debug
	return sum;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `nTips' and `nInternals' to their correct values. Uses preorder pointers to walk through tree, calling 
|	RefreshPreorder() if these pointers have never been set or have been invalidated. 
*/
void Tree::RefreshNodeCounts()
	{
	nTips		= 0;
	nInternals	= 0;
	for (TreeNode *nd = GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			++nTips;
		else
			++nInternals;
		}

	nodeCountsValid = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructs a tree from a string containing a "Newick" style (nested parentheses) tree description. If `newick' is
|	not a valid tree description, throws an XPhylogeny exception. Assumes that if any tip node names in the tree
|	description are convertible to integers, then tip node numbers should be set according to the converted node names.
|	Thus, it is considered an error to mix names that are not convertible to integers with names that are convertible,
|	and an XPhylogeny exception will be thrown if such a mixture is detected, or if a node number ends up being 
|	duplicated or out of range.
*/
void Tree::BuildFromString(const std::string &newick)
	{
	PHYCAS_ASSERT(newick.length() > 0);

	// Assume tip node names represent the 0-offset index of the taxon represented by the tip in the data matrix
	// and thus that tip node numbers should be set to the tip node name after converting it to an integer. This
	// saves having to call RectifyNumbers() for the common case of tree descriptions in which taxon names have
	// been replaced by their indices.
	//
	bool tip_numbers_equal_names = true;
	std::set<unsigned> numbers_used;

	//std::cerr << "About to build this tree: " << newick << std::endl; //POL-debug

	// Stow away any existing nodes, reset nTips and nInternals to 0, and reset data members
	//
	Clear();

	// Save first tip node and reroot at this node before returning
	TreeNode *first_tip = NULL;

	try
		{
		// At the end, we will check to make sure that if *any* edge lengths were specified then
		// *all* edge lengths were specified
		unsigned nEdgeLengths = 0;

		// Create a root node
		TreeNode *nd = GetNewNode();
		firstPreorder = nd;

		// Some flags to keep track of what we did last
		enum {
			Prev_Tok_LParen		= 0x01,	/**< previous token was a left parenthesis ('(') */
			Prev_Tok_RParen		= 0x02,	/**< previous token was a right parenthesis (')') */
			Prev_Tok_Colon		= 0x04,	/**< previous token was a colon (':') */
			Prev_Tok_Comma		= 0x08,	/**< previous token was a comma (',') */
			Prev_Tok_Name		= 0x10,	/**< previous token was a node name (e.g. '2', 'P._articulata') */
			Prev_Tok_EdgeLen	= 0x20	/**< previous token was an edge length (e.g. '0.1', '1.7e-3') */
			};
		unsigned previous = Prev_Tok_LParen;

		//std::cerr << "Prev_Tok_LParen  =" << Prev_Tok_LParen << std::endl;	// decimal 1 //POL-debug
		//std::cerr << "Prev_Tok_RParen  =" << Prev_Tok_RParen << std::endl;	// decimal 2 //POL-debug
		//std::cerr << "Prev_Tok_Colon   =" << Prev_Tok_Colon << std::endl;		// decimal 4 //POL-debug
		//std::cerr << "Prev_Tok_Comma   =" << Prev_Tok_Comma << std::endl;		// decimal 8 //POL-debug
		//std::cerr << "Prev_Tok_Name    =" << Prev_Tok_Name << std::endl;		// decimal 16 //POL-debug
		//std::cerr << "Prev_Tok_EdgeLen =" << Prev_Tok_EdgeLen << std::endl;	// decimal 32 //POL-debug

		// Some useful flag combinations 
		unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
		unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
		unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
		unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
		unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);

		// loop through the characters in newick, building up tree as we go
		std::string::const_iterator i = newick.begin();
		for (; i != newick.end(); ++i)
			{
			//std::cerr << (i == newick.end() ? "### how can we be here ###" : "### not yet at end ###") << std::endl; //POL-debug
			char ch = *i;
			//std::cerr << "[" << ch << "]" << std::endl; //POL-debug
			if (iswspace(ch))
				continue;
			switch(ch)
				{
				case ';':
					//throw XPhylogeny("premature semicolon");
					break;  

				case ')':
					// If nd is bottommost node, expecting left paren or semicolon, but not right paren
					if (nd->NoParent())
						throw XPhylogeny("too many right parentheses");
					// Expect right paren only after an edge length, a node name, or another right paren
					if (!(previous & RParen_Valid))
						throw XPhylogeny("unexpected right parenthesis");
					// Right paren means we are leaving a node for its parent, so set node number now if internal node
					if (nd->IsInternal() && nd->NumberNotYetAssigned())
						{
						nd->nodeNum = nInternals++;
						//std::cerr << "1> nInternals now " << nInternals << std::endl; //POL-debug
						}
					// Go down a level
					nd = nd->GetParent();
					if (nd->GetLeftChild()->GetRightSib() == NULL)
						throw XPhylogeny("internal node has only one child");
					previous = Prev_Tok_RParen;
					break;
				case ':':
					// Expect colon only after a node name or another right paren
					if (!(previous & Colon_Valid))
						throw XPhylogeny("unexpected colon");
					previous = Prev_Tok_Colon;
					break;

				case ',':
					// Expect comma only after an edge length, a node name, or a right paren
					if (nd->NoParent() || !(previous & Comma_Valid))
						throw XPhylogeny("unexpected comma");
					// Comma means we are leaving a node to create its sibling, so set node number now if internal
					if (nd->IsInternal() && nd->NumberNotYetAssigned())
						{
						nd->nodeNum = nInternals++;
						//std::cerr << "2> nInternals now " << nInternals << std::endl; //POL-debug
						}
					// Create the sibling
					nd->rSib = GetNewNode();
					nd->rSib->par = nd->par;
					nd = nd->rSib;
					previous = Prev_Tok_Comma;
					break;

				case '(':
					// Expect left paren only after a comma or another left paren
					if (!(previous & LParen_Valid))
						throw XPhylogeny("unexpected left parenthesis");
					// Create new_node above and to the left of the current node
					PHYCAS_ASSERT(nd->lChild == NULL);
					nd->lChild = GetNewNode();
					nd->lChild->par = nd;
					nd = nd->lChild;
					previous = Prev_Tok_LParen;
					break;

				case '\'':
					// Encountered an apostrophe, which always indicates the start of a 
					// node name (but note that node names do not have to be quoted)

					// Expect node name only after a left paren (child's name), a comma (sib's name)
					// or a right paren (parent's name)
					if (!(previous & Name_Valid))
						{
						throw XPhylogeny("1> unexpected node name");
						}

					// Get the rest of the name
					nd->nodeName.clear();
					for (++i; i != newick.end(); ++i)
						{
						ch = *i;
						if (ch == '\'')
							break;
						else if (iswspace(ch))
							nd->nodeName << ' ';
						else
							nd->nodeName << ch;
						}
					if (ch != '\'')
						{
						workspace.clear();
						workspace << "expecting single quote to mark the end of node name (";
						workspace << abbreviate(nd->GetNodeName(), 20);
						workspace << ")";
						throw XPhylogeny(workspace);
						}

					//std::cerr << "1> node name = " << nd->nodeName << std::endl; //POL-debug
					if (nd->IsTip())
						{
						if (first_tip == NULL)
							first_tip = nd;
						if (tip_numbers_equal_names)
							{
							bool ok = SetNumberFromName(nd, numbers_used);	// may throw XPhylogeny
							if (!ok)
								{
								nd->nodeNum = nTips;
								tip_numbers_equal_names = false;
								}
							}
						else
							nd->nodeNum = nTips;
						++nTips;
						nd->SetObservable(true);
						//std::cerr << "1> nTips now " << nTips << std::endl; //POL-debug
						}

					previous = Prev_Tok_Name;
					break;

				default:
					// Expecting either an edge length or an unquoted node name
					if (previous == Prev_Tok_Colon)
						{
						// Edge length expected (e.g. "235", "0.12345", "1.7e-3")
						workspace.clear();
						for (; i != newick.end(); ++i)
							{
							ch = *i;
							if (ch == ',' || ch == ')' || iswspace(ch))
								{
								--i;
								break;
								}
							bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
							if (!valid)
								throw XPhylogeny("invalid branch length character");
							workspace << ch;
							}
						//std::cerr << "edge length string = " << workspace << std::endl; //POL-debug
						double brlen = atof(workspace.c_str());
						//std::cerr << "edge length = " << brlen << std::endl; //POL-debug
						nd->SetEdgeLen(brlen);

						++nEdgeLengths;
						previous = Prev_Tok_EdgeLen;
						}
					else
						{
						// Get the node name
						nd->nodeName.clear();
						//std::cerr << "\nBuilding up node name:" << std::endl; //POL-debug
						for (; i != newick.end(); ++i)
							{
							ch = *i;
							//std::cerr << "\n  ch = " << ch << std::endl; //POL-debug
							if (ch == '(')
								throw XPhylogeny("unexpected left parenthesis inside node name");
							if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')')
								{
								--i;
								//std::cerr << "\n  backing up, *i = " << (*i) << std::endl; //POL-debug
								break;
								}
							nd->nodeName << ch;
							//std::cerr << "\n  nd->nodeName = " << nd->nodeName << std::endl; //POL-debug
							}
						//std::cerr << "2> node name = " << nd->nodeName << std::endl; //POL-debug
						//std::cerr << (i == newick.end() ? "*** at end of newick ***" : "*** not yet at end ***") << std::endl; //POL-debug

						// Expect node name only after a left paren (child's name), a comma (sib's name)
						// or a right paren (parent's name)
						if (!(previous & Name_Valid))
							{
							//std::cerr << "previous =" << previous << ", ch =" << ch << std::endl; //POL-debug
							workspace.clear();
							workspace << "unexpected node name (";
							workspace << nd->nodeName;
							workspace << ")";
							throw XPhylogeny(workspace);
							}

						if (nd->IsTip())
							{
							if (first_tip == NULL)
								first_tip = nd;
							if (tip_numbers_equal_names)
								{
								bool ok = SetNumberFromName(nd, numbers_used);	// may throw XPhylogeny
								if (!ok)
									{
									nd->nodeNum = nTips;
									tip_numbers_equal_names = false;
									}
								}
							else
								nd->nodeNum = nTips;
							++nTips;
							nd->SetObservable(true);
							//std::cerr << "2> nTips now " << nTips << std::endl; //POL-debug
							}

						previous = Prev_Tok_Name;
						}
				}
				if (i == newick.end())
					{
					//std::cerr << "@@@ At end of newick string, breaking main loop @@@" << std::endl; //POL-debug
					break;
					}

			}

		if (firstPreorder->NumberNotYetAssigned())
			{
			firstPreorder->nodeNum = nInternals++;
			//std::cerr << "3> nInternals now " << nInternals << std::endl; //POL-debug
			}
		//std::cerr << "4> nInternals now " << nInternals << std::endl; //POL-debug
		//std::cerr << "4> nTips now " << nTips << std::endl; //POL-debug
		firstPreorder->SetEdgeLen(0.0);
		if (nEdgeLengths > 0)
			{
			if (nEdgeLengths < nTips + nInternals - 1)
				throw XPhylogeny("some but not all edge lengths were specified");
			hasEdgeLens = true;
			}
		else
			hasEdgeLens = false;
		if (!nd->NoParent())
			{
			//std::cerr << "Why are we here? nd->NoParent() = " << (nd->NoParent() ? "true" : "false") << std::endl; //POL-debug
			throw XPhylogeny("too many left parentheses");
			}
		//std::cerr << "About to call RerootAtThisTip..." << std::endl; //POL-debug
		RerootAtThisTip(first_tip);
		}
	catch(XPhylogeny x)
		{
		//std::cerr << "***" << std::endl; //POL-debug
		Clear();
		throw x;
		}

	if (tip_numbers_equal_names)
		{
		PHYCAS_ASSERT(!numbers_used.empty());
//		if (numbers_used.find(0) == numbers_used.end())
//			{
//			// 0 has not been used as a tip node number, so assume that tip node numbers started at 1
			for (preorder_iterator nd = begin(); nd != end(); ++nd)
				{
				if (nd->IsTip())
					{
					nd->nodeNum = nd->nodeNum - 1;
					PHYCAS_ASSERT(nd->nodeNum >= 0);
					}
				}
//			}
		}

	// Renumber internal nodes so that they do not overlap with numbers of tip nodes
	//@POL should use for_each
	//std::cerr << "ADDING " << nTips << " to each internal node number" << std::endl; //POL temp		
// 	for (preorder_iterator nd = begin(); nd != end(); ++nd)
// 		{
// 		if (nd->IsInternal())
// 			{
// 			//std::cerr << "RENUMBERING internal node " << nd->nodeNum << " (new number is " << (nd->nodeNum + nTips) << ")" << std::endl; //POL temp		
// 			nd->nodeNum = nd->nodeNum + nTips;
// 			}
// 		}
	

	//@POL: should check to make sure tip node numbers form a sequence between 0 and ntax-1
	// If not, then it may mean that the tree descriptions left out some taxa from the data matrix
	// which is a situation that I need to deal with

	//std::cerr << "Leaving BuildFromString: nodeStorage now has " << (unsigned)nodeStorage.size() << " nodes" << std::endl; //POL-debug
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a parenthetical description of the tree (Newick format). If the tree has edge lengths, these are included in
|	the Newick description. Any nodes having non-empty node names are named in the Newick description. If a tip node
|	having node number x lacks a name, the "name" used for it in the tree description is x + 1.
*/
std::string & Tree::AppendNewick(
  std::string & s)  /**< */
	{
	if (preorderDirty)
		RefreshPreorder();

	// Create a tool for formatting floating point values
	const unsigned floatWidth = UINT_MAX;	// make output width wide enough to hold the formatted value
	const unsigned floatPrec = 5;			// number of digits following decimal point
	DoubleFormatter df(floatWidth, floatPrec); 

	bool showEdgeLens = HasEdgeLens();
	bool rooted = IsRooted();
	PHYCAS_ASSERT(!rooted); //@ this assert should go away as soon as possible

	s << '(';
	unsigned openParens = 1;

	// Start with root node, which may actually represent an extant tip (unrooted trees are 
	// rooted at one of the tips)
	TreeNode *nd = GetFirstPreorder();
	PHYCAS_ASSERT(nd->IsTipRoot());

	if (!rooted)
		{
		// If the root node has a name, use it; otherwise convert the node number to a string
		// Note that neither name nor number should be output for root node if the tree is rooted
		if (nd->nodeName.empty())
			append_unsigned(s, nd->GetNodeNumber() + 1);
		else
			s << nd->nodeName;
		}

	if (showEdgeLens && !rooted)
		{
		// note that in an unrooted tree, where the root node is actually an upside-down extant tip,
		// the root node's edge length is actually held by its only child
		s << ':';
		TreeNode *onlyChild = nd->GetLeftChild();
		PHYCAS_ASSERT(onlyChild != NULL);
		PHYCAS_ASSERT(onlyChild->rSib == NULL);
		df.FormatDouble(s, onlyChild->GetEdgeLen());
		}

	nd = nd->GetNextPreorder();
	if (nd == NULL) 
		{
		// Presumably this situation (a one-node tree) will never happen, but then it is 
		// not good to presume anything
		return s << ')'; 
		}

	s << ',';
	if (nd->GetNextPreorder() == NULL)
		{
		// Two taxon tree, presumably this won't happen, but it is easy to deal with
		if (nd->nodeName.empty())
			append_unsigned(s, nd->GetNodeNumber() + 1);
		else
			s << nd->nodeName;
		if (showEdgeLens)  
			{
			s << ":";
			// 0.0 branch length because we gave the root's child's edge length 
			// to the root because we print as a basal polytomy (not rooted at a leaf)
			df.FormatDouble(s, 0.0);
			}
		return s << ")";
		}

	// Note: no left parenthesis here; this is how we get the apparent polytomy at the root
	nd = nd->GetNextPreorder();

	// If we trip this we have a two taxon tree with 3 nodes (degree 2 node in the middle)
	PHYCAS_ASSERT(nd->GetNextPreorder());

	// Now visit all other nodes
	for (; nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			{
			// Output taxon names for tips if they are available, otherwise output tip node number plus 1
			if (nd->nodeName.empty())
				append_unsigned(s, nd->GetNodeNumber() + 1);
			else
				s << nd->nodeName;

			// Output edge length if they were provided
			if (showEdgeLens)
				{
				s << ':';
				df.FormatDouble(s, nd->GetEdgeLen());
				}

			if (nd->GetRightSib() != NULL)
				s << ',';
			else if (nd->nextPreorder != NULL)
				{
				// Go down until we find the left sib of the next preorder node,
				// outputting edge lengths as we go
				const TreeNode *tempAncNode = nd;
				do
					{
					s << ')';
					--openParens;
					PHYCAS_ASSERT(openParens > 0);
					tempAncNode = tempAncNode->par;
					PHYCAS_ASSERT(tempAncNode != NULL);

					// Output names for internal nodes if they are available, otherwise 
					// do not output anything
					if (!tempAncNode->nodeName.empty())
						s << tempAncNode->nodeName;

					if (showEdgeLens)
						{
						s << ":";
						df.FormatDouble(s, tempAncNode->GetEdgeLen());
						}
					}
				while (tempAncNode->rSib != nd->nextPreorder);
				s << ',';
				}
			}
		else
			{
			// nd is internal
			s << '(';
			++openParens;
			}
		}
	
	if (openParens > 0)
		{
		// We are at the upper right tip node of the tree now, so work back down
		// to the root, outputting branch lengths along the way
		const TreeNode *tmpNode = GetLastPreorder();
		while (openParens > 0)
			{
			s << ')';
			--openParens;
			PHYCAS_ASSERT(tmpNode != NULL);
			tmpNode = tmpNode->par;
			TreeNode *tmpNodePar = tmpNode->par;
			bool tmpNodeParIsRoot = tmpNodePar->IsTipRoot();
			if (!tmpNode->nodeName.empty())
				s << tmpNode->nodeName;
			if (showEdgeLens && !tmpNodeParIsRoot)
				{
				s << ":";
				df.FormatDouble(s, tmpNode->GetEdgeLen());
				}
			}
		}
	return s;
	}

