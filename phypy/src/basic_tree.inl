#ifndef PYPHY_BASIC_TREE_INL
#define PYPHY_BASIC_TREE_INL

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Sets firstPreorder to NULL and calls Clear() to initialize data members.
*/
inline Tree::Tree()
  : firstPreorder(NULL)
	{
	Clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls Clear() to store all nodes and delete all other allocated memory, then erases `nodeStorage' to eliminate the
|	memory allocated for TreeNode objects.
*/
inline Tree::~Tree()
	{
	Clear(); //@ not too efficient to store nodes and then immediately delete them all
	nodeStorage.erase(nodeStorage.begin(), nodeStorage.end());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Attempts to convert the name of `nd' to an integer (call it x) and then use x as the node number for `nd'. The 
|	supplied set `used' is consulted to ensure that no two nodes are given the same number. If conversion of the name to
|	the integer value x is successful, and if x is not in the set `used', then: (1) the node number of `nd' is set to x;
|	(2) x is added to the set `used'; and (3) the function returns true. If the conversion of name to number fails 
|	because name is not able to be converted to an integer, then one of two things can happen: (1) if `used' is empty, 
|	then the function returns false, which serves as a signal that this function should not be called for future nodes
|	because apparently node names do not represent integers in this tree description; (2) if `used' is not empty, then
|	an XPhylogeny exception is thrown because some node numbers have already been set according to node names and thus
|	were we to discontinue converting names to numbers, the tree would end up with a mixture of node numbers derived 
|	from node names and node numbers assigned arbitrarily. If the conversion of name to number fails because the number
|	has already been used, then an XPhylogeny exception is thrown because this is clearly an invalid tree description.
|	This function sets the `numbers_from_names' data member to true if it succeeds, and leaves it false if an exception
|	is thrown or the function returns false. Note: this function has no way of knowing whether node names (=numbers) in
|	the tree description start at 0 or 1, so `used' should be consulted after the tree is built and if 0 is not in the
|	set, a sweep of the tree should be made to decrement each tip node number.
*/
inline bool Tree::SetNumberFromName(TreeNode * nd, std::set<unsigned> & used)
	{
	PHYCAS_ASSERT(nd != NULL);
	std::string name = nd->GetNodeName();
	boost::trim(name);
	unsigned x = 0;
	try
		{
		x = boost::lexical_cast<unsigned>(name);
		}
	catch(boost::bad_lexical_cast &)
		{
		// signal that the node name could not be converted to an integer value
		x = UINT_MAX;
		}

	// This data member is set to false initially for every node examined.
	// It will thus end up being true after a tree is built only if the 
	// conversion from names to numbers succeeds for every node.
	numbers_from_names = false;

	if (x == UINT_MAX)
		{
		if (used.empty())
			{
			// First node examined, so return false to signal that this function should 
			// not be called for future node number assignments
			return false;
			}
		else
			{
			// This is not the first node examined, and previously examined nodes did have
			// names that could be converted to numbers.
			std::string s = "node name (";
			s << name;
			s << ") could not be converted to an integer to be used as the node number";
			throw XPhylogeny(s);
			}
		}

	// It is an error if x has already been used for a different node
	std::pair<std::set<unsigned>::iterator, bool> insert_result = used.insert(x);
	if (insert_result.second)
		{
		// insertion was made, so x has NOT already been used
		nd->SetNodeNum(x);
		}
	else
		{
		// insertion was not made, so set already contained x
		std::string s = "node number (";
		s << x;
		s << ") used for more than one node";
		throw XPhylogeny(s);
		}

	numbers_from_names = true;
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Used to find out where we are in Python doctest examples.
*/
inline void Tree::DebugHere(
  std::string s)	/**< is a string indicating the current location */
	{
	std::cerr << "~~~~~~~ " << s << " ~~~~~~~" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assigns all edge lengths in the tree to the supplied value `v'.
*/
inline void Tree::SetAllEdgeLens(double v)
	{
	std::for_each(begin(), end(), boost::lambda::bind(&TreeNode::SetEdgeLen, boost::lambda::_1, v));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns Tree object to just-constructed state. If preorder pointers are valid, walks tree in postorder fashion 
|	storing nodes as they are visited in `nodeStorage'. If the preorder pointers are not valid, first calls 
|	RefreshPreorder() to remedy this. If `firstPreorder' is NULL, assumes that tree has just been constructed and there
|	are thus no nodes to store.
*/
inline void Tree::Clear()
	{
	if (firstPreorder != NULL)
		{
		// Walk tree in postorder fashion, storing nodes as we go
		for (TreeNode *nd = GetLastPreorder(); nd != NULL;)
			{
			TreeNode *next = nd->GetNextPostorder();
			nd->Clear();
			nodeStorage.push_back(nd);
			nd = next;
			}

		firstPreorder = NULL;
		}

	lastPreorder		= NULL;
	nInternals			= 0;
	nTips				= 0;
	preorderDirty		= true;
	hasEdgeLens			= false;
	isRooted			= false;
	nodeCountsValid		= false;
	treeid_valid		= false;
	numbers_from_names	= false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Pulls a node out of storage, if `nodeStorage' is not empty; otherwise, allocates memory for a new TreeNode. If it
|	is known in advance how many nodes will be needed, Reserve() can be called to create all nodes needed, storing them
|	in `nodeStorage'.
*/
inline TreeNode *Tree::GetNewNode()
	{
	TreeNode *nd = NULL;
	if (nodeStorage.empty())
		{
		nd = new TreeNode();
		}
	else
		{
		nd = nodeStorage.back();
		nodeStorage.pop_back();
		}
	return nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates memory for `n' TreeNode objects and stores pointers to these in the `nodeStorage' vector. If some nodes
|	are already in `nodeStorage', only allocates enough new nodes to make the length of `nodeStorage' equal to `n'.
|	If the length of `nodeStorage' already equals or exceeds `n', returns immediately without allocating any new nodes.
*/
inline void Tree::Reserve(
  unsigned n)	/**< is the desired length of the `nodeStorage' vector */
	{
	PHYCAS_ASSERT(n > 0);
	unsigned num_existing_nodes = (unsigned)nodeStorage.size();
	if (num_existing_nodes >= n)
		return;
	
	unsigned num_nodes_needed = n - num_existing_nodes;
	for (unsigned i = 0; i < num_nodes_needed; ++i)
		{
		nodeStorage.push_back(new TreeNode());
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `firstPreorder' points to a node having no parent, no right sibling and only one child.
*/
inline bool Tree::RootValid() const
	{
	return (firstPreorder != NULL && firstPreorder->IsRoot());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `isRooted' data member.
*/
inline bool Tree::IsRooted() const
	{
	return isRooted;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `hasEdgelens' data member.
*/
inline bool Tree::HasEdgeLens() const
	{
	return hasEdgeLens;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `preorderDirty' data member.
*/
inline bool Tree::PreorderDirty() const
	{
	return preorderDirty;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of the `numbers_from_names' data member. Set member function Tree::SetNumberFromName for
|	details.
*/
inline bool Tree::TipNumbersSetUsingNames() const
	{
	return numbers_from_names;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the total number of nodes currently composing the tree (i.e. `nTips' plus `nInternals'). For unrooted 
|	trees of n taxa, one of the observable nodes serves as the root node, so this number can range from n+1 (star tree) 
|	to 2n-2 (fully-resolved). For rooted trees, the range is n+2 (star tree) to 2n-1 (fully-resolved) because the root 
|	node exists but is not identical to one of the observahble nodes. If `firstPreorder' is NULL, returns 0 because this 
|	means that the tree does not yet have any nodes. For trees that have at least one node, calls RefreshNodeCounts() if 
|	`nodeCountsValid' is false. 
*/
inline unsigned Tree::GetNNodes()
	{
	// First check to make sure that a tree really exists
	if (firstPreorder == NULL)
		return 0;

	if (!nodeCountsValid)
		RefreshNodeCounts();
	return (nTips + nInternals);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of tip nodes currently composing the tree. The data member `nTips' stores the number of degree 
|	one nodes in the tree (tip means degree one). Note that one of these degree one nodes is the root node. The number 
|	of tips equals the number of taxa for unrooted trees. For rooted trees of n taxa, there are n+1 tips because the 
|	root node is a tip node but not an observable node. If `firstPreorder' is NULL, returns 0 because this means that 
|	the tree does not yet have any nodes. For trees that have at least one node, calls RefreshNodeCounts() if 
|	`nodeCountsValid' is false. 
*/
inline unsigned Tree::GetNTips()
	{
	// First check to make sure that a tree really exists
	if (firstPreorder == NULL)
		return 0;

	if (!nodeCountsValid)
		RefreshNodeCounts();
	return nTips;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of observable tip nodes currently composing the tree. The data member `nTips' stores the number 
|	of degree one nodes in the tree (tip means degree one). Note that one of these degree one nodes is the root node. 
|	The number of observables equals the number of tips for unrooted trees. For rooted trees of n taxa, however, the 
|	number of observables equals one fewer than the number of tips because the root node is a tip but not an observable.
|	If `firstPreorder' is NULL, returns 0 because this means that the tree does not yet have any nodes. For trees that 
|	have at least one node, calls RefreshNodeCounts() if `nodeCountsValid' is false. 
*/
inline unsigned Tree::GetNObservables()
	{
	return GetNTips() - (IsRooted() ? 1 : 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of internal nodes currently composing the tree (i.e. `nInternals'). The data member `nInternals' 
|	stores the number of nodes in the tree of degree 2 or higher. Its value is 1 for a (rooted or unrooted) star tree, 
|	n-2 for a fully-resolved unrooted tree, and n-1 for a fully-resolved rooted tree. If `firstPreorder' is NULL, 
|	returns 0 because this means that the tree does not yet have any nodes. For trees that have at least one node, calls
|	RefreshNodeCounts() if `nodeCountsValid' is false. 
*/
inline unsigned Tree::GetNInternals()
	{
	// First check to make sure that a tree really exists
	if (firstPreorder == NULL)
		return 0;

	if (!nodeCountsValid)
		RefreshNodeCounts();
	return nInternals;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Return value of data member `firstPreorder'. Note: makes a (rather expensive) call to RefreshPreorder() if the
|	preorder pointers need to be updated (i.e., `preorderDirty' is true).
*/
inline TreeNode	* Tree::GetFirstPreorder()
	{
	if (preorderDirty)
		RefreshPreorder();

	return firstPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Return value of data member `lastPreorder'. Note: makes a (rather expensive) call to RefreshPreorder() if the
|	preorder pointers need to be updated (i.e., `preorderDirty' is true).
*/
inline TreeNode * Tree::GetLastPreorder()
	{
	if (preorderDirty)
		RefreshPreorder();

	return lastPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates and returns a vector of pointers to all TreeNode objects that are responsible for edge lengths. The nodes 
|	that qualify are all nodes in the tree except the root node, which is a tip node whose edge is managed by the 
|	internal node that represents its only child.
*/
inline Tree::TreeNodeVec Tree::GetNodesWithEdges()
	{
	TreeNodeVec v;

	TreeNode * nd = GetFirstPreorder();	// this is the root node, which we skip
	for (nd->GetNextPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		v.push_back(nd);
		}

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Return preorder_iterator representing beginning of the preorder sequence of nodes. Note: this function calls
|	GetFirstPreorder(), which will make the (expensive) call to RefreshPreorder() if preorder pointers need to be 
|	updated.
*/
inline preorder_iterator Tree::begin()
	{
	return preorder_iterator(GetFirstPreorder());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Return preorder_iterator representing end of the preorder sequence of nodes.
*/
inline preorder_iterator Tree::end()
	{
	return preorder_iterator();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a parenthetical description of the tree (Newick format). If the tree has edge lengths, these are included in
|	the Newick description. Any nodes having non-empty node names are named in the Newick description. Tip nodes that
|	do not have stored names are given a name that equals 1 plus their node number. Calls AppendNewick() to do the work,
|	but unlike AppendNewick() does not have a std::string reference parameter and returns a std::string by value rather
|	than by reference.
*/
inline std::string Tree::MakeNewick()
	{
	std::string s;
	AppendNewick(s);
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `treeid_valid' to false, indicating that this tree's TreeID is no longer valid because of a recent
|	change in the tree topology.
*/
inline void Tree::InvalidateTreeID()
	{
	treeid_valid = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Flags the `tree_id' data member as being in need of recalculation. Lazy recalculation is used because sometimes
|	several operations in quick succession all make changes to the topology of the tree without needing to reference
|	the `tree_id' data member. Recalculating `tree_id' for each of these would slow down the program unecessarily.
*/
inline void Tree::InvalidateID()
	{
	treeid_valid = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `nodeCountsValid' to false, indicating that the node counts reported by member functions such as
|	GetNNodes(), GetNTips(), GetNInternals(), etc., may no longer be valid due to a recent change in the tree topology.
*/
inline void Tree::InvalidateNodeCounts()
	{
	nodeCountsValid = false;
	}

}	// namespace phycas

#endif