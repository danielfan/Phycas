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

#if !defined(PHO_TREE_H)
#define PHO_TREE_H

#include <stack>
#include "phycas/src/oldphycas/tree_id.hpp"
#include "phycas/src/oldphycas/node_iterator_types.hpp"
#define ROOTING_AT_A_TIP 1

typedef TreeNode*			TNodePtr;
typedef std::stack<TreeNode *>	TreeNodeStack;

#if defined(NCL_USE_NXS_CONSOLE_OUTPUT)
#	class NxsStdOutputStream;
#endif
class FullTreeDescription;
class NxsToken;
class DrawContext;

/*----------------------------------------------------------------------------------------------------------------------
|	Exception thrown if invalid newick tree definition encountered
*/
class XBadTree
	{
	public:
		std::string msg;
		XBadTree(const char *c) : msg(c){}
	}; 
	
/*----------------------------------------------------------------------------------------------------------------------
|	Exception thrown if invalid newick tree definition encountered
*/
class XBadTreeDef 
	{
	public:
		unsigned long startpos;
		unsigned long endpos;
		std::string msg;
		XBadTreeDef(const char *m, const NxsToken &t);
		XBadTreeDef(const char *m);
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Exception thrown by DebugCheckTreeStructure if invalid navigational pointers are discovered for any node.
*/
class XBadTreeStructure
	{
	public:
		unsigned nodeNum;
		bool isLeaf;
		std::string msg;
		XBadTreeStructure(const char *m, unsigned which_nodenum, bool is_leaf);
		XBadTreeStructure(const char *m, const TreeNode*);
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the notion of a phylogenetic tree. This class has only methods necessary for representing the tree,
|	copying trees from other trees, and manipulating the tree. It does not perform any specialized activities, such 
|	as compute its own likelihood: these sorts of things are left to member functions of other classes that take trees
|	as arguments.
*/
class Tree 
	{
	public:

		// Constructors/destructor
		//
							Tree();
							Tree(const TreeID &);
							Tree(const FullTreeDescription &);
							Tree(const Tree &);
		virtual				~Tree();

		// Accessors
		//
		unsigned			GetNNodes() const;
		unsigned			GetNLeaves() const;
		unsigned			GetNInternals() const;
		unsigned			GetNStoredNodes() const;
		const TreeID		&GetID() const;
		const TreeNode		*GetRoot() const;
		const TreeNode		*GetFirstPreorder() const;
		const TreeNode		*GetLastPreorder() const;
		const TreeNode		*GetFirstPostorder() const;
		const TreeNode		*GetLastPostorder() const;

		TreeNode			*GetRoot();
		TreeNode			*GetFirstPreorder();
		TreeNode			*GetLastPreorder();
		TreeNode			*GetFirstPostorder();
		TreeNode			*GetLastPostorder();

		// Modifiers
		//
		void				InvalidateID();
		void				InvalidateNodeCounts();
		void				Clear();
		void				ChangeNInternals(signed delta);
		
		// Queries
		//
		bool				IsBinaryTree() const;
		bool				IsRooted() const;
		bool				HasEdgeLens() const;
		bool				IsValid() const;

		TreeNode		   *FindTip(unsigned tip_number);

		// Utilities
		//
		TreeNode		   *CreateTreeNode(unsigned ndN, double edgeLen, bool setSplitByNum);	
		TreeNode		   *CreateSimilarTreeNode(const TreeNode &);	//@POL Mark and Dave: ok to move this function from private to public?
		void				StoreTreeNode(TreeNode *u);
		void				DebugCheckTreeStructure(bool allow_polytomies = true, const TreeNode *nd = NULL) const;
		void				RefreshFullID(); //refreshed leaves and internals. not const because calls CreateSplit
		void				RefreshSelectionStatus(); // unselects all nodes
		void				RefreshID() const;	//refreshes internals. const because ID is mutable
		void				RefreshNodeCounts() const; // const because nLeaves and nInternals are mutable
		void				BuildTreeFromDescription(const FullTreeDescription &d);
		void				BuildTreeFromString(const char *s, bool translateNames = false, bool readAsRooted = false);
		bool				BuildTreeFromStringNoThrow(const char *s, bool translateNames = false, bool readAsRooted = false);
		void				List(int nindent = 2, TNodePtr p = NULL);
		int					Dex(TNodePtr p);
		
		bool				RerootAtTip(unsigned i, bool refresh_splits = true);

		void				Draw(DrawContext &) const; //throw XBadTree
		std::string		  & AppendNewickRepresentation(std::string & s, bool useTaxonNames = false, bool edgelens = true, bool fltlens = true, unsigned fltprec = 6) const;

		void				DebugShowPreorder(std::ostream &out);
		void				DebugShowPostorder(std::ostream &out);

		template <class FuncToCall> void PostorderTraverse(FuncToCall f);

		template <class FuncToCall> void ForAllLeaves(FuncToCall f);
		template <class FuncToCall> void ForAllLeaves(FuncToCall f) const;

#		if defined(NCL_USE_NXS_CONSOLE_OUTPUT)
			void 				DebugRerootingTest(TreeNode *nd, NxsStdOutputStream &o);
#		endif

		pre_iterator		pre_begin();
		pre_iterator		pre_end();
		
		const_pre_iterator	pre_begin() const;
		const_pre_iterator	pre_end() const;
		
		post_iterator		post_begin();
		post_iterator		post_end();
		
		const_post_iterator	post_begin() const;
		const_post_iterator	post_end() const;
		
		tip_iterator		tips_begin();
		tip_iterator		tips_end();
		
		const_tip_iterator	tips_begin() const;
		const_tip_iterator	tips_end() const;
		
		tip_iterator		leaves_begin();
		tip_iterator		leaves_end();
		
		const_tip_iterator	leaves_begin() const;
		const_tip_iterator	leaves_end() const;
	protected:
		mutable unsigned	nInternals;				/**< keeps total number of internal nodes in the tree */
		mutable unsigned	nLeaves;				/**< keeps total number of leaf nodes in the tree */
		TreeNode			*firstPreorder;			/**< pointer to the first pre-order node (equals last post-order node) */
		TreeNode			*lastPreorder;			/**< pointer to the last pre-order node (equals first post-order node) */
		bool				hasEdgelens;			/**< true if branch lengths have been specified */
		bool				isRooted;					
		mutable TreeID		id;						/**< vector of Split objects uniquely identifying the tree topology */
		mutable bool		id_valid;				/**< if false, causes GetID to recompute the id before returning it; set to true by RefreshID */
		mutable bool		node_counts_valid;		/**< if false, causes GetNLeaves, GetNNodes and GetNInternals to recompute nLeaves and nInternals */
		mutable double		max_x;					/**< computed by SetXYCoords; height (above root) of highest tip, but corresponds to horizontal length of tree since tree is shown lying on its side */
		mutable double		max_y;					/**< computed by SetXYCoords; twice number of tip nodes minus 1.0 (each tip is spaced 2.0 units apart, with first having y coordinate of 0.0 */
		bool				preordersDirty;			/**< set to false when preorder traversal pointers are set, but a function should set to true if it modifies the tree and does not leave the preorder traversal pointers valid */
		TreeNodeStack		nodeStorage;			/**< a stack for storing nodes not currently being used (e.g. an internal node detached by a death move during MCMC can be stored here) */
		
		void				TraverseTree(TreeNode *nd = NULL);

	private:
		Tree &operator=(const Tree &);
		void				CopyTreeStructure(const Tree &);
		void				DeleteSingleNode(TreeNode *);
		void				DrawTreeWithTempRoot(DrawContext &);
		TreeNode		  * BuildUnfinishedTreeToCheckString(const char *s, bool translateNames, std::set<unsigned> *);
		void 				Reroot(TreeNode *nd);
		STATIC_DATA_FUNC void FlipNodeHelper(TreeNode *nd);

		void				DeleteAllNodes();
		unsigned			GetLeafNum(const std::string &s, bool translateName);
		unsigned			GetInternalNodeNum(const std::string &s);
		//void				TranslateLine(std::string &lineBuffer, unsigned char *translation, unsigned windowWidth);
	
		friend class TreeManip;
		friend class SplitManager;
	};
typedef boost::shared_ptr<Tree> TreeShPtr;

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
	return hasEdgelens;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if and only if `firstPreorder' is not NULL.
*/
inline bool Tree::IsValid() const
	{
	return (firstPreorder != NULL);
	} 

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline TreeNode	*Tree::GetRoot()
	{
	return const_cast<TreeNode *>(const_cast<const Tree *>(this)->GetRoot());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline TreeNode	*Tree::GetFirstPreorder()
	{
	return const_cast<TreeNode *>(const_cast<const Tree *>(this)->GetFirstPreorder());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline TreeNode	*Tree::GetLastPreorder()
	{
	return const_cast<TreeNode *>(const_cast<const Tree *>(this)->GetLastPreorder());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline TreeNode	*Tree::GetFirstPostorder()
	{
	return const_cast<TreeNode *>(const_cast<const Tree *>(this)->GetFirstPostorder());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline TreeNode	*Tree::GetLastPostorder()
	{
	return const_cast<TreeNode *>(const_cast<const Tree *>(this)->GetLastPostorder());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline void Tree::ChangeNInternals(signed delta)
	{
	if (delta < 0)
		{
		unsigned u_delta = (unsigned)(-delta);
		PHYCAS_ASSERT(u_delta <= nInternals);
		nInternals -= u_delta;
		}
	else
		{
		nInternals += (unsigned)(delta);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline void Tree::InvalidateNodeCounts()
	{
	node_counts_valid = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline void Tree::InvalidateID()
	{
	id_valid = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline unsigned Tree::GetNNodes() const
	{
	if (!node_counts_valid)
		RefreshNodeCounts();
	return (nLeaves + nInternals);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline unsigned Tree::GetNLeaves() const
	{
	if (!node_counts_valid)
		RefreshNodeCounts();
	return nLeaves;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline unsigned Tree::GetNStoredNodes() const
	{
	return (unsigned)nodeStorage.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	return
*/
inline unsigned Tree::GetNInternals() const
	{
	if (!node_counts_valid)
		RefreshNodeCounts();
	return nInternals;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns a const ref to TreeID (refreshing if necessary)
*/
inline const TreeID &Tree::GetID() const
	{
	if (!id_valid)
		RefreshID();
	return id;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns a const pointer to the root of the tree (firstPreorder)
*/
inline const TreeNode *Tree::GetRoot() const
	{
	return firstPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns a const pointer to the GetFirstPreorder (firstPreorder) node in the tree
*/
inline const TreeNode *Tree::GetFirstPreorder() const
	{
	return firstPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns a const pointer to the GetLastPreorder (lastPreorder) node in the tree
*/
inline const TreeNode *Tree::GetLastPreorder() const
	{
	return lastPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns a const pointer to the GetFirstPostorder (lastPreorder) node in the tree
*/
inline const TreeNode *Tree::GetFirstPostorder() const
	{
	return lastPreorder;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns a const pointer to the GetLastPostorder (first preorder) node in the tree
*/
inline const TreeNode *Tree::GetLastPostorder() const
	{
	return firstPreorder;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Wrapper aroung BuildTreeFromString that catches exceptions.  Returns false if the tree wasn't built.
*/
inline bool Tree::BuildTreeFromStringNoThrow(const char *s, bool translateNames, bool readAsRooted)
	{
	try	
		{
		BuildTreeFromString(s, translateNames, readAsRooted);
		}
	catch (XBadTreeDef &)
		{
		return false;
		}
	return true;
	}

	
/*----------------------------------------------------------------------------------------------------------------------
|	Clears the tree and calls Tree::CopyTreeStructure
*/
inline Tree &Tree::operator=(const Tree &r)
	{
	Clear();
	CopyTreeStructure(r);
	return *this;
	}

//MTH 2-Dec-2004 added changed CreateNewickRepresentation to AppendNewickRepresentation (and changed return type to the string reference).
//POL 2-Dec-2004 added useTaxonNames to CreateNewickRepresentation
		
#endif
