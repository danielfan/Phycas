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

#ifndef PYPHY_BASIC_TREE_HPP
#define PYPHY_BASIC_TREE_HPP

#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <boost/algorithm/string.hpp>	// used by SetNumberFromName member function
#include <boost/lexical_cast.hpp>		// used by SetNumberFromName member function
#include "phycas/src/basic_tree_node.hpp"
#include "phycas/src/tree_iterators.hpp"
#include "phycas/src/xphylogeny.hpp"
#include "phycas/src/phycas_string.hpp"

//	To do next:
//		o ASCII tree drawing
//		o splits
//		o internal nodes should be numbered starting with nTips
//		o should ensure that tip node numbering sequence starts and 0 and ends at ntax-1 with no gaps
//			(see bottom of BuildFromString). (Plant to deal with valid exceptions to this rule later.)

namespace phycas{

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the notion of a phylogenetic tree. This class has only methods necessary for representing the tree, 
|	copying trees from other trees, manipulating the tree, and representing the tree graphically. It does not perform 
|	any specialized activities, such as computing its own likelihood: these sorts of activities are left to member
|	functions of other classes.
*/
class Tree
	{
	public:

		friend class TreeManip;
		friend class BushMove;

		typedef std::vector<TreeNode *> TreeNodeVec;

		// Constructors/destructor
		//								
								Tree();
		virtual					~Tree();
		// Accessors
		//
		unsigned				GetNNodes();
		unsigned				GetNTips();
		unsigned				GetNInternals();
		unsigned				GetNObservables();
		TreeNode *				GetFirstPreorder();
		TreeNode *				GetLastPreorder();
		preorder_iterator		begin();
		preorder_iterator		end();
		// Predicates
		//
		bool					RootValid() const;
		bool					IsRooted() const;
		bool					HasEdgeLens() const;
		bool					PreorderDirty() const;
		bool					TipNumbersSetUsingNames() const;
		// Utilities
		//
		void					Clear();
		void					BuildFromString(const std::string &newick); // throws XPhylogeny
		void					RectifyNumbers(std::vector<std::string> name_vector); // throws XPhylogeny
		void					RectifyNames(std::vector<std::string> name_vector); // throws XPhylogeny
		double					EdgeLenSum();
		double					calcTotalHeight();
		void					SetAllEdgeLens(double v);
		void					RerootAtThisTip(TreeNode * nd);
		void					RerootAtThisInternal(TreeNode * nd);
		void					RerootAtTip(unsigned num);
		void					RefreshPreorder(TreeNode *nd = NULL);
		std::string &			AppendNewick(std::string &);
		std::string				MakeNewick();
		bool					SetNumberFromName(TreeNode * nd, std::set<unsigned> & used);
		TreeNodeVec				GetNodesWithEdges();
		void					SelectAllNodes();
		void					UnselectAllNodes();
		// Debugging
		//
		std::string				DebugWalkTree(bool preorder = true, unsigned verbosity = 0);
		void					DebugHere(std::string s);

	protected:

		TreeNode *				FindTipNode(unsigned num);
		void					RerootHelper(TreeNode *m, TreeNode *t);
		void					GetNextNewickToken(const std::string &newick, unsigned start_pos);
		TreeNode *				GetNewNode();
		void					StoreTreeNode(TreeNode * u);
		void					Reserve(unsigned n);
		void					InvalidateID();
		void					InvalidateNodeCounts();	
		void					RefreshNodeCounts();
		void					InvalidateTreeID();
		//void					RefreshTreeID();

	protected:

		//TreeID					tree_id;			/**< A vector of splits that uniquely identify the tree topology */
		bool					treeid_valid;			/**< True if the tree_id data member is valid; if false, call RefreshTreeID to make it valid again */
		TreeNodeVec				nodeStorage;			/**< A vector of pointers to TreeNode objects */
		TreeNode *				firstPreorder;			/**< Pointer to the first preorder node (equals last postorder node) */
		TreeNode *				lastPreorder;			/**< Pointer to the last preorder node (equals first postorder node) */	
		unsigned				nTips;					/**< Total number of tip (degree = 1) nodes in the tree */
		unsigned				nInternals;				/**< Total number of internal (degree > 1) nodes in the tree */
		bool					hasEdgeLens;			/**< True if edge lengths have been specified */
		bool					isRooted;				/**< True if the tree is rooted */
		bool					preorderDirty;			/**< Set to false when preorder traversal pointers are set, but a function should set to true if it modifies the tree and does not leave the preorder traversal pointers valid */
		bool					nodeCountsValid;		/**< If false, causes functions that depend on accurate node counts, such as GetNTips(), GetNNodes(), GetNInternals() and GetNObservables(), to recompute nTips and nInternals */
		bool					numbers_from_names;		/**< True if tip node numbers were set using tip node names in tree description. */

	private:

		std::string				workspace;				/**< Used by GetNextNewickToken for storing tokens read from newick tree descriptions */
	};

typedef boost::shared_ptr<Tree> TreeShPtr;

}	// namespace phycas

#include "phycas/src/basic_tree.inl"

#endif
