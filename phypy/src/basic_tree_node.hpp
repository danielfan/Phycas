/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

#ifndef PYPHY_BASIC_TREE_NODE_HPP
#define PYPHY_BASIC_TREE_NODE_HPP

#include <string>
#include <boost/function.hpp>

namespace phycas
{

// Forward declarations: need to define these classes only if using Tree and TreeNode
// for likelihood calculations
class TipData;
class InternalData;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates a node (i.e. vertex) of a tree, with pointers to other nodes for navigation purposes, a branch length
|	data member for purposes of calculating likelihoods and drawing the tree, and pointers to unspecified TipData and
|	InternalData structures that can be used to store, for example, transition probability matrices and conditional
|	likelihood arrays.
*/
class TreeNode
	{
	public:
		friend class TreeManip;

						TreeNode();
		virtual			~TreeNode();
		
		// Typedefs
		//
		typedef boost::function< void (TipData *) >			TipDataDeleter;
		typedef boost::function< void (InternalData *) >	InternalDataDeleter;

		// Predicates
		//
		bool			HasChildren() const;
		bool			NoChildren() const;
		bool			HasParent() const;
		bool			NoParent() const;
		bool			IsTip() const;
		bool			IsObservable() const;
		bool			IsRoot() const;
		bool			IsInternal() const;
		bool			NumberNotYetAssigned() const;
		bool			EdgeLenNotYetAssigned() const;
		bool			IsSelected() const;

		// Accessors
		//
		float					GetX();
		float					GetY();
		const std::string &		GetNodeName() const;
		unsigned				GetNodeNumber() const;
		TreeNode *				GetLeftChild();
		const TreeNode *		GetLeftChildConst() const;
		TreeNode *				GetRightSib();
		const TreeNode *		GetRightSibConst() const;
		TreeNode *				GetParent();
		const TreeNode *		GetParentConst() const;
		TreeNode *				GetNextPreorder();
		TreeNode *				GetNextPostorder();
		double					GetEdgeLen() const;
		TipData *				GetTipData();
		InternalData *			GetInternalData();
		const TipData *			GetTipData() const;
		const InternalData *	GetInternalData() const;

		// Modifiers
		//
		void			SetX(float xx);
		void			SetY(float yy);
		void			SelectNode();
		void			UnselectNode();
		void			SetEdgeLen(double x);
		void			SetNodeName(std::string name);
		void			SetNodeNum(unsigned num);

		void			ResetTipData();
		void			SetTipData(TipData * d, TreeNode::TipDataDeleter f);

		void			ResetInternalData();
		void			SetInternalData(InternalData * c, TreeNode::InternalDataDeleter f);

		void			SetObservable(bool is_observable);

		void			Clear();

		// Utilities
		//
		unsigned		CountChildren() const;
		void			AppendNodeInfo(std::string &s, bool num_and_name_only = false);

	protected:

		std::string			nodeName;				/**< name of node */
		unsigned			nodeNum;				/**< for tips, this is the taxon index, ranging from 0 to ntips-1 */
		double				edgeLen;				/**< length of this node's edge */
		TreeNode *			lChild;					/**< points to leftmost child */
		TreeNode *			par;					/**< points to parent */
		TreeNode *			rSib;					/**< points to sibling on right */
		TreeNode *			nextPreorder;			/**< points to next node in preorder sequence */
		TreeNode *			prevPreorder;			/**< points to previous node in preorder sequence */
		bool				observable;				/**< true if data could be observed for this node */

		float				x;						/**< x-coordinate for purposes of drawing the tree */
		float				y;						/**< y-coordinate for purposes of drawing the tree */
		bool				selected;				/**< can be used anytime a node needs to be selected for some purpose */

		// Pointers to structures used by likelihood calculation routines. 
		//
		TipData	 *			tipData;				/**< is a pointer to a structure used to store data for tip nodes */
		TipDataDeleter		tipDataDeleter;			/**< function object used to delete memory allocated for tipData */
		InternalData *		internalData;			/**< is a pointer to a structure used to store data for internal nodes */
		InternalDataDeleter	internalDataDeleter;	/**< function object used to delete memory allocated for `internalData' */

	public:

		static const double		edgeLenEpsilon;		/**< smallest allowable edge length */
		static const double		edgeLenDefault;		/**< default edge length */
		static const unsigned	nodeNumInitValue;	/**< default number for newly-created nodes */
		static const double		edgeLenInitValue;	/**< default edge length for newly-created nodes */

		friend class Tree;	
	};

}	// namespace phycas

#include "phypy/src/basic_tree_node.inl"	// TreeNode inlined function bodies

#endif

