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

#ifndef TREE_NODE_ATTRIBUTES_H
#define TREE_NODE_ATTRIBUTES_H

class TreeNode;
class PhoMultiSiteModel;

/*----------------------------------------------------------------------------------------------------------------------
|	The info that each node of the tree should have so that the likelihood can be calculated. Nodes have an array of
|	pointers to TreeNodeAttributes, one element for each subset of the character partition. If there is no character
|	partition, the array has only the 0th. element. A particular model and data partition scheme determine the exact 
|	structure of the TreeNodeAttribute objects, which are constructed prior to calculation of likelihoods. Tree objects
|	can thus exist with NULL TreeNodeAttribute pointers; TreeNodeAttribute object represent unnecessary baggage for the
|	majority of tree operations. Because data is stored differently for terminal vs. internal nodes, some of the 
|	pointers will be NULL in any particular TreeNodeAttribute.
*/
class TreeNodeAttribute : public boost::noncopyable 
	{
	friend class TreeNode;
	friend class HKYAdHocEvaluator; //@POL can be removed - used to get star tree paper out the door

	public:
						TreeNodeAttribute(unsigned ns, unsigned nr, double *blen);
		virtual			~TreeNodeAttribute();
		
		void 			RefreshPrMatrices(const PhoMultiSiteModel *model);

	protected:
		double		*	condLike;		/**< Array of conditional likelihoods */
		TreeNode	*	lastAvoidNode;	/**< The last time the likelihood was calculated, the condlike was the conditional likelihood away from this node. Used for internal nodes only. */
		double		*	edgeLen;		/**< Pointer to the current edge length for this partition/node */
		double		***	prMatrices;		/**< Array of change probability matrices, each of which is a square matrix with dimension equal to the number of states */

		//@ perhaps replace these two bools with a char containing bits to set
		bool			pMatIsDirty;	/**< True if change probability matrices need to be recalculated. */ 
		bool			clDirty;		/**< True if conditional likelihood arrays need to be recalculated */
	};


#endif


