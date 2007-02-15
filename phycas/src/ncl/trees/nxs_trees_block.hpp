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

#ifndef NCL_NXSTREESBLOCK_H
#define NCL_NXSTREESBLOCK_H

#include "phycas/src/ncl/nxs_defs.hpp"
#include "phycas/src/ncl/nxs_basic_manager.hpp"
#include "phycas/src/ncl/nxs_block.hpp"
#include "phycas/src/ncl/taxa/nxs_taxa_manager.hpp"
#include "phycas/src/ncl/trees/newick_verifier.hpp"

class NxsTreesManager;
/*----------------------------------------------------------------------------------------------------------------------
|	Class that reads trees blocks and stores descriptions of the trees
*/
class NxsTreesBlock : public NxsCommandManagerBlock , public NxsIndexInfo, public NewickVerifier
	{
	public :
		NxsTreesBlock(NxsTaxaManager &	taxaMgr, NxsTreesManager &	treesMgr);
	
		//	Accessors
		//
		unsigned			FindIndex(const std::string &label ) const;	// returns UINT_MAX if not found
		unsigned   			FindIndexFromNamesOnly(const std::string &label ) const;
		const std::string & GetUserSuppliedLabel(unsigned i) const;
		unsigned 			GetMaxLabelLength() const;
		unsigned			GetNumTrees() const;
		unsigned   			GetSize() const	{return GetNumTrees();}
		void				AddTree(const FullTreeDescription &d);
		void				Clear();
		CmdResult 			EndEncountered();
		void				ResetCmdMgrNxsBlock();

		const VecOfTreeDescriptions & GetTreeDescriptions() const;
		unsigned 			GetDefaultTree() const;
		void				ReadNewickRepresentationAndAddToTaxa(NxsToken & token);
	private:
		
		void				AddRecognizedCommands();
		bool				ParseTree(NxsToken &);
		bool				ParseTranslate(NxsToken &);
		
		VecOfTreeDescriptions	trees;				/*	full descriptions of all known trees */
		unsigned 				defTree;			/*	index of the default tree */
		friend class SimpleNode;
		friend class NxsBasicListManager;
		
		NxsTaxaManager	  &	taxaMgr;
		NxsTreesManager	  &	treesMgr;
	};

#if defined (NCL_SUPPORT_OLD_NAMES)
	typedef NxsTreesBlock TreesBlock;
#endif

inline const VecOfTreeDescriptions &NxsTreesBlock::GetTreeDescriptions() const
	{
	return trees;
	}
	
inline unsigned NxsTreesBlock::GetDefaultTree() const
	{
	return defTree;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns the number of trees stored
*/
inline unsigned NxsTreesBlock::GetNumTrees() const 
	{
	return (unsigned) trees.size();
	}
	
inline void NxsTreesBlock::ResetCmdMgrNxsBlock()
	{
	trees.clear();
	defTree = UINT_MAX;
	translationTable.clear();
	}

inline void NxsTreesBlock::Clear()
	{
	Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns the index of the tree with the name "label", or -1 if that name is unknown.  (not case sensitive)
*/
inline unsigned NxsTreesBlock::FindIndex(
  const std::string &label) const
	{
	return GetIndexFromNexusNameGeneric(trees, label);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	returns the index of the tree with the name "label", or -1 if that name is unknown.  (not case sensitive)
*/
inline unsigned NxsTreesBlock::FindIndexFromNamesOnly(
  const std::string &label) const
	{
	return GetIndexFromNexusNameOnlyGeneric(trees, label);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the label (name) of tree i
|	Throws IndexInfo::XNotFound() if the index ind is not valid.
*/
inline const std::string & NxsTreesBlock::GetUserSuppliedLabel(unsigned i) const
	{
	return GetLabelFromVecGeneric(trees, i);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the length of the longest tree name.
*/
inline unsigned NxsTreesBlock::GetMaxLabelLength() const
	{
	return (unsigned) GetMaxNexusLabelGeneric(trees);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Adds a new tree description to the vector of trees.
*/
inline void NxsTreesBlock::AddTree(const FullTreeDescription &d)
	{
	trees.push_back(d);
	}

#endif
