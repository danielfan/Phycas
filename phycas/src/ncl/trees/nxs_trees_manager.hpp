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

#ifndef NCL_NXS_TREE_MANAGER_H
#define NCL_NXS_TREE_MANAGER_H

#include <boost/weak_ptr.hpp>
#include "phycas/src/ncl/nxs_basic_manager.hpp"
#include "phycas/src/ncl/trees/full_tree_description.hpp"
#include "phycas/src/ncl/nxs_manager_with_listeners.hpp"
#include "phycas/src/ncl/taxa/nxs_taxa_listener.hpp"
#include "phycas/src/ncl/trees/nxs_tree_listener.hpp"
class NxsTreesBlock;
class NxsTaxaManager;
typedef boost::shared_ptr<NxsTreesBlock> TreesBlockShPtr;
class NxsTreeListener;
class NxsIndexSetCmdOption;
NxsIndexSetCmdOption *InstantiateTreeSetOption(const std::string &n, NxsIndexSet *manip, const std::string &d, NxsTreesManager *tm,  bool mustBeSetExplicitly, bool  persist, CmdPermissionLevel pLevel);
class NxsTreesManager : 
  public NxsBasicListManager , 
  public ManagerWithListeners<NxsTreeListener>, 
  public NxsTaxaListener
	{
	public:
		virtual CmdResult		AddNewTrees(const std::vector<FullTreeDescription> &v, unsigned d = UINT_MAX);
		TreesBlockShPtr 		GetTreesBlockReader();
		unsigned				GetNumTrees() const;
		
		FullTreeDescription   	GetTree(unsigned i) const;
		const VecOfTreeDescriptions & GetTrees() const
			{
			return treeDescriptions;
			}
		VecOfTreeDescriptions 	GetTrees(const NxsIndexSet &) const;
		
		void					TaxaChanged(BaseTaxaManager *, NxsTaxaListener::TaxaChangeType);
	
		//  the NxsBasicListManager interface:
		//
		unsigned   				FindIndex(const std::string &label ) const;	// returns UINT_MAX if not found
		unsigned   				FindIndexFromNamesOnly(const std::string &label ) const;
		const std::string     & GetUserSuppliedLabel( unsigned index ) const;	//ASSUMES that argument is a valid index
		unsigned				GetMaxLabelLength() const;
		unsigned				GetSize() const;
		void					Clear();
	
	
		CmdResult 				NewBlockRead(NxsBlock *newTax);
	//protected:
		NxsTreesManager(NxsTaxaManager & taxaMgr);
		virtual ~NxsTreesManager();
	protected:
		void			AlertListeners(NxsTreeListener::TreeChangeType);
	
		typedef boost::shared_ptr<NxsTreesBlock> SharedTreesBlockPtr;
		typedef std::vector<SharedTreesBlockPtr> 	VecOfTreesBlockPtrs;
		
		VecOfTreesBlockPtrs	 	treeBlocks;
		VecOfTreeDescriptions	treeDescriptions;
		unsigned				defaultTree; /* most recent tree that was given an asterisk, UINT_MAX if no such trees have been read */
	public:
		NxsTaxaManager		  & taxaMgr;
	};

inline void NxsTreesManager::Clear()
	{
	treeDescriptions.clear();
	IndicesDecreased(0);
	defaultTree = UINT_MAX;
	AlertListeners(NxsTreeListener::kTreesCleared);
	}

inline unsigned NxsTreesManager::FindIndex(const std::string &label ) const
	{
	return GetIndexFromNexusNameGeneric(treeDescriptions, label);
	}	
	
inline unsigned NxsTreesManager::FindIndexFromNamesOnly(
  const std::string &label) const
	{
	return GetIndexFromNexusNameOnlyGeneric(treeDescriptions, label);
	}

inline const std::string & NxsTreesManager::GetUserSuppliedLabel( unsigned  i) const
	{
	return GetLabelFromVecGeneric(treeDescriptions, i);
	}
	
inline unsigned NxsTreesManager::GetMaxLabelLength() const
	{
	return GetMaxNexusLabelGeneric(treeDescriptions);
	}
	
inline unsigned NxsTreesManager::GetSize() const	
	{
	return GetNumTrees();
	}

inline void NxsTreesManager::TaxaChanged(
  BaseTaxaManager *,
  NxsTaxaListener::TaxaChangeType)
  	{
  	Clear();	//@should prune taxa, but we'll just clear for the time being
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns the number of trees stored
*/
inline unsigned NxsTreesManager::GetNumTrees() const 
	{
	return (unsigned) treeDescriptions.size();
	}
	
#endif

