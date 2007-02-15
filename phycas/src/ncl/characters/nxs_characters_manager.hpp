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

#ifndef NCL_CHARACTERS_MANAGER_H
#define NCL_CHARACTERS_MANAGER_H

#include <boost/weak_ptr.hpp>
#include "phycas/src/ncl/nxs_basic_manager.hpp"
#include "phycas/src/ncl/nxs_manager_with_listeners.hpp"
#include "phycas/src/ncl/taxa/nxs_taxa_listener.hpp"
#include "phycas/src/ncl/characters/nxs_char_listener.hpp"

class NxsCharactersBlock;
class NxsCharactersManager;
class NxsTaxaManager;
typedef boost::weak_ptr<NxsCharactersBlock> NxsCharsBlockAlias;
class NxsCharListener;
namespace ncl
	{
	class NxsDiscreteMatrix;
	class StoredMatrix;
	}

class NxsIndexSetCmdOption;
NxsIndexSetCmdOption * InstantiateCharSetOption(const std::string &n, NxsIndexSet *manip, const std::string &d, NxsCharactersManager *tmPtr, bool  persist, CmdPermissionLevel pLevel);

class NxsCharactersManager : 
  public NxsBasicListManager , 
  public ManagerWithListeners<NxsCharListener>, 
  public NxsTaxaListener,
  public NxsBlockManager
	{
	public:
		const NxsIndexSet & GetActiveCharacters() const	{return GetActiveSet();}
		NxsCharsBlockAlias 	GetCharsBlockReader();
		unsigned			GetNumChars() const;
		void				TaxaChanged(BaseTaxaManager *, NxsTaxaListener::TaxaChangeType);

	
		//  the NxsBasicListManager interface:
		//
		unsigned   			FindIndex(const std::string &label ) const;	// returns UINT_MAX if not found
		unsigned   			FindIndexFromNamesOnly(const std::string &label ) const;	// returns UINT_MAX if not found
		const std::string & GetUserSuppliedLabel( unsigned index ) const;	//ASSUMES that argument is a valid index
		unsigned			GetMaxLabelLength() const;
		unsigned			GetSize() const;
		void				Clear();
	
	
		unsigned GetNumMatrices() const
			{
			return (unsigned)matrices.size();
			}
			
		boost::shared_ptr<const ncl::StoredMatrix> 	GetMatrix(unsigned i) const
			{
			if (i >= GetNumMatrices())
				return boost::shared_ptr<const ncl::StoredMatrix>();
			return matrices[i];
			}
	
		CmdResult 			NewBlockRead(NxsBlock *newTax);
		NxsCharactersManager(NxsTaxaManager &);
		virtual ~NxsCharactersManager();
	protected:
		CmdResult			AddCharacters(boost::shared_ptr<ncl::NxsDiscreteMatrix> , unsigned, const NxsIndexSet *, const LabelMap *, const StateLabelMap *);
		void				AlertListeners(NxsCharListener::CharChangeType);
	
		typedef boost::shared_ptr<NxsCharactersBlock> 	SharedCharsBlockPtr;
		typedef std::vector<SharedCharsBlockPtr> 		VecCharsBlockPtrs;
		typedef boost::shared_ptr<ncl::StoredMatrix> 		StorageMatrixShPtr;
		typedef std::vector<StorageMatrixShPtr> 			VecMatrixShPtr;
		
		VecCharsBlockPtrs	 	charBlocks;
		LabelMap				charLabelMap;
		unsigned				nChars;
		VecMatrixShPtr 			matrices;
		NxsTaxaManager		  & nxsTaxaMgr;
		bool					reportNewCharacters;
	};
		
inline void NxsCharactersManager::Clear()
	{
	charLabelMap.clear();
	nChars = 0;
	IndicesDecreased(0);
	AlertListeners(NxsCharListener::kCharsCleared);
	}

inline unsigned NxsCharactersManager::FindIndex(const std::string &label ) const
	{
	return GetIndexFromNexusName(charLabelMap, label, (GetSize() > 0), true);
	}	
	
inline unsigned NxsCharactersManager::FindIndexFromNamesOnly(const std::string &label ) const
	{
	return GetIndexFromNexusName(charLabelMap, label, false);
	}	
	
inline unsigned NxsCharactersManager::GetMaxLabelLength() const
	{
	return GetMaxNexusLabelLength(charLabelMap, GetNumChars());
	}
	
inline unsigned NxsCharactersManager::GetSize() const	
	{
	return GetNumChars();
	}

inline void NxsCharactersManager::TaxaChanged(
  BaseTaxaManager *,
  NxsTaxaListener::TaxaChangeType)
  	{
  	Clear();	//@should only clear appropriate rows of data, and should fill in missing for new taxa
	}


/*----------------------------------------------------------------------------------------------------------------------
|	returns the number of characters stored
*/
inline unsigned NxsCharactersManager::GetNumChars() const 
	{
	return nChars;
	}

inline const std::string & NxsCharactersManager::GetUserSuppliedLabel( unsigned  i) const
	{
	return GetLabelFromMap(charLabelMap, i);
	}

#endif

