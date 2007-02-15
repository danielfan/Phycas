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

#ifndef TAXA_MANAGER_H
#define TAXA_MANAGER_H

#include <boost/weak_ptr.hpp>
#include "phycas/src/ncl/nxs_defs.hpp"
#include "phycas/src/ncl/nxs_basic_manager.hpp"
#include "phycas/src/ncl/nxs_manager_with_listeners.hpp"
#include "phycas/src/ncl/taxa/nxs_taxa_listener.hpp"
class BaseTaxaManager;
class BaseTaxaManager : 
  public NxsBasicListManager ,
  public ManagerWithListeners< NxsTaxaListener >
	{
	public:
		const NxsIndexSet & GetActiveTaxa() const 
			{
			return NxsBasicListManager::GetActiveSet();
			}
		virtual CmdResult   AppendTaxa(const VecString & v);
		virtual CmdResult	ReplaceTaxa(const VecString & v);
		virtual CmdResult	AppendNewTaxon(const std::string &v);
		void	Clear()
			{
			taxLabels.clear();
			NxsBasicListManager::IndicesDecreased(0); 
			AlertListeners(NxsTaxaListener::kTaxaCleared);
			}
		unsigned GetNumTaxa() const
			{
			return taxLabels.size();
			}
		//  the virtual functions from the  NxsBasicListManager interface:
		//
		
		unsigned			FindNewOrExistingIndex(const std::string &label) const;
		unsigned   			FindIndex(const std::string &label ) const;	// returns UINT_MAX if not found
		unsigned   			FindIndexFromNamesOnly(const std::string &label ) const;	// returns UINT_MAX if not found
		const std::string   GetLabel( unsigned index ) const;	//ASSUMES that argument is a valid index
		const std::string & GetUserSuppliedLabel( unsigned index ) const;	//ASSUMES that argument is a valid index
		unsigned			GetMaxLabelLength() const;
		unsigned			GetSize() const	{return GetNumTaxa();}
		
	protected:
		virtual ~BaseTaxaManager();
		
		void			AlertListeners(NxsTaxaListener::TaxaChangeType);
		
		OrderedCaseInsensitiveLabels taxLabels;
	};
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the index for the taxon label supplied (returns UINT_MAX if the label is not recognized).
*/
inline unsigned BaseTaxaManager::FindIndex(
  const std::string &label) const
	{
	return taxLabels.FindIndex(label, true, false);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the index for the taxon label supplied (returns UINT_MAX if the label is not recognized).
*/
inline unsigned BaseTaxaManager::FindNewOrExistingIndex(
  const std::string &label) const
	{
	return taxLabels.FindIndex(label, true, true);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the index for the taxon label supplied (returns UINT_MAX if the label is not recognized).
*/
inline unsigned BaseTaxaManager::FindIndexFromNamesOnly(
  const std::string &label) const
	{
	return taxLabels.FindIndex(label, false);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the label of taxon i (or throws an exception
*/
inline const std::string & BaseTaxaManager::GetUserSuppliedLabel(
  unsigned i) const
	{
	return taxLabels.GetLabel(i);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the label of taxon i (or throws an exception
*/
inline const std::string BaseTaxaManager::GetLabel(
  unsigned i) const
	{
	return GetUserSuppliedLabel(i);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of characters in the longest taxon name (or number)
*/
inline unsigned BaseTaxaManager::GetMaxLabelLength() const
	{
	return taxLabels.GetMaxLabelLength();
	}


#endif

