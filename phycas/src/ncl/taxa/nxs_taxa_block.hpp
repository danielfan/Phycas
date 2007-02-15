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

#ifndef NCL_NXSTAXABLOCK_H
#define NCL_NXSTAXABLOCK_H

#include "phycas/src/ncl/nxs_defs.hpp"
#include "phycas/src/ncl/nxs_basic_manager.hpp"
#include "phycas/src/ncl/nxs_block.hpp"
#include "phycas/src/ncl/taxa/nxs_taxa_label_reader.hpp"
class NxsTaxaManager;
/*----------------------------------------------------------------------------------------------------------------------
|	Class that reads taxa blocks and stores the list of taxon labels.
*/
class NxsTaxaBlock : public NxsCommandManagerBlock , public NxsIndexInfo, public NxsTaxaLabelReader
	{
	public:
			NxsTaxaBlock(NxsTaxaManager & taxaMgr);
	
		//	Accessors
		//
		const VecString   & GetTaxa()const;
		unsigned			FindIndex(const std::string &label ) const;
		unsigned			FindIndexFromNamesOnly(const std::string &label ) const;
		unsigned			GetMaxLabelLength() const;
		unsigned			GetNumTaxa() const;
		unsigned			GetSize() const;
		const std::string & GetUserSuppliedLabel( unsigned i ) const;
		
		//	Modifiers
		//
		void	  	AddLabel(const std::string &s);
		void  		ChangeLabel( unsigned i, const std::string &s );
		void		Clear();
		void 		ResetCmdMgrNxsBlock();
		
	protected:
				
		void		AddRecognizedCommands();
		CmdResult	EndEncountered();
		bool 		ParseDimensions(NxsToken &);	
		NxsTaxaManager & taxaMgr;
	};
	
#if defined (NCL_SUPPORT_OLD_NAMES)
	typedef NxsTaxaBlock TaxaBlock;
#endif

inline const VecString   &NxsTaxaBlock::GetTaxa() const
	{
	return taxLabels.GetLabels();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds a new taxon label (to the end of the list) and increments the number of taxa, throws NxsX_IllegalLabel
*/
inline void NxsTaxaBlock::AddLabel(
  const std::string &s)
	{
	return taxLabels.push_back(s);
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of taxa;
*/
inline unsigned NxsTaxaBlock::GetNumTaxa() const
	{
	return taxLabels.size();
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	returns the number of taxa
*/
inline unsigned NxsTaxaBlock::GetSize() const
	{
	return GetNumTaxa();
	}

				
/*----------------------------------------------------------------------------------------------------------------------
|	clears all taxa labels, sets the number of taxa to zero and calls NxsCommandManagerBlock::Reset();
*/
inline void NxsTaxaBlock::ResetCmdMgrNxsBlock()
	{
	taxLabels.clear();
	labelsRead = false;
	}

inline void NxsTaxaBlock::Clear()
	{
	Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the index for the taxon label supplied (returns UINT_MAX if the label is not recognized).
*/
inline unsigned NxsTaxaBlock::FindIndex(
  const std::string &label) const
	{
	return taxLabels.FindIndex(label, true, false);
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the index for the taxon label supplied (returns UINT_MAX if the label is not recognized).
*/
inline unsigned NxsTaxaBlock::FindIndexFromNamesOnly(
  const std::string &label) const
	{
	return taxLabels.FindIndex(label, false);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the label of taxon i (or throws an exception
*/
inline const std::string & NxsTaxaBlock::GetUserSuppliedLabel(
  unsigned i) const
	{
	return taxLabels.GetLabel(i);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of characters in the longest taxon name (or number)
*/
inline unsigned NxsTaxaBlock::GetMaxLabelLength() const
	{
	return taxLabels.GetMaxLabelLength();
	}

#endif
