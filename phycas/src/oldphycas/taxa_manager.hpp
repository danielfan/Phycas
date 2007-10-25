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

#ifndef PHO_TAXA_MANAGER_H
#define PHO_TAXA_MANAGER_H

#include "phycas/src/ncl/taxa/nxs_taxa_manager.hpp"
#include "phycas/src/oldphycas/taxstatus_settings.hpp"
class ExcludeTaxaSettings;
class IncludeTaxaSettings;
class TaxSetSettings;

/*----------------------------------------------------------------------------------------------------------------------
|	provides include/exclude and taxset and taxstatus capabilities to the NxsTaxaManager 
*/
class PhoTaxaManager : public NxsTaxaManager
	{
	public: 
		CmdResult			ReplaceTaxa(const VecString &v);
		const NxsIndexSet & GetActiveTaxa() const {return GetActiveSet();}
		const VecString	  & GetReservedSetNames() const
			{
			return reservedTaxSetNames;
			}
		CmdResult			HandleExcludeTaxa(ExcludeTaxaSettings *);
		CmdResult			HandleIncludeTaxa(IncludeTaxaSettings *);
		CmdResult			HandleTaxStatus(TaxStatusSettings *);
		CmdResult			HandleTaxSet(TaxSetSettings *);

		void				SetupExcludeTaxaCommand(NxsCommandManager *cmdMgr);
		void				SetupIncludeTaxaCommand(NxsCommandManager *cmdMgr);
		void				SetupTaxStatusCommand(NxsCommandManager *cmdMgr);
		void				SetupTaxSetCommand(NxsCommandManager *cmdMgr);
		
		PhoTaxaManager(){} //wrap_scm.cpp is the only non-singelton instantiation
	private:
	
		VecString			reservedTaxSetNames;	// was STATIC_DATA when the Manager was not a Singleton
		TaxStatusSettings 	currentTaxStatSettings; /* settings that control how the taxon status will be displayed */
	};

std::string & AppendTaxonLabel(std::string & s, unsigned int taxNumber, PhoTaxaManager & taxMgr, bool useTaxonNames);

inline std::string & AppendTaxonLabel(std::string & s, unsigned int taxNumber, PhoTaxaManager & taxMgr, bool useTaxonNames)
	{
	if (useTaxonNames) //@might need another flag for NEXUS (currently useTaxonNames will print names as single NEXUS tokens
		return s <<  GetAsNexusToken(taxMgr.GetLabel(taxNumber));
	return s << taxNumber; 
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls NxsTaxaManager::ReplaceTaxa and then displays TaxStatus to the user
*/
inline CmdResult PhoTaxaManager::ReplaceTaxa(const VecString &v)
	{
	CmdResult r =  NxsTaxaManager::ReplaceTaxa(v);
	if (r == kCmdSucceeded)
		r = HandleTaxStatus(&currentTaxStatSettings);
	return r;
	}

#endif

