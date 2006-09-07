#ifndef PHO_TAXA_MANAGER_H
#define PHO_TAXA_MANAGER_H

#include "phypy/src/ncl/taxa/nxs_taxa_manager.hpp"
#include "phypy/src/oldphycas/taxstatus_settings.hpp"
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

