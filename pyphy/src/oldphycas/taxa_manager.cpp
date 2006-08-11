//#include "phycas/force_include.h"
#include "pyphy/src/oldphycas/taxa_manager.hpp"
#include "pyphy/src/oldphycas/taxset_settings.hpp"
#include "pyphy/src/oldphycas/excludetaxa_settings.hpp"
#include "pyphy/src/oldphycas/includetaxa_settings.hpp"
#include "pyphy/src/ncl/output/nxs_output.hpp"
#include "pyphy/src/ncl/nxs_exception.hpp"
using ncl::flush;
using std::string;
inline void PrintInclExclStatus(NxsOutputStream &outStream, unsigned nActive, unsigned nExcluded);

/*----------------------------------------------------------------------------------------------------------------------
|	Prints number of active and excluded taxa.
*/
void PrintInclExclStatus(NxsOutputStream &outStream, unsigned nActive, unsigned nExcluded)
	{
	if (nActive == 1)
		outStream << "1 taxon is active, " << nExcluded << (nExcluded == 1 ? " has" : " have") << " been excluded.\n";
	else
		outStream << nActive << " taxa are active, " << nExcluded << (nExcluded == 1 ? " has" : " have") << " been excluded.\n";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	sets the currentTaxStatSettings equal to the settings in `s' and then displays a summary of the
|	taxa using these settings to control how verbose the output is.
*/	
CmdResult PhoTaxaManager::HandleTaxStatus(TaxStatusSettings *s)
	{
	currentTaxStatSettings = *s;
	NxsOutputStream * outStreamPtr = NxsOutputManager::GetInstance().GetOutputStreamPtr();
	if (outStreamPtr != NULL)
		{
		NxsOutputStream & outStream = *outStreamPtr;
		const unsigned ntax = GetNumTaxa();
		if (ntax == 1)
			outStream << "There is currently 1 taxon.\n";
		else
			outStream << "There are currently " << ntax << " taxa.\n";
		if (currentTaxStatSettings.fullDisplay)
			{
#			if defined(NCL_HAS_TABULAR_OUTPUT)
				NxsTable *table = outStream.GetTablePtr();
				table->RightJustifyColumn();
				table->AddString("#");
				table->SetLeftMargin();
				table->AddString("Name", GetMaxLabelLength());
				table->CenterColumn();
				bool someExcluded = (GetNumActive() != GetSize());
				if (currentTaxStatSettings.showExcluded && someExcluded)
					table->AddString("Status");
				table->SetRightMargin();
				table->HyphenLine();
				table->SetTopMargin();
				
				for (unsigned i = 0; i < ntax; ++i)
					{
					table->AddUInt(i+1);
					table->AddString(this->GetLabel(i));
					if (currentTaxStatSettings.showExcluded)
						{
						if (IsActive(i))
							{
							if (someExcluded)
								table->AddString(string());
							}
						else 
							table->AddString("Excluded");
						}
					}
				table->SetBottomMargin();
				outStream.PrintTable(table);
				//DisplayTaxSets(true);
#			endif
			}
		PrintInclExclStatus(outStream, GetNumActive(), GetNumInactive());
		outStream << ncl::endl;
		}
	return kCmdSucceeded;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Adds a new TaxSet.  Throws a NxsException if the name for the set is illegal
*/	
CmdResult PhoTaxaManager::HandleTaxSet(TaxSetSettings *s)
	{
	if (IsALong(s->name, NULL))
		throw NxsException("TaxSet names cannot be numbers");
	NStrCaseInsensitiveEquals compF(s->name);
	if (find_if(reservedTaxAndSetNames.begin(), reservedTaxAndSetNames.end(), compF) != reservedTaxAndSetNames.end())
		{
		string st;
		st << s->name << " is a reserved TaxSet name";
		throw NxsException(st);
		}
	if (this->FindIndex(s->name) != UINT_MAX)
		{
		string st;
		st << s->name << " cannot be used as a TaxSet name because it is a taxon label";
		throw NxsException(st);
		}
	s->taxSet.SetName(s->name);
	AddSet(s->taxSet);
	//DisplayTaxSet(s->taxSet, true);
	return kCmdSucceeded;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Flags specified taxa as excluded, issues a brief report and alerts listeners that kTaxaExcluded
*/	
CmdResult PhoTaxaManager::HandleExcludeTaxa(ExcludeTaxaSettings *s)
	{
	unsigned prevNActive = GetNumActive();
	SetActiveStatus(s->taxSet, false);
	NxsOutputStream * outStreamPtr = NxsOutputManager::GetInstance().GetOutputStreamPtr();
	if (outStreamPtr != NULL)
		{
		NxsOutputStream & outStream = *outStreamPtr;
		outStream << prevNActive - GetNumActive() << " taxa were excluded.  ";
		PrintInclExclStatus(outStream, GetNumActive(), GetNumInactive());
		outStream << flush;
		}
	AlertListeners(NxsTaxaListener::kTaxaExcluded);	
	return kCmdSucceeded;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Flags specified taxa as active, issues a brief report and alerts listeners that kTaxaIncluded
*/	
CmdResult PhoTaxaManager::HandleIncludeTaxa(IncludeTaxaSettings *s)
	{
	unsigned prevNActive = GetNumActive();
	SetActiveStatus(s->taxSet, true);
	NxsOutputStream * outStreamPtr = NxsOutputManager::GetInstance().GetOutputStreamPtr();
	if (outStreamPtr != NULL)
		{
		NxsOutputStream & outStream = *outStreamPtr;
		outStream << GetNumActive() - prevNActive << " taxa were re-included.  ";
		PrintInclExclStatus(outStream, GetNumActive(), GetNumInactive());
		outStream << flush;
		}
	AlertListeners(NxsTaxaListener::kTaxaIncluded);	
	return kCmdSucceeded;
	}
