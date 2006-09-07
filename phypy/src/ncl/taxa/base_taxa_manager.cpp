//#include "phycas/force_include.h"
#include <boost/bind.hpp>
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/taxa/base_taxa_manager.hpp"
#include "phypy/src/ncl/misc/nxs_copy.hpp"
using std::string;
BaseTaxaManager::~BaseTaxaManager()	
	{
	Clear();
	}

CmdResult BaseTaxaManager::AppendNewTaxon(const string &s)
	{
	VecString v(1, s);
	AppendTaxa(v);
	return kCmdSucceeded;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Called internally to add taxa, alert listeners and update the "all" taxset
*/
CmdResult BaseTaxaManager::AppendTaxa(const VecString &v)
	{
	if (!v.empty())
		{
		taxLabels.AppendLabels(v);
		IndicesIncreased(GetSize(), true);
		AlertListeners(NxsTaxaListener::kTaxaAdded);
		}
	return kCmdSucceeded;
	}

void BaseTaxaManager::AlertListeners(NxsTaxaListener::TaxaChangeType ct)
	{
	std::for_each(listeners.begin(), listeners.end(), boost::bind(&NxsTaxaListener::TaxaChanged, _1, this, ct));
	}
	

CmdResult BaseTaxaManager::ReplaceTaxa(const VecString &v)
	{
	if (!taxLabels.empty())
		Clear();
	return AppendTaxa(v);
	}
