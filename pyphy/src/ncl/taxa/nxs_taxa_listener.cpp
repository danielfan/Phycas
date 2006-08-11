//#include "phycas/force_include.h"
#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/taxa/base_taxa_manager.hpp"
#include "pyphy/src/ncl/taxa/nxs_taxa_listener.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	Stores alias and calls the NxsTaxaManager::AddListener
*/
NxsTaxaListener::NxsTaxaListener(BaseTaxaManager *t)
	:taxaMgr(t)
	{
	if (taxaMgr != NULL)
		taxaMgr->AddListener(this);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	removes itself from the old NxsTaxaManager list of listeners (via NxsTaxaManager::RemoveListener) and add itself to 
|	the new one	Manager's list (NxsTaxaManager::AddListener)
*/
void NxsTaxaListener::ChangeManager(BaseTaxaManager *t)
	{
	if (taxaMgr != NULL)
		taxaMgr->RemoveListener(this);
	taxaMgr = t;
	if (taxaMgr != NULL)
		taxaMgr->AddListener(this);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	removes itself from the NxsTaxaManager list of listeners (via RemoveListener)
*/
NxsTaxaListener::~NxsTaxaListener()
	{
	if (taxaMgr != NULL)
		taxaMgr->RemoveListener(this);
	}

