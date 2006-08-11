//#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/taxa/nxs_taxa_manager.hpp"	//@POL Mark, VC .NET needed this
#include "ncl/characters/nxs_characters_manager.hpp"
#include "ncl/characters/nxs_char_listener.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	Stores alias to the NxsCharactersManager and calls the NxsCharactersManager::AddListener
*/
NxsCharListener::NxsCharListener(NxsCharactersManager *t)
	:charMgr(t)
	{
	if (charMgr != NULL)
		charMgr->AddListener(this);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	removes itself from the old NxsCharactersManager list of listeners (via NxsCharactersManager::RemoveListener) and add itself to 
|	the new one	Manager's list (NxsCharactersManager::AddListener)
*/
void NxsCharListener::ChangeManager(NxsCharactersManager *t)
	{
	if (charMgr != NULL)
		charMgr->RemoveListener(this);
	charMgr = t;
	if (charMgr != NULL)
		charMgr->AddListener(this);
	}
	
	
	
/*----------------------------------------------------------------------------------------------------------------------
|	removes itself from the NxsCharactersManager list of listeners (via RemoveListener)
*/
NxsCharListener::~NxsCharListener()
	{
	if (charMgr != NULL)
		charMgr->RemoveListener(this);
	}
