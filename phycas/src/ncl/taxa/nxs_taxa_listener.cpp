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

//#include "phycas/force_include.h"
#include "phycas/src/ncl/nxs_defs.hpp"
#include "phycas/src/ncl/taxa/base_taxa_manager.hpp"
#include "phycas/src/ncl/taxa/nxs_taxa_listener.hpp"

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

