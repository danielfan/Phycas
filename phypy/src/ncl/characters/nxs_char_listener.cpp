/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/taxa/nxs_taxa_manager.hpp"	//@POL Mark, VC .NET needed this
#include "phypy/src/ncl/characters/nxs_characters_manager.hpp"
#include "phypy/src/ncl/characters/nxs_char_listener.hpp"

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
