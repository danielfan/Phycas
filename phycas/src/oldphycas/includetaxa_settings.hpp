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

/*########## includetaxa_settings.hpp ##########*/
#if !defined (PHYC_INCLUDETAXA_SETTINGS_HPP)
#define PHYC_INCLUDETAXA_SETTINGS_HPP
#include "phycas/src/ncl/misc/nxs_index_set.hpp"
class IncludeTaxaSettings
	{
	public:
		NxsIndexSet     	taxSet;
		
		IncludeTaxaSettings()
			:taxSet(std::string())
			{}
	};
#endif // if !defined (PHYC_INCLUDETAXA_SETTINGS_HPP)
/*%%%%%% /includetaxa_settings.hpp %%%%%%*/
