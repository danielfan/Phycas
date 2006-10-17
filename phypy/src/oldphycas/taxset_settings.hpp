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

/*########## taxset_settings.hpp ##########*/
#if !defined (PHYC_TAXSET_SETTINGS_HPP)
#define PHYC_TAXSET_SETTINGS_HPP
#include "phypy/src/ncl/misc/nxs_index_set.hpp"
class TaxSetSettings
	{
	public:
		std::string     	name;
		NxsIndexSet     	taxSet;
		
		TaxSetSettings()
			:name(),
			taxSet(std::string())
			{}
	};
#endif // if !defined (PHYC_TAXSET_SETTINGS_HPP)
/*%%%%%% /taxset_settings.hpp %%%%%%*/
