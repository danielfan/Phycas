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
