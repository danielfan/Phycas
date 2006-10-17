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
#include <list>
#include <ctype.h> 
#include "phypy/src/ncl/nxs_defs.hpp"
#if defined(NCL_USE_NXS_CONSOLE_OUTPUT)
#	include "phypy/src/ncl/misc/algorithm_extensions.hpp"
#	include "phypy/src/ncl/nxs_token.hpp"
#	include "phypy/src/ncl/nxs_basic_manager.hpp"
#	include "phypy/src/ncl/misc/nxs_index_set.hpp"
#	include "phypy/src/ncl/output/nxs_output_stream.hpp"
	using ncl::flush;

	void NxsStdOutputStream::DisplaySet(const NxsIndexSet &s, bool useNumbers, const NxsIndexInfo &indInfo)
		{
		*this << s.name << " = ";
		if (useNumbers)
			*this << s.GetNexusDescription(true);
		else
			{
			for (NxsIndexSet::const_iterator sIt = s.begin(); sIt != s.end(); ++sIt)
				*this << indInfo.GetLabel(*sIt) << ' ';	
			}
		*this << flush;		
		}

	NxsStdOutputStream::NxsStdOutputStream()
	  :dblFormat(UINT_MAX, 6),
	  outputToMonitor(true),
	  logFile(NULL)
		{
		message.clear();
		typist.SetOutputStream(this);
		SetWidth(kDefaultPrintWidth);
		}

#endif // defined(NCL_USE_NXS_CONSOLE_OUTPUT)
