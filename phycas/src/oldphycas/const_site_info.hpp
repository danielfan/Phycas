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

#ifndef PHO_CONST_SITE_INFO_H
#define PHO_CONST_SITE_INFO_H

#include "phycas/src/ncl/misc/nxs_data_type.hpp"
#include "phycas/src/oldphycas/compressed_2d_array.hpp"
/*--------------------------------------------------------------------------------------------------------------------------
|	simple class for holding info about the constant sites in a subset.  
|	This information is needed when efficiently scoring trees with models that include a proportion of constant sites.
|	In essence this class is packing of all of the sites that can be explained by an invariant model into groups based on
|	which states are present in every taxon.  For most DNA datasets there will be 4 groups of constant sites (one for each
|	base), but if there is substantial ambiguity there might be more (sites that could be explained by constant A or constant
|	G, for instance).
*/
class ConstSiteInfo
	{
	public:
		unsigned 				nConstStateGrps;	/* the number of constant site groups in subset that this object is associated with */
		OrdCodedArrs		 	csStates;			/* Ordinally coding of the states that could be constant */
		Compressed2DArray<int> 	patternIndices;			/* Ordinally coded list of indices of patterns for each of the constant site groups */
		
		ConstSiteInfo(unsigned nStates)	
			: nConstStateGrps(0),
			csStates(nStates)
			{}
	};

#endif
