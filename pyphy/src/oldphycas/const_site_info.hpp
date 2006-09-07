#ifndef PHO_CONST_SITE_INFO_H
#define PHO_CONST_SITE_INFO_H

#include "phypy/src/ncl/misc/nxs_data_type.hpp"
#include "phypy/src/oldphycas/compressed_2d_array.hpp"
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
