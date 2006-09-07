/*########## includetaxa_settings.hpp ##########*/
#if !defined (PHYC_INCLUDETAXA_SETTINGS_HPP)
#define PHYC_INCLUDETAXA_SETTINGS_HPP
#include "phypy/src/ncl/misc/nxs_index_set.hpp"
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
