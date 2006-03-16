/*########## taxset_settings.hpp ##########*/
#if !defined (PHYC_TAXSET_SETTINGS_HPP)
#define PHYC_TAXSET_SETTINGS_HPP
#include "ncl/misc/nxs_index_set.hpp"
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
