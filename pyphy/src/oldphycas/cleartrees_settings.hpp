/*########## cleartrees_settings.hpp ##########*/
#if !defined (PHYC_CLEARTREES_SETTINGS_HPP)
#define PHYC_CLEARTREES_SETTINGS_HPP
#include "pyphy/src/ncl/misc/nxs_index_set.hpp"
class ClearTreesSettings
	{
	public:
		NxsIndexSet     	toClear;
		
		ClearTreesSettings()
			:toClear(std::string())
			{}
	};
#endif // if !defined (PHYC_CLEARTREES_SETTINGS_HPP)
/*%%%%%% /cleartrees_settings.hpp %%%%%%*/
