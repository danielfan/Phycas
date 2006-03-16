/*########## gettrees_settings.hpp ##########*/
#if !defined (PHYC_GETTREES_SETTINGS_HPP)
#define PHYC_GETTREES_SETTINGS_HPP
#include "ncl/misc/nxs_index_set.hpp"
class GetTreesSettings
	{
	public:
		NxsIndexSet     	toSave;
		
		GetTreesSettings()
			:toSave(std::string())
			{}
	};
#endif // if !defined (PHYC_GETTREES_SETTINGS_HPP)
/*%%%%%% /gettrees_settings.hpp %%%%%%*/
