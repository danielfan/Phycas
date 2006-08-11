/*########## treeset_settings.hpp ##########*/
#if !defined (PHYC_TREESET_SETTINGS_HPP)
#define PHYC_TREESET_SETTINGS_HPP
#include "ncl/misc/nxs_index_set.hpp"
class TreeSetSettings
	{
	public:
		std::string     	name;
		NxsIndexSet     	treeSet;
		
		TreeSetSettings()
			:name(),
			treeSet(std::string())
			{}
	};
#endif // if !defined (PHYC_TREESET_SETTINGS_HPP)
/*%%%%%% /treeset_settings.hpp %%%%%%*/
