/*########## excludetaxa_settings.hpp ##########*/
#if !defined (PHYC_EXCLUDETAXA_SETTINGS_HPP)
#define PHYC_EXCLUDETAXA_SETTINGS_HPP
#include "phypy/src/ncl/misc/nxs_index_set.hpp"
class ExcludeTaxaSettings
	{
	public:
		NxsIndexSet     	taxSet;
		
		ExcludeTaxaSettings()
			:taxSet(std::string())
			{}
	};
#endif // if !defined (PHYC_EXCLUDETAXA_SETTINGS_HPP)
/*%%%%%% /excludetaxa_settings.hpp %%%%%%*/
