/*########## distribution_settings.hpp ##########*/
#if !defined (PHYC_DISTRIBUTION_SETTINGS_HPP)
#define PHYC_DISTRIBUTION_SETTINGS_HPP
#include "phycas/rand/distribution_description.hpp"
class DistributionSettings
	{
	public:
		std::string     	name;
		DistributionDescription	distrib;
		
		DistributionSettings()
			:name(),
			distrib()
			{}
	};
#endif // if !defined (PHYC_DISTRIBUTION_SETTINGS_HPP)
/*%%%%%% /distribution_settings.hpp %%%%%%*/
