/*########## phylodiversity_settings.hpp ##########*/
#if !defined (PHYC_PHYLODIVERSITY_SETTINGS_HPP)
#define PHYC_PHYLODIVERSITY_SETTINGS_HPP
#include "ncl/output/nxs_output_destination_description.hpp"
class PhylodiversitySettings
	{
	public:
		std::string     	taxonset;
		NxsOutputDestinationDescription	pdOut;
		
		PhylodiversitySettings()
			:taxonset("None"),
			pdOut()
			{}
	};
#endif // if !defined (PHYC_PHYLODIVERSITY_SETTINGS_HPP)
/*%%%%%% /phylodiversity_settings.hpp %%%%%%*/
