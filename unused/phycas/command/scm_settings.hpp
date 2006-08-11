/*########## scm_settings.hpp ##########*/
#if !defined (PHYC_SCM_SETTINGS_HPP)
#define PHYC_SCM_SETTINGS_HPP
#include "ncl/output/nxs_output_destination_description.hpp"
class SCMSettings
	{
	public:
		NxsOutputDestinationDescription	scmfile;
		
		SCMSettings()
			:scmfile()
			{}
	};
#endif // if !defined (PHYC_SCM_SETTINGS_HPP)
/*%%%%%% /scm_settings.hpp %%%%%%*/
