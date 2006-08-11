/*########## execute_settings.hpp ##########*/
#if !defined (PHYC_EXECUTE_SETTINGS_HPP)
#define PHYC_EXECUTE_SETTINGS_HPP
#include "ncl/misc/nxs_file_path.hpp"
class ExecuteSettings
	{
	public:
		NxsInFilePath   	fileToExecute;
		
		ExecuteSettings()
			:fileToExecute()
			{}
	};
#endif // if !defined (PHYC_EXECUTE_SETTINGS_HPP)
/*%%%%%% /execute_settings.hpp %%%%%%*/
