/*########## taxstatus_settings.hpp ##########*/
#if !defined (PHYC_TAXSTATUS_SETTINGS_HPP)
#define PHYC_TAXSTATUS_SETTINGS_HPP
class TaxStatusSettings
	{
	public:
		bool            	fullDisplay;
		bool            	showExcluded;
		
		TaxStatusSettings()
			:fullDisplay(true),
			showExcluded(true)
			{}
	};
#endif // if !defined (PHYC_TAXSTATUS_SETTINGS_HPP)
/*%%%%%% /taxstatus_settings.hpp %%%%%%*/
