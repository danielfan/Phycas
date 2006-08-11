/*########## alias_settings.hpp ##########*/
#if !defined (PHYC_ALIAS_SETTINGS_HPP)
#define PHYC_ALIAS_SETTINGS_HPP
class AliasSettings
	{
	public:
		std::string     	aliasKey;
		std::string     	expansion;
		
		AliasSettings()
			:aliasKey(),
			expansion()
			{}
	};
#endif // if !defined (PHYC_ALIAS_SETTINGS_HPP)
/*%%%%%% /alias_settings.hpp %%%%%%*/
