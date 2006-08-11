#ifndef NCL_NXS_CHOICE_CMD_OPTION_H
#define NCL_NXS_CHOICE_CMD_OPTION_H

#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/command/nxs_primitive_cmd_param.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	A NxsStringCmdOption that will only accept a few words.  
|	The list of legal choices can be supplied as a vector of strings, a |-separated string, or (if the list changes)
|	as a shared_ptr to a ValueHolder<VecString>
|	Note that this class takes a string and a unsigned number (the unsigned will hold the index of the choice)
*/
class NxsChoiceCmdOption : public NxsStringCmdOption
	{
	public:
		bool		AbbreviateChoices();
		std::string	GetDisplayType(bool includeIndefiniteArticle = false, bool plural = false) const;
		VecString   GetLegalChoices() const;
		std::string	GetChosenString() const
			{
			return *currentValue;
			}
		unsigned	GetValue() const
			{
			return *choiceIndex;
			}
		VecString   GetValidArgument()
			{
			RefreshChoices();
			return choices;
			}
		bool		IsCurrentlyValid();
		void SetAllowNumbers(bool n = true)
			{
			allowNumbers = n;
			}
		virtual void WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;
		
		typedef boost::function0<VecString> VecStringSource;
		
		NxsChoiceCmdOption(const std::string &n, unsigned *i, std::string *manipVal, std::string def, VecStringSource choiceProvider,   bool allowAbbreviations, bool persist, CmdPermissionLevel pLevel );
		NxsChoiceCmdOption(const std::string &n, unsigned *i, std::string *manipVal, std::string def, const VecString & choiceVec,		bool allowAbbreviations, bool persist, CmdPermissionLevel pLevel );
		NxsChoiceCmdOption(const std::string &n, unsigned *i, std::string *manipVal, std::string def, const char * choiceStr,			bool allowAbbreviations, bool persist, CmdPermissionLevel pLevel );
	
	private :
		bool 			IsValidChoice(const std::string &s);
		bool			PrepareAbbreviations() const;//const because choices is mutable (and called by many functions that should be const for CommandOptions)
		void 			RefreshChoices() const;	//const because choices is mutable (and called by many functions that should be const for CommandOptions)
		
		VecStringSource 	choiceProvider;
		unsigned		  * choiceIndex;
		mutable VecString 	choices;	//choice list may be coming from a CmdOptProvider so it is mutable
		const bool			allowAbbrevChoices;
		bool				allowNumbers; /// let the user enter a 1-based number for the choice (Warning:  numeric interpretation will be given priority in the parsing of the command)
	};
	
typedef boost::shared_ptr< NxsChoiceCmdOption > NxsChoiceCmdOptionShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	returns true if the setting is one of the legal choices
*/
inline bool NxsChoiceCmdOption::IsCurrentlyValid()
	{
	RefreshChoices();
	return IsValidChoice(*currentValue);
	}
	

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsChoiceCmdOption::AbbreviateChoices()
	{
	RefreshChoices();
	return PrepareAbbreviations();
	}

inline VecString NxsChoiceCmdOption::GetLegalChoices() const
	{
	RefreshChoices();
	return choices;
	}
	
inline void NxsChoiceCmdOption::RefreshChoices() const // logically const , choices are mutable (might need to change this it doesn't seem like this should be const)
	{
	if (choiceProvider)
		{
		choices = choiceProvider();
		PrepareAbbreviations();
		}
	}

	

#endif
