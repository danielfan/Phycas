#if !defined (NCL_NXS_FILE_CMD_OPTIONS_H)
#define NCL_NXS_FILE_CMD_OPTIONS_H

//#include "phycas/force_include.h"
#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/misc/nxs_file_path.hpp"
#include "pyphy/src/ncl/command/nxs_cmd_param.hpp"

class NxsCommand;

void AssignAllFieldsExceptAppAndReplace(NxsOutFilePath *leftOP, const NxsOutFilePath &rightOp);  //this equals operator should be used for NxsOutFilePath in the CmdOption context

#include "pyphy/src/ncl/output/nxs_output_destination_description.hpp"
class NxsOutputCmdOption : public SimpleCmdOptionInterface<NxsOutputDestinationDescription>
	{
	public:
		NxsOutputCmdOption(const std::string & name, NxsOutputDestinationDescription * manip, NxsOutputDestinationDescription def, bool persistent, CmdPermissionLevel perm)
			:SimpleCmdOptionInterface<NxsOutputDestinationDescription>(name, manip, def, true, persistent, perm)
			{
			ReturnValueToDefault();
			}
		std::string		GetCurrentValueAsString() const;
		std::string		GetDisplayType(bool includeIndefiniteArticle,  bool plural) const;
		//std::string		GetHandlersChangesToCmd() {return std::string();}	//@POL added body to see why VC reports link error
		bool			IsCurrentlyValid();
		void			StorePreviousValue();
		bool			ReadValue(NxsToken &, bool equalsAlreadyRead = false);
		void			ReturnValueToDefault();
		void			RevertValueBecauseCommandFailed();
		virtual void	WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;		
	private:
		bool InterpretRedirectString(const std::string & s);

		unsigned		previousSourceOfSetting;
	};

class NxsOutFileCmdOption : public SimpleCmdOptionInterface<NxsOutFilePath>
	{
	public:
		NxsOutFileCmdOption(	const std::string & name, 
									const std::string & prefix, 
									NxsOutFilePath *manip, 
									NxsOutFilePath def, 
									NxsCommand *cmd, 
									bool persis, 
									CmdPermissionLevel perm);
		std::string		GetCurrentValueAsString() const;
		std::string		GetDisplayType(bool includeIndefiniteArticle,  bool plural) const;
		std::string		GetHandlersChangesToCmd();
		bool			IsCurrentlyValid();
		void			StorePreviousValue();
		bool			ReadValue(NxsToken &, bool equalsAlreadyRead = false);
		void			ReturnValueToDefault();
		void			RevertValueBecauseCommandFailed();
		virtual void	WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;

	protected:
		unsigned		previousSourceOfSetting;
		std::string		prefixToOptions;
	};


class NxsInFileCmdOption : public SimpleCmdOptionInterface<NxsInFilePath>
	{
	public:
		NxsInFileCmdOption(	const std::string & n, 
							NxsInFilePath * manip, 
							NxsInFilePath def, 
							bool persis, 
							CmdPermissionLevel perm);
		std::string GetCurrentValueAsString() const;
		std::string GetDisplayType(bool includeIndefiniteArticle,  bool plural) const;
		bool		IsCurrentlyValid();
		void		StorePreviousValue();
		bool		ReadValue(NxsToken &, bool equalsAlreadyRead = false);
		void		ReturnValueToDefault();
		void 		RevertValueBecauseCommandFailed();
		void		WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;
	
	protected:
		unsigned		previousSourceOfSetting;
	};

inline std::string NxsOutFileCmdOption::GetCurrentValueAsString() const
	{
	return currentValue->GetFullName();
	}
	
inline std::string NxsInFileCmdOption::GetCurrentValueAsString() const
	{
	return currentValue->GetFullName();
	}


	
#endif

