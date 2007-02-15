/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#if !defined (NCL_NXS_FILE_CMD_OPTIONS_H)
#define NCL_NXS_FILE_CMD_OPTIONS_H

//#include "phycas/force_include.h"
#include "phycas/src/ncl/nxs_defs.hpp"
#include "phycas/src/ncl/misc/nxs_file_path.hpp"
#include "phycas/src/ncl/command/nxs_cmd_param.hpp"

class NxsCommand;

void AssignAllFieldsExceptAppAndReplace(NxsOutFilePath *leftOP, const NxsOutFilePath &rightOp);  //this equals operator should be used for NxsOutFilePath in the CmdOption context

#include "phycas/src/ncl/output/nxs_output_destination_description.hpp"
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

