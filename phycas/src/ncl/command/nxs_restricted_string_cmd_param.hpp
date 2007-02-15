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

#ifndef NCL_NXS_RESTRICTED_STRING_CMD_OPTION_H
#define NCL_NXS_RESTRICTED_STRING_CMD_OPTION_H

#include "phycas/src/ncl/nxs_defs.hpp"
#include "phycas/src/ncl/misc/nxs_index_set.hpp"
#include "phycas/src/ncl/command/nxs_cmd_param.hpp"
#include <boost/function.hpp>
/*----------------------------------------------------------------------------------------------------------------------
|	a NxsStringCmdOption that does not allow the value to be a single punctuation character (the punctuation list is 
|	(){}"'-[]/\,;:=*`+<>_	
*/
class NxsWordCmdOption : public NxsStringCmdOption
	{
	public :
		virtual bool	IsCurrentlyValid();
		NxsWordCmdOption( const std::string &n, std::string *manipVal, std::string def, bool validityIsLabile, bool persist, CmdPermissionLevel pLevel );
	};	

/*----------------------------------------------------------------------------------------------------------------------
|	A NxsWordCmdOption that will not accept some reserved words.  
|	The list of illegal words can be supplied as a vector of strings, a |-separated string, or (if the list changes)
|	as a shared_ptr to a ValueHolder<VecString>
*/
class RestrictNameCmdOption : public NxsWordCmdOption
	{
	public:
		typedef boost::function0<VecString> VecStringSource;
				RestrictNameCmdOption(const std::string &n, std::string *manipVal, std::string def, VecStringSource tabooProvider, bool persist, CmdPermissionLevel pLevel );
				RestrictNameCmdOption(const std::string &n, std::string *manipVal, std::string def, VecString tabooVec, bool persist, CmdPermissionLevel pLevel );
				RestrictNameCmdOption(const std::string &n, std::string *manipVal, std::string def, const char *tabooStr, bool persist, CmdPermissionLevel pLevel );
		bool	IsCurrentlyValid();
		void	WriteTypeInfoStateElement(NxsHiddenQueryStream & outStream) const;		
	protected:
		VecStringSource tabooListProvider;
		mutable VecString 		tabooList;
	private:
		inline void		RefreshTabooList() const;
	};

typedef boost::shared_ptr< RestrictNameCmdOption > RestrictNameCmdOptionShPtr;

inline void RestrictNameCmdOption::RefreshTabooList() const
	{
	if (tabooListProvider)
		tabooList = tabooListProvider();
	}


#endif

