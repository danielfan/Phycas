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

#if !defined(NCL_NXS_STD_USER_QUERY_H)
#define NCL_NXS_STD_USER_QUERY_H

#include "phycas/src/ncl/misc/nxs_file_path.hpp"
#include "phycas/src/ncl/nxs_exception.hpp"
namespace NxsIO
	{
	class NxsInput;
		/// thrown if NxsUserQuery functions are used but the input or output streams are invalid (and 
		/// NxsUserQuery::runWithoutInput has not been specified).
	class XNxsIOError: public NxsException
		{
		public:
			XNxsIOError(const std::string & e)
				:NxsException(e)
				{}
		};
	class UserQueryPrompter;
	}
class NxsUserQuery
	{
	public:
		NxsUserQuery(NxsWritableStream * o,  NxsIO::NxsInput * i);
		void SetOutputStream(NxsWritableStream * o)
			{
			outStream = o;
			}
		void SetInputStream(NxsIO::NxsInput * i)
			{
			inputStream = i;
			}
			
		bool CanCommunicateWithUser() const
			{
			return ((outStream != NULL) && (inputStream != NULL));
			}
			
		bool			AskUserYesNoQuery(const std::string & warningTitle, const std::string & warningMessage) const;
		NxsInFilePath   GetInputFilePath(const std::string & fileNameOrDescrip) const;
		NxsOutFilePath  GetOutputFilePath(const std::string & fileNameOrDescrip) const;
		UInt			UserChoice(const std::string & windowTitle, const std::string & query, const VecString & choices, unsigned defaultResponse, unsigned nonBlockingResponse) const;
		UInt			UserChoice(const std::string & windowTitle, const std::string & query, const std::string & choiceString, unsigned defaultResponse, unsigned nonBlockingResponse) const;
		bool			WarnUser(const std::string & warningTitle, const std::string & warningMessage) const ;
		
	private:
		template<typename T, class CMD_OPT>
		T				GenericUserQuery(const NxsIO::UserQueryPrompter & , T nonBlockingResponse, CMD_OPT & cmdOpt) const;
		NxsWritableStream * outStream;		/// alias to the outstream, used to prompt the user
		NxsIO::NxsInput   * inputStream;
		bool				runWithoutInput;/// if True, the default response will be used without querying */
	};


#endif
