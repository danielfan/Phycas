#if !defined(NCL_NXS_STD_USER_QUERY_H)
#define NCL_NXS_STD_USER_QUERY_H

#include "pyphy/src/ncl/misc/nxs_file_path.hpp"
#include "pyphy/src/ncl/nxs_exception.hpp"
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
