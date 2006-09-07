#ifndef NCL_NXSEXCEPTION_H
#define NCL_NXSEXCEPTION_H
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/misc/string_extensions.hpp"
#include <exception>
class NxsToken;
/*----------------------------------------------------------------------------------------------------------------------
|	Thrown in parts of the code that should never be reached.  If they are caught at the highest level of the program 
|	they act like asserts for release version code).
*/
class NxsFatalException: public std::exception
	{
	public: 
		std::string msg;
		NxsFatalException( const std::string &s)	: msg(s) {}
		virtual ~NxsFatalException() throw()
			{
			}
	};
	
/*----------------------------------------------------------------------------------------------------------------------
|	An exception class thrown when an error occurs in reading a Nexus file. Stores a message that will be displayed 
|	to the user and information about where in the file stream the error occurred.
*/
class NxsException: public std::exception
	{
	public:
		std::string	msg;	/*the error to display to the user */
		file_pos 	pos;	/* position in the file */
		unsigned	line;	/* line in the file */
		unsigned 	col;	/* column in the file (characters since the previous endline) */
		CmdResult	cResultCode;

		NxsException( const std::string &s, file_pos fp = 0, unsigned fl = 0U, unsigned fc = 0U, CmdResult c = kCmdFailedGenerateMessage);
		NxsException( const std::string &s, const NxsToken &t, CmdResult c = kCmdFailedGenerateMessage);
		virtual ~NxsException() throw()
			{
			}
		const std::string getMessage() const
			{
			return msg;
			}
		const char * what () const throw ()
			{
			return msg.empty() ? "Unknown Nexus Exception" : msg.c_str();
			}
	};

#if defined (NCL_SUPPORT_OLD_NAMES)
	typedef NxsException XNexus;
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	This class simply makes it easier to report unexpected end of file errors. 
*/
class NxsX_UnexpectedEOF : public NxsException
	{
	public :
		NxsX_UnexpectedEOF(NxsToken &t);
	};

#if defined (NXS_THROW_IF_EMPTY_TOKENS)
	/*----------------------------------------------------------------------------------------------------------------------
	|	This class simply makes it easier to report empty tokens ''
	*/
	class NxsX_EmptyToken : public NxsException
		{
		public :
			NxsX_EmptyToken(const NxsToken &t);
		};
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	This class simply makes it easier to report unknown commands.  The derived class is used in NxsCommandManagerBlock to 
|	skip unknown commands (if desired).  
*/
class NxsX_UnknownCommand : public NxsException
	{
	public:
		NxsX_UnknownCommand( const std::string &s, const NxsToken &t);
	};
	
#if defined (NCL_NXS_THROW_UNDEFINED)
	/*----------------------------------------------------------------------------------------------------------------------
	|	class of exceptions caused by programmer error.  Calling a function in a context (or with argurments) that results
	|	in undefined behaviour).  These do not arise from errors in a NEXUS file.   
	|	Throwing of these exceptions is wrapped in NCL_NXS_THROW_UNDEFINED so they can be turned off.
	|	Thrown using NxsX_UndefinedException("error message", __FILE__, __LINE__) so that the programmer can find the source  
	|	of the error.
	*/
	class NxsX_UndefinedException
		{
		public:
			std::string msg;
			NxsX_UndefinedException(std::string s, const char *file, int line) { msg << s << " (file " << file << ", line" << line << ')'; }
		};
#endif

#endif
