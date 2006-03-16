#if !defined MSG_EXCEPTION_H
#define MSG_EXCEPTION_H

#include "ncl/nxs_defs.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	Generic exception class that is like an assert (will only be caught at a high level leading to program termination).
|	Should not be used for conditions that are expected to occur in normal program use
*/
class MsgException
	{
	public:
		std::string msg;
		MsgException( const char * s, const char *file, int line);
		MsgException(const std::string &s, const char *file, int line);
	};

/*----------------------------------------------------------------------------------------------------------------------
|	THROW_MSG_EXCEPTION is a simple macro that adds, __FILE__ , and __LINE__ to the message
*/
#define THROW_MSG_EXCEPTION(mString) throw MsgException(mString, __FILE__, __LINE__)

#endif
