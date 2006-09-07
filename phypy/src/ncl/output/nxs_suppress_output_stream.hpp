#ifndef NCL_SUPPRESS_OUTSTREAM_H
#define NCL_SUPPRESS_OUTSTREAM_H

#include "phypy/src/ncl/nxs_defs.hpp"
#if defined (SUPPORT_TREE_OUTPUT)
	class Tree;
#endif
/*----------------------------------------------------------------------------------------------------------------------
|	Provides "put-to" operator interface that does nothing
*/
class NxsSuppressOutputStream
	{
	public :
		void		SetWidth(UInt) {}
		void		SetNDigitsAfterDecimal(UInt) {}
		unsigned	GetWidth() const	
			{
			return 80U; //bogus value here that is unlikely to break code expecting a reasonable number
			}
	};
// The NxsOutputStream Interface
//
NxsSuppressOutputStream 	& operator<<(NxsSuppressOutputStream & , int i);	
NxsSuppressOutputStream 	& operator<<(NxsSuppressOutputStream & , long l);	
NxsSuppressOutputStream 	& operator<<(NxsSuppressOutputStream & , unsigned l);	
NxsSuppressOutputStream 	& operator<<(NxsSuppressOutputStream & , unsigned long l);	
NxsSuppressOutputStream 	& operator<<(NxsSuppressOutputStream & , double d);
NxsSuppressOutputStream 	& operator<<(NxsSuppressOutputStream & , const char *c);
NxsSuppressOutputStream 	& operator<<(NxsSuppressOutputStream & , char c);
NxsSuppressOutputStream 	& operator<<(NxsSuppressOutputStream & , const std::string &s);
NxsSuppressOutputStream & operator<<(NxsSuppressOutputStream & , NxsWritableStreamManipPtr);
#if defined (SUPPORT_TREE_OUTPUT)
	NxsSuppressOutputStream & operator<<(NxsSuppressOutputStream & , const Tree &);
#endif

namespace ncl
{
inline void Prompt(NxsSuppressOutputStream & , const char * )
	{
	}
} // namespace ncl		

inline NxsSuppressOutputStream & operator<<(NxsSuppressOutputStream & outStream, NxsWritableStreamManipPtr )
	{
	return outStream;
	}

inline  NxsSuppressOutputStream & operator<<(NxsSuppressOutputStream & outStr, int) 
	{
	return outStr;
	}
	
inline  NxsSuppressOutputStream & operator<<(NxsSuppressOutputStream & outStr, long ) 
	{
	return outStr;
	}
	
inline  NxsSuppressOutputStream & operator<<(NxsSuppressOutputStream & outStr, unsigned ) 
	{
	return outStr;
	}
	
inline  NxsSuppressOutputStream & operator<<(NxsSuppressOutputStream & outStr, unsigned long ) 
	{
	return outStr;
	}
	
inline  NxsSuppressOutputStream & operator<<(NxsSuppressOutputStream & outStr, double ) 
	{
	return outStr;
	}
	
inline  NxsSuppressOutputStream & operator<<(NxsSuppressOutputStream & outStr, const char *) 
	{
	return outStr;
	}
	
inline  NxsSuppressOutputStream & operator<<(NxsSuppressOutputStream & outStr, char ) 
	{
	return outStr;
	}
	
inline  NxsSuppressOutputStream & operator<<(NxsSuppressOutputStream & outStr, const std::string &) 
	{
	return outStr;
	}

namespace ncl
{
inline  NxsSuppressOutputStream &flush(NxsSuppressOutputStream& om)
	{
	return om;
	}
	

inline  NxsSuppressOutputStream &endl(NxsSuppressOutputStream& om)
	{
	return om;
	}
} // namespace ncl

class NxsPlotStream: public NxsWritableStream
	{
	};

	
#endif
