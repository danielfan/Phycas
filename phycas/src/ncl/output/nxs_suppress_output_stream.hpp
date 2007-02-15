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

#ifndef NCL_SUPPRESS_OUTSTREAM_H
#define NCL_SUPPRESS_OUTSTREAM_H

#include "phycas/src/ncl/nxs_defs.hpp"
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
