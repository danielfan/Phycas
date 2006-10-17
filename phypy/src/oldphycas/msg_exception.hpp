/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

#if !defined MSG_EXCEPTION_H
#define MSG_EXCEPTION_H

#include "phypy/src/ncl/nxs_defs.hpp"

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
