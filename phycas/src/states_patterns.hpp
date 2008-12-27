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

#if ! defined(STATES_PATTERNS_HPP)
#define STATES_PATTERNS_HPP

#include <vector>
#include <map>

#if POLPY_NEWWAY
#	if HAVE_INTTYPES_H
#		include <inttypes.h>
#	elif HAVE_STDINT_H
#		include <stdint.h>
#	elif defined(_MSC_VER) && _MSC_VER >= 1200
#		include <basetsd.h>
		typedef   INT8 int8_t;
		typedef  UINT8 uint8_t;
		typedef  INT64 int64_t;
		typedef UINT64 uint64_t;
#	elif defined(_MSC_VER)
		typedef signed char int8_t;
		typedef unsigned char uint8_t;
		typedef long long int64_t;
		typedef unsigned long long uint64_t;
#	elif defined(_WIN32)
#		include <stdint.h>
#	endif
#else
#	include "ncl/nxscdiscretematrix.h"	    // for int8_t typedef
#endif

typedef std::vector<int8_t>				            VecStateList;
typedef std::vector<unsigned>			            VecStateListPos;

typedef double										PatternCountType;
typedef	std::map<VecStateList, PatternCountType>	PatternMapType;
typedef	std::vector<PatternCountType>				CountVectorType;

typedef std::vector< int8_t >                       StateMapping;

#endif


