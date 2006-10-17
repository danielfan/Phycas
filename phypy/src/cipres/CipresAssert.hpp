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

#if ! defined (CIPRES_ASSERT_HPP)
#define CIPRES_ASSERT_HPP

#ifdef __cplusplus
	/* C++ version */
	/* CIPRES_RELEASE_ASSERT is an assert that is not sensitive to NDEBUG */
#	include "phypy/src/cipres/compile_assert.hpp"
#	include <iostream>
#	define CIPRES_RELEASE_ASSERT(expression)  \
  	((void) ((expression) ? 0 : __TRIP_CIPRES_RELEASE_ASSERT (expression, __FILE__, __LINE__)))

#	define __TRIP_CIPRES_RELEASE_ASSERT(expression, file, lineno)  \
  	(std::cerr << file << ':' << lineno << ": failed assertion\n",	\
  	 abort (), 0)
  	 
#else //__cplusplus

	/* C version */
#	include <stdio.h>
#	define CIPRES_RELEASE_ASSERT(expression)  \
 	 ((void) ((expression) ? 0 : __TRIP_CIPRES_RELEASE_ASSERT (expression, __FILE__, __LINE__)))

#	define __TRIP_CIPRES_RELEASE_ASSERT(expression, file, lineno)  \
 	 (printf ("%s:%u: failed assertion\n", file, lineno),	\
 	  abort (), 0)

#endif //__cplusplus

#endif //#if ! defined (CIPRES_ASSERT_HPP)
