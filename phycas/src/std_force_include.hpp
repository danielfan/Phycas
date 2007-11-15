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

#ifndef STD_FORCE_INCLUDE_HPP
#define STD_FORCE_INCLUDE_HPP

// This file (or one like it) should be specified as the force include file (a header file that is 
// included first in every source code file in the project).

#define PYTHON_ONLY	// for code that only makes sense when exported to Python (most code can be used to construct a pure C++ program)
#define POL_PHYCAS	// for POL-specific code related to PHYCAS (this one needs to go)
//#define POL_PYPHY	// for POL-specific code related to boost python (this one needs to go too)
#define NO_IDL_TYPES
#define POLPY_NEWWAY  1
#define POLPY_OLDWAY  !(POLPY_NEWWAY)
//#define USING_NUMARRAY
//#define INTERFACE_WITH_CIPRES

// Uncommenting LOG_LOT_UNIFORM_CALLS macro results in a file named uniform.txt being saved that contains
// a record of every random number generated in Lot::Uniform() showing the file and line on which the call
// was made (extremely useful for tracking down where the sequence of random numbers begins to diverge)
//#define LOG_LOT_UNIFORM_CALLS
#if defined(LOG_LOT_UNIFORM_CALLS)
#   define FILE_AND_LINE __FILE__,__LINE__
#else
#   define FILE_AND_LINE
#endif

#if defined(_MSC_VER)
// The macros below allows you to insert #pragma TODO("fix this") in C++ code and have it show up in the Visual Studio IDE
// see http://mail.python.org/pipermail/python-list/2002-January/121346.html
#   define __STRINGIZE__(x) #x
#   define TODO(x) message(__FILE__"("__STRINGIZE__(__LINE__)") : " "TODO: "#x) 
#endif

#include "phycas/src/phycas_config.h"

#endif
