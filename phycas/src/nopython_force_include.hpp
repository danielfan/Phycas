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

#ifndef NOPYTHON_FORCE_INCLUDE_HPP
#define NOPYTHON_FORCE_INCLUDE_HPP

// This file (or one like it) should be specified as the force include file (a header file that is 
// included first in every source code file in the project).

#undef PYTHON_ONLY
#define POL_PHYCAS	// for POL-specific code related to PHYCAS (this one needs to go)
//#define POL_PYPHY	// for POL-specific code related to boost python (this one needs to go too)
#define NO_IDL_TYPES
#define POLPY_NEWWAY  1
#define POLPY_OLDWAY  !(POLPY_NEWWAY)

#include "phycas/src/phycas_config.h"

#endif
