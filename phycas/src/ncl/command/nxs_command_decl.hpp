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

#ifndef NCL_NXS_COMMAND_COMMON_H
#define NCL_NXS_COMMAND_COMMON_H

#include "phycas/src/ncl/nxs_defs.hpp"

class NxsCmdOption;
class NxsCommand;
class NxsTest;
typedef boost::shared_ptr<NxsTest>		NxsTestShPtr;//NOTE:this typedef is also in nxstest.h change it in both places
typedef boost::shared_ptr<NxsCommand>	NxsCommandShPtr;
typedef std::vector<NxsCommandShPtr> 	VecNxsCommandShPtr;	
typedef boost::shared_ptr<NxsCmdOption> NxsCmdOptionShPtr;
typedef std::vector<NxsCmdOptionShPtr> 	VecNxsCmdOptionShPtr;

typedef std::vector<const NxsCmdOption *> 	VecConstNxsCmdOptionPtr;


#endif
