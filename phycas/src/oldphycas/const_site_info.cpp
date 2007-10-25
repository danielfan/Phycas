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

//#include "phycas/force_include.h"
#include "phycas/src/oldphycas/const_site_info.hpp"
#include "phycas/src/ncl/misc/nxs_data_type.hpp"
using std::vector;

OrdCodedArrs::OrdCodedArrs(
  const vector<DataStorageType *> &bitCode,  
  unsigned wordLen,
  unsigned numStates)
  //POL-121803 : Compressed2DArray<NStateInt>(bitCode.size()),
  //POL-121803 nStates(numStates)
  : Compressed2DArray<NStateInt>((unsigned)bitCode.size(), numStates)
  	{
	typedef vector<DataStorageType *>::const_iterator BitCodePtrIter;
  	for (BitCodePtrIter bIt = bitCode.begin(); bIt != bitCode.end(); ++bIt)
  		{
  		vector<unsigned> ordCode = ConvertBitArrToVecUInt(*bIt, wordLen);
  		Append(ordCode);
  		}
  	}

