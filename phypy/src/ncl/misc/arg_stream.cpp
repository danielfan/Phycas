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

//#include "phycas/force_include.h"
#include "phypy/src/ncl/misc/arg_stream.hpp"
#include "phypy/src/ncl/misc/string_extensions.hpp"
#include <ctype.h>
using std::string;
ArgStream::VecToken ArgStream::Split(const string &s) const
	{
	VecToken retVec;
	string temp;
	if (!s.empty())
		{
		if (s[0] == '-' && s.length() > 2 && isgraph(s[1]) && !isdigit(s[1]) && s[1] != '-')
			{ 
			// 
			//	break "-nArg" into {"-n", "Arg"}
			//	Note that "-n Arg" indicates that the original arg was quoted, so we don't split
			retVec.push_back(string(s.begin(), s.begin() + 2));
			const string rest = string(s.begin() + 2, s.end());
			const VecToken tempVec = Split(rest);
			copy(tempVec.begin(), tempVec.end(), back_inserter(retVec)); 
			}
		else
			retVec = SplitString(s, ',');
		}
	return retVec;
	}
	
	
ArgStream::ArgStream(int argc, char *argv[])
	{
	tokenStream.reserve((unsigned) argc);
	for (int i = 1; i < argc; ++i)
		{
		string n = string(argv[i]);
		VecToken t = Split(n);
		copy(t.begin(), t.end(), back_inserter(tokenStream));
		}
	tokenIt = tokenStream.begin();
	}

