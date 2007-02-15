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

#include "phycas/src/phycas_string.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	Appends a string version of the supplied unsigned value `v' to the string `s' and returns a reference to `s'.
*/
std::string &append_unsigned(std::string &s, unsigned v)
	{
	static char tmp[128];
	sprintf(tmp, "%d", v);
	s << tmp;
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Utility function useful for reporting strings in exceptions. If supplied string `s' is longer than `max_chars', it 
|	is shortened to `max_chars' characters and an ellipsis ("...") is appended to the returned value. If `s' is less
|	than `max_chars', it does not need to be abbreviated and is thus returned in its entirety.
*/
std::string abbreviate(const std::string &s, unsigned max_chars)
	{
	if (s.length() < max_chars)
		return s;
	else
		{
		std::string ss(s.substr(0, max_chars));
		ss << "...";
		return ss;
		}
	}

