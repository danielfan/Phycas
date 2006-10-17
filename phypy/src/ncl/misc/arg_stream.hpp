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

#if !defined(ARG_STREAM_H)
#define ARG_STREAM_H
#include <string>
#include <vector>
class ArgStream
	{
	public:
		ArgStream(int argc, char *argv[]);
		bool ReadNextToken(std::string *s);
		
	private:
		typedef std::vector<std::string> VecToken;
		typedef VecToken::const_iterator VecToken_ConIt;
		
		VecToken 	Split(const std::string &s) const;

		VecToken	 	tokenStream;	
		VecToken_ConIt	tokenIt;	
	};
	
inline bool ArgStream::ReadNextToken(std::string *s)
	{
	if (tokenIt == tokenStream.end())
		return false;
	*s = *tokenIt++;
	return true;
	}
	
#endif
