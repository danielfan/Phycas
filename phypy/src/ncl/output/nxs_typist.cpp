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
#include "phypy/src/ncl/nxs_defs.hpp"
#if defined(NCL_USE_NXS_CONSOLE_OUTPUT)
#	include "phypy/src/ncl/output/nxs_output_stream.hpp"
#	include phypy/src/"ncl/output/nxs_typist.hpp"
#	if defined (C_FUNCS_IN_STD_NAMESPACE)
		using std::isgraph;
		using std::isdigit;
#	endif	
	using std::string;
	//@POL currently need to set output width to 79 to avoid blank lines on console of width 80
	// something's wrong, perhaps with typist

	/*----------------------------------------------------------------------------------------------------------------------
	| Passes the line data member to the Nxs...OutputStream object's PrintLine method for 
	| output.  Then initializes line to be a null string and sets line_len to zero.
	*/
	inline void NxsTypist::OutputLine(bool wrapping)  /* line from which function was invoked (only used in debugging) */
		{
		outStream->PrintLine(line);
		++succ_newl;
		EraseLine(wrapping);
		}

	/*----------------------------------------------------------------------------------------------------------------------
	|	Assigns kernel pointer, initializes line_len and word_len to 0 and reading_ws
	|	to false.
	*/
	NxsTypist::NxsTypist(const NxsOutputStream *o)
	  :width(50),
	  max_newl(2),
	  succ_newl(0),
	  outStream(o)
		{
		}

	/*----------------------------------------------------------------------------------------------------------------------
	|	Advances iterator through the `message'  creating lines of output as it goes, which it sends to the 
	|	the member outStream->PrintLine.
	|	Always leaves lmargin spaces 
	|	after a newline, and lmargin + margin_incr for any lines that are wrapped. 
	*/
	void NxsTypist::Run(
	  const string &message, 
	  unsigned lmargin /* = 0 */, 
	  unsigned margin_incr /* = 2 */)
		{
		BuildMarginStrings(lmargin, margin_incr);
		EraseLine(false);
		word.clear();

		const unsigned message_len = (unsigned)message.length();
		char ch = message.at(0);
		if (message.length() == 1)
			{
			if (ch != '\n')
				AppendWordToLine(message);
			OutputLine(false);
			return;
			}
		for (unsigned i = 1; i < message_len;)
			{
			if (ch == '\n') 
				{
				if (succ_newl <= max_newl) 
					OutputLine(false);
				ch = message.at(i++);
				}
			else 
				{
				if (isgraph(ch))
					{
					do 	{
						AppendCharToWord(ch);
						if (i >= message_len)
							{
							ch = ' ';
							break;
							}
						ch = message.at(i++);
						}
					while (isgraph(ch));
					if (line.length() + word.length() > width) 
						{
						if (word.length() <= width) 
							OutputLine(true);
						else
							PrintWordToShorten();
						}
					AppendWordToLine(word);
					succ_newl = 0;
					}
				else
					{
					do 	{
						AppendCharToWord(ch);
						if (i >= message_len)
							break;
						ch = message.at(i++);
						}
					while (!isgraph(ch) && ch != '\n');
					if (line.length() + word.length() > width) 
						OutputLine(true);//don't wrap extra whitespace around
					else	
						AppendWordToLine(word);
					}
				word.clear();
				}
			}
		if (ch =='\n')
			{
			if (succ_newl < max_newl) 
				OutputLine(false);
			}
		else 
			{
			if (isgraph(ch))
				{
				word = ch;
				if (line.length() + word.length() > width) 
					OutputLine(false);
				AppendWordToLine(word);
				}
			outStream->PrintLine(line);
			++succ_newl;
			}
			
		}

	/*--------------------------------------------------------------------------------------------------------------------------
	|	Used in the (presumably rare) case in which a single word is longer than the line.  This function repeatedly calls
	|	OutputLine printing out "width" characters each time and returns when < "width" characters are left in the "word" string
	|	(these characters in  word is NOT appended to the line)
	*/
	void NxsTypist::PrintWordToShorten()
		{
		string temp;
		temp.reserve(width);
		const unsigned wordLength =  (unsigned)word.length();
		unsigned i = 0;
		const unsigned firstlineWidth = width-line.length();
		if (firstlineWidth > 5)
			{
			PHYCAS_ASSERT(i + firstlineWidth < wordLength);
			temp.assign(word, i, firstlineWidth);
			AppendWordToLine(temp);
			i += firstlineWidth;
			}
		OutputLine(true);
		const unsigned otherlineWidth = width - GetIndentation();
		
		for (;; i += otherlineWidth)
			{
			if (i + otherlineWidth >= wordLength)
				{
				word.erase(0,i);
				return;
				} 
			temp.assign(word, i, otherlineWidth);
			succ_newl = 0;
			AppendWordToLine(temp);
			OutputLine(true);
			}
		}
#endif //! defined(NCL_USE_STD_OUTPUT)
