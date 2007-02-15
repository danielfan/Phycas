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

#ifndef NCL_NXSTYPIST_H
#define NCL_NXSTYPIST_H

#include "ncl/nxs_defs.hpp"
/*--------------------------------------------------------------------------------------------------------------------------
|	Used by Nxs...OutputStream (or a similar class)to wrap lines and avoid extra newlines (If more than max_newl newlines 
|	occur in succession extra newlines are ignored).  The outstream's message is const.
*/
class NxsTypist
	{
	public:
					NxsTypist(const NxsOutputStream * = NULL);
		
		void		Run(const std::string &message, unsigned lmargin = 0, unsigned margin_incr = 2);
		
	
		unsigned	GetIndentation() const;
		unsigned	GetLeftMargin() const;
		unsigned	GetMaxNewl() const;
		unsigned	GetWidth() const;
		unsigned	GetSuccessiveNewlines() const;
	
		void		SetOutputStream(const NxsOutputStream *p);
		void		SetIndentation(unsigned);
		void		SetLeftMargin(unsigned);
		void		SetMaxNewl(unsigned how_many);
		void		SetWidth(unsigned how_many_characters);
		void		SetSuccessiveNewlines(unsigned s);
	
	private:
		void		BuildMarginStrings(unsigned lmargin, unsigned margin_incr);
		void		OutputLine(bool wrapping);
		void		EraseLine(bool wrapping);
		void		AppendCharToWord(char ch);
		void		AppendWordToLine(const std::string &);
		void		PrintWordToShorten();
	
	
		unsigned	width;		/* maximum number of characters in one line */
		unsigned	max_newl;	/* maximum number of successive newlines to allow */
		unsigned 	succ_newl;	/* the number of successive newlines that have been printed */
		std::string	line;		/* workspace for storing one line's worth of information */
		std::string	word;		/* workspace for storing next word to be added to line */
		std::string	leftMargin;		/* string containing lm spaces */
		std::string	indentedMargin;		/* string containing mi spaces */
	
		const NxsOutputStream *outStream;
		
	
	};
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `succ_newl', the number of successive newlines at the end of the the most recent output.
*/
inline unsigned	NxsTypist::GetSuccessiveNewlines() const
	{
	return succ_newl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of `succ_newl', the number of successive newlines at the end of the the most recent output.
*/
inline void NxsTypist::SetSuccessiveNewlines(unsigned s)
	{
	succ_newl = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use this method to obtain the current width
|	(number of characters to display before wrapping to
|	the next line).
*/
inline unsigned NxsTypist::GetWidth() const
	{
	return width;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use this method to obtain the maximum allowed number 
|	of successive newlines (i.e. value of data member max_newl).
*/
inline unsigned NxsTypist::GetMaxNewl() const
{
	return max_newl;
}

/*----------------------------------------------------------------------------------------------------------------------
|	Use this method to set the maximum allowed number 
|	of successive newlines (i.e. value of data member max_newl).
*/
inline void NxsTypist::SetMaxNewl(
  unsigned how_many) /* number of successive newlines */
	{
	max_newl = how_many;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use this method to set the number of characters to 
|	display before wrapping to the next line.
*/
inline void NxsTypist::SetWidth(
  unsigned how_many_characters) /* number of characters */
	{
	width = how_many_characters;
	}

/*----------------------------------------------------------------------------------------------------------------------
| Sets line to equal the empty string, and sets line_len to zero.
*/
inline void NxsTypist::EraseLine(bool wrapping)
	{
	line = (wrapping ? indentedMargin : leftMargin);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets lm and mi data members and builds the corresponding strings, lm_str and mi_str, that are used to insert margin
|	spaces at the beginning of each line.
*/
inline void NxsTypist::BuildMarginStrings(unsigned lmargin, unsigned margin_incr)
	{
	leftMargin.assign(lmargin, ' ');
	indentedMargin.assign(lmargin + margin_incr, ' ');
	}

/*----------------------------------------------------------------------------------------------------------------------
| Appends character ch to end of word string 
*/
inline void NxsTypist::AppendCharToWord(
  char ch)
	{
	word.push_back(ch);
	}

/*----------------------------------------------------------------------------------------------------------------------
| Adds contents of word string to end of line string. Updates line_len by 
| adding to it word_len, then sets word to equal the empty string and sets 
| word_len to zero. Also resets succ_newl to zero.
*/
inline void NxsTypist::AppendWordToLine(const std::string &toAppend)
	{
	line.append(toAppend);
	}

inline unsigned NxsTypist::GetIndentation() const
	{
	return (unsigned)indentedMargin.size() - GetLeftMargin();
	}

inline unsigned	NxsTypist::GetLeftMargin() const
	{
	return (unsigned)leftMargin.size();
	}
	
inline void NxsTypist::SetIndentation(unsigned i) 
	{
	BuildMarginStrings(GetLeftMargin(), i);
	}

inline void NxsTypist::SetLeftMargin(unsigned i) 
	{
	const unsigned ind = GetIndentation();
	BuildMarginStrings(i, ind);
	}
		
inline void NxsTypist::SetOutputStream(const NxsOutputStream *p)
	{
	NXS_ASSERT(p != NULL);
	outStream = p;
	}

#endif
