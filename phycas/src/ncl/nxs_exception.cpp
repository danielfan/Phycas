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
#if (MWERKS_LIB_BUILD)
#	pragma export on
#endif
#include "phycas/src/ncl/nxs_defs.hpp"
#include "phycas/src/ncl/nxs_exception.hpp"
#include "phycas/src/ncl/nxs_token.hpp"
#if (MWERKS_LIB_BUILD)
#	pragma export off
#endif
using std::string;
/*----------------------------------------------------------------------------------------------------------------------
|	creates an NxsException object with the specified message, and file position information.
*/
NxsException::NxsException(
  const string	&s, /*message that describes the error */
  file_pos 		fp, /* file position where the error occurred*/
  unsigned 		fl, /* line where the error occurred*/
  unsigned 		fc, /* column where the error occurred*/
  CmdResult 	c) 	/* kCmdFailedSilent can be sent to hide throw a silent message */
	: msg(s), 
	pos(fp), 
	line(fl), 
	col(fc),
	cResultCode(c)
	{
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	creates an NxsException object with the specified message, getting file position information from the NxsToken
*/
NxsException::NxsException(
  const string 	&s, /* message that describes the error */
  const NxsToken 	&t, /* NxsToken that was supplied the last token (the token that caused the error) */
  CmdResult 	c) 		/* kCmdFailedSilent can be sent to hide throw a silent message */
  : msg(s), 
  pos(t.GetFilePosition()), 
  line(t.GetFileLine()), 
  col(t.GetFileColumn()),
  cResultCode(c)
	{
  	}
  
/*----------------------------------------------------------------------------------------------------------------------
|	Creates an NxsException exception with the message "Unexpected end of file encountered"
*/
NxsX_UnexpectedEOF::NxsX_UnexpectedEOF(
  NxsToken &t) /*NxsToken that was being read when the error occurred */
  	: NxsException( "Unexpected end of file encountered", t) 
  	{
  	}
  	
#if defined (NXS_THROW_IF_EMPTY_TOKENS)
	/*----------------------------------------------------------------------------------------------------------------------
	|	Creates an NxsException exception with the message "Illegal nexus word - pair of single quotes that are not inside a 
	|	single-quoted word"
	*/
	NxsX_EmptyToken::NxsX_EmptyToken(
	  const NxsToken &t) 
	  	: NxsException( "Illegal nexus word - pair of single quotes that are not inside a single-quoted word", t) 
	  	{
	  	}
#endif
/*----------------------------------------------------------------------------------------------------------------------
|	Calls NxsException(  const string &s, const NxsToken &t) constructor
*/
NxsX_UnknownCommand::NxsX_UnknownCommand(
  const string &s, 
  const NxsToken &t) 
  	: NxsException(s,t)
  	{
  	}
