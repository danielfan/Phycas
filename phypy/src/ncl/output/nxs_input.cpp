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
#include "phypy/src/ncl/output/nxs_input.hpp"
#include "phypy/src/ncl/output/nxs_output.hpp"
#include "phypy/src/ncl/nxs_exception.hpp"
using std::string;
using NxsIO::NxsInput;
#if defined (NCL_SOCKET_IO)
#   include "ncl/output/nxs_output.hpp"
#	include "gui/phycasGUI/network/PhycasSocket.h" //@ should get this phycas stuff out of ncl
#	include "gui/phycasGUI/network/PhycasServerSocket.h" //@ should get this phycas stuff out of ncl
#	include "gui/phycasGUI/network/PhycasSocketException.h" //@ should get this phycas stuff out of ncl
	/*----------------------------------------------------------------------------------------------------------------------
	|	appends content of a single Socket::ReceiveLine() call to *nextLine this string onto the cmdHist.
	|   \bug we need an end of message signal for user-responses in case the socket breaks the response
	*/
	void NxsIO::NxsInput::AppendNextLine(
	  string *nextLine, 		 /* on exit, contains a copy of all of the characters before the \n or \r */  
	  std::deque<string> *cmdHist) /* pointer to a deque of the recent command.  nextLine is pushed onto the front. */
		{
		Socket * s = NxsOutputManager::GetInstance().GetSpawnedSocket();
		if (s == NULL)
			throw NxsFatalException("NxsOutputManager::socket must be initialized before calling NxsInput::AppendNextLine");
		std::string n = s->ReceiveLine();
		if (n.empty())
			throw NxsFatalException("Socket disconnected in  NxsInput::AppendNextLine");
		*nextLine << n;
		if (cmdHist != NULL)
			cmdHist->push_front(n);
			
		}
	
#else
	/*----------------------------------------------------------------------------------------------------------------------
	|	uses cin.get() to append all of the characters from cin to the *nextLine until \n or \r is encountered
	|	and pushes this command onto the cmdHist)
	*/
	void NxsInput::AppendNextLine(
	  string *nextLine, 		 /* on exit, contains a copy of all of the characters before the \n or \r */  
	  std::deque<string> *cmdHist) /* pointer to a deque of the recent command.  nextLine is pushed onto the front. */
		{
		string temp;
		temp.reserve(100);
		char c = (char) std::cin.get();
		while (c != '\n' && c != '\r')
			{
			temp << c;
			c = (char) std::cin.get();
			}
		if (!temp.empty())
			{
			*nextLine << temp;
			if (cmdHist != NULL)
				cmdHist->push_front(temp);
			}
		}
#endif
