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

#include "phycas/force_include.h"
#if defined(NCL_USER_SUPPLIED_OUTPUT)
#if defined(WIN_PHOREST)
#	include <process.h>
#endif


#include "ncl/output/nxs_output.hpp"
#include "gui/phycasGUI/network/PhycasServerSocket.h" //@ temp need to find a robust socket library
using ncl::endl;
void NxsXMLSocketOutputManager::SendLine(const std::string & r) const
	{
	if (IsOpen())
		spawnedSocket->SendLine(r);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates an output stream and makes aliases to it for the outpuComment, warning, error streams and user query.
*/
NxsXMLSocketOutputManager::NxsXMLSocketOutputManager()
	:portN(4444),
	serverSocket(),
	spawnedSocket(NULL),
	errorStream("error"),
	hiddenQueryStream("hidden_query"),
	outputStream("out", "<idle/>"),
	outputCommentStream("comment"),
	plotStream("plot"),
	warningStream("warning"),
	queryOut("user_query"),
	statusStream(NULL),
	inputStream(),
	userQuery(&queryOut, &inputStream)
	{
	statusStream = &outputStream;
	returnValStream = &outputStream;
	}

bool NxsXMLSocketOutputManager::SetPortNumber(UInt newPortN)
	{
	if (newPortN != portN)
		{
		if (serverSocket)
			serverSocket = ServerSocketShPtr(); //@ need some way to terminate the session cleanly
		portN = newPortN;
		}
	return true;
	}

void NxsXMLSocketOutputManager::AcceptConnection(AnswerSocketCallBack ascb)
	{
	if (!serverSocket)
		serverSocket = ServerSocketShPtr(new ServerSocket(portN, 5));
	std::cout << "Spawning socket acceptor listening on port " << portN << std::endl;
	Socket* s = serverSocket->Accept();
#   if defined(WIN_PHOREST)
		unsigned ret;
		_beginthreadex((void *)0, 0, ascb, (void*) s, 0, &ret);
		//_CRTIMP uintptr_t __cdecl _beginthreadex(void *, unsigned, unsigned (__stdcall *) (void *), void *, unsigned, unsigned *);
#   else
		(*ascb)((void *) s);
#   endif
	}

NxsXMLSocketOutputManager::~NxsXMLSocketOutputManager()
	{
	std::cout << "Closing socket" << std::endl;
	}

NxsOutputOperationStatusID NxsXMLSocketOutputManager::StartStatusDisplay(const std::string &m, const bool /*updatesWillBePosted*/)
		{
		*statusStream << m << endl;
		return NxsOutputOperationStatusID(1); //@@
		}
		
void NxsXMLSocketOutputManager::UpdateStatusDisplay(
  const float /*proportionDone*/, 
  const NxsOutputOperationStatusID , 
  const std::string /*optionalMsg*/)
	{
	}

void NxsXMLSocketOutputManager::EndStatusDisplay(
  const NxsOutputOperationStatusID , 
  const std::string msg)
	{
	*statusStream << msg << endl;
	}
	
void NxsXMLSocketOutputManager::AlertStatusDisplay(const std::string &msg)
	{
	*statusStream << msg << endl;
	}


#endif // defined(NCL_USER_SUPPLIED_OUTPUT)
