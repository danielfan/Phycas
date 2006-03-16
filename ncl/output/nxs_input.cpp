#include "phycas/force_include.h"
#include "ncl/output/nxs_input.hpp"
#include "ncl/output/nxs_output.hpp"
#include "ncl/nxs_exception.hpp"
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
