#include "phycas/force_include.h"
#include <cstdlib>
#if defined (CONSOLE_PHOREST)
#	if defined(MONITORING_ALLOCATION)
#		define INSTANTIATE_MEMCHK
#	endif
#endif 


#include "phycas/phycas.h"
#include <sstream>
#if defined (NCL_SOCKET_IO)
#	include "gui/phycasGUI/network/PhycasSocket.h"
#	include "gui/phycasGUI/network/PhycasServerSocket.h"
#	include "gui/phycasGUI/network/PhycasSocketException.h"
#	include <string>
#endif
#include "phycas/floor.hpp"
#include "ncl/nxs_exception.hpp"
#include "ncl/misc/utilities.hpp"
#include "ncl/misc/string_extensions.hpp"
#include "ncl/nxs_token.hpp"
#include "ncl/output/nxs_input.hpp"
#include "ncl/output/nxs_input_error_wrapper.hpp"
#include "ncl/output/temp_basic_output_operators.hpp"
#include "phycas/misc/memory_check.hpp"
#include "phycas/misc/msg_exception.hpp"
using ncl::endl;
using std::string;
using std::cerr;
using std::exit;
using NxsIO::NxsInput;
#if defined (CORBA_PHYCAS)
#   include "phycas/corba/corba_wrapper.hpp"
#endif
#if defined (CORBA_SERVER_PHYCAS)
#   include "phycas/corba/tree_merge_wrapper.hpp"
#	include "CipresCommlib/CipresFacilitator.h"
#endif


#if defined (CONSOLE_PHOREST) || defined(NCL_SOCKET_IO)
 
#   if defined(PY_PHYCAS)
	#include <unistd.h>			//getcwd()
	std::string GetCWD();
	bool AddDirToPythonPath(const StrVec * = NULL);
	
	std::string GetCWD()
		{ //@ temp need to make this platform neutral
		char cwd[1000];
		::getcwd(cwd, 1000);
		return string(cwd);
		}
	
	bool AddDirToPythonPath(const StrVec * relPathsToAdd)
		{
		const std::string cwd = GetCWD();
		std::string pyAddPathsCommand;
		if (relPathsToAdd == NULL || relPathsToAdd->empty())
			StrPrintF(pyAddPathsCommand, "import sys\nsys.path.append(\'%s\')\n", cwd.c_str());
		else
			{
			pyAddPathsCommand = "import sys, os\n";
			for (StrVec::const_iterator pIt = relPathsToAdd->begin(); pIt != relPathsToAdd->end(); ++pIt)
				StrPrintF(pyAddPathsCommand, "sys.path.append(os.path.join(\'%s\', \'%s\'))\n", cwd.c_str(), pIt->c_str());
			}
		PyRun_SimpleString(pyAddPathsCommand.c_str());
		if (PyErr_Occurred())
			{
			PyErr_Print();
			return false;
			}
		return true;
		}
		

	int main(int argc, char *argv[])
		{
		if (argc < 3) 
			{
			cerr << "Usage: call pythonfile funcname [args]\n";
			return 1;
			}
		Py_Initialize();
		StrVec h(1,"testhide");
		if (!AddDirToPythonPath(&h))
			return 2;
		PyObject * pName = PyString_FromString(argv[1]);
			//@ Error checking of pName left out 
		PyObject * pModule = PyImport_Import(pName);
		Py_DECREF(pName);
		if (pModule != NULL) 
			{	// pDict and pFunc are borrowed references  (must not be Py_DECREF-ed)
			PyObject * pDict = PyModule_GetDict(pModule);
			PyObject * pFunc = PyDict_GetItemString(pDict, argv[2]);
			if (pFunc && PyCallable_Check(pFunc)) 
				{
				PyObject * pArgs = PyTuple_New(argc - 3);
				for (int i = 0; i < argc - 3; ++i) 
					{
					PyObject * pValue = PyInt_FromLong(atoi(argv[i + 3]));
					if (!pValue)
						{
						Py_DECREF(pArgs);
						Py_DECREF(pModule);
						fprintf(stderr, "Cannot convert argument\n");
						return 1;
						}
					/* pValue reference stolen here: */
					PyTuple_SetItem(pArgs, i, pValue);
					}
				PyObject * pValue = PyObject_CallObject(pFunc, pArgs);
				Py_DECREF(pArgs);
				if (pValue != NULL) 
					{
					printf("Result of call: %ld\n", PyInt_AsLong(pValue));
					Py_DECREF(pValue);
					}
				else
					{
					Py_DECREF(pModule);
					PyErr_Print();
					fprintf(stderr,"Call failed\n");
					return 1;
					}
				}
			else
				{
				if (PyErr_Occurred())
					PyErr_Print();
				fprintf(stderr, "Cannot find function \"%s\"\n", argv[2]);
				}
			Py_DECREF(pModule);
			}
		else 
			{
			PyErr_Print();
			fprintf(stderr, "Failed to load \"%s\"\n", argv[1]);
			return 1;
			}
		Py_Finalize();
		return 0;
		}

#   else //defined(PY_PHYCAS)

		int main(int argc, char *argv[])
			{
			CREATE_MEMCHK
				{
#				if ! defined (CORBA_SERVER_PHYCAS)
					char* infile_name = NULL;
#				endif
#				if defined (CORBA_PHYCAS)
					FACILITATOR->initialize(argc, argv);
					//trigger corba initialization so we crash fast, if we are going to crash
					//cipresCORBA::CorbaWrapper::GetInstance(argc, argv); 
#					if defined(NCL_SOCKET_IO)
						cerr << "calling NxOutputManager" << endll
						NxsOutputManager::GetInstance().SetPortNumber(4444);
						cerr << "back from NxOutputManager" << endll
#					endif
#				else // defined (CORBA_PHYCAS)
					if (argc > 2)
						{
						cerr << "Sorry, this program can accept at most one command line argument, which must be ";
#						if !defined(NCL_SOCKET_IO)
							cerr << "the name of a NEXUS file.";
#						else
							cerr << "the port number to server";
#						endif
						exit(1);
						}
					else if (argc > 1)
						{
#						if !defined(NCL_SOCKET_IO)
							infile_name = argv[1];
#						else
							string s(argv[1]);
							UInt portN = 4444;
							if (!IsAnUnsigned(s, &portN))
								{
								cerr << "Sorry, this program can accept at most one command line argument, which must be  the port number to serve.";
								exit(2);
								}
							NxsOutputManager::GetInstance().SetPortNumber(portN);
#						endif
						}
#				endif // defined (CORBA_PHYCAS)
				try
					{
#					if defined (CORBA_SERVER_PHYCAS)
						// TL - took PhoFloor out for now because it's ctor uses CORBAWrapper.
						cerr << "Calling RunServer" << endl;
						cipresCORBA::CorbaTreeMergeWrapper::RunServer();
						cerr << "Returned from RunServer" << endl;
						//cipresCORBA::CorbaWrapper * corbaWrapper = cipresCORBA::CorbaWrapper::GetInstance(argc, argv);
						//corbaWrapper->RunServer(floor);
#					else
						PhoFloor floor;
						floor.Run(infile_name);
#					endif
					}
				catch (...)
					{
					cerr << "Terminating as the result of an uncaught exception" << std::endl;
					exit(2);
					}
				}
#			if defined(MONITORING_ALLOCATION) && !defined(NDEBUG)
				ofstream memf("memcheck.txt");
				MEMCHK_REPORT(memf)
				memf.close();
#			endif
			return 0;
			}
#   endif //defined(PY_PHYCAS)
#endif //if defined (CONSOLE_PHOREST) || defined(NCL_SOCKET_IO)

VecString PhoFloor::GetCommandNames()
	{
	VecString cmdNames;
	const unsigned nCmds = NxsCommandManager::GetNumCommands();
	for (unsigned i = 0; i< nCmds; ++i)
		cmdNames.push_back(NxsCommandManager::GetCommandName(i));
	return cmdNames;
	}

void PhoFloor::NexusError(
  const string &msg, 
  file_pos pos, 
  unsigned line, 
  unsigned col, 
  CmdResult resCode,
  NxsBlock* currBlock)
  	{
#	if defined (HAVE_PRAGMA_UNUSED)
#		pragma unused(pos)
#	endif
	if (resCode == kCmdFailedGenerateMessage)
		{
		NxsErrorStream * errorStream = outputMgrRef.GetErrorStreamPtr();
		if (errorStream != NULL)
			{
			NxsInputErrorWrapper ies(msg, (currBlock != NULL ? currBlock->GetID() : ""), fileStack, pos, line, col);
			Emit(*errorStream, ies) << endl;
			}
		}
	if (currBlock != NULL && !fileStack.empty())
		currBlock->Reset(); 
  	}


	

#if defined (NCL_SOCKET_IO)
	
	unsigned PHYCAS__STDCALL Answer(void* a);
	PhoFloor * gSimpleApp = NULL;

	inline void PhoFloor::Prompt() 
		{
		outputMgrRef.Prompt();
		}

	void PhoFloor::Run(
	  char *)	/* the name of the NEXUS data file to execute (can be NULL) */
		{
		try
			{
			std::cout << "In Run()" << std::endl;
			gSimpleApp = this;
			while (gSimpleApp != NULL && !gSimpleApp->GetIsExiting()) 
				outputMgrRef.AcceptConnection(Answer);
			}
		catch (SocketException & e)
			{
			std::cout << "SocketException occurred:" << e.description() << "\nExiting.\n";
			throw;
			}
		}	

	unsigned PHYCAS__STDCALL Answer(void * a) 
		{
		try
			{
			NxsOutputManager & outMgrRef = NxsOutputManager::GetInstance();
			std::cout << "In Answer" << std::endl;
			Socket * s = (Socket *) a;
			outMgrRef.SetSpawnedSocket(s); //hackety-hack
			while (gSimpleApp != NULL && !gSimpleApp->GetIsExiting()) 
				{
				std::string next_command;
				next_command << ';' << s->ReceiveLine();
				if (next_command == ';') 
					break;
				string response;
				try
					{
					next_command << ";end;";
					std::cout << "About to process \"" << next_command << "\"" << std::endl;
					NxsToken token(next_command.c_str());
					gSimpleApp->ReadTokenStream(token);
					gSimpleApp->Prompt();
					}
				catch (...)
					{
					response << "Phycas is terminating because of an uncaught exception.";
					s->SendLine(response);
					break;
					}
				}
			outMgrRef.SetSpawnedSocket(NULL); //hackety-hack
			delete s;
			}
		catch (SocketException&) {}
		return 0;
		}
#else

inline void	PhoFloor::Prompt()
	{
	if (outputMgrRef.GetErrorStreamPtr())
		ncl::Prompt(*outputMgrRef.GetErrorStreamPtr(), "Phycas> ");
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Runs the command line interpreter, allowing PhoFloor to interact with user. Typically, this is the only 
|	function called in main after a PhoFloor object is created. If `infile_name' is non-NULL, the first command 
|	executed by the command interpreter will be "EXECUTE `infile_name'".
*/
void PhoFloor::Run(
  char *infile_name)	/* the name of the NEXUS data file to execute (can be NULL) */
	{
	
	quitNow = false;
	if (infile_name != NULL)
		{
		next_command.clear();
		next_command << ";execute file = " << infile_name;
		PreprocessNextCommand(&next_command);
		HandleNextCommand(next_command.c_str());
		}
		
#	if defined(DEBUGGING_FROM_RUN_NEX)
		HandleNextCommand(";execute file = '_run.nex';end;");
#	endif

	try
		{
		while (!quitNow) 
			{
#			if defined (NCL_USE_NXS_CONSOLE_OUTPUT)
				if (outputMgrRef.GetOutputStreamPtr()->GetSuccessiveNewlines() < 2)
					cerr << '\n';
#			endif
			Prompt();
#			if defined (NCL_USE_NXS_CONSOLE_OUTPUT)
				outputMgrRef.GetOutputStreamPtr()->SetSuccessiveNewlines(1);
#			endif
			next_command = ';';
			if (outputMgrRef.GetInputStreamPtr() == NULL)
				break;
			outputMgrRef.GetInputStreamPtr()->AppendNextLine(&next_command, &commandHistory);
			PreprocessNextCommand(&next_command);
			HandleNextCommand(next_command.c_str());
			}
		}
	catch (MsgException & x)
		{
		cerr << "Phorest is terminating because of an uncaught exception.\nPlease report the following message to the authors:\n";
		cerr << x.msg.c_str();
		cerr << "\n";
		}
	}
#endif
/**
 * 
 *
 * Begins with the command just entered by the user, which is stored in
 * the data member next_command, adds a semicolon (if the user failed
 * to supply one), and then adds "end;" so the whole bundle looks
 * like a very short PhoFloor block.  This is then passed to HandleNextCommand,
 * which processes it just like a real PhoFloor block in a NEXUS data file.
 */
void PhoFloor::PreprocessNextCommand(
  string *nextComm)
	{
	// If user failed to add the terminating semicolon,
	// we'll do it now. We will also remove the line feed
	// at the end and add the command "end;" to the end
	// of the line (see explanation below).
	//
	unsigned len = (unsigned) nextComm->length();
	if (len == 0)
		{
		*nextComm = ";end;";
		return;
		}
	// If character at position i-1 is a semicolon, put '\0' terminator at position i;
	// otherwise, put a semicolon at position i and terminator at i+1
	//
	if (nextComm->at(len-1) != ';') 
		*nextComm <<  ';';
	
	*nextComm << "end;";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by PhoFloor::HandleNextCommand to read a stream of NxsTokens.  Calls NxsCommandManagerBlock::Read
|	and reports errors through PhoFloor::NexusError
|	returns true if no errors occurred
*/
bool PhoFloor::ReadTokenStream(NxsToken &token)
	{
	try 	
		{
		Read(token); 
		}
	catch (NxsException & x)
		{
		NexusError(x.msg, x.pos, x.line, x.col, x.cResultCode);
		return false;
		}
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accepts a string in the form of ``;commands;end;'', packages it into a NxsToken and calls PhoFloor::ReadTokenStream
 */
bool PhoFloor::HandleNextCommand(
  const char *nextComm)
	{
	// Hold the door so messages for other applications can be processed
	//
	PhorestYield();
	NxsToken token(nextComm);
	return ReadTokenStream(token);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	NxsReader overload that indicates to the user that the beginning of a ``blockName'' block was found
*/
bool  PhoFloor::EnteringBlock(const string &blockName)
	{
	string s;
	s << "Reading " << blockName << " block";
	blockOpIDStack.push(outputMgrRef.StartStatusDisplay(s, false));
	return true;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	NxsReader overload that indicates to the user that the end of a ``blockName'' block was found
*/
void  PhoFloor::ExitingBlock(const string & blockName)
	{
	string s;
	s << "Finished with " << blockName << " block";
	const NxsOutputOperationStatusID blockOpID = (blockOpIDStack.empty() ? NxsOutputOperationStatusID(0) : blockOpIDStack.top());
	outputMgrRef.EndStatusDisplay(blockOpID, s);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	NxsBlock overload that indicates that the end of a phorest block was found
*/
CmdResult PhoFloor::EndEncountered()
	{
	return kCmdSucceeded;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overload of NxsReader function that tells the user that a disabled block is being skipped	
*/
void  PhoFloor::SkippingDisabledBlock(const string &blockName)
	{
	string s;
	s << "Skipping a " << blockName << " block (" << blockName << " blocks are currently disabled)";
	outputMgrRef.AlertStatusDisplay(s);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Overload of NxsReader function that tells the user that an unknown block is being skipped	
*/
void  PhoFloor::SkippingBlock(const string &blockName)
	{
	string s;
	s << "Skipping an unknown block (" << blockName << ')';
	outputMgrRef.AlertStatusDisplay(s);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void PhoFloor::TreesChanged(NxsTreesManager *, NxsTreeListener::TreeChangeType)
	{
	//@ dummy function
	//*outputMgrRef.GetOutputStreamPtr() << "The floor knows that trees have changed" << endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by the taxa manager when the taxon status changes.
*/
void PhoFloor::TaxaChanged(BaseTaxaManager *, NxsTaxaListener::TaxaChangeType)
	{
	//@ dummy function
	//*outputMgrRef.GetOutputStreamPtr() << "The floor knows that taxa have changed" << endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by the character manager when the char status changes.
*/
void PhoFloor::CharsChanged(NxsCharactersManager *, NxsCharListener::CharChangeType)
	{
	//@ dummy function
	//*outputMgrRef.GetOutputStreamPtr() << "The floor knows that characters have changed" << endl;
	}
