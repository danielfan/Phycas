#ifndef NCL_NXS_AUTO_COMMAND_H
#define NCL_NXS_AUTO_COMMAND_H
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/command/nxs_command_decl.hpp"
#include "phypy/src/ncl/command/nxs_command.hpp"

/*--------------------------------------------------------------------------------------------------------------------------
|	NxsManualCommand and NxsAutoCommand both demand a callable object that takes void and returns a bool to serve as 
|	the Command executor.  If a command does everything it needs to do in the parsing stage (e.g. a command that just 
|	manipulates environmental variables), NxsDummyExecuteCommandCallback can be used to simply return true when the command
|	is executed.
*/
class NxsDummyExecuteCommandCallback 
	{
		
	public :
		CmdResult operator()(){return kCmdSucceeded;}
	};

/*--------------------------------------------------------------------------------------------------------------------------
|	simple class that holds an object pointer, a pointer to a member function that takes a SettingsStruct pointer and 
|	returns true if the command succeeded and false if the command failed (and the NxsCommand should revert the settings to
|	the value before the user's illegal command).
|	Used by NxsManualCommand and NxsAutoCommand as the controlled function (takes void returns bool)
|	the overloaded () operator calls object->Function(Settings *) and hides the somewhat obscure function pointer syntax 	
|	and makes it easier to keep the object and function pointer together.
|	This is really just another name for MemFuncAccessor<class ExecutionObject bool>.  If only we could typedef templates 
*/
template<class ExecutionObject> class NxsVoidExecuteCommandCallback 
	{
	public :
		NxsVoidExecuteCommandCallback(ExecutionObject * o, CmdResult (ExecutionObject:: * h)());
		CmdResult operator()();
	protected:
		ExecutionObject * const obj;
		CmdResult (ExecutionObject:: * handlerFunction)();
	};

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<class ExecutionObject> inline NxsVoidExecuteCommandCallback<ExecutionObject>::NxsVoidExecuteCommandCallback(
  ExecutionObject * o, 			/*pointer to the object whose member function should be called */
  CmdResult (ExecutionObject:: * h)()) /*function pointer to a bool (SettingsStruct *) member function of class T */
	:obj(o),
	handlerFunction(h)
	{
	NXS_ASSERT(o != NULL);
	NXS_ASSERT(h != 0L);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	calls the contained member function.
*/
template<class ExecutionObject> inline CmdResult NxsVoidExecuteCommandCallback<ExecutionObject>::operator()()	
	{
	return (obj->*(handlerFunction))();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	simple class that holds an object pointer, a pointer to a member function that takes a SettingsStruct pointer and 
|	returns true if the command succeeded and false if the command failed (and the NxsCommand should revert the settings to
|	the value before the user's illegal command).
|	Used by NxsManualCommand and NxsAutoCommand as the controlled function (takes void returns bool)
|	the overloaded () operator calls object->Function(Settings *) and hides the somewhat obscure function pointer syntax 	
|	and makes it easier to keep the object and function pointer together.
*/
template<class ExecutionObject, class SettingsStruct> class NxsExecuteCommandCallback 
	{
		ExecutionObject * const obj;
		boost::shared_ptr<SettingsStruct> settings;
		typedef CmdResult (ExecutionObject::* CallBackPtr)(SettingsStruct *) ;
		CallBackPtr mfuncPtr;
	public :
		NxsExecuteCommandCallback(ExecutionObject *o, CallBackPtr h, boost::shared_ptr<SettingsStruct> s);
		CmdResult operator()();
	};

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<class T, class SettingsStruct> inline NxsExecuteCommandCallback<T, SettingsStruct>::NxsExecuteCommandCallback(
  T *o, 			/*pointer to the object whose member function should be called */
 CallBackPtr h, /*function pointer to a bool (SettingsStruct *) member function of class T */
  boost::shared_ptr<SettingsStruct> s)
	: obj(o),
	settings(s),
	mfuncPtr(h)
	{
	NXS_ASSERT(o != NULL);
#	if defined (MAC_PHOREST)
		NXS_ASSERT(h != NULL);	// this PHYCAS_ASSERT is causing a warning on g++ on linux.  checking != NULL on one platform should be sufficient.
#	endif
	} 
	
/*--------------------------------------------------------------------------------------------------------------------------
|	calls the contained member function.
*/
template<class T, class SettingsStruct> inline CmdResult NxsExecuteCommandCallback<T, SettingsStruct>::operator()()	
	{
	PHYCAS_ASSERT (obj != NULL );
	SettingsStruct * temp = settings.get();
	return (obj->*(mfuncPtr))(temp);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	simple class that holds an object pointer and a void-NxsToken & member function pointer.
|	for use with ManualCommands so that ParseCommand will call the function supplied in the NxsExternalCommandParser constructor.
|	the function should throw an NxsException exception if the command is illegal.
|	the overloaded () operator calls object->Function() and hides the somewhat obscure function pointer syntax and makes it 
|	keep the object and function pointer together.
*/
template<class ParsingObject> class NxsExternalCommandParser : public std::unary_function<NxsToken, bool>
	{
		ParsingObject *obj;
		typedef bool (ParsingObject::* ParseCallBackPtr)(NxsToken &) ;
		ParseCallBackPtr parserFunction;
	public :
		NxsExternalCommandParser(ParsingObject *o, ParseCallBackPtr f);
		bool operator()(NxsToken &token);
	};

/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<class ParsingObject> inline NxsExternalCommandParser<ParsingObject>::NxsExternalCommandParser(
  ParsingObject *o, 						/*pointer to the object whose member function should be called to parse a token stream*/
  ParseCallBackPtr f)	/*function pointer to a bool - void NxsToken function of class T */
	: obj(o),
	parserFunction(f)
	{
	NXS_ASSERT(obj != NULL);
	NXS_ASSERT(f != 0L);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	calls the contained member function to read the token stream
*/
template<class ParsingObject> inline bool NxsExternalCommandParser<ParsingObject>::operator()(NxsToken &token)	
	{
	return (obj->*parserFunction)(token);
	}



/*--------------------------------------------------------------------------------------------------------------------------
|	a class overrides the two pure virtual functions of NxsCommand.
|	ManualCommands should be used if you have a parsing function (that takes a reference to a NxsToken, returns void, and
|	throws NxsException exceptions in the event of errors) and an command executor (that takes void and returns bool if execution
|	is successful (NOTE returning false means the command state will be backed up to the previous state which is the same
|	behavior as when the command can't be parsed.
|
*/
template<class Parser, class Executor> class NxsManualCommand: public NxsCommand
	{
		Parser 		cmdParser;
		Executor	cmdExecutor;
		bool 				ParseCommand(NxsToken&);					
		CmdResult		    ExecuteControlledFunction();
	public :
		NxsManualCommand(const std::string &name, Parser , Executor);
	};



/*--------------------------------------------------------------------------------------------------------------------------
|	Abstract base class that both NxsAutoCommand.
|	By overriding the pure virtual ParseCommand(), this class provides automatic parsing of a NEXUS command using the CmdItem 
|	class to interpret keywords.  
|	Errors result in Parsing NxsException exceptions being thrown.
|	Still abstract because no ExecuteControlledFunction() is provided
*/
class NxsAutomaticallyParsedCommand: public NxsCommand
	{
	public :
		NxsAutomaticallyParsedCommand(const std::string &nameOfCommand);
		virtual ~NxsAutomaticallyParsedCommand(){}
		
		void 					SetReadsAsterisk(bool a);
	private:	
	
		void					ParseFailed(NxsToken&,CmdResult explain);	//throws NxsException
		bool					AllOptionsAreValid();
		void					CheckAbbreviationIsTooShort(NxsToken &);
		bool 					CheckUnnamedCmdOptIfUnread(NxsCmdOptionShPtr &optIt);
		bool 					CheckNamedCmdOptIfUnread(NxsCmdOptionShPtr &optIt);
		//void					GetErrorFromCmdOption(NxsCmdOptionShPtr, NxsToken &);
		bool					ParseCommand(NxsToken &);
		NxsCmdOptionShPtr	ParseKeyword(NxsToken &);
		void 					PrintWarnings() const;
		void 					PrepareToRead();
		void					SkipHelpRequests(NxsToken& );
	
		bool 					canReadAsterisk;
		VecString					warnings; /* command items that had non-fatal errors */
	};

/*--------------------------------------------------------------------------------------------------------------------------
|	Command class which relies on NxsAutomaticallyParsedCommand for the parsing of NxsToken, and allows the programmer
|	to supply an object and member function pointer that should be called when the command is correctly parsed.
*/
template<class CommandExecutor> class NxsAutoCommand: public NxsAutomaticallyParsedCommand
	{
		CommandExecutor cmdExecutor;
		CmdResult		ExecuteControlledFunction();
							
	public :
		NxsAutoCommand(const std::string &nameOfCommand, CommandExecutor e);
	};
	
	
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<class CommandExecutor> inline NxsAutoCommand<CommandExecutor>::NxsAutoCommand(
  const std::string &n, 	/*name of the command */
  CommandExecutor e) /* object/function pointer/argument called when the command is successfully read*/
	: NxsAutomaticallyParsedCommand(n), 
	cmdExecutor(e) 
	{
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Calls the function supplied in the constructor as the command executor
*/
template<class CommandExecutor> CmdResult NxsAutoCommand<CommandExecutor>::ExecuteControlledFunction()
	{
	CmdResult ret =  cmdExecutor();
	NotifyCmdOptionsOfExecution();
	return ret;
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Calls the function supplied in the constructor as the command executor
*/
template<class Parser, class Executor> CmdResult NxsManualCommand<Parser, Executor>::ExecuteControlledFunction()
	{
	return cmdExecutor();
	}
	
/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
template<class Parser, class Executor> inline NxsManualCommand<Parser, Executor>::NxsManualCommand(
  const std::string &n, /*name of the command */
  Parser 	parseCallback,  /*pointer to the object whose member functions provide parsing and command executions */
  Executor	executeCallback)/* command execution member function */
	: NxsCommand(n),
	cmdParser(parseCallback),
	cmdExecutor(executeCallback)
	{
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Calls the parsing function supplied in the constructor as the command executor
*/
template<class Parser, class Executor> bool NxsManualCommand<Parser, Executor>::ParseCommand(
  NxsToken &token)
	{
	return cmdParser(token);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Clears warnings from previous commands
*/
inline void NxsAutomaticallyParsedCommand::PrepareToRead()
	{
	NxsCommand::PrepareToRead();
	warnings.clear();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	tells the command whether or not it should accept * after the command name
*/
inline void NxsAutomaticallyParsedCommand::SetReadsAsterisk(
  bool a)	/* true if the command can optionally take an * as the first token */
	{
	canReadAsterisk = a;
	}



#endif
