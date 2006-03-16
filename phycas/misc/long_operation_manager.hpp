#if !defined(PHO_LONGOPERATIONMGR_H)
#define PHO_LONGOPERATIONMGR_H

#if defined(WIN_PHOREST) && defined(CONSOLE_PHOREST) // Win32 console version

//POL 18Feb205 - the following sentinel was defined to prevent windows.h from 
// including winsock.h When sockets are used in phycas (e.g. via ACE), winsock2.h
// is used rather than winsock.h, and this causes redefinition problems
//
#	define _WINSOCKAPI_	

#	include <windows.h>
#	include <wincon.h>
#	include <crtdbg.h>
	BOOL WINAPI PhoInterruptHandler(DWORD);

#elif defined(CONSOLE_PHOREST)	// Linux console version

#	include <signal.h>
	void PhoInterruptHandler(int);

#endif


/*----------------------------------------------------------------------------------------------------------------------
|	This class provides an object to manage long operations. It provides virtual functions that can be overridden by 
|	graphical interfaces to manage progress dialogs. For portable versions it handles Ctrl-C events and queries the 
|	user about what action to take upon cancelling a long operation. Only one of these objects should be instantiated
|	by the program.
*/
class LongOperationManager
	{
	public:
		enum LongOpType					/* specifies the different types of long operations possible; each can be managed separately */
			{
			LongOp_None		= 0x0000,	/* no long operations are currently in progress */
			LongOp_File		= 0x0001,	/* a file is being read and interpreted */
			LongOp_MCMC		= 0x0002,	/* an MCMC run is in progress */
			LongOp_LScore	= 0x0004,	/* finding maximum likelihood estimates of branch lengths and model parameters */
			LongOp_Search	= 0x0008,	/* a search for the optimum tree is in progress */
			LongOp_Misc		= 0x0010	/* an undefined long operation is in progress */
			};

		bool				IsRunning(LongOpType which);
		virtual void		StartLongOperation(LongOpType which);
		virtual void		StopLongOperation(LongOpType which);
		bool				CheckAbort(LongOpType which);
		void				SetInterrupted(bool wasInterrupted = true)
			{
			ctrl_c = wasInterrupted;
			}

		LongOperationManager();
		
	private:
		bool				ctrlc_handler_installed;
		unsigned			longOps;		/* on bits specify which long operations are in progress */
		bool				ctrl_c;			/* for catching interrupts - had been STATIC_DATA, but now LongOperationManager is a Singleton */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the `longOps' bit represented by `which' is currently set, and false otherwise.
*/
inline bool LongOperationManager::IsRunning(
  LongOpType which)
	{
	return ((longOps & which) != 0);
	}

#endif
