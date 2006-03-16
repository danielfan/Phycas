#include "phycas/force_include.h"
#include <cassert>
#include "phycas/misc/long_operation_manager.hpp"
#include "ncl/output/nxs_output.hpp"
#include "ncl/output/nxs_user_query.hpp"
using std::string;

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `longOps' data member to 0 and `ctrlc_handler_installed' data member to false.
*/
LongOperationManager::LongOperationManager()
	:ctrlc_handler_installed(false),
	longOps(LongOp_None),
	ctrl_c(false)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	In graphical interfaces, this function also changes the cursor to the hourglass cursor and may do things like pop
|	up a progress dialog, but here we just set the `longOps' bit corresponding to `which'.
*/
void LongOperationManager::StartLongOperation(
  LongOpType which)
	{
	assert((longOps & which) == 0);
	longOps |= which;

	if (!ctrlc_handler_installed)
		{
		ctrlc_handler_installed = true;

#		if defined(WIN_PHOREST) && defined(CONSOLE_PHOREST)
			bool installed_ok = (::SetConsoleCtrlHandler((PHANDLER_ROUTINE)PhoInterruptHandler, TRUE) == TRUE);
			assert(installed_ok);
#		elif defined(CONSOLE_PHOREST)
			::signal(SIGINT, PhoInterruptHandler);
#		endif
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	In graphical interfaces, this function also changes the cursor back to the standard arrow cursor and destroys or
|	hides any progress dialogs that have been displayed, but here we just clear the `longOps' bit corresponding to
|	`which'.
*/
void LongOperationManager::StopLongOperation(
  LongOpType which)
	{
	assert((LongOpType) (longOps & which) == which);
	longOps &= ~which;

	assert(ctrlc_handler_installed);
	if (longOps == 0)
		{
		// The fact that longOps == 0 means that the last long operation has completed, so we are free to unintall 
		// the interrupt handler now
		//
		ctrlc_handler_installed = false;

		// Nothing to do for Linux console version
		//
#		if defined(WIN_PHOREST) && defined(CONSOLE_PHOREST)
			bool uninstalled_ok = (::SetConsoleCtrlHandler((PHANDLER_ROUTINE)PhoInterruptHandler, FALSE) == TRUE);
			assert(uninstalled_ok);
#		endif
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	In portable command-line interfaces, this function detects when the user presses Ctrl-C and then prompts the user
|	for instructions about whether to continue or just stop the current long operation. Returns true if user really 
|	wants to abort the long operation, otherwise returns false. Note that this function does not call StopLongOperation.
*/
bool LongOperationManager::CheckAbort(
  LongOpType which)
	{
	NxsUserQuery * userQuery  = NxsOutput::GetUserQueryPtr();
	if (ctrl_c)
		{
		ctrl_c = false;
		if (userQuery == NULL)
			return true; // abort if we have 
		string msg				= "Abort?";
		string title			= "Ctrl-C pressed";
		if (which == LongOp_MCMC)
				{
				title.clear();
				msg = "Abort MCMC run?";
				}
		if (userQuery->UserChoice(title, msg, "No|Yes", 1, 1))
			return true;
#		if defined(CONSOLE_PHOREST) && ! defined (WIN_PHOREST)
			else
				signal(SIGINT, PhoInterruptHandler);
#		endif
		}

	return false;
	}

#if defined(WIN_PHOREST) && defined(CONSOLE_PHOREST)

BOOL WINAPI PhoInterruptHandler(DWORD dwCtrlType)
	{
	if (dwCtrlType == CTRL_C_EVENT)
		{
		LongOperationManagerSingletonHolder::Instance().SetInterrupted(true);
		return TRUE;
		}
	else
		return FALSE;
	}

#elif defined(CONSOLE_PHOREST)

// Linux console version
//
void PhoInterruptHandler(int)
	{
	LongOperationManagerSingletonHolder::Instance().SetInterrupted();
	}

#endif

