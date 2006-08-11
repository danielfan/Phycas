#if defined(WIN_PHOREST)

#include <windows.h>
#include <cassert>
#include "ncl/output/nxs_output.hpp"
#include "phycas/misc/pwd.hpp"
using ncl::endl;

PWD::PWD()
	{
	}

CmdResult PWD::HandlePWD()
	{
	char pathstr[1024];
	DWORD len = ::GetCurrentDirectory(1024, pathstr);
	if (len >= 1024)
		return kCmdFailedSilent;
	NxsOutput::GetNullCheckingStream<NxsOutput::kOutput>()  << pathstr << endl;
	return kCmdSucceeded;
	}

#endif // WIN_PHOREST

