#include "phycas/force_include.h"
#include <list>
#include <ctype.h> 
#include "ncl/nxs_defs.hpp"
#if defined(NCL_USE_NXS_CONSOLE_OUTPUT)
#	include "ncl/misc/algorithm_extensions.hpp"
#	include "ncl/nxs_token.hpp"
#	include "ncl/nxs_basic_manager.hpp"
#	include "ncl/misc/nxs_index_set.hpp"
#	include "ncl/output/nxs_output_stream.hpp"
	using ncl::flush;

	void NxsStdOutputStream::DisplaySet(const NxsIndexSet &s, bool useNumbers, const NxsIndexInfo &indInfo)
		{
		*this << s.name << " = ";
		if (useNumbers)
			*this << s.GetNexusDescription(true);
		else
			{
			for (NxsIndexSet::const_iterator sIt = s.begin(); sIt != s.end(); ++sIt)
				*this << indInfo.GetLabel(*sIt) << ' ';	
			}
		*this << flush;		
		}

	NxsStdOutputStream::NxsStdOutputStream()
	  :dblFormat(UINT_MAX, 6),
	  outputToMonitor(true),
	  logFile(NULL)
		{
		message.clear();
		typist.SetOutputStream(this);
		SetWidth(kDefaultPrintWidth);
		}

#endif // defined(NCL_USE_NXS_CONSOLE_OUTPUT)
