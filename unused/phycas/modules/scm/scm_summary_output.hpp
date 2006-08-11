#ifndef PHYCAS_SCM_SUMMARY_OUTPUT_HPP
#define PHYCAS_SCM_SUMMARY_OUTPUT_HPP

#include "phycas/modules/scm/scm_summary.hpp"

template<class OUT_STREAM>
class GenericPrinterClass<kTabSeparatedTabularLabelsOutStyle, SCMSummary *, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const SCMSummary *summary)
			{
			outStream << "supertree = " << summary->supertreeDescription;
			outStream << "\n";
			}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kTabSeparatedTabularOutStyle, SCMSummary, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM &outStream, const SCMSummary &summary)
			{
				outStream << MakeStrPrintF("%6d", summary.treeNumber + 1);
			}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kVerboseOutStyle, SCMSummary, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM &outStream, const SCMSummary &summary)
			{
			outStream << "supertree = " << summary.supertreeDescription;
			}
	};

#endif //PHYCAS_SCM_SUMMARY_OUTPUT_HPP
