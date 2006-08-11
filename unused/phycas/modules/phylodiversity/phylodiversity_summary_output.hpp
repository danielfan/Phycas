#ifndef PHYCAS_PHYLODIVERSITY_SUMMARY_OUTPUT_HPP
#define PHYCAS_PHYLODIVERSITY_SUMMARY_OUTPUT_HPP

#include "phycas/modules/phylodiversity/phylodiversity_summary.hpp"

template<class OUT_STREAM>
class GenericPrinterClass<kTabSeparatedTabularLabelsOutStyle, PhylodiversitySummary *, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const PhylodiversitySummary *)
			{
			outStream << "tree";
			outStream << "\tT";
			outStream << "\tI";
#			if defined(NEW_PD_WAY)
				outStream << "\tE";
				outStream << "\tpEI";
				outStream << "\tpET";
				outStream << "\tpIT";
#			else
				outStream << "\tEmin";
				outStream << "\tEmax";
				outStream << "\tEmax/T";
				outStream << "\tEmin/I";
				outStream << "\tpEI";
				outStream << "\tpIT";
				outStream << "\tpET";
#			endif
			outStream << "\tavgSel";
			outStream << "\tmedSel";
			outStream << "\tavgUnsel";
			outStream << "\tmedUnsel";
			outStream << "\n";
			}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kTabSeparatedTabularOutStyle, PhylodiversitySummary, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const PhylodiversitySummary & chainState)
			{
#			if defined(NEW_PD_WAY)
				outStream << MakeStrPrintF("%6d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f", 
					chainState.treeNumber,
					chainState.cumTotal, 
					chainState.cumI, 
					chainState.cumE, 
					chainState.pEI, 
					chainState.pET, 
					chainState.pIT, 
					chainState.avgSelTipLen, 
					chainState.medSelTipLen, 
					chainState.avgUnselTipLen, 
					chainState.medUnselTipLen);
#			else
				outStream << MakeStrPrintF("%6d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f", 
					chainState.treeNumber, 
					chainState.cumTotal, 
					chainState.cumEmin, 
					chainState.cumEmax, 
					chainState.cumI, 
					(chainState.cumEmax/chainState.cumTotal), 
					(chainState.cumEmin/chainState.cumI),
					chainState.pEI,
					chainState.pIT,
					chainState.pET, 
					chainState.avgSelTipLen,
					chainState.medSelTipLen,
					chainState.avgUnselTipLen,
					chainState.medUnselTipLen);
#			endif
			}
	};
template<class OUT_STREAM>
class GenericPrinterClass<kVerboseOutStyle, PhylodiversitySummary, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const PhylodiversitySummary & chainState)
			{
			outStream << "Tree " << chainState.treeNumber;
			outStream << "\n  T      = " << chainState.cumTotal;
			outStream << "\n  I      = " << chainState.cumI;
#			if defined(NEW_PD_WAY)
				outStream << "\n  E      = " << chainState.cumE;
				outStream << "\n  pEI    = " << chainState.pEI;
				outStream << "\n  pET    = " << chainState.pET;
				outStream << "\n  pIT    = " << chainState.pIT;
#			else
				outStream << "\n  Emin   = " << chainState.cumEmin;
				outStream << "\n  Emax   = " << chainState.cumEmax;
				outStream << "\n  Emax/T = " << (chainState.cumEmax/chainState.cumTotal);
				outStream << "\n  Emin/I = " << (chainState.cumEmin/chainState.cumI);
				outStream << "\n  pEI    = " << chainState.pEI;
				outStream << "\n  pIT    = " << chainState.pIT;
				outStream << "\n  pET    = " << chainState.pET;
#			endif
			outStream << "\n  Mean length of terminal branch subtending selected taxon: " << chainState.avgSelTipLen;
			outStream << "\n  Mean length of terminal branch subtending unselected taxon: " << chainState.avgUnselTipLen;
			}
	};

#endif //PHYCAS_PHYLODIVERSITY_SUMMARY_OUTPUT_HPP
