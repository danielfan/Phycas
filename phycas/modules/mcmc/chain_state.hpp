#if !defined (CHAIN_STATE_HPP) 
#define CHAIN_STATE_HPP
#include "phycas/trees/tree_output.hpp"
class ChainState
	{
	public:
		unsigned		totalSteps;
		double			lnLikelihood;
		double			cBase[4];
		double			kappa;
		double			rateVar;
		double			edgelenHyperPrior;
		Tree		  * treeAlias;
		ChainState()
			{
			clear();
			}
		
		void clear()
			{
			treeAlias = NULL;
			totalSteps = UINT_MAX;
			lnLikelihood = -DBL_MAX;
			}
	};
	
template<class OUT_STREAM>
class GenericPrinterClass<kNexusTreeCmdOutStyle, ChainState, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const ChainState & chainState)
			{
			assert(chainState.treeAlias != NULL);
			if (chainState.treeAlias != NULL)
				{
				outStream << "tree MCMC" << chainState.totalSteps << " = [&U]";
				Emit<kNewickTreeOutStyle>(outStream, *chainState.treeAlias) << ';';
				}
			}
	};
			
template<class OUT_STREAM>
class GenericPrinterClass<kTabSeparatedTabularOutStyle, ChainState, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const ChainState & chainState)
			{
			std::string tmp;

			const double cTotal = chainState.cBase[0] + chainState.cBase[1] + chainState.cBase[2] + chainState.cBase[3];
			outStream << chainState.totalSteps;

			tmp.clear();
			StrPrintF(tmp, "\t%.6f", chainState.lnLikelihood); //@POL weekly chat: only getting one or two decimal places in likelihood, but not sure this is best way to format
			outStream << tmp;

			for (unsigned i = 0; i < 4; ++i)
				{
				tmp.clear();
				StrPrintF(tmp, "\t%.6f", chainState.cBase[i]);
				outStream << tmp;
				}

			for (unsigned i = 0; i < 4; ++i)
				{
				tmp.clear();
				StrPrintF(tmp, "\t%.6f", (chainState.cBase[i]/cTotal));
				outStream << tmp;
				}

			tmp.clear();
			StrPrintF(tmp, "\t%.6f", chainState.kappa);
			outStream << tmp;

			tmp.clear();
			StrPrintF(tmp, "\t%.6f", chainState.rateVar);
			outStream << tmp;

			tmp.clear();
			StrPrintF(tmp, "\t%.6f", (1.0/chainState.rateVar));
			outStream << tmp;

			tmp.clear();
			StrPrintF(tmp, "\t%.6f", chainState.edgelenHyperPrior);
			outStream << tmp;

			//outStream << '\n';	//@POL weekly chat: why do I not need \n here but do need it for labels?
			}
	};
	
template<class OUT_STREAM>
class GenericPrinterClass<kTabSeparatedTabularLabelsOutStyle, ChainState *, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const ChainState *)
			{
			outStream << "cycle";
			outStream << "\tlnL";
			outStream << "\tcA";
			outStream << "\tcC";
			outStream << "\tcG";
			outStream << "\tcT";
			outStream << "\tfreqA";
			outStream << "\tfreqC";
			outStream << "\tfreqG";
			outStream << "\tfreqT";
			outStream << "\tkappa";
			outStream << "\tvariance";
			outStream << "\tshape";
			outStream << "\tmeanedgelen";
			outStream << '\n';	//@POL weekly chat: why doesn't endl work here? all output was going on one line
			}
	};

#endif
