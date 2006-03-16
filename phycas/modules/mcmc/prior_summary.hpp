#if !defined (PRIOR_SUMMARY_HPP) 
#define PRIOR_SUMMARY_HPP
#include "ncl/output/temp_basic_output_operators.hpp"
#include "phycas/command/mcmc_settings.hpp"
#include "phycas/trees/tree_counter.hpp"

using std::fabs; // VKJ 11/08 - added so would compile

class Tree;
class PriorSummary
	{
	public:
		PriorSummary(const MCMCSettings & inSettings, 
		  UInt						inNTax)
			:settings(inSettings),
			nTax(inNTax)
			{
			RecalcTopologyPriors();
			}
			
		PriorSummary() 
			:nTax(0)
			{
			}
		void RecalcConstTopologyPriors(double polytomyLnPriorRatio);
		void RecalcFlatTopologyPriors();
		void RecalcGeometricTopologyPriors(double pTopoPrior);
		/*----------------------------------------------------------------------------------------------------------------------
		|	Recalculates `vTopoPriorByResClass', `vExpectedResClassPrior' and `vInducedPolytomyPrior' vectors. Uses one of three
		|	functions: RecalcGeometricTopologyPrior (if geometric topology prior specified), RecalcConstTopologyPrior (if the
		|	constant factor approach to computing priors was specified), or RecalcFlatTopologyPrior (if a flat prior over all
		|	topologies was specified.
		*/
		void	RecalcTopologyPriors()
			{		//@POL, I removed the conditional compilation here because we now have a user setting.  OK?  mth
			if (nTax > 0)
				{
				if (settings.resClassPrior)
					RecalcGeometricTopologyPriors(settings.pTopoPrior);
				else
					RecalcConstTopologyPriors(settings.polytomyLnPriorRatio);
				}
			}
		double GetLnTopologyPrior(const Tree & tree);
		double GetLnInducedPrior(UInt m) const;
		
		MCMCSettings		settings;
		UInt				nTax;
		VecDbl				vInducedResClassPrior;
		VecDbl				vExpectedResClassPrior;
		VecDbl				vTopoPriorByResClass;
		mutable TreeCounter treeCounter;					/* calculates and stores vector of numbers of trees of all possible resolutions for a given number of taxa */
	};

template<int FORMAT_HINT, class OUT_STREAM>
class GenericPrinterClass<FORMAT_HINT, PriorSummary, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const PriorSummary & ps)
			{
#			if defined(C_FUNCS_IN_STD_NAMESPACE)
				using std::exp;
#			endif
			if (!ps.settings.resClassPrior)
				{
				if (ps.settings.topoPriorFlat || fabs(ps.settings.polytomyLnPriorRatio) < 1.0e-6)
					{
					outStream << "\n\nTopology prior: all tree topologies equally probable";
					outStream << "\n  (i.e. star tree has same prior probability as a fully-resolved tree)";
					}
				else
					{
					outStream << "\n\nTopology prior: log prior ratio in favor of less resolved topology equal to " << ps.settings.polytomyLnPriorRatio;
					outStream << "\n  (tree with one fewer internal node favored by a factor of " << exp(ps.settings.polytomyLnPriorRatio) << ")\n";
					}
				}
			else
				{
				outStream << "\n\nTopology prior: geometric prior with p = " << ps.settings.pTopoPrior;
				if (ps.settings.pTopoPrior < 0.25)
					outStream << "\n  (favors less-resolved trees slightly more than resolved ones)";
				else if (ps.settings.pTopoPrior < 0.75)
					outStream << "\n  (favors less-resolved trees somewhat more than resolved ones)";
				else
					outStream << "\n  (favors less-resolved trees much more than resolved ones)";
				}
			if (ps.settings.bushMoveWeight == 0.0)
				outStream << "Note:  less-than-fully-resolved trees are not considered in this sampling because the bushMoveWeight is set to 0)\n";
			if (ps.nTax < 3)
				return;
			unsigned m;
			double unresolved_induced_prior = 0.0;
			double unresolved_expected_prior = 0.0;
			double resolved_induced_prior =  exp(ps.vInducedResClassPrior[ps.nTax - 3]);
			double resolved_expected_prior = exp(ps.vExpectedResClassPrior[ps.nTax - 3]);
			std::string tmp;
			outStream << "\n      Internal nodes      Log(no. trees)       Induced prior       Nominal prior";
			outStream << "\n--------------------------------------------------------------------------------";
			for (m = 0; m < ps.nTax - 2; ++m)
				{
				outStream << MakeStrPrintF("\n%20d%20.1f%20.5f%20.5f", m + 1, log(ps.treeCounter.GetUnrootedCount(ps.nTax, m+1)), exp(ps.vInducedResClassPrior[m]), exp(ps.vExpectedResClassPrior[m]));
				if (m < ps.nTax - 3)
					{
					unresolved_induced_prior += exp(ps.vInducedResClassPrior[m]);
					unresolved_expected_prior += exp(ps.vExpectedResClassPrior[m]);
					}
				}
			outStream << "\n--------------------------------------------------------------------------------" << endl;
			outStream << "\nInduced prior summary:";
			outStream << "\n    proportion due to unresolved classes:   " << unresolved_induced_prior;
			outStream << "\n    proportion due to fully-resolved class: " << resolved_induced_prior;
			outStream << "\n    unresolved/fully-resolved prior ratio:  " << (unresolved_induced_prior/resolved_induced_prior);
			outStream << "\nExpected prior summary:";
			outStream << "\n    proportion due to unresolved classes:   " << unresolved_expected_prior;
			outStream << "\n    proportion due to fully-resolved class: " << resolved_expected_prior;
			outStream << "\n    unresolved/fully-resolved prior ratio:  " << (unresolved_expected_prior/resolved_expected_prior) << '\n';
			}
	};

#endif

