#ifndef TOPO_PRIOR_CALCULATOR_HPP
#define TOPO_PRIOR_CALCULATOR_HPP

#include <boost/shared_ptr.hpp>

namespace phycas
{

class Tree;
typedef boost::shared_ptr<Tree> TreeShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Computes topological priors used by BushMove to handle polytomous trees in MCMC analyses. Also provides several
|	utility functions for computing the number of tree topologies with varying degrees of resolution. This class 
|	combines functions from two classes used in the former non-Python version of Phycas, TreeCounter and PriorSummary.
*/
class TopoPriorCalculator
	{
	public:
										TopoPriorCalculator();
										~TopoPriorCalculator();

		bool							IsResolutionClassPrior() const;
		bool							IsPolytomyPrior() const;

		bool							IsRooted() const;
		bool							IsUnrooted() const;

		void							SetNTax(unsigned n);
		unsigned						GetNTax() const;

		void							ChooseRooted();
		void							ChooseUnrooted();

		double							GetCount(unsigned n, unsigned m);
		double							GetSaturatedCount(unsigned n);
		double							GetTotalCount(unsigned n);

		std::vector<double>				GetCountsVect();

		void							ChooseResolutionClassPrior();
		void							ChoosePolytomyPrior();

		void							SetC(double c);
		double							GetC() const;

		double							GetLnTopologyPrior(unsigned m);
		double							GetLnNormalizedTopologyPrior(unsigned m);
		double							GetLnNormConstant();

		std::vector<double>				GetTopoPriorVect();
		std::vector<double>				GetRealizedResClassPriorsVect();

		void							Reset();

	private:

		void							RecalcCountsAndPriorsImpl(unsigned n); // clears topo_priors_dirty
		void							AddNextCount(unsigned k);
		void							RecalcPriorsImpl();

		unsigned						ntax;						/**< TopoPriorCalculator is currently holding counts for this number of taxa */
		bool							is_rooted;					/**< If false, ntax is number of tips in an unrooted tree; if true, ntax is number of tips in a rooted tree */
		bool							is_resolution_class_prior;	/**< True for the resolution class prior, false for the polytomy prior */
		double							C;							/**< Determines strength of prior */

		bool							topo_priors_dirty;			/**< True if `ntax', `C', `is_rooted', or `is_resolution_class_prior' has changed since Reset was last called. Causes Reset to be called if any function is called that requires access to either `counts' or `polytomy_prior' */

		std::vector<double> 			counts;						/**< Vector of length `ntax' holding counts for trees having all possible numbers of internal nodes */
		std::vector<double>				topology_prior;				/**< Vector the mth element of which holds the unnormalized prior probability of a tree with m internal nodes; the 0th element holds the normalizing constant */
	};

typedef boost::shared_ptr<TopoPriorCalculator> TopoPriorCalculatorShPtr;

}	// namespace phycas

#include "pyphy/src/topo_prior_calculator.inl"

#endif
