#ifndef PHO_GIBBSHYPERPARAMETER_H
#define PHO_GIBBSHYPERPARAMETER_H

#include <boost/shared_ptr.hpp>
#include "phycas/modules/mcmc/gibbs_param.hpp"

#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::log;
	using std::pow;
	using std::exp;
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	A GibbsParameter that affects edge length priors and thus needs to include the recalculated edge length priors in
|	its full conditional distribution.
*/
class GibbsHyperParameter : public GibbsParameter
	{
	public:
								GibbsHyperParameter(double v, ProbabilityDistributionShPtr p, EdgeLenPriorCalculatorShPtr c, Tree *t);
		
		void					RefreshEdgeLenPrior(double v = 0.0);
		double					CalcLnProb(double v);	// overrides GibbsParameter virtual function

	private:

		ProbabilityDistributionShPtr	edgelenPrior;
		EdgeLenPriorCalculatorShPtr		edgeLenPriorCalc;
		Tree							*tree;
	};

typedef boost::shared_ptr<GibbsHyperParameter> GibbsHyperParameterShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Calls base class (SliceSampler) constructor, initializes `lnL' to 0.0 and `val' to the supplied value.
*/
inline GibbsHyperParameter::GibbsHyperParameter(
  double v,							/**< is the numerical value with which to initialize this parameter */
  ProbabilityDistributionShPtr p,	/**< points to the ProbabilityDistribution object that represents the edge length prior */
  EdgeLenPriorCalculatorShPtr c,	/**< points to the object that recomputes the joint edge length prior */
  Tree *t)							/**< is the tree containing the edges for which priors are needed */
  : GibbsParameter(v)
	{
	edgelenPrior		= p;
	edgeLenPriorCalc	= c;
	tree				= t;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log of the probability density for the value `v' of the hyperparameter. Uses the base class version of
|	the function to do most of the work, but turns off calculation of the log-likelihood, as the likelihood is, by 
|	definition, simply a constant in the full conditional distribution of a hyperparameter.
*/
inline double GibbsHyperParameter::CalcLnProb(
  double v)	/**< is the numerical value for which the log-probability density is to be calculated */
	{
	bool prev_ignoring_data = ignoring_data;
	ignoring_data = true;
	double ln_hyperprior = GibbsParameter::CalcLnProb(v);
	ignoring_data = prev_ignoring_data;

#if defined(HYPERDOOF)
	char tmps[256];
	std::ofstream hyperdoof;
	hyperdoof.open("hyperdoof.txt", std::ios::out | std::ios::app);
	hyperdoof << "\nCalcLnProb(";
	sprintf(tmps, "%.6f", v);
	hyperdoof << tmps << ')' << std::endl;
	hyperdoof.close();
#endif

	double ln_edgelen_prior = 0.0;
	if (ln_hyperprior > -DBL_MAX)
		{
		RefreshEdgeLenPrior(v);
		ln_edgelen_prior = edgeLenPriorCalc->GetLnEdgeLenPrior(*tree);
		}
	double ln_prior = ln_hyperprior + ln_edgelen_prior;

#if defined(HYPERDOOF)
	hyperdoof.open("hyperdoof.txt", std::ios::out | std::ios::app);
	if (ln_hyperprior > -DBL_MAX)
		{
		sprintf(tmps, "  ln_prior = %.6f", ln_prior);
		hyperdoof << tmps << std::endl;
		sprintf(tmps, "    ln_hyperprior = %.6f", ln_hyperprior);
		hyperdoof << tmps << std::endl;
		sprintf(tmps, "    ln_edgelen_prior = %.6f", ln_edgelen_prior);
		hyperdoof << tmps << std::endl;
		}
	else
		{
		hyperdoof << "   hyperprior out of bounds" << std::endl;
		}
	hyperdoof.close();
#endif

	return ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new edge length prior ProbabilityDistribution object with mean equal to `val'.
*/
inline void GibbsHyperParameter::RefreshEdgeLenPrior(double v)
	{
	double mean = v;
	if (mean == 0.0)
		mean = val;
	assert(mean > 0.0);
	//@POL better to create new object or add SetMean member function to ProbabilityDistribution?
	edgelenPrior->SetMeanAndVariance(mean, pow(mean, 2.0));
	}

#endif
