#include "phycas/force_include.h"
#include "phycas/modules/mcmc/gibbs_param.hpp"
#include "phycas/trees/hky_ad_hoc.hpp" // to be deprecated

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides SliceSampler pure virtual function. Informs model that parameter value has changed, then asks model to
|	calculate the log-likelihood and log-prior, returning the log of the posterior density.
*/
double GibbsParameter::CalcLnProb(double v)
	{
	assert(model != NULL);

	// Base class version increments func_evals counter
	//
	SliceSampler::CalcLnProb(val);

	// Give model heads-up that a parameter value has changed (ParameterChanged currently invalidates all 
	// conditional likelihood arrays)
	//
	val = v;
	model->ParameterChanged();

	lnL = (ignoring_data ? 0.0 : model->CalcLogLikelihood());
	double log_prior = priorDistr->GetRelativeLnPDF(val);
	double log_posterior = -DBL_MAX;
	if ((lnL > -DBL_MAX) && (log_prior > -DBL_MAX))
		{
		log_posterior = lnL + log_prior;
		}

	return log_posterior;
	}
