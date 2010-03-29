/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/basic_tree_node.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/dirichlet_move.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/tree_manip.hpp"
#include "phycas/src/gtr_model.hpp"

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor sets tuning parameter `psi' and `max_psi' to 300.0 and calls reset().
*/
DirichletMove::DirichletMove() : MCMCUpdater()
	{
	is_move = true;
	dim = 0;
	boldness = 0.0;
	max_psi = 1.0;
	max_psi = 300.0;
	psi = max_psi;
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|   Resets the 'dir_forward' and 'dir_reverse' shared pointers.
*/
void DirichletMove::reset()
	{
    dir_forward.reset();
    dir_reverse.reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'psi', which is the tuning parameter for this move.
*/
void DirichletMove::setPsi(
  double x) /* is the new value for psi */
	{
	psi = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'min_psi', which is the minimum value of the tuning parameter for this move.
|   The minimum value is used whenever boldness comes into play, for example during a path sampling analysis when the
|   target distribution changes from posterior to prior and boldness is adjusted accordingly.
*/
void DirichletMove::setPriorTuningParam(
  double x) /* is the new value for max_psi */
	{
	min_psi = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'max_psi', which is the maximum value of the tuning parameter for this move.
|   The maximum value is the default value, used for exploring the posterior distribution.
*/
void DirichletMove::setPosteriorTuningParam(
  double x) /* is the new value for max_psi */
	{
	max_psi = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'dim', which is the number of parameters updated jointly by this move.
*/
void DirichletMove::setDimension(
  unsigned d) /* is the new value for psi */
	{
	dim = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member 'psi', which is the tuning parameter for this move, based on a boldness value
|	that ranges from 0 (least bold) to 100 (most bold). The actual formula used is psi = 3*(100 - x), with minimum of 
|	0 and maximum of `max_psi' enforced after the calculation. Current parameter values are multiplied by the value `psi'
|	to obtain the parameters of a new Dirichlet distribution used for sampling the proposed new values. The boldest
|	distribution that can be used for sampling is a flat dirichlet, however, so if any of the parameters are less than
|	1 after multiplication by `psi', they are simply set to 1.
*/
void DirichletMove::setBoldness(
  double x) /* is the new boldness value */
	{
	boldness = x;
	if (boldness < 0.0)
		boldness = 0.0;
	else if (boldness > 100.0)
		boldness = 100.0;

    //  Assume:
    //    min_psi = 5
    //    max_psi = 300
    //
    //  psi = min_psi + (max_psi - min_psi)*(100 - boldness)/100
    //
    //  boldness   psi calculation
    //      0      5 + (300 - 5)*(100 - 0)/100   = 300
    //    100      5 + (300 - 5)*(100 - 100)/100 = 5
    //
	psi = min_psi + (max_psi - min_psi)*(100.0-boldness)/100.0;
	//std::cerr << boost::str(boost::format("####### DirichletMove::setBoldness(%.5f), boldness = %.5f, psi = %.5f, max_psi = %.5f, min_psi = %.5f") % x % boldness % psi % max_psi % min_psi) << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member 'dim', which is the number of parameters updated jointly by this move.
*/
unsigned DirichletMove::getDimension() const
	{
	return dim;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member 'psi', which is the tuning parameter for this move.
*/
double DirichletMove::getPsi() const
	{
	return psi;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|   Chooses a new vector of state frequencies using a sharp Dirichlet distribution centered at the original frequencies.
*/
void DirichletMove::proposeNewState()
	{
    // copy the current parameters from the model to the data member orig_params
    getParams();
    
    // create vector of Dirichlet parameters for selecting new parameter values (by multiplying
    // each current parameter value by the value `psi')
    unsigned sz = (unsigned)orig_params.size();
    c_forward.resize(sz);
	std::transform(orig_params.begin(), orig_params.end(), c_forward.begin(), 1.0 + boost::lambda::_1*psi);

	//std::cerr << "c_forward:" << std::endl;
	//std::copy(c_forward.begin(), c_forward.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << std::endl;
	
	// create Dirichlet distribution and sample from it to get proposed frequencies
    dir_forward = DirichletShPtr(new DirichletDistribution(c_forward));
    dir_forward->SetLot(getLot().get());
    new_params = dir_forward->Sample();

	//std::cerr << "new_params:" << std::endl;
	//std::copy(new_params.begin(), new_params.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << std::endl;
	
    // create vector of Dirichlet parameters for selecting old frequencies (needed for Hasting ratio calculation)
    c_reverse.resize(new_params.size());
	std::transform(new_params.begin(), new_params.end(), c_reverse.begin(), 1.0 + boost::lambda::_1*psi);
    dir_reverse = DirichletShPtr(new DirichletDistribution(c_reverse));
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move. The Hastings ratio is:
|>
|	Pr(new state -> old state)   Dir(c'[i]*psi) density evaluated at c
|	-------------------------- = -------------------------------------
|	Pr(old state -> new state)   Dir(c[i]*psi) density evaluated at c'
|>
|   where c[i] is the ith current parameter and c'[i] is the ith proposed new parameter.
*/
double DirichletMove::getLnHastingsRatio() const
	{
    return (dir_reverse->GetLnPDF(orig_params) - dir_forward->GetLnPDF(new_params));
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	This move does not change the model dimension, so the Jacobian is irrelevant.
*/
double DirichletMove::getLnJacobian() const
	{
	return 0.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is accepted. Simply calls the reset() function.
*/
void DirichletMove::accept()
	{
	MCMCUpdater::accept();
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Calls the reset() function after reinstating the original state frequencies
|   and ensuring that all conditional likelihood arrays will be recalculated when the likelihood is next calculated.
*/
void DirichletMove::revert()
	{
	MCMCUpdater::revert();
    setParams(orig_params);

    // invalidate all CLAs
    likelihood->useAsLikelihoodRoot(NULL);	

	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool DirichletMove::update()
	{
	if (is_fixed)
		return false;
		
	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);

	proposeNewState();

	// Note that it may seem incorrect to use a Dirichlet prior here because transformed versions of the
	// orig_params values are used in the likelihood function. However, the transformation only adds a constant
	// (i.e., -log(nsubsets)) to the log prior and thus this constant cancels out in the ln_accept_ratio.
	// More care should be taken when reporting the log prior or when estimating marginal likelihoods.
    double prev_ln_prior		= (p->getSAMCLikelihoodOnly() ? 0.0 : mv_prior->GetLnPDF(orig_params));
	double prev_ln_like			= p->getLastLnLike();
	PHYCAS_ASSERT(!use_working_prior || mv_working_prior);
	double prev_ln_working_prior = (use_working_prior ? mv_working_prior->GetLnPDF(orig_params) : 0.0);

    // replace current parameter values with new ones
    setParams(new_params);
    likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs

	double curr_ln_prior		= (p->getSAMCLikelihoodOnly() ? 0.0 : mv_prior->GetLnPDF(new_params));
	double curr_ln_like			= (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
	double curr_ln_working_prior = (use_working_prior ? mv_working_prior->GetLnPDF(new_params) : 0.0);

    double prev_posterior = 0.0;
	double curr_posterior = 0.0;
	unsigned prev_samc_index = UINT_MAX;
	unsigned curr_samc_index = UINT_MAX;
	
	if (p->doingSAMC())
		{
		prev_posterior		= prev_ln_like + prev_ln_prior;
		prev_samc_index		= p->getSAMCEnergyLevel(prev_posterior);
		double prev_theta	= p->getSAMCWeight(prev_samc_index);
		prev_posterior		-= prev_theta;
		
		curr_posterior		= curr_ln_like + curr_ln_prior;
		curr_samc_index		= p->getSAMCEnergyLevel(curr_posterior);
		double curr_theta	= p->getSAMCWeight(curr_samc_index);
		curr_posterior		-= curr_theta;
		}
    else
		{
		if (is_standard_heating)
			{
			prev_posterior = heating_power*(prev_ln_like + prev_ln_prior);
			curr_posterior = heating_power*(curr_ln_like + curr_ln_prior);
			if (use_working_prior)
				{
				prev_posterior += (1.0 - heating_power)*prev_ln_working_prior;
				curr_posterior += (1.0 - heating_power)*curr_ln_working_prior;
				}
			}
		else
			{
			prev_posterior = heating_power*prev_ln_like + prev_ln_prior;
			curr_posterior = heating_power*curr_ln_like + curr_ln_prior;
			}
		}
	double ln_hastings			= getLnHastingsRatio();
	double ln_accept_ratio		= curr_posterior - prev_posterior + ln_hastings;

    double lnu = std::log(rng->Uniform(FILE_AND_LINE));
    
    //     std::cerr << "\nDirichletMove::update" << std::endl;
    //     std::cerr << boost::str(boost::format("  boldness = %.1f, psi = %.5f\n") % boldness % psi);
    //     std::vector<double>::const_iterator origit = orig_params.begin();
    //     std::vector<double>::const_iterator newit = new_params.begin();
    //     double min_new = 1.0;
    //     double max_new = 0.0;
    //     for (unsigned i = 0; i < 61; ++i, ++origit, ++newit)
    //         {
    //         std::cerr << boost::str(boost::format("%20.6f %20.6f") % (*origit) % (*newit)) << std::endl;
    //         if (*newit < min_new)
    //             min_new = *newit;
    //         if (*newit > max_new)
    //             max_new = *newit;
    //         }
    //     std::cerr << "  min_new         = " << min_new << std::endl;
    //     std::cerr << "  max_new         = " << max_new << std::endl;
    //     std::cerr << "  curr_posterior  = " << curr_posterior << std::endl;
    //     std::cerr << "  prev_posterior  = " << prev_posterior << std::endl;
    //     std::cerr << "  ln_posterior_ratio = " << (curr_posterior - prev_posterior) << std::endl;
    //     std::cerr << "  ln_hastings        = " << ln_hastings << std::endl;
    //     std::cerr << "  ln_accept_ratio    = " << ln_accept_ratio << std::endl;
    //     std::cerr << "  lnu                = " << lnu << std::endl;
        
	if (ln_accept_ratio >= 0.0 || lnu <= ln_accept_ratio)
		{
	    if (save_debug_info)
    	    {
			debug_info = boost::str(boost::format("ACCEPT, prev_ln_working_prior = %.5f, curr_ln_working_prior = %.5f, prev_ln_prior = %.5f, curr_ln_prior = %.5f, prev_ln_like = %.5f, curr_ln_like = %.5f, lnu = %.5f, ln_accept_ratio = %.5f\norig_params: ") % prev_ln_working_prior % curr_ln_working_prior % prev_ln_prior % curr_ln_prior % prev_ln_like % curr_ln_like % lnu % ln_accept_ratio);
			for (std::vector<double>::const_iterator it = orig_params.begin(); it != orig_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\nnew_params:  ";
			for (std::vector<double>::const_iterator it = new_params.begin(); it != new_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\nc_forward:  ";
			for (std::vector<double>::const_iterator it = c_forward.begin(); it != c_forward.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\ndir_forward:  ";
			debug_info += dir_forward->GetDistributionDescription();
			}
		p->setLastLnPrior(curr_ln_prior);
		p->setLastLnLike(curr_ln_like);
		if (p->doingSAMC())
			p->updateSAMCWeights(curr_samc_index);
		accept();
		return true;
		}
	else
		{
	    if (save_debug_info)
    	    {
			debug_info = boost::str(boost::format("REJECT, prev_ln_working_prior = %.5f, curr_ln_working_prior = %.5f, prev_ln_prior = %.5f, curr_ln_prior = %.5f, prev_ln_like = %.5f, curr_ln_like = %.5f, lnu = %.5f, ln_accept_ratio = %.5f\norig_params: ") % prev_ln_working_prior % curr_ln_working_prior % prev_ln_prior % curr_ln_prior % prev_ln_like % curr_ln_like % lnu % ln_accept_ratio);
			for (std::vector<double>::const_iterator it = orig_params.begin(); it != orig_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\nnew_params:  ";
			for (std::vector<double>::const_iterator it = new_params.begin(); it != new_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\nc_forward:  ";
			for (std::vector<double>::const_iterator it = c_forward.begin(); it != c_forward.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\ndir_forward:  ";
			debug_info += dir_forward->GetDistributionDescription();
			}
		curr_ln_like	= p->getLastLnLike();
		curr_ln_prior	= p->getLastLnPrior();
		if (p->doingSAMC())
			p->updateSAMCWeights(prev_samc_index);
		revert();
		return false;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version simply returns an empty vector. Override this function in derived classes to return a vector
|	of parameter values.
*/
double_vect_t DirichletMove::listCurrValuesFromModel()
	{
	double_vect_t v;
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override this base class version to add the current parameter vector to the data already stored in 
|	`mv_fitting_sample'.
*/
void DirichletMove::educateWorkingPrior()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use samples in `fitting_sample' to parameterize `working_prior'. This function is called during s-cubed style
|	steppingstone sampling after the initial phase of sampling from the posterior so that the working prior can be
|	used for the remaining phases. Assumes `fitting_sample' has more than 1 element. Assigns a DirichletDistribution 
|	object to `working_prior'.
*/
void DirichletMove::finalizeWorkingPrior()
	{
	if (!isFixed())
		{
		PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
		PHYCAS_ASSERT(mv_fitting_sample.size() > 1);

		// Let a, b, c, d be the parameters of a Dirichlet(a,b,c,d).
		// Let phi = a + b + c + d
		// Let mu_1, mu_2, mu_3, and mu_4 be sample means
		// Let s_1^2, s_2^2, s_3^2, and s_4^2 be sample variances
		// Noting that mu_1 = a/phi, mu_2 = b/phi, mu_3 = c/phi and mu_4 = d/phi
		// and noting that s_1^2 = a*(b + c + d)/[phi^2*(phi + 1)], etc., and
		// letting z = s_1^2/mu_1 + s_2^2/mu_2 + s_3^2/mu_3 + s_4^2/mu_4,
		// phi can be estimated as (3/z) - 1, where the 3 is really k-1 
		// and k is the number of Dirichlet parameters. Now, 
		// a = mu_1*phi, b = mu_2*phi, c = mu_3*phi and d = mu_4*phi
		
		// First compute the sample means and variances using the data stored in mv_working_prior
		double_vect_t sums(dim, 0.0);
		double_vect_t ss(dim, 0.0);
		double n = 0.0;
		for (double_vect_vect_t::iterator i = mv_fitting_sample.begin(); i != mv_fitting_sample.end(); ++i)
			{
			unsigned k = 0;
			n += 1.0;
			for (double_vect_t::iterator j = (*i).begin(); j != (*i).end(); ++j)
				{
				double v = (*j);
				sums[k] += v;
				ss[k++] += v*v;
				}
			}
			
		double_vect_t means(dim, 0.0);
		double_vect_t variances(dim, 0.0);
		for (unsigned i = 0; i < dim; ++i)
			{
			double mean = sums[i]/n;
			means[i] = mean;
			variances[i] = (ss[i] - n*mean*mean)/(n - 1.0);
			}
		
		// Now compute the Dirichlet parameters
		double z = 0.0;
		for (unsigned i = 0; i < dim; ++i)
			{
			z += variances[i]/means[i];
			}
			
		double phi = (double)(dim - 1)/z - 1.0;
		double_vect_t params;
		for (unsigned i = 0; i < dim; ++i)
			{
			params.push_back(phi*means[i]);
			}
			
		mv_working_prior = MultivarProbDistShPtr(new DirichletDistribution(params));
		//	if (params.size() == 4)
		//		std::cerr << boost::str(boost::format("@@@@@@@@@ working prior is Dirichlet(%g,%g,%g,%g) for updater %s") % params[0] % params[1] % params[2] % params[3] % getName()) << std::endl;
		//	else if (params.size() == 6)
		//		std::cerr << boost::str(boost::format("@@@@@@@@@ working prior is Dirichlet(%g,%g,%g,%g,%g,%g) for updater %s") % params[0] % params[1] % params[2] % params[3] % params[4] % params[5] % getName()) << std::endl;
		//	else
		//		std::cerr << boost::str(boost::format("@@@@@@@@@ working prior is Dirichlet(not 4 or 6 params) for updater %s") % getName()) << std::endl;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor simply calls the base class (DirichletMove) constructor.
*/
StateFreqMove::StateFreqMove() : DirichletMove()
	{
	dim = 4;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class version adds the current vector of state frequencies to the data already stored in 
|	`mv_fitting_sample'.
*/
void StateFreqMove::educateWorkingPrior()
	{
	if (!isFixed())
		{
		double_vect_t rfreqs;
		getCurrValuesFromModel(rfreqs);
		mv_fitting_sample.push_back(rfreqs);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the state frequencies of the associated HKY or GTR model to those in the supplied vector `v'.
*/
void StateFreqMove::sendCurrValuesToModel(const double_vect_t & v)
	{
	PHYCAS_ASSERT(dim == v.size());
	model->setStateFreqsUnnorm(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current state frequencies from the model, storing them in the supplied vector `v'.
*/
void StateFreqMove::getCurrValuesFromModel(double_vect_t & v) const
	{
	PHYCAS_ASSERT(dim > 0);
	if (model)
		{
    	const std::vector<double> & rfreqs = model->getStateFreqs();
    	v.resize(rfreqs.size());
		PHYCAS_ASSERT(dim == rfreqs.size());
    	std::copy(rfreqs.begin(), rfreqs.end(), v.begin());
    	}
    else
    	{
    	v.assign(dim, 1.0/(double)dim);
    	}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current state frequencies from the model, returning them as an anonymous vector.
*/
double_vect_t StateFreqMove::listCurrValuesFromModel()
	{
	PHYCAS_ASSERT(dim > 0);
	double_vect_t v(dim);
	if (model)
		{
    	const std::vector<double> & rfreqs = model->getStateFreqs();
		PHYCAS_ASSERT(dim == rfreqs.size());
    	std::copy(rfreqs.begin(), rfreqs.end(), v.begin());
    	}
    else
    	{
    	v.assign(dim, 1.0/(double)dim);
    	}
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current state frequencies from the model, storing them in the data member `orig_params'.
*/
void StateFreqMove::getParams()
	{
	getCurrValuesFromModel(orig_params);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the state frequencies in the model with those supplied in the vector `v'.
*/
void StateFreqMove::setParams(
  const std::vector<double> & v)    /*< is the vector of parameter values to send to the model */
	{
    model->setStateFreqsUnnorm(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor simply calls the base class (DirichletMove) constructor.
*/
RelRatesMove::RelRatesMove() : DirichletMove()
	{
	dim = 6;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class version adds the current vector of relative rates to the data already stored in 
|	`mv_fitting_sample'.
*/
void RelRatesMove::educateWorkingPrior()
	{
	if (!isFixed())
		{
		double_vect_t rrates;
		getCurrValuesFromModel(rrates);
		mv_fitting_sample.push_back(rrates);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the relative rates of the associated GTR model to those in the supplied vector `v'.
*/
void RelRatesMove::sendCurrValuesToModel(const double_vect_t & v)
	{
	PHYCAS_ASSERT(dim == v.size());
	GTR * gtr_model = dynamic_cast<GTR *>(model.get());
	PHYCAS_ASSERT(gtr_model);
	gtr_model->setRelRates(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, storing them in the supplied vector `v'.
*/
void RelRatesMove::getCurrValuesFromModel(double_vect_t & v) const
	{
	PHYCAS_ASSERT(dim > 0);
	GTR * gtr_model = dynamic_cast<GTR *>(model.get());
	if (gtr_model)
		{
		const std::vector<double> & rrates = gtr_model->getRelRates();
		PHYCAS_ASSERT(dim == rrates.size());
		v.resize(rrates.size());
		std::copy(rrates.begin(), rrates.end(), v.begin());
    	}
    else
    	{
    	v.assign(dim, 1.0/(double)dim);
    	}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, returning them as an anonymous vector.
*/
double_vect_t RelRatesMove::listCurrValuesFromModel()
	{
	PHYCAS_ASSERT(dim > 0);
	double_vect_t v(dim);
	GTR * gtr_model = dynamic_cast<GTR *>(model.get());
	if (gtr_model)
		{
		const double_vect_t & rrates = gtr_model->getRelRates();
		PHYCAS_ASSERT(dim == rrates.size());
		std::copy(rrates.begin(), rrates.end(), v.begin());
    	}
    else
    	{
    	v.assign(dim, 1.0/(double)dim);
    	}
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, storing them in the data member `orig_params'.
*/
void RelRatesMove::getParams()
	{
	getCurrValuesFromModel(orig_params);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the relative rates in the model with those supplied in the vector `v'.
*/
void RelRatesMove::setParams(
  const std::vector<double> & v)    /*< is the vector of parameter values to send to the model */
	{
    GTR * gtr_model = dynamic_cast<GTR *>(model.get());
    gtr_model->setRelRates(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor simply calls the base class (DirichletMove) constructor.
*/
SubsetRelRatesMove::SubsetRelRatesMove() : DirichletMove()
	{
	dim = 1;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class version adds the current vector of subset relative rates to the data already stored in 
|	`mv_fitting_sample'.
*/
void SubsetRelRatesMove::educateWorkingPrior()
	{
	if (!isFixed())
		{
		double_vect_t rrates;
		getCurrValuesFromModel(rrates);
		std::transform(rrates.begin(), rrates.end(), rrates.begin(), boost::lambda::_1/(double)dim);
		mv_fitting_sample.push_back(rrates);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the relative rates of the associated partition model to those in the supplied vector `v'.
*/
void SubsetRelRatesMove::sendCurrValuesToModel(const double_vect_t & v)
	{
	PHYCAS_ASSERT(dim == v.size());
	partition_model->setSubsetRelRatesVect(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, storing them in the supplied vector `v'.
*/
void SubsetRelRatesMove::getCurrValuesFromModel(double_vect_t & v) const
	{
	const std::vector<double> & rrates = partition_model->getSubsetRelRatesVect();
	PHYCAS_ASSERT(dim == rrates.size());
	v.resize(rrates.size());
	std::copy(rrates.begin(), rrates.end(), v.begin());
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, returning them as an anonymous vector.
*/
double_vect_t SubsetRelRatesMove::listCurrValuesFromModel()
	{
	double_vect_t v(dim);
	const double_vect_t & rrates = partition_model->getSubsetRelRatesVect();
	PHYCAS_ASSERT(dim == rrates.size());
	std::copy(rrates.begin(), rrates.end(), v.begin());
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, storing them in the data member `orig_params'. Because subset
|	relative rates have mean 1 rather than sum 1, and since the base class DirichletMove expects the parameters to 
|	sum to 1, each value in `orig_params' is transformed by dividing by `dim' (the length of the `orig_params' vector).
*/
void SubsetRelRatesMove::getParams()
	{
	getCurrValuesFromModel(orig_params);

	//std::cerr << "orig_params before transforming:" << std::endl;
	//std::copy(orig_params.begin(), orig_params.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << std::endl;

	std::transform(orig_params.begin(), orig_params.end(), orig_params.begin(), boost::lambda::_1/(double)dim);

	//std::cerr << "orig_params after transforming by dividing by 2:" << std::endl;
	//std::copy(orig_params.begin(), orig_params.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the relative rates in the model with those supplied in the vector `v'. The values supplied are multiplied
|	by `dim' (the length of the vector) before sending to the model. See explanation in documentation for 
|	SubsetRelRatesMove::getParams().
*/
void SubsetRelRatesMove::setParams(
  const std::vector<double> & v)    /*< is the vector of parameter values to send to the model */
	{
	double_vect_t vcopy(dim);
	std::transform(v.begin(), v.end(), vcopy.begin(), boost::lambda::_1*(double)dim);

	//std::cerr << "v:" << std::endl;
	//std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << std::endl;
	//std::cerr << "vcopy (equals v*2):" << std::endl;
	//std::copy(vcopy.begin(), vcopy.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << std::endl;

    sendCurrValuesToModel(vcopy);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Sets `partition_model' to the supplied `m'.
*/
void SubsetRelRatesMove::setPartitionModel(
  PartitionModelShPtr m)    /*< is a shared pointer to the new partition model */
	{
	partition_model = m;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Retrieves current value for the parameter being managed by this updater, then returns log of the working prior 
|	probability density at that value. If this updater is not a prior steward, simply returns 0.0.
*/
double SubsetRelRatesMove::recalcWorkingPrior() const
	{
	if (isFixed() || !isPriorSteward())
		return 0.0;
		
	double lnwp = 0.0;
	double_vect_t values;
	getCurrValuesFromModel(values);	
	std::transform(values.begin(), values.end(), values.begin(), boost::lambda::_1/(double)dim);
	try 
		{
		PHYCAS_ASSERT(mv_working_prior);
		lnwp = mv_working_prior->GetLnPDF(values);
		}
	catch(XProbDist &)
		{
		PHYCAS_ASSERT(0);
		}
	double n = (double)dim;
	lnwp -= (n - 1.0)*log(n);
	return lnwp;
	}

