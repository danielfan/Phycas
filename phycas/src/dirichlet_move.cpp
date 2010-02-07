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
	
	// create Dirichlet distribution and sample from it to get proposed frequencies
    dir_forward = DirichletShPtr(new DirichletDistribution(c_forward));
    dir_forward->SetLot(getLot().get());
    new_params = dir_forward->Sample();

	//temporary!!
    //     std::cerr << boost::str(boost::format(">>>>>>>>>> boldness = %.1f, psi = %.5f\n") % boldness % psi);
    //     std::cerr << ">>>>>>>>>> orig_params: ";
    //     std::copy(orig_params.begin(), orig_params.end(), std::ostream_iterator<double>(std::cerr," "));
    //     std::cerr << "\n";
    //     std::cerr << ">>>>>>>>>> c_forward: ";
    //     std::copy(c_forward.begin(), c_forward.end(), std::ostream_iterator<double>(std::cerr," "));
    //     std::cerr << "\n";
    //     std::cerr << ">>>>>>>>>> new_params: ";
    //     std::copy(new_params.begin(), new_params.end(), std::ostream_iterator<double>(std::cerr," "));
    //     std::cerr << "\n" << std::endl;

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
	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);

	proposeNewState();

	// Note that it is not really correct to use a Dirichlet prior here because transformed versions of the
	// orig_params values are used in the likelihood function. However, the transformation only adds a constant
	// (i.e., -log(nsubsets)) to the log prior and thus this constant cancels out in the ln_accept_ratio.
	// More care should be taken when reporting the log prior or when estimating marginal likelihoods.
    double prev_ln_prior		= mvprior->GetRelativeLnPDF(orig_params);
	double prev_ln_like			= p->getLastLnLike();

    // replace current parameter values with new ones
    setParams(new_params);
    likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs

	double curr_ln_prior		= mvprior->GetRelativeLnPDF(new_params);
	double curr_ln_like			= likelihood->calcLnL(tree);

    double prev_posterior = 0.0;
	double curr_posterior = 0.0;
    if (is_standard_heating)
        {
        prev_posterior = heating_power*(prev_ln_like + prev_ln_prior);
	    curr_posterior = heating_power*(curr_ln_like + curr_ln_prior);
        }
    else
        {
        prev_posterior = heating_power*prev_ln_like + prev_ln_prior;
	    curr_posterior = heating_power*curr_ln_like + curr_ln_prior;
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
			debug_info = boost::str(boost::format("ACCEPT, prev_ln_prior = %.5f, curr_ln_prior = %.5f, prev_ln_like = %.5f, curr_ln_like = %.5f, lnu = %.5f, ln_accept_ratio = %.5f\norig_params: ") % prev_ln_prior % curr_ln_prior % prev_ln_like % curr_ln_like % lnu % ln_accept_ratio);
			for (std::vector<double>::const_iterator it = orig_params.begin(); it != orig_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\nnew_params:  ";
			for (std::vector<double>::const_iterator it = new_params.begin(); it != new_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			}
		p->setLastLnPrior(curr_ln_prior);
		p->setLastLnLike(curr_ln_like);
		accept();
		return true;
		}
	else
		{
	    if (save_debug_info)
    	    {
			debug_info = boost::str(boost::format("REJECT, prev_ln_prior = %.5f, curr_ln_prior = %.5f, prev_ln_like = %.5f, curr_ln_like = %.5f, lnu = %.5f, ln_accept_ratio = %.5f\norig_params: ") % prev_ln_prior % curr_ln_prior % prev_ln_like % curr_ln_like % lnu % ln_accept_ratio);
			for (std::vector<double>::const_iterator it = orig_params.begin(); it != orig_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			debug_info += "\nnew_params:  ";
			for (std::vector<double>::const_iterator it = new_params.begin(); it != new_params.end(); ++it)
				debug_info += boost::str(boost::format("%15.8f ") % (*it));
			}
		curr_ln_like	= p->getLastLnLike();
		curr_ln_prior	= p->getLastLnPrior();
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
|	The constructor simply calls the base class (DirichletMove) constructor.
*/
StateFreqMove::StateFreqMove() : DirichletMove()
	{
	dim = 4;
	}

#if POLPY_NEWWAY
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
void StateFreqMove::getCurrValuesFromModel(double_vect_t & v)
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
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current state frequencies from the model, storing them in the data member `orig_params'.
*/
void StateFreqMove::getParams()
	{
#if POLPY_NEWWAY
	getCurrValuesFromModel(orig_params);
#else //old way
	PHYCAS_ASSERT(dim > 0);
	if (model)
		{
    	const std::vector<double> & rfreqs = model->getStateFreqs();
    	orig_params.resize(rfreqs.size());
		PHYCAS_ASSERT(dim == rfreqs.size());
    	std::copy(rfreqs.begin(), rfreqs.end(), orig_params.begin());
    	}
    else
    	{
    	orig_params.assign(dim, 1.0/(double)dim);
    	}
#endif
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

#if POLPY_NEWWAY
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
void RelRatesMove::getCurrValuesFromModel(double_vect_t & v)
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
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, storing them in the data member `orig_params'.
*/
void RelRatesMove::getParams()
	{
#if POLPY_NEWWAY
	getCurrValuesFromModel(orig_params);
#else //old way
	PHYCAS_ASSERT(dim > 0);
	GTR * gtr_model = dynamic_cast<GTR *>(model.get());
	if (gtr_model)
		{
		const std::vector<double> & rrates = gtr_model->getRelRates();
		PHYCAS_ASSERT(dim == rrates.size());
		orig_params.resize(rrates.size());
		std::copy(rrates.begin(), rrates.end(), orig_params.begin());
    	}
    else
    	{
    	orig_params.assign(dim, 1.0/(double)dim);
    	}
#endif
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

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	The constructor simply calls the base class (DirichletMove) constructor.
*/
SubsetRelRatesMove::SubsetRelRatesMove() : DirichletMove()
	{
	dim = 1;
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
void SubsetRelRatesMove::getCurrValuesFromModel(double_vect_t & v)
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
|	Obtains the current relative rates from the model, storing them in the data member `orig_params'.
*/
void SubsetRelRatesMove::getParams()
	{
	getCurrValuesFromModel(orig_params);
	
	// need to modify the values obtained from the model because the model maintains relative rates (mean 1.0)
	// but for Dirichlet move we need to work with parameters that sum to 1.0
	std::transform(orig_params.begin(), orig_params.end(), orig_params.begin(), boost::lambda::_1/(double)dim);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the relative rates in the model with those supplied in the vector `v'.
*/
void SubsetRelRatesMove::setParams(
  const std::vector<double> & v)    /*< is the vector of parameter values to send to the model */
	{
	// need to modify the values in v because the model maintains relative rates (mean 1.0)
	// but for Dirichlet move we work with parameters that sum to 1.0
	double_vect_t vcopy(v.size());
	std::transform(v.begin(), v.end(), vcopy.begin(), boost::lambda::_1*(double)dim);
	
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

#endif

