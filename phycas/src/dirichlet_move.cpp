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

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor sets tuning parameter `psi' and `maxpsi' to 300.0 and calls reset().
*/
DirichletMove::DirichletMove() : MCMCUpdater()
	{
	dim = 0;
	boldness = 0.0;
	minpsi = 1.0;
	maxpsi = 300.0;
	psi = maxpsi;
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
|	0 and maximum of `maxpsi' enforced after the calculation. Current parameter values are multiplied by the value `psi'
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

    // compute psi from boldness value
    //  
    //  minpsi = 5
    //  maxpsi = 300
    //  % represents result in column to the left
    //
    //  boldness     100-%    (maxpsi-minpsi)*%/100   minpsi+%  
    //      0         100              295               300
    //    100           0                0                 5
    //
    //  psi = minpsi + (maxpsi - minpsi)*(100-boldness)/100
    //
	psi = minpsi + (maxpsi - minpsi)*(100.0-boldness)/100.0;
	std::cerr << boost::str(boost::format("####### x = %.5f, boldness = %.5f, psi = %.5f") % x % boldness % psi) << std::endl;
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
	//std::cerr << boost::str(boost::format(">>>>>>>>>> boldness = %.1f, psi = %.5f\n") % boldness % psi);
	//std::cerr << ">>>>>>>>>> orig_params: ";
	//std::copy(orig_params.begin(), orig_params.end(), std::ostream_iterator<double>(std::cerr," "));
	//std::cerr << "\n";
	//std::cerr << ">>>>>>>>>> c_forward: ";
	//std::copy(c_forward.begin(), c_forward.end(), std::ostream_iterator<double>(std::cerr," "));
	//std::cerr << "\n";
	//std::cerr << ">>>>>>>>>> new_params: ";
	//std::copy(new_params.begin(), new_params.end(), std::ostream_iterator<double>(std::cerr," "));
	//std::cerr << "\n" << std::endl;

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
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Calls the reset() function after reinstating the original state frequencies
|   and ensuring that all conditional likelihood arrays will be recalculated when the likelihood is next calculated.
*/
void DirichletMove::revert()
	{
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
	if (ln_accept_ratio >= 0.0 || lnu <= ln_accept_ratio)
		{
        //std::cerr << "ACCEPTED\n" << std::endl;
		p->setLastLnPrior(curr_ln_prior);
		p->setLastLnLike(curr_ln_like);
		accept();
		return true;
		}
	else
		{
        //std::cerr << "rejected\n" << std::endl;
		curr_ln_like	= p->getLastLnLike();
		curr_ln_prior	= p->getLastLnPrior();
		revert();
		return false;
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
|	Obtains the current state frequencies from the model, storing them in the data member `orig_params'.
*/
void StateFreqMove::getParams()
	{
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
    	orig_params.assign(dim, 1.0);
    	}
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
|	Obtains the current relative rates from the model, storing them in the data member `orig_params'.
*/
void RelRatesMove::getParams()
	{
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
    	orig_params.assign(dim, 1.0);
    	}
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

