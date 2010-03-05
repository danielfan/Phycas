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

#include "phycas/src/mcmc_param.hpp"
#include "phycas/src/codon_model.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor and initializes its HKY pointer to NULL. Also sets the
|	`curr_value' data member to 4.0 and refreshes `curr_ln_prior' accordingly.
*/
KappaParam::KappaParam()
  : MCMCUpdater(), hky(NULL), codon(NULL)
	{
	curr_value = 4.0;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor.
*/
KappaParam::~KappaParam() 
	{
	//std::cerr << "\n>>>>> KappaParam dying..." << std::endl;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool KappaParam::update()
	{
	//@POL should probably put next two lines in Model::update and have Model::update call a virtual Model::updateImpl
	// because it will be hard to remember to put these lines in every overloaded update function
	if (is_fixed)
		return false;
	slice_sampler->Sample();

	ChainManagerShPtr p = chain_mgr.lock();
	if (p->doingSAMC())
		{
		double logf = slice_sampler->GetLastSampledYValue();
		unsigned i = p->getSAMCEnergyLevel(logf);
		p->updateSAMCWeights(i);
		}
	
    if (save_debug_info)
        {
        debug_info = str(boost::format("KappaParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class version adds the current kappa parameter value to the data already stored in 
|	`fitting_sample'.
*/
void KappaParam::educateWorkingPrior()
	{
	if (!isFixed())
		{
		PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
		double kappa = getCurrValueFromModel();
		fitting_sample.push_back(kappa);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use samples in `fitting_sample' to parameterize `working_prior'. This function is called during s-cubed style
|	steppingstone sampling after the initial phase of sampling from the posterior so that the working prior can be
|	used for the remaining phases. Assumes `fitting_sample' has more than 1 element. Assigns a GammaDistribution object
|	to `working_prior'.
*/
void KappaParam::finalizeWorkingPrior()
	{
	PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
	fitGammaWorkingPrior();
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the kappa value in either `hky' or `codon' (whichever is relevant) to `v'.
*/
void KappaParam::sendCurrValueToModel(double v)
	{
	if (hky != NULL)
		hky->setKappa(v);
	else
		codon->setKappa(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the kappa value in either `hky' or `codon' (whichever is relevant).
*/
double KappaParam::getCurrValueFromModel() const
	{
	if (hky != NULL)
        return hky->getKappa();
    else
        return codon->getKappa();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	KappaParam is a functor whose operator() returns a value proportional to the full-conditional posterior probability 
|	density for a particular value of kappa, the transition/transversion rate ratio. If the supplied kappa value `k' is
|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative
|	infinity).
*/
double KappaParam::operator()(
  double k)	/**< is a new value for the parameter kappa */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (k > 0.0)
		{
    	
		sendCurrValueToModel(k);
		recalcPrior();

		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

		if (p->doingSAMC())
			{
			double log_posterior = curr_ln_like + curr_ln_prior;
			unsigned i = p->getSAMCEnergyLevel(log_posterior);
			double curr_theta = p->getSAMCWeight(i);
			return (log_posterior - curr_theta);
			}
		else 
			{
			if (is_standard_heating)
				if (use_working_prior)
					{
					PHYCAS_ASSERT(working_prior);
					double curr_ln_working_prior = working_prior->GetLnPDF(k);
					return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
					}
				else 
					return heating_power*(curr_ln_like + curr_ln_prior);
			else
				return heating_power*curr_ln_like + curr_ln_prior;
			}
		}
    else
        return ln_zero;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	This override of the base class virtual function uses dynamic_cast to set the `hky' data member from the supplied
|	model shared pointer. Having an HKY pointer allows this parameter to call HKY-specific member functions that are not
|	present in the base class Model (such as setKappa).
*/
void KappaParam::setModel(ModelShPtr m)
	{
	MCMCUpdater::setModel(m);
	Model * p = m.get();
	if (m->isCodonModel())
    	codon = dynamic_cast<Codon *>(p);	// forces inclusion of "phycas/src/likelihood_models.hpp"
    else
    	hky = dynamic_cast<HKY *>(p);	// forces inclusion of "phycas/src/likelihood_models.hpp"

	//POL tried unsuccessfully to get this to compile as an inlined function, but VC gave me this 
	// error (which makes sense):
	//
	// ...mcmc_param.inl(43) : error C2680: 'phycas::HKY *' : 
	//  invalid target type for dynamic_cast
	//  'HKY' : class must be defined before using in a dynamic_cast
	//
	// Decided to add this comment because otherwise I will forget this and be tempted to move the 
	// function back to mcmc_param.inl
	}

}
