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

#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/gtr_model.hpp"
#include "phycas/src/codon_model.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/mcmc_param.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/basic_tree.hpp"

namespace phycas
{

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor and initializes its HKY pointer to NULL. Also sets the
|	`curr_value' data member to 4.0 and refreshes `curr_ln_prior' accordingly.
*/
EdgeLenParam::EdgeLenParam()
  : MCMCUpdater(), my_node(NULL)
	{
	curr_value = 0.01;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor.
*/
EdgeLenParam::~EdgeLenParam() 
	{
	//std::cerr << "\n>>>>> EdgeLenParam dying..." << std::endl;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool EdgeLenParam::update()
	{
	if (is_fixed)
		return false;
		
	if (my_node->IsInternal())
		likelihood->useAsLikelihoodRoot(my_node);
	else 
		likelihood->useAsLikelihoodRoot(my_node->GetParent());
	likelihood->invalidateAwayFromNode(*my_node);

	double current_brlen = getCurrValueFromModel();
	//double last_sampled_brlen = slice_sampler->GetLastSampledXValue();
	//if (fabs(current_brlen - last_sampled_brlen) > 0.00001)
	//	{
	//	std::cerr << "BAD BAD! BAD!!" << std::endl;
	//	}
	slice_sampler->SetXValue(current_brlen);
	slice_sampler->Sample();
	
	ChainManagerShPtr p = chain_mgr.lock();
	
#if defined(SAMC_TWO)
	if (p->doingSAMC())
		{
		double logf = slice_sampler->GetLastSampledYValue();
		unsigned i = p->getSAMCEnergyLevel(logf);
		p->updateSAMCWeights(i);
		}
#endif
	
    if (save_debug_info)
        {
        debug_info = str(boost::format("EdgeLenParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Override of base class version adds the current edge length to the data already stored in `fitting_sample'.
*/
void EdgeLenParam::educateWorkingPrior()
	{
	if (!isFixed())
		{
		PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
		double edgelen = getCurrValueFromModel();
		fitting_sample.push_back(edgelen);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use samples in `fitting_sample' to parameterize `working_prior'. This function is called during s-cubed style
|	steppingstone sampling after the initial phase of sampling from the posterior so that the working prior can be
|	used for the remaining phases. Assumes `fitting_sample' has more than 1 element. Assigns a GammaDistribution object
|	to `working_prior'.
*/
void EdgeLenParam::finalizeWorkingPrior()
	{
	PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
	fitGammaWorkingPrior();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the edge length to `v'.
*/
void EdgeLenParam::sendCurrValueToModel(double v)
	{
	PHYCAS_ASSERT(my_node != NULL);
	my_node->SetEdgeLen(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the edge length currently stored in the tree.
*/
double EdgeLenParam::getCurrValueFromModel() const
	{
	PHYCAS_ASSERT(my_node != NULL);
	return my_node->GetEdgeLen();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets private data member `my_node' to point to the supplied TreeNode `nd'.
*/
void EdgeLenParam::setTreeNode(TreeNode & nd)
	{
	my_node = &nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	EdgeLenParam is a functor whose operator() returns a value proportional to the full-conditional posterior probability 
|	density for a particular value of the edge length managed by this object. If the supplied kappa value `v' is
|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative
|	infinity).
*/
double EdgeLenParam::operator()(
  double v)	/**< is a new value for the edge length parameter */
	{
	PHYCAS_ASSERT(my_node != NULL);
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (v > 0.0)
		{
		sendCurrValueToModel(v);
		recalcPrior();

		// specifying NULL as likelihood root invalidates all CLAs
		// but here we want to specify my_node in order to minimize recalculation of CLAs
		//likelihood->storeAllCLAs(tree);
		//if (my_node->IsInternal())
		//	likelihood->useAsLikelihoodRoot(my_node);
		//else 
		//	likelihood->useAsLikelihoodRoot(my_node->GetParent());
		//likelihood->invalidateAwayFromNode(*my_node);
		//likelihood->invalidateBothEnds(my_node);
		//likelihood->useAsLikelihoodRoot(NULL);//temp
		
		//std::cerr << "@@@@@@@@@@ Updating edge length for " << my_node->GetNodeNumber() << std::endl;
			
		//likelihood->startTreeViewer(tree, boost::str(boost::format("Before calcLnL: new edge length = %.6f") % v));
		
		curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
		
		//likelihood->startTreeViewer(tree, boost::str(boost::format("After calcLnL: curr_ln_like = %.6f") % curr_ln_like));
		
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

#if defined(SAMC_TWO)
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
					double curr_ln_working_prior = working_prior->GetLnPDF(v);
					return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
					}
				else 
					return heating_power*(curr_ln_like + curr_ln_prior);
			else
				return heating_power*curr_ln_like + curr_ln_prior;
			}
#else	//not SAMC_TWO
		if (is_standard_heating)
			if (use_working_prior)
				{
				PHYCAS_ASSERT(working_prior);
				double curr_ln_working_prior = working_prior->GetLnPDF(v);
				return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
				}
			else 
				return heating_power*(curr_ln_like + curr_ln_prior);
		else
			return heating_power*curr_ln_like + curr_ln_prior;
#endif	//SAMC_TWO
		}
    else
        return ln_zero;
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	This override of the base class virtual function uses dynamic_cast to set the `gtr' data member from the supplied
|	model shared pointer. Having a GTR pointer allows this parameter to call GTR-specific member functions not present 
|	in the base class Model (such as setRelRates).
*/
void GTRRateParam::setModel(ModelShPtr m)
	{
	MCMCUpdater::setModel(m);
	Model * p = m.get();
	gtr = dynamic_cast<GTR *>(p);	// forces inclusion of "likelihood_models.hpp"

	// If tempted to move this to mcmc_param.inl, see comment in KappaParam::setModel function
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the appropriate relative rate in `model' to the value `v'.
*/
void GTRRateParam::sendCurrValueToModel(double v)
	{
	gtr->setRelRateUnnorm(which, v);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the appropriate relative rate value in `model'.
*/
double GTRRateParam::getCurrValueFromModel() const
	{
    return gtr->getRelRateUnnorm(which);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	GTRRateParam is a functor whose operator() returns a value proportional to the full-conditional posterior 
|	probability density for a particular value of one of the six relative rates. If the supplied rate value `r' is
|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative
|	infinity).
*/
double GTRRateParam::operator()(
  double r)	/**< is a new value for the relative rate parameter */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

    if (r > 0.0)
		{
		PHYCAS_ASSERT(which < 6);
		PHYCAS_ASSERT(gtr);
		sendCurrValueToModel(r);
		recalcPrior();
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
        curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

#if defined(SAMC_TWO)
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
					double curr_ln_working_prior = working_prior->GetLnPDF(r);
					return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
					}
				else 
					return heating_power*(curr_ln_like + curr_ln_prior);
			else
				return heating_power*curr_ln_like + curr_ln_prior;
			}				
#else	//not SAMC_TWO
		if (is_standard_heating)
			if (use_working_prior)
				{
				PHYCAS_ASSERT(working_prior);
				double curr_ln_working_prior = working_prior->GetLnPDF(r);
				return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
				}
			else 
				return heating_power*(curr_ln_like + curr_ln_prior);
		else
			return heating_power*curr_ln_like + curr_ln_prior;
#endif	//SAMC_
		}
    else
        return ln_zero;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This override of the base class virtual function uses dynamic_cast to set the `codon' data member from the supplied
|	model shared pointer. Having a Codon pointer allows this parameter to call Codon-specific member functions that are
|	not present in the base class Model (such as setOmega).
*/
void OmegaParam::setModel(ModelShPtr m)
	{
	MCMCUpdater::setModel(m);
	Model * p = m.get();
	codon = dynamic_cast<Codon *>(p);	// forces inclusion of "phycas/src/likelihood_models.hpp"

	//POL tried unsuccessfully to get this to compile as an inlined function, but VC gave me this 
	// error (which makes sense):
	//
	// ...mcmc_param.inl(43) : error C2680: 'phycas::Codon *' : 
	//  invalid target type for dynamic_cast
	//  'Codon' : class must be defined before using in a dynamic_cast
	//
	// Decided to add this comment because otherwise I will forget this and be tempted to move the 
	// function back to mcmc_param.inl
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the value of `omega' in `model' to `v'.
*/
void OmegaParam::sendCurrValueToModel(double v)
	{
	codon->setOmega(v);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the value of `omega' in `model'.
*/
double OmegaParam::getCurrValueFromModel() const
	{
    return codon->getOmega();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	OmegaParam is a functor whose operator() returns a value proportional to the full-conditional posterior probability 
|	density for a particular value of omega, the nonsynonymous/synonymous rate ratio. If the supplied omega value `w' is
|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative
|	infinity).
*/
double OmegaParam::operator()(
  double w)	/**< is a new value for the parameter omega */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (w > 0.0)
		{
		PHYCAS_ASSERT(codon);
		sendCurrValueToModel(w);
		recalcPrior();

		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

#if defined(SAMC_TWO)
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
					double curr_ln_working_prior = working_prior->GetLnPDF(w);
					return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
					}
				else 
					return heating_power*(curr_ln_like + curr_ln_prior);
			else
				return heating_power*curr_ln_like + curr_ln_prior;
			}
#else	//not SAMC_TWO
		if (is_standard_heating)
			if (use_working_prior)
				{
				PHYCAS_ASSERT(working_prior);
				double curr_ln_working_prior = working_prior->GetLnPDF(w);
				return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
				}
			else 
				return heating_power*(curr_ln_like + curr_ln_prior);
		else
			return heating_power*curr_ln_like + curr_ln_prior;
#endif	//SAMC_TWO
		}
    else
        return ln_zero;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the value of `pinvar' in `model' to `v'.
*/
void PinvarParam::sendCurrValueToModel(double v)
	{
	model->setPinvar(v);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the value of `pinvar' in `model'.
*/
double PinvarParam::getCurrValueFromModel() const
	{
    return model->getPinvar();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	PinvarParam is a functor whose operator() returns a value proportional to the full-conditional posterior
|	probability density for a particular value of the proportion of invariable sites parameter. If the supplied 
|	proportion of invariable sites `pinv' is out of bounds (i.e. < 0.0 or >= 1.0), the return value is `ln_zero' 
|	(closest we can come to a log posterior equal to negative infinity).
*/
double PinvarParam::operator()(
  double pinv)	/**< is a new value for the proportion of invariable sites parameter */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (pinv >= 0.0 && pinv < 1.0)
		{
		sendCurrValueToModel(pinv);
		recalcPrior(); // base class function that recomputes curr_ln_prior for the value curr_value
		likelihood->recalcRelativeRates();	// must do this whenever model's rate heterogeneity status changes
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

#if defined(SAMC_TWO)
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
					double curr_ln_working_prior = working_prior->GetLnPDF(pinv);
					return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
					}
				else 
					return heating_power*(curr_ln_like + curr_ln_prior);
			else
				return heating_power*curr_ln_like + curr_ln_prior;
			}
#else	//not SAMC_TWO
		if (is_standard_heating)
			if (use_working_prior)
				{
				PHYCAS_ASSERT(working_prior);
				double curr_ln_working_prior = working_prior->GetLnPDF(pinv);
				return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
				}
			else 
				return heating_power*(curr_ln_like + curr_ln_prior);
		else
			return heating_power*curr_ln_like + curr_ln_prior;
#endif	//SAMC_TWO
		}
    else
        return ln_zero;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the appropriate unnormalized state frequency in `model' to `v'.
*/
void StateFreqParam::sendCurrValueToModel(double v)
	{
	model->setStateFreqUnnorm(which, v);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the appropriate unnormalized state frequency in `model'.
*/
double StateFreqParam::getCurrValueFromModel() const
	{
    return model->getStateFreqUnnorm(which);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	StateFreqParam is a functor whose operator() returns a value proportional to the full-conditional posterior 
|	probability density for a particular value of f, the quantity governing a particular relative base frequency (but 
|	f is not a relative base frequency itself, as it is not constrained to lie between 0.0 and 1.0, nor is it 
|	constrained by any of the analogous quantities governing the other three base frequencies. If the supplied base 
|	frequency parameter value `f' is out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to 
|	a log posterior equal to negative infinity).
*/
double StateFreqParam::operator()(
  double f) /**< is a new value for the base frequency parameter */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

    if (f > 0.0)
		{
		PHYCAS_ASSERT((model->isCodonModel() && which < 61) || which < 4);
		//state_freq.freqs[which] = f;
        //model.normalizeFreqs(state_freq.freqs);
		sendCurrValueToModel(f);
		recalcPrior();
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
        curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

#if defined(SAMC_TWO)
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
					double curr_ln_working_prior = working_prior->GetLnPDF(f);
					return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
					}
				else 
					return heating_power*(curr_ln_like + curr_ln_prior);
			else
				return heating_power*curr_ln_like + curr_ln_prior;
			}
#else	//not SAMC_TWO
		if (is_standard_heating)
			if (use_working_prior)
				{
				PHYCAS_ASSERT(working_prior);
				double curr_ln_working_prior = working_prior->GetLnPDF(f);
				return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
				}
			else 
				return heating_power*(curr_ln_like + curr_ln_prior);
		else
			return heating_power*curr_ln_like + curr_ln_prior;
#endif	//SAMC_TWO
		}
    else
        return ln_zero;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set value of hyperparameter in `model' to `v'.
*/
void HyperPriorParam::sendCurrValueToModel(double v)
	{
	if (external_edges)
		model->setExternalEdgelenHyperparam(v);
	else
		model->setInternalEdgelenHyperparam(v);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return value of hyperparameter from model.
*/
double HyperPriorParam::getCurrValueFromModel() const
	{
	if (external_edges)
		return model->getExternalEdgelenHyperparam();
	else 
		return model->getInternalEdgelenHyperparam();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	HyperPriorParam is a functor whose operator() returns a value proportional to the full-conditional posterior 
|	probability density for a particular value of mu, the mean of the edge length prior distribution. If the supplied 
|	value `mu' is out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior 
|	equal to negative infinity).
*/
double HyperPriorParam::operator()(
  double mu) 	/**< is a new parameter value for the prior affecting all edge lengths */
	{
	double edgeLensLnPrior = 0.0;
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (mu > 0.0)
		{
		sendCurrValueToModel(mu);
		recalcPrior();

		// no need to invalidate all CLAs because only the prior is changing
		//@POL the likelihood is not even needed here because mu does not appear in the likelihood function anywhere
		curr_ln_like = likelihood->calcLnL(tree);

		edgelen_master_param->setPriorMeanAndVariance(mu, mu*mu); //@POL note that this implicitly assumes an exponential edge length prior
        edgeLensLnPrior = edgelen_master_param->recalcPrior();

        // Store current log-likelihood value (if next updater is a move, it will then not have to calculate it)
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

#if defined(SAMC_TWO)
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
					double curr_ln_working_prior = working_prior->GetLnPDF(mu);
					return heating_power*(curr_ln_like + edgeLensLnPrior + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
					}
				else 
					return heating_power*(curr_ln_like + edgeLensLnPrior + curr_ln_prior);
			else
				return heating_power*curr_ln_like + edgeLensLnPrior + curr_ln_prior;
			}
#else	//not SAMC_TWO
		if (is_standard_heating)
			if (use_working_prior)
				{
				PHYCAS_ASSERT(working_prior);
				double curr_ln_working_prior = working_prior->GetLnPDF(mu);
				return heating_power*(curr_ln_like + edgeLensLnPrior + curr_ln_prior) + (1.0 - heating_power)*curr_ln_working_prior;
				}
			else 
				return heating_power*(curr_ln_like + edgeLensLnPrior + curr_ln_prior);
		else
			return heating_power*curr_ln_like + edgeLensLnPrior + curr_ln_prior;
#endif	//SAMC_TWO
		}
    else
        return ln_zero;
	}

}	// namespace phycas
