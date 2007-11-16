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

//#include "phycas/force_include.h"
#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/mcmc_param.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/basic_tree.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	This override of the base class virtual function uses dynamic_cast to set the `hky' data member from the supplied
|	model shared pointer. Having an HKY pointer allows this parameter to call HKY-specific member functions that are not
|	present in the base class Model (such as setKappa).
*/
void KappaParam::setModel(ModelShPtr m)
	{
	MCMCUpdater::setModel(m);
	Model * p = m.get();
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

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set `curr_value' to the kappa value in `model'.
*/
void KappaParam::setCurrValueFromModel()
	{
    curr_value = hky->getKappa();
	}
#endif

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
		PHYCAS_ASSERT(hky);
		hky->setKappa(k);
		curr_value = k;
		recalcPrior();

		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = likelihood->calcLnL(tree);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        if (is_standard_heating)
            return heating_power*(curr_ln_like + curr_ln_prior);
        else
            return heating_power*curr_ln_like + curr_ln_prior;
		}
    else
        return ln_zero;
	}

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

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set `curr_value' to the corresponding relative rate value in `model'.
*/
void GTRRateParam::setCurrValueFromModel()
	{
    curr_value = gtr->getRelRateUnnorm(which);
	}
#endif

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
		gtr->setRelRateUnnorm(which, r);
        curr_value = r;
		recalcPrior();
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
        curr_ln_like = likelihood->calcLnL(tree);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        if (is_standard_heating)
            return heating_power*(curr_ln_like + curr_ln_prior);
        else
            return heating_power*curr_ln_like + curr_ln_prior;
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

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set `curr_value' to the value of `omega' in `model'.
*/
void OmegaParam::setCurrValueFromModel()
	{
    curr_value = codon->getOmega();
	}
#endif

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
		codon->setOmega(w);
		curr_value = w;
		recalcPrior();

		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = likelihood->calcLnL(tree);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        if (is_standard_heating)
            return heating_power*(curr_ln_like + curr_ln_prior);
        else
            return heating_power*curr_ln_like + curr_ln_prior;
		}
    else
        return ln_zero;
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set `curr_value' to the value of `gamma_shape' in `model'.
*/
void DiscreteGammaShapeParam::setCurrValueFromModel()
	{
    curr_value = model->getShape();
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	DiscreteGammaShapeParam is a functor whose operator() returns a value proportional to the full-conditional posterior
|	probability density for a particular value of the gamma shape parameter. If the supplied gamma shape value `a' is
|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative
|	infinity).
*/
double DiscreteGammaShapeParam::operator()(
  double a)	/**< is a new value for the gamma shape parameter */
	{
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (a > 0.0)
		{
		curr_value = a;
		recalcPrior(); // base class function that recomputes curr_ln_prior for the value curr_value
		if (shape_inverted)
			model->setShape(1.0/a);	// change the gamma shape parameter in the model
		else
			model->setShape(a);	// change the gamma shape parameter in the model
		likelihood->recalcRelativeRates();	// must do this whenever model's shape parameter changes
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = likelihood->calcLnL(tree);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        if (is_standard_heating)
            return heating_power*(curr_ln_like + curr_ln_prior);
        else
            return heating_power*curr_ln_like + curr_ln_prior;
		}
    else
        return ln_zero;
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set `curr_value' to the value of `pinvar' in `model'.
*/
void PinvarParam::setCurrValueFromModel()
	{
    curr_value = model->getPinvar();
	}
#endif

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
		curr_value = pinv;
		recalcPrior(); // base class function that recomputes curr_ln_prior for the value curr_value
		model->setPinvar(pinv);	// change the proportion of invariable sites parameter in the model
		likelihood->recalcRelativeRates();	// must do this whenever model's rate heterogeneity status changes
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
		curr_ln_like = likelihood->calcLnL(tree);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        if (is_standard_heating)
            return heating_power*(curr_ln_like + curr_ln_prior);
        else
            return heating_power*curr_ln_like + curr_ln_prior;
		}
    else
        return ln_zero;
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set `curr_value' to the corresponding value of the `pinvar' in `model'.
*/
void StateFreqParam::setCurrValueFromModel()
	{
    curr_value = model->getStateFreqUnnorm(which);
	}
#endif

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
		PHYCAS_ASSERT(which < 4);
		model->setStateFreqUnnorm(which, f);
		//state_freq.freqs[which] = f;
        //model.normalizeFreqs(state_freq.freqs);
        curr_value = f;
		recalcPrior();
		likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
        curr_ln_like = likelihood->calcLnL(tree);
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        if (is_standard_heating)
            return heating_power*(curr_ln_like + curr_ln_prior);
        else
            return heating_power*curr_ln_like + curr_ln_prior;
		}
    else
        return ln_zero;
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set `curr_value' to the corresponding value of the `pinvar' in `model'.
*/
void HyperPriorParam::setCurrValueFromModel()
	{
    //@POL hyperparameters should be part of the model, but currently they are not
	}
#endif

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
		// no need to invalidate all CLAs because only the prior is changing
		curr_ln_like = likelihood->calcLnL(tree);
		curr_value = mu;
		recalcPrior();

		edgelen_master_param->setPriorMeanAndVariance(mu, mu*mu); //@POL note that this implicitly assumes an exponential edge length prior
        edgeLensLnPrior = edgelen_master_param->recalcPrior();

        // Store current log-likelihood value (if next updater is a move, it will then not have to calculate it)
		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        if (is_standard_heating)
            return heating_power*(curr_ln_like + edgeLensLnPrior + curr_ln_prior);
        else
            return heating_power*curr_ln_like + edgeLensLnPrior + curr_ln_prior;
		}
    else
        return ln_zero;
	}

}	// namespace phycas
