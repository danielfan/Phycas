#include "phycas/force_include.h"
#include "CipresCommlib/CipresDataMatrixHelper.h"
#include "pyphy/likelihood/likelihood_models.hpp"
#include "pyphy/likelihood/tree_likelihood.hpp"
#include "phycas/rand/probability_distribution.hpp"
#include "pyphy/likelihood/mcmc_param.hpp"
#include "pyphy/likelihood/mcmc_chain_manager.hpp"
#include "pyphy/phylogeny/basic_tree.hpp"

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
	hky = dynamic_cast<HKY *>(p);	// forces inclusion of "pyphy/likelihood/likelihood_models.hpp"

	//POL tried unsuccessfully to get this to compile as an inlined function, but VC gave me this 
	// error (which makes sense):
	//
	// C:\Synchronized\Projects\phycasdev\pyphy\likelihood\mcmc_param.inl(43) : error C2680: 'phycas::HKY *' : 
	//  invalid target type for dynamic_cast	//  'HKY' : class must be defined before using in a dynamic_cast	//	// Decided to add this comment because otherwise I will forget this and be tempted to move the 	// function back to mcmc_param.inl	}/*----------------------------------------------------------------------------------------------------------------------|	KappaParam is a functor whose operator() returns a value proportional to the full-conditional posterior probability |	density for a particular value of kappa, the transition/transversion rate ratio. If the supplied kappa value `k' is|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative|	infinity).*/double KappaParam::operator()(  double k)	/**< is a new value for the parameter kappa */	{	curr_ln_like = ln_zero;	curr_ln_prior = 0.0;	if (k > 0.0)
		{
		assert(hky);		hky->setKappa(k);		curr_value = k;		recalcPrior();
		curr_ln_like = likelihood->calcLnL(tree);		ChainManagerShPtr p = chain_mgr.lock();		assert(p);		p->setLastLnLike(curr_ln_like);		}	return curr_ln_like + curr_ln_prior;	}/*----------------------------------------------------------------------------------------------------------------------|	This override of the base class virtual function uses dynamic_cast to set the `gtr' data member from the supplied|	model shared pointer. Having a GTR pointer allows this parameter to call GTR-specific member functions not present |	in the base class Model (such as setRelRates).*/void GTRRateParam::setModel(ModelShPtr m)	{	MCMCUpdater::setModel(m);	Model * p = m.get();	gtr = dynamic_cast<GTR *>(p);	// forces inclusion of "pyphy/likelihood/likelihood_models.hpp"	// If tempted to move this to mcmc_param.inl, see comment in KappaParam::setModel function	}/*----------------------------------------------------------------------------------------------------------------------|	GTRRateParam is a functor whose operator() returns a value proportional to the full-conditional posterior |	probability density for a particular value of one of the six relative rates. If the supplied rate value `r' is|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative|	infinity).*/double GTRRateParam::operator()(  double r)	/**< is a new value for the relative rate parameter */	{	curr_ln_like = ln_zero;	curr_ln_prior = 0.0;    if (r > 0.0)
		{
		assert(which < 6);
		assert(gtr);
		gtr->setRelRateUnnorm(which, r);
        curr_value = r;		recalcPrior();
        curr_ln_like = likelihood->calcLnL(tree);		ChainManagerShPtr p = chain_mgr.lock();		assert(p);		p->setLastLnLike(curr_ln_like);		}	return curr_ln_like + curr_ln_prior;	}/*----------------------------------------------------------------------------------------------------------------------|	DiscreteGammaShapeParam is a functor whose operator() returns a value proportional to the full-conditional posterior|	probability density for a particular value of the gamma shape parameter. If the supplied gamma shape value `a' is|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative|	infinity).*/double DiscreteGammaShapeParam::operator()(  double a)	/**< is a new value for the gamma shape parameter */	{	curr_ln_like = ln_zero;	curr_ln_prior = 0.0;	if (a > 0.0)
		{
		curr_value = a;		recalcPrior(); // base class function that recomputes curr_ln_prior for the value curr_value
#if POLPY_NEWWAY
		if (invert_shape)
			model->setShape(1.0/a);	// change the gamma shape parameter in the model
		else
			model->setShape(a);	// change the gamma shape parameter in the model
#else
		model->setShape(a);	// change the gamma shape parameter in the model
#endif
		likelihood->recalcRelativeRates();	// must do this whenever model's shape parameter changes
		curr_ln_like = likelihood->calcLnL(tree);		ChainManagerShPtr p = chain_mgr.lock();		assert(p);		p->setLastLnLike(curr_ln_like);		}	return curr_ln_like + curr_ln_prior;	}#if POLPY_NEWWAY/*----------------------------------------------------------------------------------------------------------------------|	PatternSpecificRateParam is a functor whose operator() returns a value proportional to the full-conditional |	posterior probability density for a particular value of a pattern-specific rate parameter. If the supplied pattern|	specific rate value `r' is out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log |	posterior equal to negative infinity).*/double PatternSpecificRateParam::operator()(  double r)	/**< is a new value for the pattern-specific rate parameter */	{	curr_ln_like = ln_zero;	curr_ln_prior = 0.0;	if (r > 0.0)
		{
		curr_value = r;		recalcPrior(); // base class function that recomputes curr_ln_prior for the value curr_value
		model->setPatternSpecificRate(r, which_pattern);	// change the pattern-specific rate parameter in the model
		curr_ln_like = likelihood->calcLnL(tree); //@POL only *need* to recalculate likelihood for one pattern, but currently the setPatternSpecificRate function renormalizes all the rates, which forces us to recalculate everything (any way we can avoid having to normalize while updating one of these parameters?)		ChainManagerShPtr p = chain_mgr.lock();		assert(p);		p->setLastLnLike(curr_ln_like);		}	return curr_ln_like + curr_ln_prior;	}#endif/*----------------------------------------------------------------------------------------------------------------------|	PinvarParam is a functor whose operator() returns a value proportional to the full-conditional posterior|	probability density for a particular value of the proportion of invariable sites parameter. If the supplied |	proportion of invariable sites `pinv' is out of bounds (i.e. < 0.0 or >= 1.0), the return value is `ln_zero' |	(closest we can come to a log posterior equal to negative infinity).*/double PinvarParam::operator()(  double pinv)	/**< is a new value for the proportion of invariable sites parameter */	{	curr_ln_like = ln_zero;	curr_ln_prior = 0.0;	if (pinv >= 0.0 && pinv < 1.0)
		{
		curr_value = pinv;		recalcPrior(); // base class function that recomputes curr_ln_prior for the value curr_value
		model->setPinvar(pinv);	// change the proportion of invariable sites parameter in the model
		likelihood->recalcRelativeRates();	// must do this whenever model's rate heterogeneity status changes
		curr_ln_like = likelihood->calcLnL(tree);		ChainManagerShPtr p = chain_mgr.lock();		assert(p);		p->setLastLnLike(curr_ln_like);		}	return curr_ln_like + curr_ln_prior;	}/*----------------------------------------------------------------------------------------------------------------------|	BaseFreqParam is a functor whose operator() returns a value proportional to the full-conditional posterior |	probability density for a particular value of f, the quantity governing a particular relative base frequency (but |	f is not a relative base frequency itself, as it is not constrained to lie between 0.0 and 1.0, nor is it |	constrained by any of the analogous quantities governing the other three base frequencies. If the supplied base |	frequency parameter value `f' is out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to |	a log posterior equal to negative infinity).*/double BaseFreqParam::operator()(  double f) /**< is a new value for the base frequency parameter */	{	curr_ln_like = ln_zero;	curr_ln_prior = 0.0;    if (f > 0.0)
		{
		assert(which < 4);
		model->setStateFreqUnnorm(which, f);
		//base_freq.freqs[which] = f;
        //model.normalizeFreqs(base_freq.freqs);        curr_value = f;		recalcPrior();
        curr_ln_like = likelihood->calcLnL(tree);		ChainManagerShPtr p = chain_mgr.lock();		assert(p);		p->setLastLnLike(curr_ln_like);		}	return curr_ln_like + curr_ln_prior;	}/*----------------------------------------------------------------------------------------------------------------------|	HyperPriorParam is a functor whose operator() returns a value proportional to the full-conditional posterior |	probability density for a particular value of mu, the mean of the edge length prior distribution. If the supplied |	value `mu' is out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior |	equal to negative infinity).*/double HyperPriorParam::operator()(  double mu) 	/**< is a new parameter value for the prior affecting all edge lengths */	{	double edgeLensLnPrior = 0.0;	curr_ln_like = ln_zero;	curr_ln_prior = 0.0;	if (mu > 0.0)		{		curr_ln_like = likelihood->calcLnL(tree);		curr_value = mu;		recalcPrior();
		// MCMCChainManager::getLnPrior() just returns the last prior value calculated		// by each edge length parameter, so we need to recalculate all of these using 		// the new edge length prior parameter value mu so that the call to 		// chain_mgr->calcJointLnPrior() will return the correct joint prior		ChainManagerShPtr p = chain_mgr.lock();		assert(p);		//std::cerr << "  calling MCMCChainManager::recalcEdgeLenPriors" << std::endl;		edgeLensLnPrior = p->recalcEdgeLenPriors(mu, mu*mu);	//@POL implicit assumption that edge length prior is exponential		p->setLastLnLike(curr_ln_like);		}	return curr_ln_like + edgeLensLnPrior + curr_ln_prior;	}/*----------------------------------------------------------------------------------------------------------------------|	EdgeLenParam is a functor whose operator() returns a value proportional to the full-conditional posterior |	probability density for a particular value `x', which is an edge length for the node managed by this EdgeLenParam.|	If the supplied value `x' is out of bounds (i.e. <= 0.0), the return value is `ln_zero' (closest we can come to a|	log posterior equal to negative infinity).*/double EdgeLenParam::operator()(  double x) /**< is a new edge length for the node being managed by this object */	{	curr_ln_like = ln_zero;	curr_ln_prior = 0.0;	if (x > 0.0)
		{
		if (nd->IsRoot())
			{
			assert(nd->GetLeftChild());
			nd->GetLeftChild()->SetEdgeLen(x);
			}
		else
			{
			nd->SetEdgeLen(x);
			}
		curr_ln_like = likelihood->calcLnL(tree);		curr_value = x;		recalcPrior();
		ChainManagerShPtr p = chain_mgr.lock();		assert(p);		p->setLastLnLike(curr_ln_like);		}	return curr_ln_like + curr_ln_prior;	}}	// namespace phycas