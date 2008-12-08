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

#if ! defined(MCMC_UPDATER_INL)
#define MCMC_UPDATER_INL

#include <limits>									// for std::numeric_limits
#include "phycas/src/slice_sampler.hpp"	// see MCMCUpdater::createSliceSampler
#include "phycas/src/probability_distribution.hpp"	// see MCMCUpdater::setPriorMeanAndVariance

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor for the base class of all move and parameter classes used for MCMC. MCMCUpdater-derived objects manage 
|	either Metropolis-Hastings proposals (moves) or model quantities updated via slice sampling (parameters) in a 
|	Bayesian analysis.
*/
inline MCMCUpdater::MCMCUpdater()
  : curr_value(0.1), 
  curr_ln_prior(0.0), 
  curr_ln_like(0.0), 
  has_slice_sampler(false),
  is_move(false),
  is_master_param(false),
  is_hyper_param(false), 
  is_fixed(false),
  slice_max_units(UINT_MAX),
  heating_power(1.0),
  is_standard_heating(true),
  save_debug_info(false)
	{
	//ln_zero = std::log(std::numeric_limits<double>::denorm_min()); // this doesn't work, lnL can get much lower than the log of dnorm_min!
	ln_zero = -DBL_MAX;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `is_fixed'. While `is_fixed' is true, then the update function always returns 
|	immediately without ever modifying the parameter value.
*/
inline bool MCMCUpdater::isFixed() const
	{
	return is_fixed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Fixes the value of this parameter to the current value. Subsequent calls to the update member function will have no
|	effect until freeParameter is called. Sets the value of the data member `is_fixed' to true.
*/
inline void MCMCUpdater::fixParameter()
	{
	is_fixed = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `is_fixed' to the default value (false). Subsequent calls to the update member
|	function will update the value of the parameter.
*/
inline void MCMCUpdater::freeParameter()
	{
	is_fixed = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Performs an update. For moves this involves calling proposeNewState, computing the posterior density of the new 
|	state, and deciding whether to accept or reject the proposed move. For parameters, this involves choosing the next
|	parameter value by slice sampling. This base class version does nothing. In derived classes, update should return 
|   true if it changed anything, false if whatever was proposed was rejected and thus nothing was changed.
*/
inline bool MCMCUpdater::update()
	{
    if (save_debug_info)
        {
        debug_info = "In MCMCUpdater::update(): this base class method does nothing";
        }
    return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the value of the data member `is_move' is false, and vice versa.
*/
inline bool	MCMCUpdater::isParameter() const
	{
	return !is_move;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `is_master_param', which is true if this is a master parameter (a parameter
|	object that computes the joint prior density for several model parameters but does not update any of them. This 
|	function can be queried if it is important to know whether the data member 'curr_value' is meaningful (it is not
|	in the case of a master parameter).
*/
inline bool	MCMCUpdater::isMasterParameter() const
	{
	return is_master_param;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `is_hyper_param', which is true if this represents a hyperparameter (a model 
|	parameter that is part of the prior specification but does not appear in the likelihood function).
*/
inline bool	MCMCUpdater::isHyperParameter() const
	{
	return is_hyper_param;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `has_slice_sampler'.
*/
inline bool	MCMCUpdater::hasSliceSampler() const
	{
	return has_slice_sampler;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `is_move'.
*/
inline bool	MCMCUpdater::isMove() const
	{
	return is_move;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing, just returns 0.0. Derived classes representing MCMC moves should override 
|	this to return the log of the Hastings ratio (only necessary if the proposal is not symmetric).
*/
inline double MCMCUpdater::getLnHastingsRatio() const
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing, just returns 0.0. Derived classes representing MCMC moves should override 
|	this to return the log of the Jacobian (only necessary if the proposal changes the model dimension).
*/
inline double MCMCUpdater::getLnJacobian() const
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing. Derived classes representing MCMC moves should override this to implement the
|	actual proposed move. This function should not do the accept/revert component, however. That should be done in
|	update() using the functions accept() and revert().
*/
inline void MCMCUpdater::proposeNewState()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing. Derived classes representing MCMC moves should override this to do anything
|	necessary to accept the state created by the last call to proposeNewState().
*/
inline void MCMCUpdater::accept()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing. Derived classes representing MCMC moves should override this to do anything
|	necessary to revert back to the state that existed before the last call to proposeNewState().
*/
inline void MCMCUpdater::revert()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Gets the value of the `weight' data member used to determine the number of times this updater's update() function
|	should be called in each update cycle.
*/
inline unsigned MCMCUpdater::getWeight() const
	{
	return weight;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of the `weight' data member used to determine the number of times this updater's update() function
|	should be called in each update cycle.
*/
inline void MCMCUpdater::setWeight(unsigned w)
	{
	weight = w;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the boldness of proposals. The meaning of boldness is interpreted individually by derived classes.
*/
inline void MCMCUpdater::setBoldness(double b)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value of data member `slice_max_units', which is the maximum number of units used for slice sampling 
|	updates. If a slice sampler has been created, calls the setMaxUnits function of `slice_sampler'.
*/
inline void MCMCUpdater::setMaxUnits(unsigned max_units)
	{
	slice_max_units = max_units;
	if (slice_sampler)
		{
		slice_sampler->SetMaxUnits(slice_max_units);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This base class version does nothing. Derived classes should override this to compute the log posterior density of
|	the supplied value conditional on all other parameter values.
*/
inline double MCMCUpdater::operator()(double)
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This creates a SliceSampler object and assigns it to the `slice_sampler' data member, which is a shared_ptr that 
|	points to nothing initially. The SliceSampler constructor takes a boost::shared_ptr<Lot> (which we have available)  
|	FuncToSampleShPtr (this object). Thus, the SliceSampler cannot be created in the constructor because "this" does
|	not yet fully exist, hence the need for a createSliceSampler() member function, which needs to be called at some 
|	point after the MCMCUpdater-derived object is created and of course before its `slice_sampler' data member starts 
|	being used.
*/
inline void MCMCUpdater::createSliceSampler() 
	{
	PHYCAS_ASSERT(!slice_sampler);	// don't want to do this more than once
	slice_sampler.reset(new SliceSampler(rng, shared_from_this())); // forces inclusion of "phycas/src/slice_sampler.hpp"
	slice_sampler->SetMaxUnits(slice_max_units);
	slice_sampler->SetXValue(curr_value);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `slice_sampler' data member, which may or may not point to a SliceSampler object. Some MCMCUpdater-
|	derived classes do not use slice sampling to update their parameter's value, so these have an empty `slice_sampler'
|	member. Also, MCMCUpdater-derived classes that represent moves rather than parameters will generally not use their 
|	`slice_sampler' data members.
*/
inline SliceSamplerShPtr MCMCUpdater::getSliceSampler()
	{
	return slice_sampler;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   If this function is called with the value true, the string data member `debug_info' will be cleared and filled upon
|   each update to reflect the last value of the parameter or move. This is useful for debugging for observing exactly
|   where the MCMC analysis begins to diverge from its behavior prior to a change in the code. Call the public member
|   function getDebugInfo to retrieve the `debug_info' string. Calling this function with the value false stops the 
|   saving of debug information, thus increasing the efficiency of the MCMC analysis.
*/
inline void MCMCUpdater::setSaveDebugInfo(
  bool save_info)       /**> is either true (to begin saving debug info) or false (to stop saving debug info) */
	{
	save_debug_info = save_info;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a copy of `debug_info', which holds information about the last update if `save_debug_info' is true. Call
|   the public member function setSaveDebugInfo to change the state of `save_debug_info'. Note that this function leaves
|   the `debug_info' string in its current state, so subsequent calls before the next update will return the same 
|   information.
*/
inline std::string MCMCUpdater::getDebugInfo() const
	{
	return debug_info;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `prior' data member, which points to a ProbabilityDistribution object. This accessor will primarily be
|	of use to moves, as they may use the master edge length updater's prior distribution to directly compute the 
|	relative log-prior over only the edge lengths they modify during an update.
*/
inline ProbDistShPtr MCMCUpdater::getPriorDist()
	{
	return prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The MCMCChainManager::finalize function calls this function for each MCMCUpdater it knows about. This provides a way 
|	for each updater to call back the MCMCChainManager when it needs the joint prior over all parameters.
*/
inline void MCMCUpdater::setChainManager(
 ChainManagerWkPtr p)		/**< is a pointer to the MCMCChainManager containing this updater */
	{
	chain_mgr = p;
	if (has_slice_sampler)
		createSliceSampler();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Useful primarily for debugging, this function uses the supplied value `x' rather than `curr_value' to set the value
|   of `curr_ln_prior' for this updater.
*/
inline double MCMCUpdater::setCurrLnPrior(
  double x) /**< is the value to use, temporarily, to set the prior */
    {
	try 
		{
		curr_ln_prior = prior->GetLnPDF(x);
		}
	catch(XProbDist &)
		{
		curr_ln_prior = prior->GetRelativeLnPDF(x);
		}
	return curr_ln_prior;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Calls `prior'->GetLnPDF(`curr_value') to recalculate `curr_ln_prior'. This function is important because the 
|	tempting getLnPrior() member function only returns the value of `curr_ln_prior' (it does not recalculate anything).
*/
inline double MCMCUpdater::recalcPrior()
	{
	try 
		{
		curr_ln_prior = prior->GetLnPDF(curr_value);
		}
	catch(XProbDist &)
		{
		curr_ln_prior = prior->GetRelativeLnPDF(curr_value);
		}
	return curr_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the log likelihood just after this updater's update() member function was last called.
*/
inline double MCMCUpdater::getLnLike() const
	{
	return curr_ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a const reference to the `name' data member.
*/
inline const std::string & MCMCUpdater::getName() const
	{
	return name;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a const reference to the string returned by the GetDistributionDescription() function of the `prior' data 
|	member.
*/
inline std::string MCMCUpdater::getPriorDescr() const
	{
	return prior->GetDistributionDescription();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of the updater (the value of the parameter just after its operator() was last called).
*/
inline double MCMCUpdater::getCurrValue() const
	{
	return curr_value;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `name' data member of this move to the supplied string.
*/
inline void MCMCUpdater::setName(const std::string & s)
	{
	name = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `rng' data member, which is a shared pointer to the pseudorandom number generator (Lot object).
*/
inline LotShPtr MCMCUpdater::getLot()
	{
	return rng;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `rng' data member to the supplied Lot shared pointer.
*/
inline void MCMCUpdater::setLot(LotShPtr r)
	{
	rng = r;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `tree' data member to the supplied Tree shared pointer. Also calls the setTree() member function of the
|	`tree_manipulator' data member. 
*/
inline void MCMCUpdater::setTree(TreeShPtr p)
	{
	tree = p;
	tree_manipulator.setTree(p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `likelihood' data member to the supplied TreeLikelihood shared pointer.
*/
inline void MCMCUpdater::setTreeLikelihood(TreeLikeShPtr L)
	{
	likelihood = L;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `model' data member to the supplied Model shared pointer.
*/
inline void MCMCUpdater::setModel(ModelShPtr m)
	{
	model = m;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `prior' data member to the supplied ProbabilityDistribution shared pointer.
*/
inline void MCMCUpdater::setPrior(ProbDistShPtr p)
	{
	prior = p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `mvprior' data member to the supplied MultivariateProbabilityDistribution shared pointer.
*/
inline void MCMCUpdater::setMultivarPrior(MultivarProbDistShPtr p)
	{
	mvprior = p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the current value to the supplied value `x'. The value is stored in the data member `curr_value' and passed to
|   the slice sampler if and when it is created in the member function MCMCUpdater::createSliceSampler().
*/
inline void MCMCUpdater::setStartingValue(double x)
	{
    curr_value = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If data member `prior' actually points to a ProbabilityDistribution object, this function calls the 
|	ProbabilityDistribution::SetMeanAndVariance function to change the mean and variance of the prior distribution. If
|	`prior' does not currently point to anything, then the function does nothing.
*/
inline void MCMCUpdater::setPriorMeanAndVariance(double m, double v)
	{
	if (prior)
		{
		//std::cerr << "    setting prior mean, variance to " << m << ", " << v << " for " << getName() << std::endl;
		prior->SetMeanAndVariance(m, v);	// forces inclusion of "phycas/src/probability_distribution.hpp"
		}
	}

} // namespace phycas
#endif
