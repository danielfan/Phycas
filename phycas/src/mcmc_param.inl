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

#if !defined(MCMC_PARAM_INL)
#define MCMC_PARAM_INL

#include <limits>										// for std::numeric_limits
#include <numeric>										// for std::accumulate
#include <boost/lambda/lambda.hpp>						// for boost::lambda::_1, boost::lambda::_2
#include <boost/lambda/bind.hpp>						// for boost::lambda::bind
#include "phycas/src/basic_tree.hpp"				// for Tree::begin() and Tree::end()

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor and initializes its HKY pointer to NULL. Also sets the
|	`curr_value' data member to 4.0 and refreshes `curr_ln_prior' accordingly.
*/
inline KappaParam::KappaParam()
  : MCMCUpdater(), hky(NULL)
	{
	curr_value = 4.0;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
inline void KappaParam::update()
	{
	//@POL should probably put next two lines in Model::update and have Model::update call a virtual Model::updateImpl
	// because it will be hard to remember to put these lines in every overloaded update function
	if (is_fixed)
		return;
	slice_sampler->Sample();

    if (save_debug_info)
        {
        debug_info = str(boost::format("KappaParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor and initializes its Codon pointer to NULL. Also sets the
|	`curr_value' data member to 1.0 and refreshes `curr_ln_prior' accordingly.
*/
inline OmegaParam::OmegaParam()
  : MCMCUpdater(), codon(NULL)
	{
	curr_value = 1.0;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
inline void OmegaParam::update()
	{
	//@POL should probably put next two lines in Model::update and have Model::update call a virtual Model::updateImpl
	// because it will be hard to remember to put these lines in every overloaded update function
	if (is_fixed)
		return;
	slice_sampler->Sample();
    
    if (save_debug_info)
        {
        debug_info = str(boost::format("OmegaParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor and initializes its GTR pointer to NULL. Also sets the
|	`curr_value' data member to 1.0 and refreshes `curr_ln_prior' accordingly. Assumes `w' is greater than or equal to
|	zero and less than 6.
*/
inline GTRRateParam::GTRRateParam(
  unsigned w)		/**< The 0-based index of the relative rate being managed by this object (0=AC, 1=AG, 2=AT, 3=CG, 4=CT and 5=GT) */
  : MCMCUpdater(), gtr(NULL), which(w)
	{
	PHYCAS_ASSERT(w >= 0);
	PHYCAS_ASSERT(w < 6);
	curr_value = 1.0;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
inline void GTRRateParam::update()
	{
	//@POL should probably put next two lines in Model::update and have Model::update call a virtual Model::updateImpl
	// because it will be hard to remember to put these lines in every overloaded update function
	if (is_fixed)
		return;
	slice_sampler->Sample();
    
    if (save_debug_info)
        {
        debug_info = str(boost::format("GTRRateParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor. Also sets the `curr_value' data member to 0.5 (or 2.0
|	if `invert' is true. Sets `invert_shape' to the value of `invert'.
*/
inline DiscreteGammaShapeParam::DiscreteGammaShapeParam(bool invert)
  : MCMCUpdater(), invert_shape(invert)
	{
	curr_value = (invert_shape ? 2.0 : 0.5);
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
inline void DiscreteGammaShapeParam::update()
	{
	if (is_fixed)
		{
		return;
		}
	slice_sampler->Sample();
    
    if (save_debug_info)
        {
        debug_info = str(boost::format("DiscreteGammaShapeParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor. Also sets the `curr_value' data member to 0.5.
*/
inline PinvarParam::PinvarParam()
  : MCMCUpdater()
	{
	curr_value = 0.5;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
inline void PinvarParam::update()
	{
	if (is_fixed)
		return;
	slice_sampler->Sample();
    
    if (save_debug_info)
        {
        debug_info = str(boost::format("PinvarParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The FlexRateParam constructor requires a reference to an unsigned variable that holds the number of spacers (fake
|	rates that serve to determine the strength of the prior for relative rate parameters). It also requires a reference
|	to a double variable that holds the maximum possible value for a relative rate. Finally, a reference to the vector
|	of relative rates is required because in order to update a particular rate parameter, one needs to know the rates 
|	on either side of it.
*/
inline FlexRateParam::FlexRateParam(
  unsigned & s,					/**< is the number of spacers between each rate in the order statistics prior */
  double & ub,					/**< is the maximum possible value for a rate */
  std::vector<double> & rr)		/**< is the relative rate parameter vector (this parameter needs to know the value on either side when updating) */
  : MCMCUpdater(),  nspacers(s), upper_bound(ub), rel_rates(rr), left_value(0.0), right_value(1.0), which(0)
	{
	curr_value = rel_rates[which];
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The MCMCChainManager::finalize function calls this function for each MCMCUpdater it knows about. This provides a way 
|	for each updater to call back the MCMCChainManager when it needs the joint prior over all parameters.
*/
inline void FlexRateParam::setChainManager(
 ChainManagerWkPtr p)		/**< is a pointer to the MCMCChainManager containing this updater */
	{
	MCMCUpdater::setChainManager(p);
	//refreshLeftRightValues();
	//double starting_x = (left_value + right_value)/2.0;
	//slice_sampler->SetXValue(starting_x);
	//slice_sampler->SetXValue(curr_value);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Determines the quantities `curr_value', `left_value' and `right_value'. The smallest this parameter can become is 
|	`left_value', which should equal `rel_rates'[`which' - 1] (unless `which' is 0, in which case `left_value' is 0.0).
|	The largest this parameter can become is `right_value', which should equal `rel_rates'[`which' + 1] (unless `which'
|	is `num_rates' - 1, in which case `right_value' is `upper_bound'). The `curr_value' of this parameter is equal to
|	`rel_rates'[`which'].
*/
inline void FlexRateParam::refreshLeftRightValues()
	{
	PHYCAS_ASSERT(rel_rates.size() > 0);
	unsigned last = (unsigned)rel_rates.size() - 1;
	if (which == 0 && which == last)
		{
		// only one rate (i.e. no rate heterogeneity)
		left_value  = 0.0;
		right_value = upper_bound;
		}
	else if (which == 0)
		{
		left_value = 0.0;
		right_value = rel_rates[which + 1];
		}
	else if (which == last)
		{
		left_value = rel_rates[which - 1];
		right_value = upper_bound;
		}
	else
		{
		left_value = rel_rates[which - 1];
		right_value = rel_rates[which + 1];
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Order statistics prior" (user gets no choice of prior in this case).
*/
inline std::string FlexRateParam::getPriorDescr() const
	{
	return std::string("Order statistics prior");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the order statistics prior used for relative rates in the FLEX model, and sets `curr_ln_prior'. The order 
|	statistics prior is largest when `curr_value' is between `left_value' and `right_value' and becomes tiny if 
|	`curr_value' gets close to either `left_value' or `right_value'. Between each pair of rates lies one or more 
|   imaginary spacers. The more spacers lying between each adjacent pair of rates, the harder it is for rates to get 
|   close together. The number of spacers between each adjacent pair is specified by the data member `nspacers', which 
|   is set in the constructor.
*/
inline double FlexRateParam::recalcPrior()
	{
	refreshLeftRightValues();
	PHYCAS_ASSERT(curr_value == rel_rates[which]);
	curr_ln_prior = (double)nspacers*(std::log(curr_value - left_value) + std::log(right_value - curr_value));
	return curr_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
inline void FlexRateParam::update()
	{
	if (is_fixed)
		return;

    // Ok, it looks like the following cannot work, but it does! The key is the assignment to `which'
    // and the fact that SliceSampler::Sample() calls FlexRateParam::operator() repeatedly. Inside
    // the operator() function, `which' is used to determine which of the relative rates is 
    // modified. See FlexRateParam::operator() in mcmc_flexcat_param.cpp for details.
	unsigned nrates = rel_rates.size();
	for (unsigned i = 0; i < nrates; ++i)
		{
		which = i;
		curr_value = rel_rates[which];
		slice_sampler->SetXValue(curr_value);
		slice_sampler->Sample();
		}

    if (save_debug_info)
        {
        debug_info = "FlexRateParam ";
    	for (unsigned i = 0; i < nrates; ++i)
            debug_info += str(boost::format("\n  %d %f") % i % rel_rates[i]);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The FlexProbParam constructor requires a reference to the vector of category probabilities so that it can deposit
|	a changed value in the appropriate place after an update.
*/
inline FlexProbParam::FlexProbParam(
  std::vector<double> & rp)		/**< is the category probability parameter vector */
  : MCMCUpdater(), which(0), rate_probs(rp)
	{
	curr_value = 1.0;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

inline void FlexProbParam::setLot(LotShPtr r)
	{
	rng = r;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member to update each of the probabilities in the 
|   vector `rate_probs'.
*/
inline void FlexProbParam::update()
	{
	if (is_fixed)
		return;

    // Ok, it looks like the following cannot work, but it does! The key is the assignment to `which'
    // and the fact that SliceSampler::Sample() calls FlexProbParam::operator() repeatedly. Inside
    // the operator() function, `which' is used to determine which of the rate probabilities is 
    // modified. See FlexProbParam::operator() in mcmc_flexcat_param.cpp for details.
	unsigned nrates = rate_probs.size();
	for (unsigned i = 0; i < nrates; ++i)
		{
		which = i;
		curr_value = rate_probs[which];
		slice_sampler->SetXValue(curr_value);
		slice_sampler->Sample();
		}

    if (save_debug_info)
        {
        debug_info = "FlexProbParam ";
    	for (unsigned i = 0; i < nrates; ++i)
            debug_info += str(boost::format("\n  %d %f") % i % rate_probs[i]);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The StateFreqParam constructor requires the caller to specify a value for `which'. `which' is 0, 1, 2, or 3 for 
|	nucleotide models and  determines the particular state frequency (e.g. A, C, G, or T, respectively) this object 
|	represents. It calls the base class (MCMCUpdater) constructor. Also sets the `curr_value' data member to 1.0 and 
|	refreshes `curr_ln_prior' accordingly.
*/
inline StateFreqParam::StateFreqParam(
  unsigned w)		/**< The 0-based index of the base frequency being managed by this object (0=A, 1=C, 2=G and 3=T) */
  : MCMCUpdater(),  which(w)
	{
	curr_value = 1.0;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
inline void StateFreqParam::update()
	{
	if (is_fixed)
		return;
	slice_sampler->Sample();

    if (save_debug_info)
        {
        debug_info = str(boost::format("StateFreqParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor of HyperPriorParam simply calls the base class (MCMCUpdater) constructor. Also sets the `curr_value'
|	data member to 0.1 and refreshes `curr_ln_prior' accordingly.
*/
inline HyperPriorParam::HyperPriorParam()
  : MCMCUpdater()
	{
	curr_value = 0.1;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
inline void HyperPriorParam::update()
	{
	if (is_fixed)
		return;
	if (slice_sampler)
		{
		slice_sampler->Sample();
		}

    if (save_debug_info)
        {
        if (slice_sampler)
            debug_info = str(boost::format("HyperPriorParam %f") % (slice_sampler->GetLastSampledXValue()));
        else
            debug_info = "HyperPriorParam (no slice sampler)";
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor does nothing except set both `has_slice_sampler' and `is_move' to false.
*/
inline EdgeLenMasterParam::EdgeLenMasterParam()
  : MCMCUpdater()
	{
	has_slice_sampler = false;
	is_move = false;
	is_master_param = true;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Member function that exists only to facilitate using boost::lambda::bind to be used in the getLnPrior() function. It
|	returns the log of the probability density evaluated at the current edge length associated with `nd'.
*/
inline double EdgeLenMasterParam::lnPriorOneEdge(TreeNode & nd) const
	{
	if (nd.IsTipRoot())
		return 0.0;

	double v = nd.GetEdgeLen();

	double retval = 0.0;
	try 
		{
		retval = prior->GetLnPDF(v);
		}
	catch(XProbDist &)
		{
		retval = prior->GetRelativeLnPDF(v);
		}

	return retval;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the joint log prior over all edges in the associated tree and sets `curr_ln_prior'.
*/
inline double EdgeLenMasterParam::recalcPrior()
	{
	curr_ln_prior = std::accumulate(tree->begin(), tree->end(), 0.0,
		boost::lambda::_1 += boost::lambda::bind(&EdgeLenParam::lnPriorOneEdge, this, boost::lambda::_2));

	return curr_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If using a hyperprior model in which the mean of the edge length prior is itself a parameter in the model, when the 
|	value of the hyperparameter changes, this change must be propagated to all EdgeLenParam objects so that they can
|	calculate their prior correctly. This function provides the means for changing the mean and variance of the prior
|	distribution. It simply calls prior->SetMeanAndVariance(), and then recalculates `curr_ln_prior' using the new prior 
|	distribution.
*/
inline void EdgeLenMasterParam::setPriorMeanAndVariance(double m, double v)
	{
	MCMCUpdater::setPriorMeanAndVariance(m, v);
	recalcPrior();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The parameter `node' is the node in the tree associated with the edge whose length is managed by this 
|	MCMCUpdater-derived object.
*/
inline EdgeLenParam::EdgeLenParam(
  TreeNode * node)		/**< is the node in the tree associated with the edge whose length is managed by this object */
	: EdgeLenMasterParam(), nd(node)
	{
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
inline void EdgeLenParam::update()
	{
	if (is_fixed)
		return;
	slice_sampler->Sample();

    if (save_debug_info)
        {
        debug_info = str(boost::format("EdgeLenParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	}

} // namespace phycas

#endif
