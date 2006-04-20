#if !defined(MCMC_PARAM_INL)
#define MCMC_PARAM_INL

#include <limits>										// for std::numeric_limits
#include <numeric>										// for std::accumulate
#include <boost/lambda/lambda.hpp>						// for boost::lambda::_1, boost::lambda::_2
#include <boost/lambda/bind.hpp>						// for boost::lambda::bind
#include "pyphy/phylogeny/basic_tree.hpp"				// for Tree::begin() and Tree::end()
//#include "phycas/rand/probability_distribution.hpp"		// see EdgeLenParam::setPriorMeanAndVariance

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
	assert(w >= 0);
	assert(w < 6);
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
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor. Also sets the `curr_value' data member to 0.5 (or 2.0
|	if `invert' is true. Sets `invert_shape' to the value of `invert'.
*/
#if POLPY_NEWWAY
inline DiscreteGammaShapeParam::DiscreteGammaShapeParam(bool invert)
  : MCMCUpdater(), invert_shape(invert)
#else
inline DiscreteGammaShapeParam::DiscreteGammaShapeParam()
  : MCMCUpdater()
#endif
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
	}

#if POLPY_NEWWAY	// pattern-specific rates
/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor. Also sets the `curr_value' data member to 0.5.
*/
inline PatternSpecificRateParam::PatternSpecificRateParam(unsigned which)
  : which_pattern(which), MCMCUpdater()
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
inline void PatternSpecificRateParam::update()
	{
	if (is_fixed)
		{
		return;
		}
	slice_sampler->Sample();
	}
#endif

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
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The BaseFreqParam constructor requires the caller to specify a value for `which'. `which' is 0, 1, 2 or 3 and 
|	determines the particular base frequency (A, C, G or T, respectively) this object represents. It calls the base
|	class (MCMCUpdater) constructor. Also sets the `curr_value' data member to 1.0 and refreshes `curr_ln_prior' 
|	accordingly.
*/
inline BaseFreqParam::BaseFreqParam(
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
inline void BaseFreqParam::update()
	{
	if (is_fixed)
		return;
	slice_sampler->Sample();
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
		//std::cerr << "In HyperPriorParam::update" << std::endl;
		slice_sampler->Sample();
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
	if (nd.IsRoot())
		return 0.0;

	double v = nd.GetEdgeLen();
	double retval = prior->GetLnPDF(v);

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
	}

} // namespace phycas

#endif
