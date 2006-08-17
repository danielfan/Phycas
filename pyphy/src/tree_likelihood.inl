#if !defined(TREE_LIKELIHOOD_INL)
#define TREE_LIKELIHOOD_INL

#include "pyphy/src/internal_data.hpp"
#include "pyphy/src/tip_data.hpp"
#include "pyphy/src/edge_endpoints.hpp"

namespace phycas
{

// **************************************************************************************
// ***** TreeLikelihood inlines *********************************************************
// **************************************************************************************

/*----------------------------------------------------------------------------------------------------------------------
|	TreeLikelihood constructor.
*/
inline TreeLikelihood::TreeLikelihood(
  ModelShPtr mod)		/**< is the substitution model */
  :
  likelihood_root(0),
  store_site_likes(true),
  no_data(false),
  nTaxa(0),
  num_patterns(0),
  num_states(mod->getNStates()),
  num_rates(mod->getNRatesTotal()),
  model(mod), 
  rate_means(mod->getNRatesTotal(), 1.0), 
  rate_probs(mod->getNRatesTotal(), 1.0), 
  nevals(0)
	{
	mod->recalcRatesAndProbs(rate_means, rate_probs);
	underflow_policy.setTriggerSensitivity(50);
	underflow_policy.setCorrectToValue(10000.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the setTriggerSensitivity function of the data member `underflow_policy' to set the number of edges that must
|	be traversed before taking action to prevent underflow.
*/
inline void TreeLikelihood::setUFNumEdges(

  unsigned nedges)	/**< is the number of edges to traverse before taking action to prevent underflow */

	{

	underflow_policy.setTriggerSensitivity(nedges);

	}


/*----------------------------------------------------------------------------------------------------------------------

|	Returns number of bytes allocated for each CLA. This equals sizeof(LikeFltType) times the product of the number of

|	patterns, number of rates and number of states. Calls corresponding function of data member `cla_pool' to get the

|	value returned.

*/

inline unsigned TreeLikelihood::bytesPerCLA() const

	{

	return cla_pool.bytesPerCLA();

	}



/*----------------------------------------------------------------------------------------------------------------------

|	Returns the number of CondLikelihood objects created since the `cla_pool' data member was constructed, or since the

|	last call to the function clearStack of `cla_pool', which resets the value to zero. Calls corresponding function of

|	data member `cla_pool' to get the value returned.

*/

inline unsigned TreeLikelihood::numCLAsCreated() const

	{

	return cla_pool.numCLAsCreated();

	}



/*----------------------------------------------------------------------------------------------------------------------

|	Returns the current number of CondLikelihood objects stored in `cla_pool'. The total number of CLAs currently

|	checked out to the tree can be obtained as TreeLikelihood::numCLAsCreated() minus TreeLikelihood::numCLAsStored().

|	Calls corresponding function of data member `cla_pool' to get the value returned.

*/

inline unsigned TreeLikelihood::numCLAsStored() const

	{

	return cla_pool.numCLAsStored();

	}



/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of `likelihood_root'. See TreeLikelihood::useAsLikelihoodRoot for more information about
|	the meaning of the likelihood root.
*/
inline TreeNode * TreeLikelihood::getLikelihoodRoot()

	{

	return likelihood_root;

	}



/*----------------------------------------------------------------------------------------------------------------------
|	Specifies the (internal) node to use as the likelihood root (it will be stored in the `likelihood_root' data member). 
|	The likelihood root is separate from the actual root of the tree, and specifies the node used when harvesting the 
|	log-likelihood from surrounding conditional likelihood arrays (CLAs). It behooves one to set the likelihood root to 
|	that node requiring the fewest CLA recalculations. Specifying NULL for the likelihood root will result in the 
|	unconditional recalculation of all CLAs in the entire tree, and subsequently the likelihood root will be set to the
|	subroot node of the tree (the only child of the root).
*/
inline void	TreeLikelihood::useAsLikelihoodRoot(

  TreeNode * nd)

	{

	likelihood_root = nd;

	}



/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of the node currently serving as the likelihood root. If `likelihood_root' is NULL, returns -1
|	instead to indicate that no node is currently designated as the likelihood root. This function was written primarily
|	for use by the TreeViewer.py application for debugging purposes.
*/
inline int TreeLikelihood::getLikelihoodRootNodeNum() const

	{

	return (likelihood_root ? (int)likelihood_root->GetNodeNumber() : -1);

	}



/*----------------------------------------------------------------------------------------------------------------------
|	Returns a shared pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a 
|	particular state being assigned to `focal_nd' when `avoid' is the likelihood root. This function obtains the correct

|	shared pointer, but if that pointer does not point to an object, it goes to `cla_pool' to get another one.
*/
inline CondLikelihoodShPtr getCondLikePtr(
  TreeNode * focal_nd,	/**< is the focal node */
  TreeNode * avoid) 	/**< is the focal node neighbor (node closer to likelihood root) */
	{
	EdgeEndpoints e(focal_nd, avoid);
	return getCondLikePtr(e);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a shared pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a 
|	particular state being assigned to `focal_nd' when `avoid' is the likelihood root. 
*/
inline ConstCondLikelihoodShPtr getValidCondLikePtr(
  const TreeNode * focal_nd,	/**< is the focal node */
  const TreeNode * avoid) 		/**< is the focal node neighbor (node closer to likelihood root) */
	{
	ConstEdgeEndpoints e(focal_nd, avoid);
	return getValidCondLikePtr(e);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a shared pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a
|	particular state being assigned to the focal node of `edge' when `edge' focal neighbor is the likelihood root. This

|	function obtain the correct shared pointer, but if that pointer does not point to a CondLikelihood object, it goes 
|	to `cla_pool' to get one.
*/
inline CondLikelihoodShPtr getCondLikePtr(
  EdgeEndpoints edge) /**< is the edge specifying the focal node and the focal node neighbor */
	{
	TreeNode * actual_child = edge.getActualChild();
	if (actual_child == edge.getFocalNode())
		{
		// focal node F is a child of N, the focal neighbor (likelihood root is somewhere below N)
		//
		//      \   /
		//       \ /
		//        F
		//   \   / <-- filial CLA of focal node is returned
		//    \ /
		//     N
		//    / 
		//
		PHYCAS_ASSERT(actual_child->IsInternal());
		InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getChildCondLikePtr();
		}

	// focal neighbor N is a child of the focal node F
	PHYCAS_ASSERT(actual_child == edge.getFocalNeighbor());
	if (actual_child->IsInternal())
		{
		// focal neighbor N is an internal node (likelihood root is somewhere above N)
		//
		//      \   /
		//       \ /
		//        N
		//   \   /
		//    \ / <-- parental CLA of focal neighbor is returned
		//     F
		//    / 
		//
		InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getParentalCondLikePtr();		
		}

	// focal neighbor N is a tip node (N equals the likelihood root in this case)
	//
	//        N
	//   \   /
	//    \ / <-- parental CLA of focal neighbor is returned
	//     F
	//    / 
	//
	TipData * child_tip_data = actual_child->GetTipData();
	PHYCAS_ASSERT(child_tip_data != NULL);
	return child_tip_data->getParentalCondLikePtr();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a pointer to the CondLikelihood object that would be used to compute the likelihood conditional on a
|	particular state being assigned to the focal node of `edge' when `edge' focal neighbor is the likelihood root.
|	This is a const version of the corresponding getCondLikePtr function, to be used when it can be assumed that the
|	conditional likelihood arrays are up-to-date.
*/
inline ConstCondLikelihoodShPtr getValidCondLikePtr(
  ConstEdgeEndpoints edge) /**< is the edge comprising the focal node and focal node neighbor */
	{
#if 1
	const TreeNode * actual_child = edge.getActualChild();
	if (actual_child == edge.getFocalNode())
		{
		// focal node F is a child of N, the focal neighbor (likelihood root is somewhere below N)
		//
		//      \   /
		//       \ /
		//        F
		//   \   / <-- filial CLA of focal node is returned
		//    \ /
		//     N
		//    / 
		//
		PHYCAS_ASSERT(actual_child->IsInternal());
		const InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getValidChildCondLikePtr();
		}

	// focal neighbor N is a child of the focal node F
	PHYCAS_ASSERT(actual_child == edge.getFocalNeighbor());
	if (actual_child->IsInternal())
		{
		// focal neighbor N is an internal node (likelihood root is somewhere above N)
		//
		//      \   /
		//       \ /
		//        N
		//   \   /
		//    \ / <-- parental CLA of focal neighbor is returned
		//     F
		//    / 
		//
		const InternalData * child_internal_data = actual_child->GetInternalData();
		PHYCAS_ASSERT(child_internal_data != NULL);
		return child_internal_data->getValidParentalCondLikePtr();		
		}

	// focal neighbor N is a tip node (N equals the likelihood root in this case)
	//
	//        N
	//   \   /
	//    \ / <-- parental CLA of focal neighbor is returned
	//     F
	//    / 
	//
	const TipData * child_tip_data = actual_child->GetTipData();
	PHYCAS_ASSERT(child_tip_data != NULL);
	return child_tip_data->getValidParentalCondLikePtr();
#else
	const TreeNode * c = edge.getActualChild();
	if (edge.getFocalNode() == c)
		{
		PHYCAS_ASSERT(c->IsInternal());
		const InternalData * childInternalData = c->GetInternalData();
		PHYCAS_ASSERT(childInternalData != NULL);
		return childInternalData->getValidChildCondLike();
		}
	// moving up the tree in calculations (root to leaves).
	PHYCAS_ASSERT(c == edge.getFocalNeighbor());
	if (c->IsInternal())
		{
		const InternalData * childInternalData = c->GetInternalData();
		PHYCAS_ASSERT(childInternalData != NULL);
		return childInternalData->getValidParentalCondLike();		
		}
	const TipData * childTipData = c->GetTipData();
	PHYCAS_ASSERT(childTipData != NULL);
	return childTipData->getValidParentalCondLike();
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the value of the `no_data' data member to true.
*/
inline void TreeLikelihood::setNoData()
	{
	no_data = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the value of the `no_data' data member to false.
*/
inline void TreeLikelihood::setHaveData()
	{
	no_data = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the `num_patterns' data member.
*/
inline void TreeLikelihood::setNPatterns(
  unsigned nPatterns)	/**< is the number of patterns */
	{
	num_patterns = nPatterns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the `model' data member, replacing the model defined in the constructor.
*/
inline void TreeLikelihood::replaceModel(
  ModelShPtr m)
	{
	model = m;
	recalcRelativeRates();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the values of the `num_states' and `num_rates' data members according to the model, then calls the 
|	recalcRatesAndProbs function of the model to force recalculation of the `rate_means' and `rate_probs' vectors. 
|	Should be called after changing the number of rate categories, the gamma shape parameter, or the pinvar parameter of
|	the model. Note that trees on which likelihoods need to be calculated also need to be re-equipped by calling 
|	prepareForLikelihood if the number of rate categories changes (not done by this function). 
*/
inline void TreeLikelihood::recalcRelativeRates()
	{
	num_states = model->getNStates();
	num_rates = model->getNRatesTotal();
	model->recalcRatesAndProbs(rate_means, rate_probs);
	likelihood_rate_site.resize(num_rates*num_patterns, 0.0);
	if (!no_data)
		cla_pool.setCondLikeDimensions(num_patterns, num_rates, num_states);
	underflow_policy.setDimensions(num_patterns, num_rates, num_states);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `cla_pool' data member.
*/
inline const CondLikelihoodStorage & TreeLikelihood::getCLAStorage() const
	{
	return cla_pool;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `nTaxa' data member.
*/
inline unsigned TreeLikelihood::getNTaxa() const
	{
	return nTaxa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `num_patterns' data member.
*/
inline unsigned TreeLikelihood::getNPatterns() const
	{
	return num_patterns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `num_rates' data member. Assumes that the length of the `rate_means' vector 
|	equals the value of the `num_rates' data member.
*/
inline unsigned TreeLikelihood::getNRatesTotal() const
	{
	PHYCAS_ASSERT(rate_means.size() == num_rates);
	return num_rates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the `num_states' data member.
*/
inline unsigned TreeLikelihood::getNStates() const
	{
	return num_states;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns a copy of the (shared_ptr) data member `model'.
*/
inline ModelShPtr TreeLikelihood::getModel() const
	{
	return model;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `state_list'.
*/
inline const VecStateList & TreeLikelihood::getStateList() const
	{
	return state_list;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `state_list_pos'.
*/
inline const VecStateListPos & TreeLikelihood::getStateListPos() const
	{
	return state_list_pos;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `rate_means'.
*/
inline const std::vector<double> & TreeLikelihood::getRateMeans() const
	{
	return rate_means;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `rate_probs'.
*/
inline const std::vector<double> & TreeLikelihood::getRateProbs() const
	{
	return rate_probs;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `category_boundaries'.
*/
inline std::vector<double> TreeLikelihood::getCategoryLowerBoundaries() const
	{
	std::vector<double> tmp_means;
	std::vector<double> returned_boundaries;
	model->recalcGammaRatesAndBoundaries(tmp_means, returned_boundaries);
	return returned_boundaries;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls TreeLikelihood::simulateImpl specifying that the transition probabilities should be recalculated (or 
|	calculated the first time) before beginning simulations. After TreeLikelihood::simulateFirst is called once, 
|	TreeLikelihood::simulate can be called many more times to generate more data sets using the same transition 
|	probabilities.
*/
inline void TreeLikelihood::simulateFirst(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar)
	{
	simulateImpl(sim_data, t, rng, nchar, true);	// true means recalculate transition probabilities
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls TreeLikelihood::simulateImpl specifying that the transition probabilities should NOT be recalculated. Only
|	call this function after calling TreeLikelihood::simulateFirst at least once, otherwise the transition probabilities
|	will contain garbage.
*/
inline void TreeLikelihood::simulate(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar)
	{
	simulateImpl(sim_data, t, rng, nchar, false);	// false means do not recalculate transition probabilities
	}

} // namespace phycas

#endif
