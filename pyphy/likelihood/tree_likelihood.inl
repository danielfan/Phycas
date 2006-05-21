#if !defined(TREE_LIKELIHOOD_INL)
#define TREE_LIKELIHOOD_INL

#include "pyphy/likelihood/internal_data.hpp"
#include "pyphy/likelihood/tip_data.hpp"
#include "pyphy/phylogeny/edge_endpoints.hpp"

namespace phycas
{

// **************************************************************************************
// ***** TreeLikelihood inlines *********************************************************
// **************************************************************************************

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline TreeLikelihood::TreeLikelihood(
  ModelShPtr mod)		/**< is the substitution model */
  :
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
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Modifier function that sets the value of the `no_data' data member to true.
*/
inline void TreeLikelihood::setNoData()
	{
	no_data = true;
	}

inline CondLikelihood * getCondLike(TreeNode *focalNd, TreeNode *avoid) 
	{
	EdgeEndpoints e(focalNd, avoid);
	return getCondLike(e);
	}

inline const CondLikelihood * getCondLike(const TreeNode *focalNd, const TreeNode *avoid) 
	{
	ConstEdgeEndpoints e(focalNd, avoid);
	return getCondLike(e);
	}


inline const CondLikelihood * getCondLike(ConstEdgeEndpoints edge) 
	{
	const TreeNode * c = edge.getActualChild();
	if (edge.getFocalNode() == c)
		{
		assert(c->IsInternal());
		const InternalData * childInternalData = c->GetInternalData();
		assert(childInternalData != NULL);
		return childInternalData->getChildCondLike();
		}
	/// moving up the tree in calculations (root to leaves).
	assert(c == edge.getFocalNeighbor());
	if (c->IsInternal())
		{
		const InternalData * childInternalData = c->GetInternalData();
		assert(childInternalData != NULL);
		return childInternalData->getParentalCondLike();		
		}
	else
		{
		const TipData * childTipData = c->GetTipData();
		assert(childTipData != NULL);
		return childTipData->getParentalCondLike();
		}
	}
	
inline CondLikelihood * getCondLike(EdgeEndpoints edge) 
	{
	ConstEdgeEndpoints c(edge.first, edge.second);
	const CondLikelihood * cl = getCondLike(c);
	return const_cast<CondLikelihood *>(cl);
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
	assert(rate_means.size() == num_rates);
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
