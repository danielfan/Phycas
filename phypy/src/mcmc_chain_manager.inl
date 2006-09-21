#if ! defined(MCMC_CHAIN_MANAGER_INL)
#define MCMC_CHAIN_MANAGER_INL

#include <algorithm>				// for std::for_each	
#include <numeric>					// for std::accumulate	
#include <boost/lambda/lambda.hpp>	// for boost::lambda::_1, boost::lambda::_2
#include <boost/lambda/bind.hpp>	// for boost::lambda::bind
//#include <boost/lambda/if.hpp>		// for boost::lambda::if_then
//#include <boost/lambda/casts.hpp>	// for boost::lambda::ll_dynamic_cast

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	The MCMCChainManager constructor does nothing.
*/
inline MCMCChainManager::MCMCChainManager() 
  : last_ln_like(0.0), last_ln_prior(0.0), dirty(true)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `all_updaters' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through all updaters using a construct such as this: for x in chainMgr.getAllUpdaters(): ...
*/
inline const MCMCUpdaterVect & MCMCChainManager::getAllUpdaters() const 
	{
	return all_updaters;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `moves' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through the stored moves using a construct such as this: for x in chainMgr.getMoves(): ...
*/
inline const MCMCUpdaterVect & MCMCChainManager::getMoves() const 
	{
	return moves;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `model_params' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through the model-specific parameters using a construct such as this: for x in chainMgr.getModelParams(): ...
*/
inline const MCMCUpdaterVect & MCMCChainManager::getModelParams() const 
	{
	return model_params;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `edge_len_params' data member as a const MCMCUpdaterVect reference. Used so that Python programs can iterate 
|	through the edge length parameters using a construct such as this: for x in chainMgr.getEdgeLenParams(): ...
*/
inline const MCMCUpdaterVect & MCMCChainManager::getEdgeLenParams() const 
	{
	return edge_len_params;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the `edge_len_hyperparam' data member, which is a shared pointer to a HyperPriorParam.
*/
inline MCMCUpdaterShPtr MCMCChainManager::getEdgeLenHyperparam()
	{
	return edge_len_hyperparam;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds a pointer to an MCMCUpdater-derived object to the `moves' vector. The `moves' vector is a vector of shared 
|	pointers to ensure that each move object is eventually destroyed.
*/
inline void	MCMCChainManager::addMove(MCMCUpdaterShPtr p) 
	{
	moves.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds a pointer to an MCMCUpdater-derived object to the `model_params' vector. The `model_params' vector is a 
|	vector of shared pointers to ensure that each parameter object is eventually destroyed.
*/
inline void	MCMCChainManager::addModelParam(MCMCUpdaterShPtr p) 
	{
	model_params.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds an EdgeLenParam shared pointer to the `edge_len_params' vector.
*/
inline void	MCMCChainManager::addEdgeLenParam(MCMCUpdaterShPtr p) 
	{
	edge_len_params.push_back(p);
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Specifies a pointer to a HyperPriorParam to be stored in the `edge_len_hyperparam' data member.
*/
inline void	MCMCChainManager::addEdgeLenHyperparam(MCMCUpdaterShPtr p) 
	{
	edge_len_hyperparam = p;
	dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current log likelihood. Each call to an update() member function of one of the updaters refreshes the
|	`last_ln_like' data member.
*/
inline double MCMCChainManager::getLastLnLike() const 
	{
	return last_ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Used by an MCMCUpdater-derived object to set the current log prior just before it returns from its update()
|	function. Ordinarily, the calcJointLnPrior() would be used instead to set `last_ln_prior', but for moves that need
|	to revert this function is used because the previous log prior was saved and can thus easily be restored.
*/
inline void MCMCChainManager::setLastLnPrior(double ln_prior) 
	{
	last_ln_prior = ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Used by an MCMCUpdater-derived object to set the current log likelihood just before it returns from its update()
|	function.
*/
inline void MCMCChainManager::setLastLnLike(double ln_like) 
	{
	last_ln_like = ln_like;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current log joint prior. Each call to the calcJointLnPrior() member function refreshes the 
|	`last_ln_prior' data member.
*/
inline double MCMCChainManager::getLastLnPrior() const 
	{
	return last_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Clears all vectors (all_updaters, moves, model_params, edge_len_params) and resets edge_len_hyperparam.
*/
inline void MCMCChainManager::clear()
	{
	all_updaters.clear();
	moves.clear();
	model_params.clear();
	edge_len_params.clear();
	edge_len_hyperparam.reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `chain_mgr' data member of all parameters so that they can call MCMCChainManager::getLnPrior when any one 
|	of them needs to compute the joint prior over all parameters (or wishes to refresh the MCMCChainManager's 
|	`last_ln_like' data member. The MCMCUpdater::`chain_mgr' data member is a weak_ptr (if a shared_ptr were used, we 
|	would have a cycle: MCMCChainManager holds a shared_ptr to the MCMCUpdater object, and MCMCUpdater would hold a 
|	shared_ptr to MCMCChainManager). The weak_ptr that is passed to the setChainManager member function of each 
|	MCMCUpdater object is created from a shared_ptr to this, which we can get because MCMCChainManager was derived from
|	boost::enable_shared_from_this<MCMCChainManager>. Using the shared_from_this() function presumes that there exists 
|	a shared_ptr already managing this. Thus, a MCMCChainManager object should not be created within C++ without using 
|	shared_ptr to manage it.
*/
inline void MCMCChainManager::finalize()
	{
	// Determine how big the all_updaters vector needs to be
	unsigned sz = (unsigned)(moves.size() + edge_len_params.size() + model_params.size()) + (edge_len_hyperparam ? 1 : 0);
	if (sz == 0)
		{
		throw XLikelihood("must add at least one updater object to the chain manager before calling finalize()");
		}

	// Clear all_updaters vector
	all_updaters.clear();
	all_updaters.resize(sz);

	// Add the moves first
	moves_end = std::copy(moves.begin(), moves.end(), all_updaters.begin());
	moves_begin = all_updaters.begin();

	// Next add the normal edge length parameters
	edgelens_end = std::copy(edge_len_params.begin(), edge_len_params.end(), moves_end);
	edgelens_begin = moves_end;

	// Next add the edge length hyperparameter (if used)
	hyperparam_end = edgelens_end;
	if (edge_len_hyperparam)
		*hyperparam_end++ = edge_len_hyperparam;
	hyperparam_begin = edgelens_end;

	// Next add the remaining model parameters
	model_params_end = std::copy(model_params.begin(), model_params.end(), hyperparam_end);
	model_params_begin = hyperparam_end;

	params_begin = edgelens_begin;
	params_end = model_params_end;

	// Call each updater's setChainManager member function, supplying a shared pointer to this object
	// This allows each updater to ask the MCMCChainManager to calculate the joint prior over all
	// parameters when needed for computing the posterior density during slice sampling
	std::for_each(all_updaters.begin(), all_updaters.end(), 
		boost::lambda::bind(&MCMCUpdater::setChainManager, *boost::lambda::_1, ChainManagerWkPtr(shared_from_this())));

	// Call each parameter's recalcPrior member function to make sure their curr_ln_prior data members are
	// up to date
	std::for_each(params_begin, params_end, 
		boost::lambda::bind(&MCMCUpdater::recalcPrior, *boost::lambda::_1));

	//@POL Wondering why we need all these vectors? It is convenient to have one vector of updaters, but not all of 
	// these updaters are necessarily created in MCMCChainManager::addMCMCUpdaters. When finalize() is called, we 
	// can be sure that all the updaters that are going to be defined have been defined, and at that point we can 
	// quickly build the vector all_updaters. In the future, if it turns out we never need moves, model_params, or
	// edge_len_params after this point, we can just clear them at this point because we have iterators now set to
	// the beginning and end of each section of the all_updaters vector anyway.

	dirty = false;
	}

// Using USE_PARTIAL_SUM_FUNCTOR is slightly faster, but makes the assumption that the edge length prior
// distribution is Exponential, which is not good. This assumption is not necessary, but without it
// speed suffers and the alternative is preferable.
//#define USE_PARTIAL_SUM_FUNCTOR

#if defined(USE_PARTIAL_SUM_FUNCTOR)
struct PartialPriorSum : public std::unary_function<double, void>
	{
	double sum;
	ExponentialDistribution & distr;

	PartialPriorSum(ExponentialDistribution & p) : sum(0.0), distr(p) 
		{}

	void operator()(const double & nxt)
		{
		sum += distr.GetLnPDF(nxt);
		}

	double result()
		{
		return sum;
		}
	};
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Uses first edge length parameter in the `all_updaters' vector to compute the log-prior for each value in the 
|	supplied vector of edge lengths, returning the sum. Useful for efficiency purposes for moves that modify only a few 
|	edge lengths.
*/
inline double MCMCChainManager::partialEdgeLenPrior(const std::vector<double> & edge_len_vect) const
	{
	if (edgelens_end - edgelens_begin == 0)
		{
		throw XLikelihood("partialEdgeLenPrior requires a master edge length prior to have been added to the chain manager");
		}

	MCMCUpdaterShPtr u = *edgelens_begin;
	ProbDistShPtr p = u->getPriorDist();

#if defined(USE_PARTIAL_SUM_FUNCTOR)
	ExponentialDistribution * d = dynamic_cast<ExponentialDistribution *>(p.get());
	PHYCAS_ASSERT(d);
	//@POL assuming edge length distribution is Exponential for now (allows GetLnPDF to be inlined)
	return std::for_each(edge_len_vect.begin(), edge_len_vect.end(), PartialPriorSum(*d)).result();
#else
	double partial_prior_sum = 0.0;
	for (std::vector<double>::const_iterator it = edge_len_vect.begin(); it != edge_len_vect.end(); ++it)
		{
		partial_prior_sum += p->GetLnPDF(*it);
		}
	return partial_prior_sum;
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Visits each parameter and calls its getLnPrior function to compute the joint log prior over all parameters. The
|	accumulate algorithm requires its fourth argument to be a binary functor, for which an unnamed boost lambda function
|	object is used. This function refreshes the value of the data member `last_ln_prior'.
*/
inline double MCMCChainManager::calcJointLnPrior()
	{
	if (all_updaters.empty())
		{
		throw XLikelihood("should not call MCMCChainManager::calcJointLnPrior before calling MCMCChainManager::finalize");
		}

	last_ln_prior = std::accumulate(params_begin, params_end, 0.0, 
		boost::lambda::_1 += boost::lambda::bind(&MCMCUpdater::getLnPrior, *boost::lambda::_2));

	return last_ln_prior;
	}

struct EdgeLenPriorUpdater
	{
	double mu;
	double var;
	double ln_prior;

	EdgeLenPriorUpdater(double m, double v) : mu(m), var(v), ln_prior(0.0) {}
	void operator()(MCMCUpdaterShPtr u)
		{
		u->setPriorMeanAndVariance(mu, var);
		ln_prior += u->recalcPrior();
		}
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Visit all EdgeLenParam objects in the `edge_len_params' list and tell each of them to use the supplied new value 
|	`mu' as the mean of their prior distribution.
*/
inline double MCMCChainManager::recalcEdgeLenPriors(
  double mu,	/**< is a new value for the edge length prior mean */
  double var)	/**< is a new value for the edge length prior variance */
	{
	if (dirty)
		{
		throw XLikelihood("cannot call MCMCChainManager::recalcEdgeLenPriors before calling MCMCChainManager::finalize");
		}

	EdgeLenPriorUpdater x = std::for_each(edgelens_begin, edgelens_end, EdgeLenPriorUpdater(mu, var));
	return x.ln_prior;

	// old way (rel_0_1_0)
	//std::for_each(edgelens_begin, edgelens_end,
	//	boost::lambda::bind(&MCMCUpdater::setPriorMeanAndVariance, *boost::lambda::_1, mu, var));
	}

} // namespace phycas

#endif