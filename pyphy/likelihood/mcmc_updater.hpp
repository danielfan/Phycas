#if ! defined(MCMC_UPDATER_HPP)
#define MCMC_UPDATER_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "pyphy/phylogeny/tree_manip.hpp"			// for TreeManip
#include "phycas/modules/mcmc/slice_sampler.hpp"
#include "phycas/rand/probability_distribution.hpp"

namespace phycas
{

class Tree;
typedef boost::shared_ptr<Tree>					TreeShPtr;

class TreeLikelihood;
typedef boost::shared_ptr<TreeLikelihood>		TreeLikeShPtr;

class Model;
typedef boost::shared_ptr<Model>				ModelShPtr;

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>		ChainManagerWkPtr;

class MCMCUpdater;
typedef boost::shared_ptr<MCMCUpdater>			MCMCUpdaterShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Class encapsulating both MCMC moves (i.e. Metropolis-Hastings updates) and model parameters (slice sampling 
|	updates). Provides member functions for specifying a name to be used for this updater, as well as objects needed to
|	compute the poseterior density: a tree, a substitution model, a pseudo-random number generator, and a likelihood 
|	calculator. These should all be set before adding the MCMCUpdater to the MCMCChainManager. Inside 
|	MCMCChainManager::finalize(), the MCMCChainManager calls the setChainManager() itself for all MCMCUpdater objects
|	that have been added to it. The key virtual member function of this interface is update(). MCMC move objects 
|	override update() to perform a Metropolis-Hastings proposal and either accept or reject it. Parameter objects 
|	override update() to perform a single slice sampling update of the model parameter they are managing. Below is a 
|	list of some of the differences between moves and parameters:
|	
|	Moves have these properties:
|	1) do NOT use the createSliceSampler() or getSliceSampler() member functions and thus leave their `slice_sampler' 
|		data member empty
|	2) do NOT use the `curr_value' data member or the getCurrValue() member function
|	3) do NOT use the setPrior() and setPriorMeanAndVariance() member functions
|	4) do NOT use the operator() member function
|	5) DO use the getLnHastingsRatio(), getLnJacobian(), proposeNewState(), accept() and revert() member functions
|	6) some use the `tree_manipulator' data member to effect topological rearrangements
|	7) the `has_slice_sampler' data member is always false
|	8) the `is_move' data member is always true
|	
|	Parameters have these properties:
|	1) do NOT use the createSliceSampler() or getSliceSampler() member functions and thus leave their `slice_sampler' 
|		data member empty
|	2) DO use the `curr_value' data member or the getCurrValue() member function
|	3) DO use the setPrior() and setPriorMeanAndVariance() member functions
|	4) DO usethe operator() member function
|	5) do NOT use the getLnHastingsRatio(), getLnJacobian(), proposeNewState(), accept() and revert() member functions
|	6) never use the `tree_manipulator' data member
|	7) the `has_slice_sampler' data member is usually true (exception is EdgeLenMasterParam)
|	8) the `is_move' data member is always false
|	
|	The `weight' data member determines how many times in a single update cycle the update() function of this 
|	particular MCMCUpdater will be called. Typically, parameters have weight 1 whereas moves (e.g. LargetSimonMove) 
|	have a higher weight so that their proposals will be attempted several times for every parameter update.
*/
class MCMCUpdater : public AdHocDensity, public boost::enable_shared_from_this<MCMCUpdater>
	{
	public:
								MCMCUpdater();
								virtual	~MCMCUpdater() 
									{
									//std::cerr << "MCMCUpdater dying..." << std::endl;
									}

		// Predicates
		bool					isParameter() const;
		bool					isMasterParameter() const;
		bool					isHyperParameter() const;
		bool					hasSliceSampler() const;
		bool					isMove() const;

		// Accessors
		const std::string &		getName() const;
		unsigned				getWeight() const;
		double					getLnLike() const;
		double					getLnPrior() const;
		virtual std::string 	getPriorDescr() const;
		LotShPtr				getLot();

		// Accessors used only by parameters
		SliceSamplerShPtr		getSliceSampler();
		double					getCurrValue() const;

		// Accessors used only by moves
		ProbDistShPtr			getPriorDist();

		// Modifiers
		virtual void			setName(const std::string & s);
		virtual void			setWeight(unsigned w);
		void					setMaxUnits(unsigned max_units);
		virtual void			setLot(LotShPtr p);
		virtual void			setTree(TreeShPtr p);
		virtual void			setTreeLikelihood(TreeLikeShPtr p);
		virtual void			setModel(ModelShPtr p);
		virtual void			setChainManager(ChainManagerWkPtr p);

		// Modifiers used only by parameters
		virtual void			setPrior(ProbDistShPtr p);
		virtual void			setPriorMeanAndVariance(double m, double v);
		virtual double			recalcPrior();

		// Utilities
		virtual void			update();
		virtual double			recalcLike();

		// Utilities used only by parameters
		bool					isFixed() const;
		void					fixParameter();
		void					freeParameter();
		void					createSliceSampler();
		virtual double			operator()(double);

		// Note: some member functions could be made pure virtuals were it not for a bug in the 
		// boost::lambda library that causes compiles to fail if attempting to use boost::lambda::bind 
		// with a shared_ptr to an abstract base class. Even if the shared_ptr points to a derived class, 
		// the boost::lambda library attempts (without success) to create an object of the abstract base 
		// class. We have two choices: give up using boost::lambda::bind (painful) or avoid abstract base 
		// classes that are going to be used with shared_ptr (also painful, but less so). Hopefully, the 
		// bug will eventually get fixed and we can have our cake and eat it too. For a discussion of 
		// this bug by the author of the boost::lambda library (Jaakko Jarvi), see 
		// http://lists.boost.org/Archives/boost/2003/10/55139.php

	protected:

		// Virtual functions needed only by moves
		//
		virtual double			getLnHastingsRatio() const;
		virtual double			getLnJacobian() const;
		virtual void 			proposeNewState();
		virtual void			accept();
		virtual void			revert();

	protected:

		std::string				name;					/**< Storage for the name of this parameter to be used in reporting to the user */
		unsigned				weight;					/**< The number of times this updater's update() function should be called in each update cycle */
		TreeShPtr				tree;					/**< The tree on which the likelihood will be calculated */
		TreeManip				tree_manipulator;		/**< The object that facilitates topological rearrangements */
		ModelShPtr				model;					/**< The substitution model to be used in computing the likelihood */
		TreeLikeShPtr			likelihood;				/**< The object that knows how to calculate the likelihood */
		LotShPtr				rng;					/**< The pseudorandom number generator object used in updating parameter value */
		ProbDistShPtr			prior;					/**< The probability distribution serving as the prior for a parameter (not used by moves) */
		SliceSamplerShPtr		slice_sampler;			/**< The slice sampler used by parameters for updating (not used by moves) */
		ChainManagerWkPtr		chain_mgr;				/**< The object that knows how to compute the joint log prior density */
		double					curr_value;				/**< The current value of this parameter */
		double					curr_ln_prior;			/**< The value of log prior after most recent call to this object's update() */
		double					curr_ln_like;			/**< The value of log likelihood after most recent call to this object's update() */
		double					ln_zero;				/**< The value to return for the log posterior density if the parameter value is out of range */
		bool					has_slice_sampler;		/**< If true, then createSliceSampler() will be called inside setChainManager, ensuring that all parameters that need a slice sampler have one before they begin updating */
		bool					is_move;				/**< True if this updater is a move (i.e. encapsulates a Metropolis-Hastings proposal); if false, then this updater is considered a parameter */
		bool					is_master_param;		/**< True if this updater is a master parameter (does not update any model parameters but can compute the joint prior density for several model parameters) */
		bool					is_hyper_param;			/**< True if this updater represents a hyperparameter (a model parameter that is part of the prior specification but not the likelihood function) */
		bool					is_fixed;				/**< If true, update returns immediately so parameter is never updated */
		unsigned				slice_max_units;		/**< Maximum number of units used by `slice_sampler' */
	};

typedef std::vector<MCMCUpdaterShPtr>		MCMCUpdaterVect;
typedef MCMCUpdaterVect::iterator			MCMCUpdaterIter;

} // namespace phycas

#include "pyphy/likelihood/mcmc_updater.inl"

#endif
