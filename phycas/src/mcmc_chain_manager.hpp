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

#if ! defined(MCMC_CHAIN_MANAGER_HPP)
#define MCMC_CHAIN_MANAGER_HPP

#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "phycas/src/mcmc_updater.hpp"		// for MCMCUpdaterShPtr (plus LotShPtr, ModelShPtr, TreeShPtr, TreeLikeShPtr, etc.)
#include "phycas/src/mcmc_param.hpp"			// for MCMCUpdaterShPtr (plus LotShPtr, ModelShPtr, TreeShPtr, TreeLikeShPtr, etc.)

struct CIPRES_Matrix;

namespace phycas
{

class MCMCChainManager;
typedef boost::enable_shared_from_this<MCMCChainManager> MCMCChainManagerThisShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Manages a vector of pointers to MCMCUpdater objects representing MCMC moves or model parameters. A very useful
|	member function is getAllUpdaters(), which returns a reference to the `all_updaters' vector. One can iterate through
|	this vector, calling the update() member function of each updater visited, which provides the core of any MCMC
|	algorithm. It also keeps track of the current log likelihood, log prior and log posterior density, so it is useful
|	for both reporting these and as a means for an MCMC move to know the log posterior of its starting state. Because 
|	it knows about all of the parameters being used in an MCMC analysis, this class can compute the joint prior 
|	probability density over all the parameters it manages.
|	
|	To add updaters, call addMove(), addModelParam(), addEdgeLenParam() or addEdgeLenHyperparam(). Updaters can be added
|	in any order because the finalize() member function (which must be called after all updaters have been added) creates
|	the master `all_updaters' vector with updaters laid out in a logical sequence. The finalize() member function also 
|	sets the chain_mgr callback pointer of each updater being managed, and initializes various iterators that make it
|	easy to iterate through all edge length parameters, all non-edge-length parameters, or all moves. 
|
|	In the case of MCMCMC (Metropolis-coupled Markov chain Monte Carlo), this class encapsulates one of several coupled
|	Markov chains that run independently except when swaps occur. This functionality has not yet been implemented, 
|	however.
|	
|	In the case of partitioned analyses, it is anticipated that parameters and moves will be knowledgeable about the
|	particular data subset to which they pertain. Thus, in partitioned analyses, there will simply be more parameter
|	updaters (and perhaps more moves) added to the `all_updaters' vector. This functionality has also not yet been
|	implemented.
|	
|	The boost::enable_shared_from_this base class adds the capability to obtain a shared_ptr from this (see 
|	MCMCChainManager::finalize), which is needed because we do not have access to the shared_ptr that Boost Python uses 
|	to manage the MCMCChainManager object created on the Python end.
*/
class MCMCChainManager : public MCMCChainManagerThisShPtr
	{
	public:
								MCMCChainManager();
								~MCMCChainManager();
								
		void 					releaseUpdaters();


		void					addMove(MCMCUpdaterShPtr p);
		void					addModelParam(MCMCUpdaterShPtr p);
		void					addEdgeLenParam(MCMCUpdaterShPtr p);
		void					addEdgeLenHyperparam(MCMCUpdaterShPtr p);
        void                    setEdgeLenHyperparam(unsigned i, double v);

		void					finalize();

		double					recalcEdgeLenPriors(double mu, double var);
		void					refreshLastLnLike();
		void					refreshLastLnPrior();
		
#if USING_EDGE_SPECIFIC_WORKING_PRIORS
        double                  calcExternalEdgeLenWorkingPrior(const TreeNode & nd, double v) const;
        double                  calcInternalEdgeLenWorkingPrior(const TreeNode & nd, double v) const;
#else
        double                  calcExternalEdgeLenWorkingPrior(double v) const;
        double                  calcInternalEdgeLenWorkingPrior(double v) const;
#endif
		
		double 					recalcLnWorkingPrior() const;

		double					getLastLnPrior() const;
		void					setLastLnPrior(double ln_prior);

        double                  calcExternalEdgeLenPriorUnnorm(double v) const;
        double                  calcInternalEdgeLenPriorUnnorm(double v) const;
		//double					partialEdgeLenPrior(const std::vector<double> & edge_len_vect) const;
		double					calcJointLnPrior();

		double					getLastLnLike() const;
		void					setLastLnLike(double ln_like);

		void					clear();
		
		void					updateAllUpdaters();

		const MCMCUpdaterVect &	getAllUpdaters() const;

		const MCMCUpdaterVect &	getMoves() const;
		const MCMCUpdaterVect &	getModelParams() const;
		const MCMCUpdaterVect &	getEdgeLenParams() const;

        MCMCUpdaterVect         getEdgeLenHyperparams() const;

        void					addMCMCUpdaters(ModelShPtr m, TreeShPtr t, TreeLikeShPtr like, LotShPtr r, unsigned max_units, unsigned weight, int subset_pos);
        
        void					debugUpdaterReport(std::string s);

#if POLPY_NEWWAY
		void 					initSAMC(const double_vect_t & elevels, double t0, double eta);
		unsigned 				calcRFDistance(TreeShPtr ref_tree) const;
		unsigned				getSAMCRobinsonFouldsBest() const;
		unsigned 				getSAMCEnergyLevel(double logf)	const;
		double 					getSAMCWeight(unsigned i) const;
		void 					updateSAMCGain(unsigned t);
		void 					updateSAMCWeights(unsigned i);
		bool					doingSAMC() const;
		double					getSAMCGain() const;
		double					getSAMCBest() const;
		double 					getSAMCPiRMSE() const;
		bool 					getSAMCLikelihoodOnly() const;
		void					setSAMCBest(double new_best);
		void 					setSAMCRefTree(TreeShPtr ref_tree);
		void 					setSAMCLikelihoodOnly(bool likelihood_only_status);
		uint_vect_t				getSAMCCounts() const;
		double_vect_t			getSAMCLogLevels() const;
		double_vect_t 			getSAMCWeights() const;
		double_vect_t 			getSAMCPi() const;
#endif

	public:
		// magic workaround for avoiding bogus error when compiling with darwin GCC 3.3 build 1640: 
		// error: cannot dynamic_cast `p' (of type `class phycas::MCMCChainManager*') to type `void*' (source type is not polymorphic)
		// see "Does Boost.Python work with Mac OS X?" entry in Boost Python FAQ http://www.boost.org/libs/python/doc/v2/faq.html
#		if defined(__MACH__) && defined(__APPLE_CC__) && (__APPLE_CC__ >= 1493 && __APPLE_CC__ <= 1666)
			bool dummy_;
#		endif

	protected:

		double					last_ln_like;			/**< The current log likelihood of the chain */
		double					last_ln_prior;			/**< The current log joint prior of the chain */

	private:

#if POLPY_NEWWAY
		// SAMC-related 
		bool 					doing_samc;				/**< If true, SAMC analysis will be performed */
		bool 					samc_likelihood_only;	/**< If true, SAMC analysis will be based on the normalized likelihood (improper prior) rather than the usual posterior (proper prior) */
		double 					samc_best;				/**< The best log posterior that SAMC has seen so far */
		unsigned				samc_rfbest;			/**< The smallest Robinson-Foulds distance (between samc_ref_tree and a tree that arises during MCMC) that SAMC has seen so far */
		TreeShPtr				samc_ref_tree;			/**< The best tree known for the data set being analyzed (SAMC should find this tree if it is working correctly). Ok to leave set to nothing. */
		uint_vect_t 			samc_count;				/**< Vector of counts of the number of times each energy level has been sampled in a SAMC analysis */
		double_vect_t			samc_loglevels;			/**< Vector of energy level boundaries (log scale) sorted from smallest (element 0) to largest (last element) */
		double_vect_t			samc_theta;				/**< Vector of weights associated with each energy level in a SAMC analysis (these weights affect the acceptance probability of Metropolis moves) */
		double_vect_t			samc_pi;				/**< Vector of frequencies associated with each energy level (SAMC analysis is expected to visit level i in proportion to samc_pi[i]) */
		double					samc_gain;				/**< The factor determining the magnitude of the adjustment made to weights in samc_theta after each update */
		double					samc_t0;				/**< The gain at time t is samc_gain = (samc_t0/max(samc_t0,t))^samc_eta */
		double					samc_eta;				/**< The gain at time t is samc_gain = (samc_t0/max(samc_t0,t))^samc_eta */
#endif
		bool					dirty;					/**< If true, means just constructed or at least one updater has been added since finalize() was last called; false means ready to use the all_updaters vector */

		MCMCUpdaterIter			moves_begin;			/**< Iterator positioned at the first move */
		MCMCUpdaterIter			moves_end;				/**< Iterator positioned just beyond the last move */

		MCMCUpdaterIter			edgelens_begin;			/**< Iterator positioned at the first edge length parameter */
		MCMCUpdaterIter			edgelens_end;			/**< Iterator positioned just beyond the last edge length parameter */

		MCMCUpdaterIter			hyperparams_begin;		/**< Iterator positioned at the first edge length hyperparameter (equals hyperparams_end() otherwise) */
		MCMCUpdaterIter			hyperparams_end;		/**< Iterator positioned just beyond the last edge length hyperparameter if any exist (equals hyperparams_begin otherwise) */

		MCMCUpdaterIter			model_params_begin;		/**< Iterator positioned at the first non-edge-length model parameter */
		MCMCUpdaterIter			model_params_end;		/**< Iterator positioned just beyond the last non-edge-length model parameter */

		MCMCUpdaterIter			params_begin;			/**< Iterator positioned at the first parameter (useful for iterating over all non-move updaters) */
		MCMCUpdaterIter			params_end;				/**< Iterator positioned just beyond the last parameter (useful for iterating over all non-move updaters) */

		MCMCUpdaterVect			moves;					/**< Vector of moves */
		MCMCUpdaterVect			model_params;			/**< Vector of model-specific parameters (e.g. kappa) */
		MCMCUpdaterVect			edge_len_params;		/**< Vector of edge length parameters */
		MCMCUpdaterVect		    edge_len_hyperparams;	/**< Vector of edge length hyperparameters */
		MCMCUpdaterVect			all_updaters;			/**< Vector of all updaters (includes moves, model_params, edge_len_params and edge_len_hyperparams); empty until finalize() is called, after which it contains a copy of everything in moves, model_params, edge_len_params and edge_len_hyperparam */
	};

typedef boost::shared_ptr<MCMCChainManager>		ChainManagerShPtr;
typedef boost::weak_ptr<MCMCChainManager>		ChainManagerWkPtr;

double calcEdgeLenLnPrior(const TreeNode &x, double edge_len, ChainManagerShPtr & chain_mgr);


} // namespace phycas

//#include "phycas/src/mcmc_chain_manager.inl"

#endif
