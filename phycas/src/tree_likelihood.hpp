/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
|																			  |
|  This program is free software; you can redistribute it and/or modify		  |
|  it under the terms of the GNU General Public License as published by		  |
|  the Free Software Foundation; either version 2 of the License, or		  |
|  (at your option) any later version.										  |
|																			  |
|  This program is distributed in the hope that it will be useful,			  |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of			  |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the			  |
|  GNU General Public License for more details.								  |
|																			  |
|  You should have received a copy of the GNU General Public License along	  |
|  with this program; if not, write to the Free Software Foundation, Inc.,	  |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.				  |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#if ! defined(TREE_LIKELIHOOD_HPP)
#define TREE_LIKELIHOOD_HPP

#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include "ncl/nxscxxdiscretematrix.h"

#include "phycas/src/states_patterns.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/cond_likelihood.hpp"
#include "phycas/src/cond_likelihood_storage.hpp"
#include "phycas/src/underflow_manager.hpp"
#include "phycas/src/univent_prob_mgr.hpp"

#if POLPY_NEWWAY
#include "phycas/src/partition_model.hpp"
#endif

namespace phycas
{
unsigned ** getNodeSMat(TreeNode * nd);

typedef const double * const * const * ConstPMatrices;
typedef std::vector<unsigned int> StateListPos;
class CondLikelihood;
class Tree;

template<typename T> class GenericEdgeEndpoints;
typedef GenericEdgeEndpoints<TreeNode *> EdgeEndpoints;
typedef GenericEdgeEndpoints<const TreeNode *> ConstEdgeEndpoints;

class TipData;
class InternalData;

class SimData;
typedef boost::shared_ptr<SimData>	SimDataShPtr;

class Lot;
typedef boost::shared_ptr<Lot>	LotShPtr;

class Univents;
Univents & getUniventsRef(TreeNode &);
const Univents & getUniventsConstRef(const TreeNode &);

/*----------------------------------------------------------------------------------------------------------------------
|	Used for computing the likelihood on a tree.
*/
class TreeLikelihood
	{
	public:

#if POLPY_NEWWAY
										TreeLikelihood(PartitionModelShPtr);
#else
										TreeLikelihood(ModelShPtr);
#endif
								virtual ~TreeLikelihood(); //needed as long as startTreeViewer is virtual

		// Accessors
		unsigned						getNTaxa() const;
#if POLPY_NEWWAY
		PartitionModelShPtr				getPartitionModel() const;
		const double_vect_t &			getRateMeans(unsigned i) const;
		const double_vect_t &			getRateProbs(unsigned i) const;
		unsigned						getNRatesTotal(unsigned i) const;
		double_vect_t					getCategoryLowerBoundaries(unsigned i) const;
		unsigned						getNStates(unsigned i) const;
		const state_list_vect_t &		getStateList() const;
		const state_list_pos_vect_t &	getStateListPos() const;
		void							replacePartitionModel(PartitionModelShPtr);
		const count_vect_t &			getPatternCounts() const;
#else
		unsigned						getNPatterns() const;
		const std::vector<double> &		getRateMeans() const;
		const std::vector<double> &		getRateProbs() const;
		unsigned						getNRatesTotal() const;
		std::vector<double>				getCategoryLowerBoundaries() const;
		unsigned						getNStates() const;
		ModelShPtr						getModel() const;
		const VecStateList &			getStateList() const;
		const StateListPos &			getStateListPos() const;
		void							replaceModel(ModelShPtr);
		const CountVectorType &			getPatternCounts() const;
#endif
		void							recalcRelativeRates();
		const std::vector<unsigned> &	getListOfAllMissingSites() const;
		const std::vector<double> &		getSiteLikelihoods() const;
		const std::vector<double> &		getSiteUF() const;
		bool							storingSiteLikelihoods() const;
		const std::vector<unsigned> &	getCharIndexToPatternIndex() const;

		// Modifiers
#if POLPY_NEWWAY
#else // old way
		void							setNPatterns(unsigned npatterns);
#endif
		void							setNoData();
		void							setHaveData();
		void							storeSiteLikelihoods(bool yes);

		// Utilities
#if POLPY_NEWWAY
		void							releasePartitionModel();
#else
		void							releaseModel();
#endif
		double							calcLogLikeAtSubstitutionSaturation() const;

		void							addOrphanTip(TreeShPtr t, unsigned row, std::string name);
		void							addDecoratedInternalNode(TreeShPtr t, unsigned num = UINT_MAX);
		void							prepareForSimulation(TreeShPtr);
		void							prepareForLikelihood(TreeShPtr);
		void							prepareInternalNodeForLikelihood(TreeNode * nd);

#if POLPY_NEWWAY    // partitioning data
		void							copyDataFromDiscreteMatrix(const NxsCXXDiscreteMatrix &, const std::vector<unsigned> & partition_info);
#else
		void							copyDataFromDiscreteMatrix(const NxsCXXDiscreteMatrix &, const std::vector<unsigned> & partition_info);
		//obsolete
		//void							copyDataFromDiscreteMatrix(const NxsCXXDiscreteMatrix &);
#endif
		void							copyDataFromSimData(SimDataShPtr sim_data);

		bool							invalidateNode(TreeNode * ref_nd, TreeNode * neighbor_closer_to_likelihood_root);
		bool							invalidateBothEnds(TreeNode * ref_nd, TreeNode * unused = NULL);
		bool							invalidateBothEndsDiscardCache(TreeNode * ref_nd, TreeNode * unused = NULL);
		void							invalidateAwayFromNode(TreeNode & focalNode);

		bool							restoreFromCacheNode(TreeNode * ref_nd, TreeNode * neighbor_closer_to_likelihood_root);
		bool							restoreFromCacheBothEnds(TreeNode * ref_nd, TreeNode * unused = NULL);
		bool							restoreFromCacheParentalOnly(TreeNode * ref_nd, TreeNode * unused = NULL);
		bool							discardCacheBothEnds(TreeNode * ref_nd, TreeNode * unused = NULL);
		void							restoreFromCacheAwayFromNode(TreeNode & focalNode);
		void							discardCacheAwayFromNode(TreeNode & focalNode);

		const CondLikelihoodStorageShPtr	getCLAStorage() const;
		
		unsigned						bytesPerCLA() const;
		unsigned						numCLAsCreated() const;
		unsigned						numCLAsStored() const;

		TreeNode *						storeAllCLAs(TreeShPtr t);
		bool							debugCheckCLAsRemainInTree(TreeShPtr t) const;
		bool							debugCheckForUncachedCLAs(TreeShPtr t) const;
		bool							debugCheckCLAsRemainInNode(TreeNode * nd) const;

		bool							isValid(const TreeNode *focal, const TreeNode *avoidNd);
		void							refreshCLA(TreeNode & nd, const TreeNode * avoid);
		double							calcLnLFromNode(TreeNode & focal_node);
		double							calcLnL(TreeShPtr);

		std::string						listPatterns(bool translate);
#if POLPY_NEWWAY
		std::string						getStateStr(unsigned i, state_code_t state) const;
#else // old way
		std::string						getStateStr(int8_t state) const;
#endif

		void							setDebug(bool on);

		void							simulateFirst(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar);
		void							simulate(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar);

		const UniventProbMgr &          GetUniventProbMgrConstRef() const {return univentProbMgr;}
		UniventProbMgr &                GetUniventProbMgrRef() {return univentProbMgr;}
		void							useUnimap(bool yes_or_no = true);
		bool							isUsingUnimap();
		void							fullRemapping(TreeShPtr t, LotShPtr rng, bool doSampleUnivents);
		void							debugCheckSMatrix(TreeShPtr t);
		std::string						debugShowSMatrix();
		void							slideNode(double fraction, TreeNode * slider, TreeNode * other);
		void							swapInternalDataAndEdgeLen(TreeNode * nd1, TreeNode * nd2);
		void							addDataTo(SimData & other);

		unsigned						getNumLikelihoodEvals() const;
		void							incrementNumLikelihoodEvals();
		void							resetNumLikelihoodEvals();

		TreeNode *						getLikelihoodRoot();
		void							useAsLikelihoodRoot(TreeNode * nd);

		int								getLikelihoodRootNodeNum() const;
		void							debugSaveCLAs(TreeShPtr t, std::string fn, bool overwrite);
		virtual int						startTreeViewer(TreeShPtr, std::string, unsigned site = 0) const {return 0;}

		void							setUFNumEdges(unsigned nedges);

#if POLPY_NEWWAY
		void							calcPMatTranspose(unsigned i, double * * * transPMats, const uint_vect_t & stateListPosVec, double edgeLength);
		void							calcPMat(unsigned i, double * * * p, double edgeLength); //

		void							calcCLATwoTips(CondLikelihood & condLike, const TipData & leftTip, const TipData & rightTip);
		void							calcCLAOneTip(CondLikelihood & condLike, const TipData & leftChild, const InternalData & rightChild, const CondLikelihood & rightCondLike);
		void							calcCLANoTips(CondLikelihood & condLike, const InternalData & leftChild, const CondLikelihood & leftCondLike, const InternalData & rightChild, const CondLikelihood & rightCondLike);
  
		void							conditionOnAdditionalTip(CondLikelihood & condLike, const TipData & tipData);
		void							conditionOnAdditionalInternal(CondLikelihood & condLike, const InternalData & child, const CondLikelihood & childCondLike);
#else
		void							calcPMatTranspose(double * * *, const StateListPos &, double);
		void							calcPMat(double * * *, double); //

		void							calcCLATwoTips(CondLikelihood &, const TipData &, const TipData &);
		void							calcCLAOneTip(CondLikelihood &, const TipData &, ConstPMatrices, const CondLikelihood &);
		void							calcCLANoTips(CondLikelihood &, ConstPMatrices, const CondLikelihood &, ConstPMatrices, const CondLikelihood &);
  
		void							conditionOnAdditionalTip(CondLikelihood &, const TipData &);
		void							conditionOnAdditionalInternal(CondLikelihood &, ConstPMatrices , const CondLikelihood &);
#endif

		double							harvestLnL(EdgeEndpoints & focalEdge);
		double							harvestLnLFromValidEdge(ConstEdgeEndpoints & focalEdge);

		unsigned						buildConstantStatesVector();

		CondLikelihoodStorageShPtr		getCondLikelihoodStorage();
		void							flagNodeWithInvalidUnivents(TreeNode *nd)
		{
			invalidUniventMappingNodes.insert(nd);
		}
		void                            setLot(LotShPtr r)
		{
		    localRng = r;
		}

		void							remapUniventsForNode(TreeShPtr, TreeNode *);

	protected:

		UnderflowManager				underflow_manager;		/**< The object that takes care of underflow correction when computing likelihood for large trees */

		TreeNode *						likelihood_root;		/**< If not NULL< calcLnL will use this node as the likelihood root, then reset it to NULL before returning */
		CondLikelihoodStorageShPtr		cla_pool;

		bool							store_site_likes;		/**< If true, calcLnL always stores the site likelihoods in the `site_likelihood' data member; if false, the `site_likelihood' data member is not updated by calcLnL */
		bool							no_data;				/**< If true, calcLnL always returns 0.0 (useful for allowing MCMC to explore the prior) */

		unsigned						nTaxa;					/**< The number of taxa */
#if POLPY_NEWWAY
		PartitionModelShPtr				partition_model;		/**< The object that holds information about the model applied to each subset of the data partition */
		
		//@POL these next four should logically be inside the PartitionModel class
		double_vect_vect_t				rate_means;				/**< Vector of relative rate vectors (one rate vector for each partition subset) */ 
		double_vect_vect_t				rate_probs;				/**< Vector of relative rate probability vectors  (one probability vector for each partition subset) */
		state_list_vect_t				state_list;				/**< `state_list[i]' provides a vector of state code definitions for subset i */
		state_list_pos_vect_t			state_list_pos;			/**< `state_list_pos[i]' is a vector of positions of states in `state_list[i]' */
#else
		unsigned						num_patterns;			/**< The number of site patterns */
		unsigned						num_states;				/**< The number of states */
		unsigned						num_rates;				/**< The number of relative rate categories */
		ModelShPtr						model;					/**< The substitution model */
		std::vector<double>				rate_means;				/**< Vector of relative rates */ //POL_BOOKMARK rate_means declaration
		std::vector<double>				rate_probs;				/**< Vector of relative rate probabilities */
		VecStateList					state_list;				/**< The global lookup table for decoding coded states */
		StateListPos					state_list_pos;			/**< The vector of positions of states in `state_list' */
#endif

		bool							debugging_now;			/**< For debugging, indicates whether user wants to see debugging output */

	protected:

		// Utilities
		void							buildPatternReprVector(std::vector<std::string> &, TreeShPtr);
		TipData *						allocateTipData(unsigned);
		InternalData *					allocateInternalData();

#if POLPY_NEWWAY    // partitioning data
		void	 						storePattern(pattern_map_t & pattern_map, pattern_to_sites_map_t & pattern_to_site_map, const std::vector<int8_t> & pattern, const unsigned pattern_index, const pattern_count_t weight);
		unsigned						compressDataMatrix(const NxsCXXDiscreteMatrix &, const std::vector<unsigned> & partition_info);
		void							calcPMatCommon(unsigned i, double * * * pMatrices, double edgeLength);
#else
		void	 						storePattern(const std::vector<int8_t> & pattern, const unsigned pattern_index, const PatternCountType weight);
		unsigned						compressDataMatrix(const NxsCXXDiscreteMatrix &);
		void							calcPMatCommon(double * * *, double);
#endif

		void							calcTMatForSim(TipData &, double);
		void							simulateImpl(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar, bool refresh_probs);

	private:
#if POLPY_NEWWAY
#else
		typedef std::list<unsigned>					IndexList;
		typedef std::map<VecStateList, IndexList>	PatternToIndex;
		PatternToIndex 					patternToIndex;			/**< map that associates a list of character indices with each pattern; for example, if some pattern is found at sites 0, 15, and 167, then patternToIndex.first is the pattern and patternToIndex.second is the list<unsigned> [0, 15, 167] */
#endif

	protected:

		bool							using_unimap;			/**< if true, uniformized mapping likelihoods will be used; if false, Felsenstein-style integrated likelihoods will be used */
		UniventProbMgr					univentProbMgr;
		unsigned * *					treeSMat;					/**< sMat[i][j] is total number of univents in which state i changes to state j */
		unsigned						nunivents;				/**< total number of univents over all edges and all sites */
		std::vector<unsigned>			obs_state_counts;
		std::set<TreeNode *>			invalidUniventMappingNodes;
        LotShPtr                        localRng;
		unsigned						nevals;					/**< For debugging, records the number of times the likelihood is calculated */
        
	public: //@POL these should be protected rather than public

#if POLPY_NEWWAY
		count_vect_t					pattern_counts;				/**< vector of pattern counts */
		uint_vect_t						subset_offset;				/**< `subset_offset'[i] holds the index into `pattern_vect' where the patterns from subset i begin. The number of elements is one greater than the number of subsets (the last element holds num_patterns to make it easy to find the end of any subset, including the last subset). */
		pattern_vect_t					pattern_vect;				/**< vector of patterns (all patterns for a given partition subset are contiguous, but within a subset, order may differ from data file order) */
		double_vect_t					site_likelihood;			/**< site_likelihood[pat] stores the site likelihood for pattern pat, but only if `store_site_likes' is true */
		pattern_to_sites_t				pattern_to_sites;			/**< vector of lists that provides a list of character indices for each pattern in `pattern_vect'. For example, if pattern j is found at sites 0, 15, and 167, then pattern_to_sites[j] is the list [0, 15, 167] */
		uint_vect_t						charIndexToPatternIndex; 	/**< maps original character index to the position of the corresponding element in `pattern_vect' */
		uint_vect_t						constant_states;			/**< keeps track of the states for potentially constant sites. See TreeLikelihood::buildConstantStatesVector for description of the structure of this vector. */
		uint_vect_t						all_missing;				/**< keeps track of sites excluded automatically because they have missing data for all taxa. */
		double_vect_t					site_uf;					/**< site_uf[pat] stores the underflow correction factor used for pattern pat, but only if `store_site_likes' is true */
#else
		//unused? std::vector<double>	likelihood_rate_site;		/**< Vector of likelihoods for each rate/site combination */
		CountVectorType					pattern_counts;				/**< vector of pattern counts */
		PatternMapType					pattern_map;				/**< keys are patterns, values are pattern counts */
		std::vector<double>				site_likelihood;			/**< site_likelihood[pat] stores the site likelihood for pattern pat, but only if `store_site_likes' is true */
		std::vector<unsigned>			charIndexToPatternIndex; 	/**< maps original character index to the position of the corresponding element in `pattern_map' */
		std::vector<unsigned>			constant_states;			/**< keeps track of the states for potentially constant sites. See TreeLikelihood::buildConstantStatesVector for description of the structure of this vector. */
		std::vector<unsigned>			all_missing;				/**< keeps track of sites excluded automatically because they have missing data for all taxa. */
		std::vector<double>				site_uf;					/**< site_uf[pat] stores the underflow correction factor used for pattern pat, but only if `store_site_likes' is true */
#endif


	};

/// used to get access to a CLA to write it
CondLikelihoodShPtr getCondLikePtr(EdgeEndpoints edge);
CondLikelihoodShPtr getCondLikePtr(TreeNode * focal_nd, TreeNode * avoid);

/// read-only access to a CLA (must be valid already)
ConstCondLikelihoodShPtr getValidCondLikePtr(ConstEdgeEndpoints edge);
ConstCondLikelihoodShPtr getValidCondLikePtr(const TreeNode * focal_nd, const TreeNode * avoid);

} // namespace phycas

//#include "tree_likelihood.inl"

#endif
