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

#if ! defined(TREE_LIKELIHOOD_HPP)
#define TREE_LIKELIHOOD_HPP

#include "phycas/src/states_patterns.hpp"
#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include "phycas/src/cipres/AllocateMatrix.hpp"
#include "phycas/src/cipres/ConfigDependentHeaders.h"	// int8_t typedef
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/cond_likelihood.hpp"
#include "phycas/src/cond_likelihood_storage.hpp"
#include "phycas/src/underflow_policy.hpp"

#if POLPY_NEWWAY
#include "phycas/src/square_matrix.hpp"
#endif

struct CIPRES_Matrix;

namespace CipresNative
{
class DiscreteMatrix;
}
	
namespace phycas
{
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

#if defined(INTERFACE_WITH_CIPRES)
class PhycasIDLishMatrix
	{
	friend class TreeLikelihood;

	public:
		typedef std::vector<CIPR_State_t> CIPRStateVect;

		PhycasIDLishMatrix(unsigned nt, unsigned nc, std::string s, int dt, unsigned ns)
		  : ntaxa(nt), nchar(nc), symbols(s), datatype(dt), nstates(ns)
			{
			state_lookup.resize(symbols.length());
			matrix.resize(ntaxa);
			}

		void setCharStateLookup(unsigned n, const CIPRStateVect & state_codes)
			{
			if (n >= state_lookup.size())
				state_lookup.resize(n + 1);
			state_lookup[n] = state_codes;
			}

		void replaceRow(unsigned n, const CIPRStateVect & row)
			{
			PHYCAS_ASSERT(row.size() == nchar);
			PHYCAS_ASSERT(n < matrix.size());
			matrix[n] = row;
			}

	private:
		unsigned ntaxa;
		unsigned nchar;
		std::string symbols;
		int datatype;
		unsigned nstates;
		std::vector<CIPRStateVect> matrix;
		std::vector<CIPRStateVect> state_lookup;
	};
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Used for computing the likelihood on a tree.
*/
class TreeLikelihood
	{
	public:

										TreeLikelihood(ModelShPtr);
								virtual ~TreeLikelihood(); //needed as long as startTreeViewer is virtual

		// Accessors
		unsigned						getNTaxa() const;
		unsigned						getNPatterns() const;
		unsigned						getNRatesTotal() const;
		unsigned						getNStates() const;
		ModelShPtr						getModel() const;
		const VecStateList &			getStateList() const;
		const VecStateListPos &			getStateListPos() const;
		const std::vector<double> &		getRateMeans() const;
		const std::vector<double> &		getRateProbs() const;
		void							recalcRelativeRates();
		std::vector<double>				getCategoryLowerBoundaries() const;

		// Modifiers
		void							setNPatterns(unsigned npatterns);
		void							replaceModel(ModelShPtr);
		void							setNoData();
		void							setHaveData();

		// Utilities
        void                            addOrphanTip(TreeShPtr t, unsigned row, std::string name);
        void                            addDecoratedInternalNode(TreeShPtr t, unsigned num = UINT_MAX);
        void							prepareForSimulation(TreeShPtr);
		void							prepareForLikelihood(TreeShPtr);
		void							prepareInternalNodeForLikelihood(TreeNode * nd);

#if defined(INTERFACE_WITH_CIPRES)
		void							copyDataFromIDLMatrix(const PhycasIDLishMatrix & m);
#endif
		void							copyDataFromDiscreteMatrix(const CipresNative::DiscreteMatrix &);
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

		const CondLikelihoodStorage &	getCLAStorage() const;
		unsigned						bytesPerCLA() const;
		unsigned						numCLAsCreated() const;
		unsigned						numCLAsStored() const;

		TreeNode *						storeAllCLAs(TreeShPtr t);
		bool							debugCheckCLAsRemainInTree(TreeShPtr t) const;
        bool                            debugCheckForUncachedCLAs(TreeShPtr t) const;
        bool                            debugCheckCLAsRemainInNode(TreeNode * nd) const;

		bool							isValid(const TreeNode *focal, const TreeNode *avoidNd);
		void 							refreshCLA(TreeNode & nd, const TreeNode * avoid);
		double							calcLnLFromNode(TreeNode & focal_node);
		double							calcLnL(TreeShPtr);

		std::string						listPatterns(bool translate);
		std::string						getStateStr(int8_t state) const;

        void                            setDebug(bool on);

		void							simulateFirst(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar);
		void							simulate(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar);

#if POLPY_NEWWAY
		void							useUnimap(bool yes_or_no = true);
        bool                            isUsingUnimap();
        unsigned                        sampleM(int8_t start_state, int8_t end_state, double transition_prob, double edgelen, LotShPtr rng);
        void                            unimapEdgeOneSite(StateTimeList & state_time_vect, int8_t start_state, int8_t end_state, double transition_prob, double edgelen, LotShPtr rng);
        void                            nielsenMapping(TreeShPtr t, LotShPtr rng);
        double                          calcUnimapLnL(TreeShPtr t);
        void                            recalcUMatVect();
        void                            recalcSMatrix(TreeShPtr t);
        std::string                     debugShowSMatrix();
        void                            slideNode(double fraction, TreeNode * slider, TreeNode * other);
        void                            copyStateTimeListVect(TreeNode * nd, StateTimeListVect & stcopy);
        void                            revertStateTimeListVect(TreeNode * nd, StateTimeListVect & stcopy);
        void                            swapInternalDataAndEdgeLen(TreeNode * nd1, TreeNode * nd2);
#endif
		void							addDataTo(SimData & other);

		unsigned						getNEvals();
		void							resetNEvals();

		TreeNode * 						getLikelihoodRoot();
		void							useAsLikelihoodRoot(TreeNode * nd);

		int								getLikelihoodRootNodeNum() const;
		void							debugSaveCLAs(TreeShPtr t, std::string fn, bool overwrite);
		virtual int						startTreeViewer(TreeShPtr, std::string, unsigned site = 0) const {return 0;}

		void							setUFNumEdges(unsigned nedges);

		void							calcPMatTranspose(double * * *, const StateListPos &, double);
		void							calcPMat(double * * *, double); //

		void							calcCLATwoTips(CondLikelihood &, const TipData &, const TipData &);
		void							calcCLAOneTip(CondLikelihood &, const TipData &, ConstPMatrices, const CondLikelihood &);
		void							calcCLANoTips(CondLikelihood &, ConstPMatrices, const CondLikelihood &, ConstPMatrices, const CondLikelihood &);
  
		void							conditionOnAdditionalTip(CondLikelihood &, const TipData &);
		void							conditionOnAdditionalInternal(CondLikelihood &, ConstPMatrices , const CondLikelihood &);

		double							harvestLnL(EdgeEndpoints & focalEdge);
		double							harvestLnLFromValidEdge(ConstEdgeEndpoints & focalEdge);

		CondLikelihoodStorage & getCondLikelihoodStorage();
	protected:

		//NaiveUnderflowPolicy			underflow_policy;		/**< The object that takes care of underflow correction when computing likelihood for large trees */
		SimpleUnderflowPolicy			underflow_policy;		/**< The object that takes care of underflow correction when computing likelihood for large trees */

		TreeNode *						likelihood_root;		/**< If not NULL< calcLnL will use this node as the likelihood root, then reset it to NULL before returning */
		CondLikelihoodStorage			cla_pool;				/**< Stores currently unused CondLikelihood objects */

		bool							store_site_likes;		/**< If true, calcLnL always stores the site likelihoods in the `site_likelihood' data member; if false, the `site_likelihood' data member is not updated by calcLnL */
		bool							no_data;				/**< If true, calcLnL always returns 0.0 (useful for allowing MCMC to explore the prior) */

		unsigned						nTaxa;					/**< The number of taxa (i.e. the number of elements in each pattern stored in pattern_map) */
		unsigned						num_patterns;			/**< The number of site patterns */
		unsigned						num_states;				/**< The number of states */
		unsigned						num_rates;				/**< The number of relative rate categories */
		ModelShPtr						model;					/**< The substitution model */
		VecStateList					state_list;				/**< The global lookup table for decoding coded states */
		VecStateListPos					state_list_pos;			/**< The vector of positions of states in `state_list' */

        //POL_BOOKMARK
		std::vector<double>				rate_means;				/**< Vector of relative rates */
		std::vector<double>				rate_probs;				/**< Vector of relative rate probabilities */

		unsigned						nevals;					/**< For debugging, records the number of times calcLnL() is called */

        bool                            debugging_now;          /**< For debugging, indicates whether user wants to see debugging output */

	protected:

		// Utilities
		void							buildPatternReprVector(std::vector<std::string> &, TreeShPtr);
		TipData *						allocateTipData(unsigned);
		InternalData *					allocateInternalData();

#if defined(INTERFACE_WITH_CIPRES)
		unsigned						compressDataMatrix(unsigned ntax, unsigned nchar, const CipresNative::DiscreteMatrix &);
		unsigned						compressIDLMatrix(unsigned ntax, unsigned nchar, const CIPR_StateSet_t * const * matrix);
		unsigned						buildPatternMapFromRawMatrix(unsigned ntax, unsigned nchar, const CIPR_StateSet_t * const * mat);
#else
		unsigned						compressDataMatrix(const CipresNative::DiscreteMatrix &);
#endif

		void							calcTMatForSim(TipData &, double);
		void							simulateImpl(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar, bool refresh_probs);
		void							calcPMatCommon(double * * *, double);


    protected:

#if POLPY_NEWWAY
        bool                            using_unimap;           /**< if true, uniformized mapping likelihoods will be used; if false, Felsenstein-style integrated likelihoods will be used */
        std::vector<SquareMatrix>       uMatVect;               /**< uMat[m][i][j] is the (i,j)th element of the matrix obtained by raising the uniformized transition matrix to the power m */               
        double * *                      uMat;                   /**< uMat[i][j] is the log of the corresponding element of the uniformized transition matrix element */
        double                          lambda;                 /**< rate at which uniformization events, or univents, occur */
        unsigned * *                    sMat;                   /**< sMat[i][j] is total number of univents in which state i changes to state j */
        unsigned                        nunivents;              /**< total number of univents over all edges and all sites */
        unsigned                        maxm;                   /**< the maximum number of univents on any edge for any site (determines length of logmfact vector) */
        std::vector<double>             logmfact;               /**< contains the log of the factorial of m for m = 0..maxm */
        std::vector<unsigned>           obs_state_freqs;        /**< vector of length num_states holding the frequency with which each state was observed in the subroot node (used in computing the unimap likelihood) */
#endif

	public: //@POL these should be protected rather than public

		std::vector<double>				likelihood_rate_site;	/**< Vector of likelihoods for each rate/site combination */
		CountVectorType					pattern_counts;			/**< vector of pattern counts */
		PatternMapType					pattern_map;			/**< keys are patterns, values are pattern counts */
		std::vector<double>				site_likelihood;		/**< site_likelihood[pat] stores the site likelihood for pattern pat, but only if `store_site_likes' is true */
		std::map<unsigned, unsigned>	charIndexToPatternIndex; /**< maps original character index to the index in compressed pattern "matrix" */
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
