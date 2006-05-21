#if ! defined(TREE_LIKELIHOOD_HPP)
#define TREE_LIKELIHOOD_HPP

#include "pyphy/common/states_patterns.hpp"
#include <vector>
#include <cassert>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include "CipresCommlib/AllocateMatrix.hpp"
#include "CipresCommlib/ConfigDependentHeaders.h"	// int8_t typedef
#include "pyphy/likelihood/likelihood_models.hpp"
#include "pyphy/likelihood/cond_likelihood.hpp"

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

/*----------------------------------------------------------------------------------------------------------------------
|	Used for computing the likelihood on a tree.
*/
class TreeLikelihood
	{
	//friend class CalcTransitionMatrixForOneNode;

	public:

										TreeLikelihood(ModelShPtr);
										//~TreeLikelihood() 
										//	{
										//	std::cerr << "TreeLikelihood dying..." << std::endl;
										//	}

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
		void							prepareForSimulation(TreeShPtr);
		void							prepareForLikelihood(TreeShPtr);
		void							prepareInternalNodeForLikelihood(TreeNode * nd);
		void							copyDataFromDiscreteMatrix(const CipresNative::DiscreteMatrix &);
		void							copyDataFromSimData(SimDataShPtr sim_data);

#if POLPY_NEWWAY
		bool							invalidateNode(const TreeNode * refNd, const TreeNode * neighborCloserToEffectiveRoot);
		void							invalidateAwayFromNode(TreeNode & focalNode);
		bool							isValid(const TreeNode *focal, const TreeNode *avoidNd);
		void 							refreshCLA(TreeNode & nd, const TreeNode * avoid);
		double							calcLnLFromNode(TreeNode & focal_node);
		double							calcLnL(TreeShPtr);
#else
		double							calcLnL(TreeShPtr);
#endif
		std::string						listPatterns(bool translate);
		std::string						getStateStr(int8_t state) const;

		void							simulateFirst(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar);
		void							simulate(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar);

		void							addDataTo(SimData & other);

		unsigned						getNEvals();
		void							resetNEvals();

		
	protected:

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

		std::vector<double>				rate_means;				/**< Vector of relative rates */
		std::vector<double>				rate_probs;				/**< Vector of relative rate probabilities */

		unsigned						nevals;					/**> For debugging, records the number of times calcLnL() is called */

	protected:

		// Utilities
		void							buildPatternReprVector(std::vector<std::string> &, TreeShPtr);
		TipData *						allocateTipData(unsigned);
		InternalData *					allocateInternalData();
		unsigned						compressDataMatrix(const CipresNative::DiscreteMatrix &);

		void							calcTMatForSim(TipData &, double);
		void							simulateImpl(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar, bool refresh_probs);
		void							calcPMatTranspose(double * * *, const StateListPos &, double);
		void							calcPMat(double * * *, double); //
		void							calcPMatCommon(double * * *, double);

		void							calcCLATwoTips(CondLikelihood &, const TipData &, const TipData &);
		void							calcCLAOneTip(CondLikelihood &, const TipData &, ConstPMatrices, const CondLikelihood &);
		void							calcCLANoTips(CondLikelihood &, ConstPMatrices, const CondLikelihood &, ConstPMatrices, const CondLikelihood &);
  
		void							conditionOnAdditionalTip(CondLikelihood &, const TipData &);
		void							conditionOnAdditionalInternal(CondLikelihood &, ConstPMatrices , const CondLikelihood &);

		double							harvestLnL(EdgeEndpoints & focalEdge);
		double							harvestLnLFromValidEdge(ConstEdgeEndpoints & focalEdge);

	public: //@POL these should be protected rather than public

		std::vector<double>				likelihood_rate_site;	/**< Vector of likelihoods for each rate/site combination */
		CountVectorType					pattern_counts;			/**< vector of pattern counts */
		PatternMapType					pattern_map;			/**< keys are patterns, values are pattern counts */
		std::vector<double>				site_likelihood;		/**< site_likelihood[pat] stores the site likelihood for pattern pat, but only if `store_site_likes' is true */
		std::map<unsigned, unsigned>	charIndexToPatternIndex; /**< maps original character index to the index in compressed pattern "matrix" */
	};

CondLikelihood * getCondLike(EdgeEndpoints edge);
CondLikelihood * getCondLike(TreeNode *focalNd, TreeNode *avoid);
const CondLikelihood * getCondLike(ConstEdgeEndpoints edge);
const CondLikelihood * getCondLike(const TreeNode *focalNd, const TreeNode *avoid);

} // namespace phycas

#include "tree_likelihood.inl"

#endif
