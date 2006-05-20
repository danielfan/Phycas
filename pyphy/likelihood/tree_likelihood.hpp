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

struct CIPRES_Matrix;

namespace CipresNative
{
class DiscreteMatrix;
}
	
namespace phycas
{

class Tree;
class EdgeEndpoints;
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
		void							calcPMatTranspose(TipData &, double);
		void							calcPMat(InternalData &, double);
		void							calcPMatCommon(double * * *, double);

		void							calcCLATwoTips(InternalData &, const TipData &, const TipData &);
		void							calcCLAOneTip(InternalData &, const TipData &, const InternalData &);
		void							calcCLANoTips(InternalData &, const InternalData &, const InternalData &);

		void							conditionOnAdditionalTip(InternalData &, const TipData &);
		void							conditionOnAdditionalInternal(InternalData &, const InternalData &);

		double							harvestLnL(EdgeEndpoints & focalEdge);
		double							harvestLnLFromValidEdge(EdgeEndpoints & focalEdge);

	public: //@POL these should be protected rather than public

		std::vector<double>				likelihood_rate_site;	/**< Vector of likelihoods for each rate/site combination */
		CountVectorType					pattern_counts;			/**< vector of pattern counts */
		PatternMapType					pattern_map;			/**< keys are patterns, values are pattern counts */
		std::map<unsigned, unsigned>	charIndexToPatternIndex; /*maps original character index to the index in compressed pattern "matrix" */
	};

} // namespace phycas

#include "tree_likelihood.inl"

#endif
