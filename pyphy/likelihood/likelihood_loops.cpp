#include "phycas/force_include.h"
#include "CipresCommlib/CipresDataMatrixHelper.h"
#include "pyphy/phylogeny/basic_tree.hpp"
#include "pyphy/likelihood/likelihood_models.hpp"
#include "pyphy/likelihood/tree_likelihood.hpp"
#include "pyphy/likelihood/tip_data.hpp"
#include "pyphy/likelihood/internal_data.hpp"
#include "CipresCommlib/util_copy.hpp"
#include <numeric>

#include <iostream>
using std::cerr;
using std::endl;

#include <cmath>
using std::log;

using std::accumulate;
template<typename T>
std::vector<T> scaleVector(T scaler, const std::vector<T> & orig);
template<typename T>
void transpose(T * * mat, unsigned dim);

using std::vector;

template<typename T>
vector<T> scaleVector(T scaler, const vector<T> & orig)
	{
	const unsigned origLen = (const unsigned)orig.size();
	vector<T> scaled(origLen, scaler);
	cip_imult(orig.begin(), orig.end(), scaled.begin());
	return scaled;
	}

template<typename T>
void transpose(T * * mat, unsigned dim)
	{
	for (unsigned i = 0; i < dim; ++i)
		for (unsigned j = i + 1; j < dim; ++j)
			std::swap(mat[i][j], mat[j][i]);
	}

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates a transition matrix for every rate category given the supplied `edgeLength' and the current `model'. This
|	function calculates the basic square transition matrix in which rows represent "from" states and columns represent
|	"to" states. For tips, the transition matrix generated by this function is transposed by both 
|	TreeLikelihood::calcPMatTranspose and TreeLikelihood::calcTMatForSim, and in the case of 
|	TreeLikelihood::calcPMatTranspose also augmented by adding rows corresponding to tip-specific ambiguities observed.
*/
void TreeLikelihood::calcPMatCommon(
  double * * *	pMatrices, 
  double		edgeLength)
	{
	assert(num_rates > 0);
	const vector<double> scaledEdges = scaleVector(edgeLength, rate_means);
	model->calcPMatrices(pMatrices, &scaledEdges[0], num_rates); //PELIGROSO

#if 0
	// debugging code
	if (use_pattern_specific_rates)
		{
		std::cerr << "\n\nHere are all " << rate_means.size() << " rate_means used for computing pMatrices:" << std::endl;
		for (std::vector<double>::const_iterator it = rate_means.begin(); it != rate_means.end(); ++it)
			{
			std::cerr << (*it) << std::endl;
			}
		}
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls calcPMatCommon to compute the transition matrix for an internal node given the value of `edgeLength'.
*/
void TreeLikelihood::calcPMat(
  InternalData &	cla, 
  double			edgeLength)
	{
	calcPMatCommon(cla.getPMatrices(), edgeLength);
	}

#if 0
void DebugShowNuclTransMatrix(double * * m, const char * title)
	{
	std::cerr << title << std::endl;
	std::cerr << str(boost::format(" %12.5f %12.5f %12.5f %12.5f") % m[0][0] % m[0][1] % m[0][2] % m[0][3]) << std::endl;
	std::cerr << str(boost::format(" %12.5f %12.5f %12.5f %12.5f") % m[1][0] % m[1][1] % m[1][2] % m[1][3]) << std::endl;
	std::cerr << str(boost::format(" %12.5f %12.5f %12.5f %12.5f") % m[2][0] % m[2][1] % m[2][2] % m[2][3]) << std::endl;
	std::cerr << str(boost::format(" %12.5f %12.5f %12.5f %12.5f") % m[3][0] % m[3][1] % m[3][2] % m[3][3]) << std::endl;
	std::cerr << std::endl;
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates tip-specific transition matrices for purposes of simulation. These are transposed but not augmented
|	transition matrices. They are not augmented because no ambiguities are generated when simulating data and thus no
|	additional rows are necessary.
*/
void TreeLikelihood::calcTMatForSim(
  TipData &	tipData, 
  double	edgeLength)
	{
	double * * * transPMats = tipData.getTransposedPMatrices();
	calcPMatCommon(transPMats,  edgeLength);

	// For each rate category, transpose the num_states x num_states portion of the matrices
	for (unsigned rate = 0; rate < num_rates; ++rate)
		{
		double * * pMat = transPMats[rate];
		
		//std::cerr << "Edge length: " << edgeLength << std::endl;
		//DebugShowNuclTransMatrix(pMat, "Untransposed");

		transpose(pMat, num_states);

		//DebugShowNuclTransMatrix(pMat, "Transposed");
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates tip-specific transition matrices for purposes of computing likelihoods. The resulting matrices are
|	T matrices: transposed and augmented transition matrices.
*/
void TreeLikelihood::calcPMatTranspose(
  TipData &	tipData, 
  double	edgeLength)
	{
	double * * * transPMats = tipData.getTransposedPMatrices();
	calcPMatCommon(transPMats,  edgeLength);

	// For each rate category, transpose the num_states x num_states portion of the matrices
	// and fill in the ambiguity codes by summing columns
	const std::vector<unsigned int> stateListPosVec(tipData.getConstStateListPos());
	const unsigned nPartialAmbigs = (unsigned)stateListPosVec.size();
	const unsigned int * const stateListPosArr = &stateListPosVec[0]; //PELIGROSO
	const VecStateList stateListVec = state_list;
	const int8_t * const stateListArr = &stateListVec[0]; //PELIGROSO
	for (unsigned rate = 0; rate < num_rates; ++rate)
		{
		double * * pMat = transPMats[rate];
		transpose(pMat, num_states);
		unsigned currPMatRowIndex = num_states;
		std::fill(pMat[currPMatRowIndex], pMat[currPMatRowIndex] + num_states, 1.0);
		++currPMatRowIndex;
		for (unsigned ambigCode = 0; ambigCode < nPartialAmbigs; ++ambigCode, ++currPMatRowIndex)
			{
			unsigned indexIntoStateList = stateListPosArr[ambigCode];
			const unsigned nObservedStates = stateListArr[indexIntoStateList++];
			unsigned currObservedState = stateListArr[indexIntoStateList++];
			cip_copy(pMat[currObservedState], pMat[currObservedState] + num_states, pMat[currPMatRowIndex]);
			for (unsigned m = 1; m < nObservedStates; ++m)
				{
				currObservedState = stateListArr[indexIntoStateList++];
				cip_iadd(pMat[currObservedState], pMat[currObservedState] + num_states, pMat[currPMatRowIndex]);
				}
			}
		}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Computes the conditional likelihood arrays at an internal node subtending two tips.
*/
void TreeLikelihood::calcCLATwoTips(
  InternalData &	internalData, 
  const TipData &	leftTip, 
  const TipData &	rightTip)
	{
	const double * const * const * leftPMatricesTrans = leftTip.getConstTransposedPMatrices();
	const int8_t * leftStateCodes = leftTip.getConstStateCodes();
	const double * const * const * rightPMatricesTrans = rightTip.getConstTransposedPMatrices();
	const int8_t * rightStateCodes = rightTip.getConstStateCodes();
	double * cla = internalData.getCLA(); //PELIGROSO

	for (unsigned r = 0; r < num_rates; ++r)
		{
		const double * const * leftPMatT = leftPMatricesTrans[r];
		const double * const * const rightPMatT = rightPMatricesTrans[r];
		for (unsigned p = 0; p < num_patterns; ++p, cla += num_states)
			{
			const double * leftPMatTRow = leftPMatT[leftStateCodes[p]];
			const double * rightPMatTRow = rightPMatT[rightStateCodes[p]];
			for (unsigned s = 0; s < num_states; ++s)
				cla[s] = leftPMatTRow[s]*rightPMatTRow[s];
			}
		}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Computes the conditional likelihood arrays at an internal node having a tip node as its left child and an internal
|	node for a right child.
*/
void TreeLikelihood::calcCLAOneTip(
  InternalData &		ndInfo, 
  const TipData &		leftChild, 
  const InternalData &	rightChild)
	{
	// Get transition probability matrices for the left and right child nodes of this node
	// These are 3D because there are potentially several 2D transition matrices, one
	// for each relative rate
	const double * const * const * leftPMatricesTrans = leftChild.getConstTransposedPMatrices();
	const double * const * const * rightPMatrices = rightChild.getConstPMatrices();

	// cla is the conditional likelihood array we are updating
	double * cla = ndInfo.getCLA(); //PELIGROSO

	// These are the conditional likelihood arrays of the left and right child nodes of this node
	const int8_t * leftStateCodes = leftChild.getConstStateCodes();
	const double * rightCLA = rightChild.getConstCLA();

	// conditional likelihood arrays are laid out as follows for DNA data:
	//
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	// |   pattern 1   |               |               |               |               |
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	// | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T |
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	//
	for (unsigned r = 0; r < num_rates; ++r)
		{
		const double * const * const leftPMatrixT = leftPMatricesTrans[r];
		const double * const * rightPMatrix = rightPMatrices[r];
		for (unsigned pat = 0; pat < num_patterns; ++pat, rightCLA += num_states)
			{
			const double * leftPMatT_pat = leftPMatrixT[leftStateCodes[pat]];
			for (unsigned i = 0; i < num_states; ++i)
				{
				double right_side = 0.0;
				const double * rightP_i = rightPMatrix[i];
				for (unsigned j = 0; j < num_states; ++j)
					right_side += rightP_i[j]*rightCLA[j];
				*cla++ = (leftPMatT_pat[i]*right_side);
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the conditional likelihood arrays at an internal node having two internal nodes for children.
*/
void TreeLikelihood::calcCLANoTips(
  InternalData &		ndInfo, 
  const InternalData &	leftNode, 
  const InternalData &	rightNode)
	{
	// This function updates the conditional likelihood array of a node assuming that the
	// conditional likelihood arrays of its left and right children have already been 
	// updated

	// Get transition probability matrices for the left and right child nodes of this node
	// These are 3D because there are potentially several 2x2 transition matrices, one
	// for each relative rate
	const double * const * const * leftPMatrices  = leftNode.getConstPMatrices();
	const double * const * const * rightPMatrices = rightNode.getConstPMatrices();

	// cla is the conditional likelihood array we are updating
	double * cla = ndInfo.getCLA(); //PELIGROSO

	// These are the conditional likelihood arrays of the left and right child nodes of this node
	const double * leftCLA  = leftNode.getConstCLA();
	const double * rightCLA = rightNode.getConstCLA();

	// conditional likelihood arrays are laid out as follows for DNA data:
	//
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	// |                            rate 1                             | ...
	// +---+---+---+---+---+---+---+---+---------------+---+---+---+---+
	// |   pattern 1   |   pattern 2   |      ...      |   pattern n   | ...
	// +---+---+---+---+---+---+---+---+---------------+---+---+---+---+
	// | A | C | G | T | A | C | G | T |      ...      | A | C | G | T | ...
	// +---+---+---+---+---+---+---+---+---------------+---+---+---+---+
	//
	for (unsigned r = 0; r < num_rates; ++r)
		{
		const double * const * leftPMatrix  = leftPMatrices[r];
		const double * const * rightPMatrix = rightPMatrices[r];
		for (unsigned pat = 0; pat < num_patterns; ++pat, leftCLA += num_states, rightCLA += num_states)
			{
			for (unsigned i = 0; i < num_states; ++i)
				{
				double left_side  = 0.0;
				double right_side = 0.0;
				const double * leftP_i  = leftPMatrix[i];
				const double * rightP_i = rightPMatrix[i];
				for (unsigned j = 0; j < num_states; ++j)
					{
					left_side  += (leftP_i[j]*leftCLA[j]);
					right_side += (rightP_i[j]*rightCLA[j]);
					}
				*cla++ = (left_side*right_side);
				}
			}
		}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Modifies the conditional likelihood arrays for an internal node having more than two children. This function 
|	handles situations in which an additional child (beyond the first two) is a tip. The conditional likelihood arrays
|	should already have been calculated for the first two children using TreeLikelihood::calcCLATwoTips,
|	TreeLikelihood::calcCLAOneTip or TreeLikelihood::calcCLANoTips.
*/
void TreeLikelihood::conditionOnAdditionalTip(
  InternalData &	ndInfo,
  const TipData &	tipData)
	{
	double * cla = ndInfo.getCLA(); //PELIGROSO
	const double * const * const * tipPMatricesTrans = tipData.getConstTransposedPMatrices();
	const int8_t * tipStateCodes = tipData.getConstStateCodes();
	
	// conditional likelihood arrays are laid out as follows for DNA data:
	//
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	// |   pattern 1   |               |               |               |               |
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	// | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T |
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	//
	for (unsigned r = 0; r < num_rates; ++r)
		{
		const double * const * const tipPMatrixT = tipPMatricesTrans[r];
		for (unsigned pat = 0; pat < num_patterns; ++pat)
			{
			const double * tipPMatT_pat = tipPMatrixT[tipStateCodes[pat]];
			for (unsigned i = 0; i < num_states; ++i)
				*cla++ *= tipPMatT_pat[i];
			}
		}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Modifies the conditional likelihood arrays for an internal node having more than two children. This function 
|	handles situations in which an additional child (beyond the first two) is an internal node. The conditional 
|	likelihood arrays should already have been calculated for the first two children using 
|	TreeLikelihood::calcCLATwoTips, TreeLikelihood::calcCLAOneTip or TreeLikelihood::calcCLANoTips.
*/
void TreeLikelihood::conditionOnAdditionalInternal(
  InternalData &		ndInfo,
  const InternalData &	childInfo)
	{
	// cla is the conditional likelihood array we are updating
	double * cla = ndInfo.getCLA(); //PELIGROSO

	const double * const * const * childPMatrices = childInfo.getConstPMatrices();
	const double * childCLA = childInfo.getConstCLA();

	// conditional likelihood arrays are laid out as follows for DNA data:
	//
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	// |   pattern 1   |               |               |               |               |
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	// | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T |
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	//
	for (unsigned r = 0; r < num_rates; ++r)
		{
		const double * const * childPMatrix = childPMatrices[r];
		for (unsigned pat = 0; pat < num_patterns; ++pat, childCLA += num_states)
			{
			for (unsigned i = 0; i < num_states; ++i)
				{
				double childLike = 0.0;
				const double * childP_i = childPMatrix[i];
				for (unsigned j = 0; j < num_states; ++j)
					childLike += childP_i[j]*childCLA[j];
				*cla++ *= childLike;
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called after conditional likelihood arrays are brought up-to-date, this function calculates the likelihood for each
|	data pattern, as well as the total log-likelihood.
|	
|	Conditional likelihood arrays are laid out as follows for DNA data:
|>
|	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
|	|                            rate 1                             | ...
|	+---+---+---+---+---+---+---+---+---------------+---+---+---+---+
|	|   pattern 1   |   pattern 2   |      ...      |   pattern n   | ...
|	+---+---+---+---+---+---+---+---+---------------+---+---+---+---+
|	| A | C | G | T | A | C | G | T |      ...      | A | C | G | T | ...
|	+---+---+---+---+---+---+---+---+---------------+---+---+---+---+
|>	
*/
double TreeLikelihood::harvestLnL(
  const InternalData & rootCLA)
	{
	// Get pointer to start of array where likelihoods for each site and rate combination will be stored
	//@POL why have likelihood_rate_site? I think it is only used in this function. Try replacing it 
	// everywhere with a simple double value and see if it makes a difference
	unsigned sz = likelihood_rate_site.size();
	if (sz != num_rates*num_patterns)
		{
		std::cerr << "sz           = " << sz << std::endl;
		std::cerr << "num_rates    = " << num_rates << std::endl;
		std::cerr << "num_patterns = " << num_patterns << std::endl;
		}
	assert(likelihood_rate_site.size() == num_rates*num_patterns);
	double * like_rate_site = (double *)(&likelihood_rate_site[0]); //PELIGROSO

	// Get pointer to start of array where site likelihoods will be stored
	assert(site_likelihoods.size() == num_patterns);
	double * siteLike = (double *)(&site_likelihoods[0]); //PELIGROSO

	// Get pointer to start of array holding pattern counts
	assert(pattern_counts.size() == num_patterns);
	PatternCountType * counts = (PatternCountType *)(&pattern_counts[0]); //PELIGROSO

	// Get state frequencies from model and alias rate category probability array for speed
	const double * stateFreq = &model->getStateFreqs()[0]; //PELIGROSO
	const double * rateCatProbArray = &rate_probs[0]; //PELIGROSO
	double rateCatProb = rateCatProbArray[0];

	// Get conditional likelihood arrays for the root node
	assert(rootCLA.getCLASize() == num_rates*num_patterns*num_states);
	const double * cla = rootCLA.getConstCLA();

	double lnLikelihood = 0.0;
	if (num_rates == 1)
		{
		// Compute lnLikelihood assuming rate homogeneity
		// The like_rate_site vector has one element for every pattern in this case
		// For each pattern p, compute current site likelihood, storing this value in
		// both like_rate_site[p] and siteLike[p]. Add this value to lnLikelihood.
		for (unsigned p = 0; p < num_patterns; ++p)
			{
			*like_rate_site = 0.0;
			const double * sf = stateFreq;
			for (unsigned j = 0; j < num_states; ++j)
				{
				double state_freq = *sf++;
				double cond_like = *cla++;
				*like_rate_site += (state_freq*cond_like);
				}
			double count = (double)(*counts++);
			double site_likelihood = *like_rate_site++;
			siteLike[p] = site_likelihood;
			lnLikelihood += count*log(site_likelihood);
			}
		}
	else
		{
		// Compute lnLikelihood assuming rate heterogeneity
		// Compute array of site likelihoods
		memset(siteLike, 0, num_patterns*sizeof(double));
		for (unsigned r = 0; r < num_rates; ++r)
			{
			rateCatProb = rateCatProbArray[r];
		
			for (unsigned p = 0; p < num_patterns; ++p, cla += num_states, ++like_rate_site)
				{
				*like_rate_site = 0.0;
				for (unsigned j = 0; j < num_states; ++j)
					{
					*like_rate_site += stateFreq[j]*cla[j];
					}
				siteLike[p] += rateCatProb*(*like_rate_site);
				}
			}

		// Compute log-likelihood
		//@POL note that an additional loop is needed because rates are not nested within patterns in the cla
		for (unsigned k = 0; k < num_patterns; ++k)
			{
			assert(siteLike[k] > 0.0);
			lnLikelihood += counts[k]*log(siteLike[k]);
			}
		}

	return lnLikelihood;
	}

} //namespace phycas
