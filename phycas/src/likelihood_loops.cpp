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

#define DO_UNDERFLOW_POLICY	
//#if defined(_MSC_VER)
//#	pragma warning(error : 4710) // reports functions not inlined as errors
//#endif

//#include "phycas/force_include.h"
#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/cond_likelihood.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/tip_data.hpp"
#include "phycas/src/internal_data.hpp"
#include "phycas/src/edge_endpoints.hpp"
#include "phycas/src/cipres/util_copy.hpp"
#include <numeric>

//#define LOG_SITELIKES

#include <iostream>
using std::cerr;
using std::endl;

#include <cmath>
using std::log;

using std::accumulate;

template<typename T>
void scaleVector(std::vector<T> & result, const std::vector<T> & orig, T scaler);

template<typename T>
void transpose(T * * mat, unsigned dim);

using std::vector;

template<typename T>
void scaleVector(std::vector<T> & scaled, const vector<T> & orig, T scaler)
	{
	const unsigned origLen = (const unsigned)orig.size();
	scaled.resize(origLen);
	scaled.assign(origLen, scaler);
	cip_imult(orig.begin(), orig.end(), scaled.begin());
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
	PHYCAS_ASSERT(num_rates > 0);
	vector<double> scaledEdges;
	scaleVector(scaledEdges, rate_means, edgeLength);
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
|	
*/
void TreeLikelihood::calcPMat(
  double * * *	p, 
  double		edgeLength)
	{
	calcPMatCommon(p, edgeLength);
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
  double * * * 			transPMats,
  const StateListPos &	stateListPosVec, 
  double				edgeLength)
	{
	calcPMatCommon(transPMats,  edgeLength);

#if 0 // POLPY_NEWWAY    // Rota bug
    if (debugging_now)
        {
        std::cerr << "*** After calcPMatCommon, edgeLength = " << edgeLength << " ***" << std::endl;
        for (unsigned r = 0; r < num_rates; ++r)
            {
            std::cerr << "rate " << r << std::endl;
            std::cerr << str(boost::format("%12.5f %12.5f %12.5f %12.5f") % transPMats[r][0][0] % transPMats[r][0][1] % transPMats[r][0][2] % transPMats[r][0][3]) << std::endl;
            std::cerr << str(boost::format("%12.5f %12.5f %12.5f %12.5f") % transPMats[r][1][0] % transPMats[r][1][1] % transPMats[r][1][2] % transPMats[r][1][3]) << std::endl;
            std::cerr << str(boost::format("%12.5f %12.5f %12.5f %12.5f") % transPMats[r][2][0] % transPMats[r][2][1] % transPMats[r][2][2] % transPMats[r][2][3]) << std::endl;
            std::cerr << str(boost::format("%12.5f %12.5f %12.5f %12.5f") % transPMats[r][3][0] % transPMats[r][3][1] % transPMats[r][3][2] % transPMats[r][3][3]) << std::endl;
            }
        }
#endif

	// For each rate category, transpose the num_states x num_states portion of the matrices
	// and fill in the ambiguity codes by summing columns
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
  CondLikelihood & 		condLike,
  const TipData &		leftTip, 
  const TipData &		rightTip)
	{
	const double * const * const * leftPMatricesTrans = leftTip.getConstTransposedPMatrices();
	const int8_t * leftStateCodes = leftTip.getConstStateCodes();
	const double * const * const * rightPMatricesTrans = rightTip.getConstTransposedPMatrices();
	const int8_t * rightStateCodes = rightTip.getConstStateCodes();
	LikeFltType * cla = condLike.getCLA();
	for (unsigned r = 0; r < num_rates; ++r)
		{
		const double * const * leftPMatT = leftPMatricesTrans[r];
		const double * const * const rightPMatT = rightPMatricesTrans[r];
		for (unsigned p = 0; p < num_patterns; ++p, cla += num_states)
			{
			const double * leftPMatTRow = leftPMatT[leftStateCodes[p]];
			const double * rightPMatTRow = rightPMatT[rightStateCodes[p]];
            for (unsigned s = 0; s < num_states; ++s)
#if 1 // POLPY_NEWWAY    // Rota bug
				cla[s] = leftPMatTRow[s]*rightPMatTRow[s];
#else
                {
                //if (debugging_now && r == 0 && p == 162 && s == 1)
                //    {
                //    std::cerr << "*** rate = " << r << ", pattern = " << p << ", state = " << s << " ***" << std::endl;
                //    }
				cla[s] = leftPMatTRow[s]*rightPMatTRow[s];
                if (cla[s] < 0.0)
                    {
                    std::cerr << "*** Doh! cla negative in calcCLATwoTips for rate = " << r << ", pattern = " << p << ", state = " << s << " ***" << std::endl;
                    std::exit(0);
                    }
                }
#endif
			}
		}
#if defined(DO_UNDERFLOW_POLICY)
	underflow_policy.twoTips(condLike);
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the conditional likelihood arrays at an internal node having a tip node as its left child and an internal
|	node for a right child.
*/
void TreeLikelihood::calcCLAOneTip(
  CondLikelihood &			condLike,
  const TipData &			leftChild,
  ConstPMatrices			rightPMatrices,
  const CondLikelihood &	rightCondLike)
	{
	// Get transition probability matrices for the left and right child nodes of this node
	// These are 3D because there are potentially several 2D transition matrices, one
	// for each relative rate
	ConstPMatrices leftPMatricesTrans = leftChild.getConstTransposedPMatrices();

	// cla is the conditional likelihood array we are updating
	LikeFltType * cla = condLike.getCLA();

	// These are the conditional likelihood arrays of the left and right child nodes of this node
	const int8_t * leftStateCodes = leftChild.getConstStateCodes();
	const double * rightCLA = rightCondLike.getCLA();
	
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
		if (num_states == 4)
			{
			//POL 16-June-2006 Unrolling the nested loops across states here as well as in TreeLikelihood::calcCLANoTips
			// resulted in a 26.6% speedup on Windows using green.nex and SVN version 97
			for (unsigned pat = 0; pat < num_patterns; ++pat)
				{
				const double * leftPMatT_pat = leftPMatrixT[leftStateCodes[pat]];
				double rightCLA0 = *rightCLA++;
				double rightCLA1 = *rightCLA++;
				double rightCLA2 = *rightCLA++;
				double rightCLA3 = *rightCLA++;

				// *** from state 0
				const double * rightP_i = rightPMatrix[0];
				double right0 = (*rightP_i++)*rightCLA0;
				double right1 = (*rightP_i++)*rightCLA1;
				double right2 = (*rightP_i++)*rightCLA2;
				double right3 = (*rightP_i++)*rightCLA3;
				double right_side = right0 + right1 + right2 + right3;
				*cla++ = ((*leftPMatT_pat++)*right_side);

				// *** from state 1
				rightP_i = rightPMatrix[1];
				right0 = (*rightP_i++)*rightCLA0;
				right1 = (*rightP_i++)*rightCLA1;
				right2 = (*rightP_i++)*rightCLA2;
				right3 = (*rightP_i++)*rightCLA3;
				right_side = right0 + right1 + right2 + right3;
				*cla++ = ((*leftPMatT_pat++)*right_side);

				// *** from state 2
				rightP_i = rightPMatrix[2];
				right0 = (*rightP_i++)*rightCLA0;
				right1 = (*rightP_i++)*rightCLA1;
				right2 = (*rightP_i++)*rightCLA2;
				right3 = (*rightP_i++)*rightCLA3;
				right_side = right0 + right1 + right2 + right3;
				*cla++ = ((*leftPMatT_pat++)*right_side);

				// *** from state 3
				rightP_i = rightPMatrix[3];
				right0 = (*rightP_i++)*rightCLA0;
				right1 = (*rightP_i++)*rightCLA1;
				right2 = (*rightP_i++)*rightCLA2;
				right3 = (*rightP_i++)*rightCLA3;
				right_side = right0 + right1 + right2 + right3;
				*cla++ = ((*leftPMatT_pat++)*right_side);
				}
			}
		else	// if (num_states == 4)
			{
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
		} // loop over rates

#if defined(DO_UNDERFLOW_POLICY)
	underflow_policy.oneTip(condLike, rightCondLike, pattern_counts);
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the conditional likelihood arrays at an internal node having two internal nodes for children.
*/
void TreeLikelihood::calcCLANoTips(
  CondLikelihood &			condLike,
  ConstPMatrices 			leftPMatrices,
  const CondLikelihood &	leftCondLike, 
  ConstPMatrices 			rightPMatrices,
  const CondLikelihood &	rightCondLike)
	{
	
	// This function updates the conditional likelihood array of a node assuming that the
	// conditional likelihood arrays of its left and right children have already been 
	// updated

	// cla is the conditional likelihood array we are updating
	LikeFltType * cla = condLike.getCLA();

	// These are the conditional likelihood arrays of the left and right child nodes of this node
	const LikeFltType * leftCLA  = leftCondLike.getCLA();
	const LikeFltType * rightCLA = rightCondLike.getCLA();

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
		double const * const * leftPMatrix  = leftPMatrices[r];
		double const * const * rightPMatrix = rightPMatrices[r];
		if (num_states == 4)
			{
			//POL 16-June-2006 Unrolling the nested loops across states here as well as in TreeLikelihood::calcCLAOneTip
			// resulted in a 26.6% speedup on Windows using green.nex and SVN version 97
			for (unsigned pat = 0; pat < num_patterns; ++pat)
				{
				double leftCLA0 = *leftCLA++;
				double leftCLA1 = *leftCLA++;
				double leftCLA2 = *leftCLA++;
				double leftCLA3 = *leftCLA++;
				double rightCLA0 = *rightCLA++;
				double rightCLA1 = *rightCLA++;
				double rightCLA2 = *rightCLA++;
				double rightCLA3 = *rightCLA++;

				// *** from state 0 ***
				const double * leftP_i  = leftPMatrix[0];
				const double * rightP_i = rightPMatrix[0];

				// to state 0
				double left0 = (*leftP_i++)*leftCLA0;
				double right0 = (*rightP_i++)*rightCLA0;

				// to state 1
				double left1 = (*leftP_i++)*leftCLA1;
				double right1 = (*rightP_i++)*rightCLA1;

				// to state 2
				double left2 = (*leftP_i++)*leftCLA2;
				double right2 = (*rightP_i++)*rightCLA2;

				// to state 3
				double left3 = (*leftP_i++)*leftCLA3;
				double right3 = (*rightP_i++)*rightCLA3;

				double left_side  = left0 + left1 + left2 + left3;
				double right_side = right0 + right1 + right2 + right3;
				*cla++ = (left_side*right_side);

				// *** from state 1 ***
				leftP_i  = leftPMatrix[1];
				rightP_i = rightPMatrix[1];

				// to state 0
				left0 = (*leftP_i++)*leftCLA0;
				right0 = (*rightP_i++)*rightCLA0;

				// to state 1
				left1 = (*leftP_i++)*leftCLA1;
				right1 = (*rightP_i++)*rightCLA1;

				// to state 2
				left2 = (*leftP_i++)*leftCLA2;
				right2 = (*rightP_i++)*rightCLA2;

				// to state 3
				left3 = (*leftP_i++)*leftCLA3;
				right3 = (*rightP_i++)*rightCLA3;

				left_side  = left0 + left1 + left2 + left3;
				right_side = right0 + right1 + right2 + right3;
				*cla++ = (left_side*right_side);

				// *** from state 2 ***
				leftP_i  = leftPMatrix[2];
				rightP_i = rightPMatrix[2];

				// to state 0
				left0 = (*leftP_i++)*leftCLA0;
				right0 = (*rightP_i++)*rightCLA0;

				// to state 1
				left1 = (*leftP_i++)*leftCLA1;
				right1 = (*rightP_i++)*rightCLA1;

				// to state 2
				left2 = (*leftP_i++)*leftCLA2;
				right2 = (*rightP_i++)*rightCLA2;

				// to state 3
				left3 = (*leftP_i++)*leftCLA3;
				right3 = (*rightP_i++)*rightCLA3;

				left_side  = left0 + left1 + left2 + left3;
				right_side = right0 + right1 + right2 + right3;
				*cla++ = (left_side*right_side);

				// *** from state 3 ***
				leftP_i  = leftPMatrix[3];
				rightP_i = rightPMatrix[3];

				// to state 0
				left0 = (*leftP_i++)*leftCLA0;
				right0 = (*rightP_i++)*rightCLA0;

				// to state 1
				left1 = (*leftP_i++)*leftCLA1;
				right1 = (*rightP_i++)*rightCLA1;

				// to state 2
				left2 = (*leftP_i++)*leftCLA2;
				right2 = (*rightP_i++)*rightCLA2;

				// to state 3
				left3 = (*leftP_i++)*leftCLA3;
				right3 = (*rightP_i++)*rightCLA3;

				left_side  = left0 + left1 + left2 + left3;
				right_side = right0 + right1 + right2 + right3;
				*cla++ = (left_side*right_side);
				}
			}
		else	// if (num_states == 4)
			{
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
		}	// loop over rates

#if defined(DO_UNDERFLOW_POLICY)
	underflow_policy.noTips(condLike, leftCondLike, rightCondLike, pattern_counts);
#endif
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Modifies the conditional likelihood arrays for an internal node having more than two children. This function 
|	handles situations in which an additional child (beyond the first two) is a tip. The conditional likelihood arrays
|	should already have been calculated for the first two children using TreeLikelihood::calcCLATwoTips,
|	TreeLikelihood::calcCLAOneTip or TreeLikelihood::calcCLANoTips.
*/
void TreeLikelihood::conditionOnAdditionalTip(
  CondLikelihood &	condLike,
  const TipData &	tipData)
	{
	LikeFltType * cla = condLike.getCLA();
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
	//@POL need to deal with underflow here
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Modifies the conditional likelihood arrays for an internal node having more than two children. This function 
|	handles situations in which an additional child (beyond the first two) is an internal node. The conditional 
|	likelihood arrays should already have been calculated for the first two children using 
|	TreeLikelihood::calcCLATwoTips, TreeLikelihood::calcCLAOneTip or TreeLikelihood::calcCLANoTips.
*/
void TreeLikelihood::conditionOnAdditionalInternal(
  CondLikelihood &			condLike,
  ConstPMatrices			childPMatrices,
  const CondLikelihood &	childCondLike)
	{
	double * cla = condLike.getCLA();
	const double * childCLA = childCondLike.getCLA();

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
	//@POL need to deal with underflow here
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called after neighboring conditional likelihood arrays are brought up-to-date.
|	This function calculates the likelihood for each data pattern, as well as the total log-likelihood.
|	The site-likelihoods are not stored.
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
   EdgeEndpoints & focal_edge)
	{
	// Get the focal node
	TreeNode * focal_node = focal_edge.getFocalNode();
	PHYCAS_ASSERT(focal_node != NULL);

	// If the focal neighbor is NULL, let the parent of the focal node be the focal neighbor
	TreeNode * focal_neighbor = focal_edge.getFocalNeighbor();
	if (focal_neighbor == NULL)
		{
		focal_neighbor = focal_node->GetParent();
		focal_edge.setFocalNeighbor(focal_neighbor);
		}
	PHYCAS_ASSERT(focal_neighbor != NULL);

	// Recompute the conditional likelihood array of the focal node
	refreshCLA(*focal_node, focal_neighbor);
	ConstEdgeEndpoints c(focal_node, focal_neighbor);
	return harvestLnLFromValidEdge(c);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates log-likelihood using the focal node of `focal_edge' as the likelihood root. Assumes that all the 
|	necessary conditional likelihood arrays that are needed are up-to-date.
*/
double TreeLikelihood::harvestLnLFromValidEdge(
   ConstEdgeEndpoints & focal_edge)	/**< is the edge containing the focal node that will serve as the likelihood root */	
	{
	PHYCAS_ASSERT(focal_edge.getFocalNode() != NULL);
	PHYCAS_ASSERT(focal_edge.getFocalNeighbor() != NULL);
	const TreeNode * focalNeighbor = focal_edge.getFocalNeighbor();
	const TreeNode * focalNode = focal_edge.getFocalNode();
	PHYCAS_ASSERT(focalNode->IsInternal());
	const TreeNode * actualChild = focal_edge.getActualChild();
	const double focalEdgeLen = actualChild->GetEdgeLen();
	
	ConstCondLikelihoodShPtr focalCondLike = getValidCondLikePtr(focal_edge);
	PHYCAS_ASSERT(focalCondLike);
	const LikeFltType * focalNodeCLA = focalCondLike->getCLA(); //PELIGROSO
	PHYCAS_ASSERT(focalNodeCLA != NULL);
	const unsigned singleRateCLALength = num_patterns*num_states;

	// Get pointer to start of array holding pattern counts
	PHYCAS_ASSERT(pattern_counts.size() == num_patterns);
	PatternCountType * counts = (PatternCountType *)(&pattern_counts[0]); //PELIGROSO

	// Get state frequencies from model and alias rate category probability array for speed
	const double * stateFreq = &model->getStateFreqs()[0]; //PELIGROSO
	const double * rateCatProbArray = &rate_probs[0]; //PELIGROSO

	if (store_site_likes)
		site_likelihood.clear();

#if defined(LOG_SITELIKES)
	std::ofstream tmpf("doof.txt", std::ios::out | std::ios::app);
	tmpf << "pat\tfreqA\tCLAA\tPrA\tfreqC\tCLAC\tPrC\tfreqG\tCLAG\tPrG\tfreqT\tCLAT\tPrT\tuncorr\tcorr\tcount\tlnL" << std::endl;
#endif

	double lnLikelihood = 0.0;
	if (focalNeighbor->IsTip())
		{
		const TipData & tipData = *focalNeighbor->GetTipData();
		double * * * p = tipData.getMutableTransposedPMatrices();
		calcPMatTranspose(p, tipData.getConstStateListPos(),  focalEdgeLen);
		const double * const * const * tipPMatricesTrans = tipData.getConstTransposedPMatrices();
		const int8_t * tipStateCodes = tipData.getConstStateCodes();
		std::vector<const double *> focalNdCLAPtr(num_rates);
		for (unsigned i = 0; i < num_rates; ++i)
			focalNdCLAPtr[i] = focalNodeCLA + singleRateCLALength*i;
			
		for (unsigned pat = 0; pat < num_patterns; ++pat)
			{
#if defined(LOG_SITELIKES)
			tmpf << pat << '\t';
#endif
			double siteLike = 0.0;
			for (unsigned r = 0; r < num_rates; ++r)
				{
				const double * const * const tipPMatrixT = tipPMatricesTrans[r];
				const double * tipPMatT_pat = tipPMatrixT[tipStateCodes[pat]];
				double siteRateLike = 0.0;
                for (unsigned i = 0; i < num_states; ++i)
                    {
#if defined(LOG_SITELIKES)
        			tmpf << stateFreq[i] << '\t';
		        	tmpf << focalNdCLAPtr[r][i] << '\t';
			        tmpf << tipPMatT_pat[i] << '\t';
			        tmpf << (stateFreq[i]*focalNdCLAPtr[r][i]*tipPMatT_pat[i]) << '\t';
#endif
					siteRateLike += stateFreq[i] * focalNdCLAPtr[r][i] * tipPMatT_pat[i];
                    }
				siteLike += rateCatProbArray[r] * siteRateLike;
				focalNdCLAPtr[r] += num_states;
				}
			double site_lnL = std::log(siteLike);

#if defined(LOG_SITELIKES)
			tmpf << site_lnL << '\t';
#endif

#if defined(DO_UNDERFLOW_POLICY)
			underflow_policy.correctSiteLike(site_lnL, pat, focalCondLike);
#endif

			if (store_site_likes)
				site_likelihood.push_back(site_lnL);
			lnLikelihood += counts[pat]*site_lnL;

#if defined(LOG_SITELIKES)
			tmpf << site_lnL << '\t';
			tmpf << counts[pat] << '\t';
			tmpf << lnLikelihood << std::endl;
#endif
			}		
#if defined(DO_UNDERFLOW_POLICY)
		underflow_policy.correctLnLike(lnLikelihood, focalCondLike);
#endif

#if defined(LOG_SITELIKES)
		tmpf << "1\t\t\t\t\t" << lnLikelihood << std::endl;
#endif
		}
	else
		{
		const InternalData * neighborID = focalNeighbor->GetInternalData();
		calcPMat(neighborID->getMutablePMatrices(), focalEdgeLen);
		const double * const * const * childPMatrices = neighborID->getConstPMatrices();
		
		ConstCondLikelihoodShPtr neighborCondLike = getValidCondLikePtr(focalNeighbor, focalNode); // 

		const double * focalNeighborCLA = neighborCondLike->getCLA(); //PELIGROSO
		std::vector<const double *> focalNdCLAPtr(num_rates);
		std::vector<const double *> focalNeighborCLAPtr(num_rates);
		
		for (unsigned i = 0; i < num_rates; ++i)
			{
			focalNdCLAPtr[i] = focalNodeCLA + singleRateCLALength*i;
			focalNeighborCLAPtr[i] = focalNeighborCLA + singleRateCLALength*i;
			}

		for (unsigned pat = 0; pat < num_patterns; ++pat)
			{
			double siteLike = 0.0;
			for (unsigned r = 0; r < num_rates; ++r)
				{
				const double * const * childPMatrix = childPMatrices[r];
				const double * neigborCLAForRate = focalNeighborCLAPtr[r];
				double siteRateLike = 0.0;
				for (unsigned i = 0; i < num_states; ++i)
					{
					double neigborLike = 0.0;
					const double * childP_i = childPMatrix[i];
					for (unsigned j = 0; j < num_states; ++j)
						neigborLike += childP_i[j]*neigborCLAForRate[j];
					siteRateLike += stateFreq[i]*focalNdCLAPtr[r][i]*neigborLike;
					}
				siteLike += rateCatProbArray[r]*siteRateLike;
				focalNdCLAPtr[r] += num_states;
				focalNeighborCLAPtr[r] += num_states;
				}
			double site_lnL = std::log(siteLike);

#if defined(LOG_SITELIKES)
			tmpf << '\t';
			tmpf << pat << '\t';
			tmpf << site_lnL << '\t';
#endif

#if defined(DO_UNDERFLOW_POLICY)
			underflow_policy.correctSiteLike(site_lnL, pat, focalCondLike);
			underflow_policy.correctSiteLike(site_lnL, pat, neighborCondLike);
#endif
			if (store_site_likes)
				site_likelihood.push_back(site_lnL);
			lnLikelihood += counts[pat]*site_lnL;

#if defined(LOG_SITELIKES)
			tmpf << site_lnL << '\t';
			tmpf << counts[pat] << '\t';
			tmpf << lnLikelihood << std::endl;
#endif
			}
#if defined(DO_UNDERFLOW_POLICY)
		underflow_policy.correctLnLike(lnLikelihood, focalCondLike);
		underflow_policy.correctLnLike(lnLikelihood, neighborCondLike);
#endif

#if defined(LOG_SITELIKES)
		tmpf << "0\t\t\t\t\t" << lnLikelihood << std::endl;
#endif
		}

#if defined(LOG_SITELIKES)
	tmpf.close();
#endif

	//std::ofstream doof("counts_patterns.txt");
	//std::vector<double>::iterator sit = site_likelihood.begin();
	//for (PatternMapType::iterator pit = pattern_map.begin(); pit != pattern_map.end(); ++pit, ++sit)
	//	{
	//	PHYCAS_ASSERT(sit != site_likelihood.end());
	//	doof << str(boost::format("%12.5f\t%.1f\t|\t") % (*sit) % pit->second);
	//	std::copy(pit->first.begin(), pit->first.end(), std::ostream_iterator<int>(doof, "\t")); 
	//	doof << std::endl;
	//	}
	//doof.close();

	return lnLikelihood;
	}

} //namespace phycas
