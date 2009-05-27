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

#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/basic_tree_node.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/unimap_sample_ambig_move.hpp"
#include "phycas/src/topo_prior_calculator.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/tree_manip.hpp"
#include "phycas/src/tip_data.hpp"

//#define DEBUG_LOG
#if defined(DEBUG_LOG)
#   include <fstream>
#endif

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor sets `origEdgelen' to 0.0, `origNode' to NULL, `likeRoot' to NULL, and `lambda' to 1.0.
*/
UnimapSampleAmbigMove::UnimapSampleAmbigMove(
  TreeLikeShPtr treeLike,
  TreeShPtr t,
  ModelShPtr model,
  unsigned weight)
  	{
  	this->setTreeLikelihood(treeLike);
  	this->setTree(t);
  	this->setModel(model);
  	this->setWeight(weight);
  	const unsigned numPatterns = treeLike->getNPatterns();
  	const unsigned numStates = treeLike->getNStates();
	for (TreeNode * nd = tree->GetFirstPreorder(); nd ; nd = nd->GetNextPreorder())
		{
		if (nd->IsTip())
			{
			std::list<unsigned> ambigIndList;
			TipData * td = nd->GetTipData();
			PHYCAS_ASSERT(td);
			const int8_t * scArray = td->getConstStateCodes();
	        for (unsigned i = 0; i < numPatterns; ++i)
	        	{
	        	const int_state_code_t sc = scArray[i];
				if (sc < 0 || sc >= (int_state_code_t)numStates)
        			ambigIndList.push_back(i);
	        	}
	        if (!ambigIndList.empty())
	        	{
	        	std::vector<unsigned> ambigIndVec(ambigIndList.begin(), ambigIndList.end());
	        	ambigTipToAmbigCol[nd] = ambigIndVec;
				}
			}
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets `origNode' to NULL and `origEdgelen' to 0.0. These variables are used to save quantities required for reverting
|	a proposed move that was not accepted, so this function is called at the end of both accept() and revert() to reset the
|	object for the next proposal. It is also called by the constructor to initialize these variables.
*/
void UnimapSampleAmbigMove::reset()
	{
	}


/*--------------------------------------------------------------------------------------------------------------------------
*/
double UnimapSampleAmbigMove::getLnHastingsRatio() const
	{
	PHYCAS_ASSERT(false);
	return 0.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	This move does not change the model dimension, so the Jacobian is irrelevant.
*/
double UnimapSampleAmbigMove::getLnJacobian() const
	{
	PHYCAS_ASSERT(false);
	return 0.0;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is accepted.
*/
void UnimapSampleAmbigMove::accept()
	{
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	not called because this move draws states from the posterior
*/
void UnimapSampleAmbigMove::revert()
	{
	PHYCAS_ASSERT(false);
	}

void fillEndStateProb(std::vector<double> & esProb, const double * const * pMat, const unsigned numStates, const int8_t * const stateListArr, const unsigned neighborSC, unsigned indexIntoStateList, const bool samplingTipRoot);

inline void fillEndStateProb(std::vector<double> & esProb, const double * const * pMat, const unsigned numStates, const int8_t * const stateListArr, const unsigned neighborSC, unsigned indexIntoStateList, const bool samplingTipRoot)
	{
	PHYCAS_ASSERT(neighborSC >= 0 && neighborSC < numStates);
	esProb.assign(numStates, 0.0);
	//std::cerr << "esProb.size() == " << esProb.size() << '\n';
	double totalProb = 0.0;
	if (samplingTipRoot)
		{
		const unsigned nObservedStates = stateListArr[indexIntoStateList++];
		for (unsigned m = 0; m < nObservedStates; ++m)
			{
			unsigned currObservedState = stateListArr[indexIntoStateList++];
			const double * pMatRow = pMat[currObservedState];
			const double el = pMatRow[neighborSC];
			esProb[currObservedState] = el;
			totalProb += el;
			}
		}
	else
		{
		const unsigned nObservedStates = stateListArr[indexIntoStateList++];
		const double * pMatRow = pMat[neighborSC];
		for (unsigned m = 0; m < nObservedStates; ++m)
			{
			unsigned currObservedState = stateListArr[indexIntoStateList++];
			const double el = pMatRow[currObservedState];
			esProb[currObservedState] = el;
			totalProb += el;
			}
		}	
	// normalize (in the case of non-? states, totalProb will be a sub probability.
	for (unsigned i = 0; i < numStates; ++i)
		esProb[i] /= totalProb;
	}
	
	
void UnimapSampleAmbigMove::sampleTipsAsDisconnected()
	{
	AmbigTipMap::iterator ndIt = ambigTipToAmbigCol.begin();
	for (; ndIt != ambigTipToAmbigCol.end(); ++ndIt)
		{
		TreeNode * nd = ndIt->first;
		const std::vector<unsigned> & ambigInds = ndIt->second;
		sampleNewStateArrayForNodeAsDisconnected(nd, ambigInds);
		}
	}


	
/*--------------------------------------------------------------------------------------------------------------------------
| Updates a single node.
*/
void UnimapSampleAmbigMove::sampleNewStateArrayForNodeAsDisconnected(TreeNode * nd, const std::vector<unsigned> & ambigInds)
	{
	PHYCAS_ASSERT(nd->IsTip());
	TipData * td = nd->GetTipData();
	PHYCAS_ASSERT(td);	
	Univents & univents = td->getUniventsRef();
	std::vector<int_state_code_t> & scVec = univents.getEndStatesVecRef();
	const int_state_code_t * scArray = td->getConstStateCodes();

	const std::vector<double> & sf = model->getStateFreqs();
	const double * iniSF = & sf[0];
  	const unsigned numStates = likelihood->getNStates();
	std::vector<const double *> fakePMat(numStates, iniSF);
	const double * const * pMat = &fakePMat[0];

	const StateListPos & stateListPosVec = td->getConstStateListPos();
	const unsigned nPartialAmbigs = (unsigned)stateListPosVec.size();
	const VecStateList & stateListVec = likelihood->getStateList();
	const int8_t * const stateListArr = &stateListVec[0];
	const unsigned int * const stateListPosArr = &stateListPosVec[0];
	const unsigned numStateCodes = numStates + 1 + nPartialAmbigs;
	std::vector<double> emptyCell;
	std::vector< std::vector<double> > endStateProb(numStateCodes, emptyCell);
	//std::cerr << "nPartialAmbigs = " << nPartialAmbigs << '\n';
	
	for (std::vector<unsigned>::const_iterator iIt = ambigInds.begin(); iIt != ambigInds.end(); ++iIt)
		{
		unsigned i = *iIt;
		const int_state_code_t ambigSC = scArray[i];
		std::vector<double> & esProb = endStateProb[ambigSC];
		//std::cerr << "ambigSC = " << (int)ambigSC << " endStateProb.size() = " << endStateProb.size() <<'\n';

		if (esProb.empty())
			{
			if (ambigSC > (int_state_code_t) numStates)
				{
				unsigned indexIntoStateList = stateListPosArr[ambigSC];
				fillEndStateProb(esProb, pMat, numStates, stateListArr, 0, indexIntoStateList, false);
				}
 			else
				{
				PHYCAS_ASSERT(ambigSC < 0 || ambigSC == (int) numStates);
				esProb.resize(numStates);
				const double * pRow = *pMat;
				for (unsigned j = 0 ; j < numStates ; ++j)
					esProb[j] = pRow[j];
				}
	 	 	}
//		std::cerr << "in sampleNewStateArrayForNodeAsDisconnected esProb.size() == " << esProb.size() << " numStates = " << numStates << '\n';

		const int_state_code_t chosenStateCode = (int_state_code_t) rng->MultinomialDraw(&esProb[0], numStates);
		PHYCAS_ASSERT(chosenStateCode >= 0 && chosenStateCode < (int_state_code_t)numStates);
		scVec[i] = chosenStateCode;
		}


	likelihood->flagNodeWithInvalidUnivents(nd);
	univents.setValid(false);
	}


	
/*--------------------------------------------------------------------------------------------------------------------------
| Updates a single node.
*/
void UnimapSampleAmbigMove::proposeNewStateArrayForNode(TreeNode * nd, const std::vector<unsigned> & ambigInds)
	{
	PHYCAS_ASSERT(nd->IsTip());
	TipData * td = nd->GetTipData();
	PHYCAS_ASSERT(td);
	const bool samplingTipRoot = nd->IsTipRoot();
	const TreeNode * neighbor = (samplingTipRoot ? nd->GetLeftChildConst() : nd->GetParentConst());
	PHYCAS_ASSERT(neighbor);
	const Univents & neighborUnivents = getUniventsConstRef(*neighbor);
	const std::vector<int_state_code_t> & neighborStateCodes = neighborUnivents.getEndStatesVecConstRef();
	
	Univents & univents = td->getUniventsRef();
	std::vector<int_state_code_t> & scVec = univents.getEndStatesVecRef();
	const int_state_code_t * scArray = td->getConstStateCodes();
	double * * * pMats = td->getTransposedPMatrices();

	const double edgeLen = (samplingTipRoot ? neighbor->GetEdgeLen() : nd->GetEdgeLen());
	likelihood->calcPMat(pMats, edgeLen);


	const StateListPos & stateListPosVec = td->getConstStateListPos();
	const unsigned nPartialAmbigs = (unsigned)stateListPosVec.size();
	const VecStateList & stateListVec = likelihood->getStateList();
	const int8_t * const stateListArr = &stateListVec[0];
	PHYCAS_ASSERT(likelihood->getNRatesTotal() == 1);
	const double * const * pMat = const_cast<const double * const *>(pMats[0]); //@ assumes no rate het
	const unsigned int * const stateListPosArr = &stateListPosVec[0];
  	const unsigned numStates = likelihood->getNStates();
	const unsigned numStateCodes = numStates + 1 + nPartialAmbigs;

	std::vector<double> emptyCell;
	std::vector< std::vector<double> > emptyRow(numStateCodes, emptyCell);
	std::vector<std::vector< std::vector<double> > > endStateProb(numStates, emptyRow);
	
	
	for (std::vector<unsigned>::const_iterator iIt = ambigInds.begin(); iIt != ambigInds.end(); ++iIt)
		{
		unsigned i = *iIt;
		const int_state_code_t neighborSC = neighborStateCodes[i];
		const int_state_code_t ambigSC = scArray[i];
		std::vector<double> & esProb = endStateProb[neighborSC][ambigSC];
		if (esProb.empty())
			{
			if (ambigSC > (int)numStates)
				{
				unsigned indexIntoStateList = stateListPosArr[ambigSC];
				fillEndStateProb(esProb, pMat, numStates, stateListArr, neighborSC, indexIntoStateList, samplingTipRoot);
				}
			else
				{
				PHYCAS_ASSERT(ambigSC < 0 || ambigSC == (int) numStates);
				esProb.resize(numStates);
				const double * pRow = pMat[neighborSC];
				for (unsigned j = 0 ; j < numStates ; ++j)
					esProb[j] = pRow[j];
				}
 	 	 	}
		const int_state_code_t chosenStateCode = (int_state_code_t) rng->MultinomialDraw(&esProb[0], numStates);
		PHYCAS_ASSERT(chosenStateCode >= 0 && chosenStateCode < (int_state_code_t)numStates);
		scVec[i] = chosenStateCode;
		}

	
	// this next calcPMatTranspose is overkill, but I don't want to leave the node's data
	//	in an inconsistent state. Above we call calcPMatCommon which may flag the pmat as "clean"
	//	but it may not have all of the partial ambiguity columns, and it will not be 
	//	the transpose of the pMat (which is what it usually stored at a tip).
	likelihood->calcPMatTranspose(pMats, stateListPosVec, nd->GetEdgeLen());

	likelihood->flagNodeWithInvalidUnivents(nd);
	univents.setValid(false);
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Updates all nodes in the tree.
*/
void UnimapSampleAmbigMove::proposeNewState()
	{
	AmbigTipMap::iterator ndIt = ambigTipToAmbigCol.begin();
	for (; ndIt != ambigTipToAmbigCol.end(); ++ndIt)
		{
		TreeNode * nd = ndIt->first;
		const std::vector<unsigned> & ambigInds = ndIt->second;
		proposeNewStateArrayForNode(nd, ambigInds);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool UnimapSampleAmbigMove::update()
	{
	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
	proposeNewState();
	accept();
	reset();
	return true;
	}

