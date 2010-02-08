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

//#include "phycas/force_include.h"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/tip_data.hpp"
#include "phycas/src/cond_likelihood.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor for TipData objects that will only be used for simulating data and thus needs only to allocate the 
|	`state' and `pMatrixTranspose' data members.
*/
TipData::TipData(
  unsigned nRates,		/**< is the number of relative rate categories */
  unsigned nStates,		/**< is the number of states in the model */
  CondLikelihoodStorageShPtr cla_storage)
	:
	state(-1), 
	pMatrixTranspose(NULL),
	cla_pool(cla_storage),
	sMat(0L)
	{
#if DISABLED_UNTIL_SIMULATION_WORKING_WITH_PARTITIONING
	const unsigned nToStates = nStates;	// simulated data never has ambiguities, so num. rows in T matrix is just nStates
	ownedPMatrices.Initialize(nRates, nToStates, nStates);
	pMatrixTranspose = ownedPMatrices.ptr;
	sMat =  NewTwoDArray<unsigned>(nStates, nStates); //@ we should make this only true in unimap mode!!!
	for (unsigned i = 0; i < nStates*nStates ; ++i)
		sMat[0][i] = 0;
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor for TipData objects that will be used in likelihood calculations and thus needs to store the observed
|	data for the tip as well as information about local state codes.
*/
TipData::TipData(
  bool                              			using_unimap,       /**< is true if tips are to be prepared for uniformized mapping likelihood; it is false if tips are to be prepared for Felsenstein-style integrated likelihoods */
  unsigned				            			nPatterns,			/**< is the number of site patterns */
  PartitionModelShPtr							partition,			/**< is the PartitionModel object containing information about the number of states and rates for each subset */
  const state_list_pos_vect_t &					positions,			/**< is the vector of vectors of positions of each state into the `state_codes' vector */
  const state_list_vect_t &						states,				/**< is the vector of vectors of `state_codes' */ 
  CondLikelihoodStorageShPtr 					cla_storage)		/**< is the pool of available conditional likelihood arrays */
	:
    unimap(using_unimap),
	state(-1), 
	cla_pool(cla_storage)
	{
	const unsigned num_subsets = partition->getNumSubsets();
	state_list_pos	= state_list_pos_vect_t(num_subsets);
	state_codes 	= state_list_vect_t(num_subsets);
	
	pMatrixTranspose.resize(num_subsets);	
	for (unsigned i = 0; i < num_subsets; ++i)
		{
		// copy state list positions for subset i
		state_list_pos[i].resize(positions[i].size());
		std::copy(positions[i].begin(), positions[i].end(), state_list_pos[i].begin());
		
		// copy state codes for subset i
		state_codes[i].resize(states[i].size());
		std::copy(states[i].begin(), states[i].end(), state_codes[i].begin());
		
		// allocate memory for the ScopedThreeDMatrix of transition probabilities for subset i
		const unsigned num_obs_states	= partition->subset_num_states[i] + 1 + positions[i].size();
		const unsigned num_rates		= partition->subset_num_rates[i];
		const unsigned num_states		= partition->subset_num_states[i];
		pMatrixTranspose[i].Initialize(num_rates, num_obs_states, num_states);
		}

#if DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
    if (using_unimap)
        {
        univents.resize(nPatterns);
        univents.setEndStates(state_codes.get());
        std::vector<state_code_t> & scVec = univents.getEndStatesVecRef();
        for (unsigned i = 0; i < nPatterns; ++i)
        	{
        	const state_code_t sc = state_codes[i];
        	////////////////////////////////////////////////////////////////////
        	// @@@@
        	// This is really dangerous. 
        	// In this loop, we alter the ambiguous data to have state 0.
        	// This should work because we are calling sampleTipsAsDisconnected before MCMC (in MCMCManager.py)
        	// It is hard to guarantee this (we could add a boolean flag to the univents class
        	//		so that we could at least flag this univent as temporarily bogus).
        	// It would be hard to do the correct sampling here because we would need the 
        	//		neighboring nodes and a Lot instance.
        	////////////////////////////////////////////////////////////////////
        	if (sc < 0 || (int)sc >= (int)nStates)
        		scVec[i] = 0; 
        	// throw XLikelihood("Sorry, we currently do not support data sets with ambiguity or gaps when you are using uniformization-based methods");
        	}
		sMat =  NewTwoDArray<unsigned>(nStates, nStates); //@ we should make this only true in unimap mode!!!
		for (unsigned i = 0; i < nStates*nStates ; ++i)
			sMat[0][i] = 0;
        }
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor ensures that all CLA structures are returned to `cla_pool'. The `ownedPMatrices' and `univents' data 
|   members takes care of deleting themselves when they go out of scope.
*/
TipData::~TipData()
	{
#if DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
	DeleteTwoDArray<unsigned>(sMat);
#endif
	// Invalidate the parental CLAs if they exist
	if (parWorkingCLA)
		{
		cla_pool->putCondLikelihood(parWorkingCLA);
		parWorkingCLA.reset();
		}

	// Remove cached parental CLAs if they exist
	if (parCachedCLA)
		{
		cla_pool->putCondLikelihood(parCachedCLA);
		parCachedCLA.reset();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `parWorkingCLA' data member. If `parWorkingCLA' does not currently point to anything, a CondLikelihood 
|	object is first retrieved from `cla_pool', so this function always returns a shared pointer that actually points to
|	something.
*/
CondLikelihoodShPtr TipData::getParentalCondLikePtr()
	{
	if (!parWorkingCLA)
		parWorkingCLA = cla_pool->getCondLikelihood();
	return parWorkingCLA;
	}

}	// namespace phycas
