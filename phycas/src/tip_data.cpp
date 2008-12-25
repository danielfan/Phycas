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
#if POLPY_NEWWAY	//CLAShPtr
  CondLikelihoodStorageShPtr cla_storage)
#else
  CondLikelihoodStorage & cla_storage)
#endif
	:
	//parCLAValid(false)
	//parWorkingCLA(NULL),
	//parCachedCLA(NULL),
	state(-1), 
	pMatrixTranspose(NULL),
	cla_pool(cla_storage)
	{
	const unsigned nToStates = nStates;	// simulated data never has ambiguities, so num. rows in T matrix is just nStates
	ownedPMatrices.Initialize(nRates, nToStates, nStates);
	pMatrixTranspose = ownedPMatrices.ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor for TipData objects that will be used in likelihood calcuations and thus need to store the observed
|	data for the tip as well as information about local state codes.
*/
TipData::TipData(
  bool                              using_unimap,       /**< is true if tips are to be prepared for uniformized mapping likelihood; it is false if tips are to be prepared for Felsenstein-style integrated likelihoods */
  unsigned				            nPatterns,			/**< is the number of site patterns */
  const std::vector<unsigned int> &	stateListPosVec,	/**< is the vector of positions of each state into the `state_codes' vector */
  boost::shared_array<const int8_t>	stateCodesShPtr,	/**< is the `state_codes' vector */ 
  unsigned							nRates,				/**< is the number of relative rate categories */
  unsigned							nStates,			/**< is the number of states in the model */
  double * * *						pMatTranspose,		/**< is an alias to the rates by states by states pMatrix array, may be NULL */
  bool								managePMatrices, 	/**< if true, a 3D matrix will be allocated (if pMat is also NULL, the pMatrices will alias ownedPMatrices.ptr) */ 
#if POLPY_NEWWAY	//CLAShPtr
  CondLikelihoodStorageShPtr 		cla_storage)
#else
  CondLikelihoodStorage & 			cla_storage)
#endif
	:
	//parCLAValid(false),
	//parWorkingCLA(NULL),
	//parCachedCLA(NULL),
    unimap(using_unimap),
	state(-1), 
	state_list_pos(stateListPosVec), 
	state_codes(stateCodesShPtr), 
	pMatrixTranspose(pMatTranspose),
	cla_pool(cla_storage)
	{
	const unsigned nObservedStates = nStates + 1 + state_list_pos.size();
    if (using_unimap)
        {
        univents.resize(nPatterns);
        univents.setEndStates(state_codes.get());
        }
	if (managePMatrices)
		{
		ownedPMatrices.Initialize(nRates, nObservedStates, nStates);
		//@POL 31-Oct-2005 what if managePMatrices is true, but pMatrixTranspose is not NULL?
		if (pMatrixTranspose == NULL)
			pMatrixTranspose = ownedPMatrices.ptr;
		}
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor ensures that all CLA structures are returned to `cla_pool'. The `ownedPMatrices' and `univents' data 
|   members takes care of deleting themselves when they go out of scope.
*/
TipData::~TipData()
	{
	// Invalidate the parental CLAs if they exist
	if (parWorkingCLA)
		{
#if POLPY_NEWWAY	//CLAShPtr
		cla_pool->putCondLikelihood(parWorkingCLA);
#else
		cla_pool.putCondLikelihood(parWorkingCLA);
#endif
		parWorkingCLA.reset();
		}

	// Remove cached parental CLAs if they exist
	if (parCachedCLA)
		{
#if POLPY_NEWWAY	//CLAShPtr
		cla_pool->putCondLikelihood(parCachedCLA);
#else
		cla_pool.putCondLikelihood(parCachedCLA);
#endif
		parCachedCLA.reset();
		}
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `parWorkingCLA' data member. If `parWorkingCLA' does not currently point to anything, a CondLikelihood 
|	object is first retrieved from `cla_pool', so this function always returns a shared pointer that actually points to
|	something.
*/
CondLikelihoodShPtr TipData::getParentalCondLikePtr()
	{
#if POLPY_NEWWAY	//CLAShPtr
	if (!parWorkingCLA)
		parWorkingCLA = cla_pool->getCondLikelihood();
#else
	if (!parWorkingCLA)
		parWorkingCLA = cla_pool.getCondLikelihood();
#endif
	return parWorkingCLA;
	}

}	// namespace phycas
