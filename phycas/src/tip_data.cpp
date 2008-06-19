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
#include "phycas/src/cipres/CipresDataMatrixHelper.h"
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
  CondLikelihoodStorage & cla_storage)
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

#if POLPY_NEWWAY
StateTimeListVect * TipData::getStateTimeListVect()
	{
	return &state_time;
	}
#endif
/*----------------------------------------------------------------------------------------------------------------------
|	Constructor for TipData objects that will be used in likelihood calcuations and thus need to store the observed
|	data for the tip as well as information about local state codes.
*/
TipData::TipData(
#if POLPY_NEWWAY
  bool                              using_unimap,       /**< is true if tips are to be prepared for uniformized mapping likelihood; it is false if tips are to be prepared for Felsenstein-style integrated likelihoods */
  unsigned				            nPatterns,			/**< is the number of site patterns */
#endif
  const std::vector<unsigned int> &	stateListPosVec,	/**< is the vector of positions of each state into the `state_codes' vector */
  boost::shared_array<const int8_t>	stateCodesShPtr,	/**< is the `state_codes' vector */ 
  unsigned							nRates,				/**< is the number of relative rate categories */
  unsigned							nStates,			/**< is the number of states in the model */
  double * * *						pMatTranspose,		/**< is an alias to the rates by states by states pMatrix array, may be NULL */
  bool								managePMatrices, 	/**< if true, a 3D matrix will be allocated (if pMat is also NULL, the pMatrices will alias ownedPMatrices.ptr) */ 
  CondLikelihoodStorage & 			cla_storage)
	:
	//parCLAValid(false),
	//parWorkingCLA(NULL),
	//parCachedCLA(NULL),
	state(-1), 
	state_list_pos(stateListPosVec), 
	state_codes(stateCodesShPtr), 
	pMatrixTranspose(pMatTranspose),
#if POLPY_NEWWAY
    unimap(using_unimap),
    mdot(0),
#endif
	cla_pool(cla_storage)
	{
	const unsigned nObservedStates = nStates + 1 + state_list_pos.size();
#if POLPY_NEWWAY
    if (using_unimap)
        {
        state_time.resize(nPatterns);
        }
#endif
	if (managePMatrices)
		{
		ownedPMatrices.Initialize(nRates, nObservedStates, nStates);
		//@POL 31-Oct-2005 what if managePMatrices is true, but pMatrixTranspose is not NULL?
		if (pMatrixTranspose == NULL)
			pMatrixTranspose = ownedPMatrices.ptr;
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
		parWorkingCLA = cla_pool.getCondLikelihood();
	return parWorkingCLA;
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of univents for site `site'. Assumes `unimap' is true and `size' is less than the length of the
|   `state_time' vector.
*/
unsigned TipData::getNumUnivents(
  unsigned site) const      /**< is the site of interest */
    {
    if (!unimap || site >= (unsigned)state_time.size())
        return 0;
    return (unsigned)state_time[site].size();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of univent states for site `site'. Assumes `unimap' is true and `size' is less than the length of 
|   the `state_time' vector. This function is not particularly efficient, and it intended primarily for transferring
|   univent states to Python code for debugging purposes.
*/
std::vector<unsigned> TipData::getUniventStates(
  unsigned site) const      /**< is the site of interest */
    {
    std::vector<unsigned> v;
    if (!unimap || site >= (unsigned)state_time.size())
        return v;
    v.resize(state_time[site].size());
    unsigned i = 0;
    for (StateTimeList::const_iterator it = state_time[site].begin(); it != state_time[site].end(); ++it, ++i)
        {
        v[i] = (unsigned)it->first;
        }
    return v;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of univent times for site `site'. Assumes `unimap' is true and `size' is less than the length of 
|   the `state_time' vector. This function is not particularly efficient, and it intended primarily for transferring
|   univent times to Python code for debugging purposes.
*/
std::vector<double> TipData::getUniventTimes(
  unsigned site) const      /**< is the site of interest */
    {
    std::vector<double> v;
    if (!unimap || site >= (unsigned)state_time.size())
        return v;
    v.resize(state_time[site].size());
    unsigned i = 0;
    for (StateTimeList::const_iterator it = state_time[site].begin(); it != state_time[site].end(); ++it, ++i)
        {
        v[i] = (double)it->second;
        }
    return v;
    }

#endif

}	// namespace phycas
