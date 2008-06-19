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
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/internal_data.hpp"

namespace phycas
{

#if POLPY_NEWWAY
StateTimeListVect * InternalData::getStateTimeListVect()
	{
	return &state_time;
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor initializes the conditional likelihood vectors using `nRates' and `nStates', and initializes the 
|	`pMatrices' data member using `pMat'. The parameter `managePMatrices' determines whether to allocate space for the 
|	transition probability matrices.
*/
InternalData::InternalData(
#if POLPY_NEWWAY
  bool                  using_unimap,       /**< is true if internal nodes are to be prepared for uniformized mapping likelihood; it is false if internal nodes are to be prepared for Felsenstein-style integrated likelihoods */
#endif
  unsigned				nPatterns,			/**< is the number of site patterns */
  unsigned				nRates,				/**< is the number of relative rate categories */
  unsigned				nStates,			/**< is the number of states in the model */
  double * * *			pMat,				/**< is an alias to the rates by states by states `pMatrix' array, and may be NULL */
  bool 					managePMatrices, 	/**< if true, a 3D matrix will be allocated (if `pMat' is also NULL, `pMatrices' will alias `ownedPMatrices.ptr') */ 
  CondLikelihoodStorage & cla_storage)
	:
	//parCLAValid(false),
	//parWorkingCLA(NULL),	
	//parCachedCLA(NULL),
	//childCLAValid(false),
	//childWorkingCLA(NULL),
	//childCachedCLA(NULL),
	state(-1), 
	pMatrices(pMat),
#if POLPY_NEWWAY
    unimap(using_unimap),
    mdot(0),
#endif
	cla_pool(cla_storage)
	{
#if POLPY_NEWWAY
    if (using_unimap)
        {
        state_time.resize(nPatterns);
        }
#endif
	if (managePMatrices)
		{
		ownedPMatrices.Initialize(nRates, nStates, nStates);
		if (pMatrices == NULL)
			pMatrices = ownedPMatrices.ptr;
		}
	}

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|   Calls the swap method of the `state_time' data member, supplying as the argument other's state_time data member.
*/
void InternalData::swapStateTime(
  InternalData * other) /**< is a pointer to the other InternalData structure involved in the swap  */
    {
    state_time.swap(other->state_time);
    }
#endif

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of univents for site `site'. Assumes `unimap' is true and `size' is less than the length of the
|   `state_time' vector.
*/
unsigned InternalData::getNumUnivents(
  unsigned site) const      /**< is the site of interest */
    {
    if (!unimap || site >= (unsigned)state_time.size())
        return 0;
    return (unsigned)state_time[site].size();
    }
#endif

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of univent states for site `site'. Assumes `unimap' is true and `size' is less than the length of 
|   the `state_time' vector. This function is not particularly efficient, and it intended primarily for transferring
|   univent states to Python code for debugging purposes.
*/
std::vector<unsigned> InternalData::getUniventStates(
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
#endif

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of univent times for site `site'. Assumes `unimap' is true and `size' is less than the length of 
|   the `state_time' vector. This function is not particularly efficient, and it intended primarily for transferring
|   univent times to Python code for debugging purposes.
*/
std::vector<double> InternalData::getUniventTimes(
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
