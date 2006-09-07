//#include "phycas/force_include.h"
#include "phypy/src/cipres/CipresDataMatrixHelper.h"
#include "phypy/src/basic_tree.hpp"
#include "phypy/src/likelihood_models.hpp"
#include "phypy/src/tree_likelihood.hpp"
#include "phypy/src/tip_data.hpp"
#include "phypy/src/cond_likelihood.hpp"

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

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor for TipData objects that will be used in likelihood calcuations and thus need to store the observed
|	data for the tip as well as information about local state codes.
*/
TipData::TipData(
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
	cla_pool(cla_storage)
	{
	const unsigned nObservedStates = nStates + 1 + state_list_pos.size();
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

}	// namespace phycas
