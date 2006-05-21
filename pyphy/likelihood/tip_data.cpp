#include "phycas/force_include.h"
#include "CipresCommlib/CipresDataMatrixHelper.h"
#include "pyphy/phylogeny/basic_tree.hpp"
#include "pyphy/likelihood/likelihood_models.hpp"
#include "pyphy/likelihood/tree_likelihood.hpp"
#include "pyphy/likelihood/tip_data.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor for TipData objects that will only be used for simulating data and thus needs only to allocate the 
|	`state' and `pMatrixTranspose' data members.
*/
TipData::TipData(
  unsigned nRates,		/**< is the number of relative rate categories */
  unsigned nStates)		/**< is the number of states in the model */
	:
	parValidCLA(NULL),
	parCachedCLA(NULL),
	state(-1), 
	pMatrixTranspose(NULL)
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
  bool								managePMatrices) 	/**< if true, a 3D matrix will be allocated (if pMat is also NULL, the pMatrices will alias ownedPMatrices.ptr) */ 
	:
	parValidCLA(NULL),
	parCachedCLA(NULL),
	state(-1), 
	state_list_pos(stateListPosVec), 
	state_codes(stateCodesShPtr), 
	pMatrixTranspose(pMatTranspose)
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

}	// namespace phycas
