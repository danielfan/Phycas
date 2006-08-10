//#include "phycas/force_include.h"
#include "pyphy/src/cipres/CipresDataMatrixHelper.h"
#include "pyphy/src/likelihood_models.hpp"
#include "pyphy/src/tree_likelihood.hpp"
#include "pyphy/src/internal_data.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor initializes the conditional likelihood vectors using `nRates' and `nStates', and initializes the 
|	`pMatrices' data member using `pMat'. The parameter `managePMatrices' determines whether to allocate space for the 
|	transition probability matrices.
*/
InternalData::InternalData(
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
	cla_pool(cla_storage)
	{
	if (managePMatrices)
		{
		ownedPMatrices.Initialize(nRates, nStates, nStates);
		if (pMatrices == NULL)
			pMatrices = ownedPMatrices.ptr;
		}
	}

}	// namespace phycas
