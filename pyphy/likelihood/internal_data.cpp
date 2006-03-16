#include "phycas/force_include.h"
#include "CipresCommlib/CipresDataMatrixHelper.h"
#include "pyphy/likelihood/likelihood_models.hpp"
#include "pyphy/likelihood/tree_likelihood.hpp"
#include "pyphy/likelihood/internal_data.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor initializes the conditional likelihood vectors using `nRates' and `nStates', and initializes the 
|	`pMattrices' data member using `pMat'. The parameter `managePMatrices' determines whether to allocate space for the 
|	transition probability matrices.
*/
InternalData::InternalData(
  unsigned				nPatterns,			/**< is the number of site patterns */
  unsigned				nRates,				/**< is the number of relative rate categories */
  unsigned				nStates,			/**< is the number of states in the model */
  double * * *			pMat,				/**< is an alias to the rates by states by states `pMatrix' array, and may be NULL */
  bool 					managePMatrices) 	/**< if true, a 3D matrix will be allocated (if `pMat' is also NULL, `pMatrices' will alias `ownedPMatrices.ptr') */ 
	: state(-1), cla(nPatterns*nRates*nStates), pMatrices(pMat)
	{
	if (managePMatrices)
		{
		ownedPMatrices.Initialize(nRates, nStates, nStates);
		if (pMatrices == NULL)
			pMatrices = ownedPMatrices.ptr;
		}
	}

}	// namespace phycas
