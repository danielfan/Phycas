#include "phycas/force_include.h"
#include "phycas/likelihood/tree_node_attributes.hpp"
/*----------------------------------------------------------------------------------------------------------------------
|	Creates an attribute to be used by HKYAdHocEvaluator class when computing likelihoods of trees. Allocates `condLike'
|	array of length 4*`nsites'*`nrates'. Allocates `prMatrices' and initializes `prMatrices[0]' (no partitioning 
|	allowed at the moment).
*/
TreeNodeAttribute::TreeNodeAttribute(
  unsigned nsites,		/**< is the number of patterns */
  unsigned nrates,		/**< is the number of rate categories to be used with discrete gamma rate heterogeneity */
  double *blen)			/**< is a pointer to an edge length */
	{
	assert(blen == NULL);
#	if defined(HAVE_PRAGMA_UNUSED) && defined (NDEBUG)
#		pragma unused(blen)
#	endif
	assert(nsites > 0);
	assert(nrates > 0);

	lastAvoidNode = NULL;
	pMatIsDirty = true;
	clDirty = true;

	// Allocate memory for and initialize condLike array
	//
	unsigned total_len = 4*nsites*nrates;
	condLike = new double[total_len];
	for (unsigned i = 0; i < total_len; ++i)
		condLike[i] = 1.0;

	//@ not using edgeLen
	edgeLen = NULL;

	// Allocate memory for and initialize transition probability matrices
	//
	//@POL Allocating one prMatrix for each rate category, but not accommodating character partitioning yet
	prMatrices = new double**[nrates];
	for (unsigned r = 0; r < nrates; ++r)
		{
		prMatrices[r] = new double*[4];

		prMatrices[r][0] = new double[4];
		prMatrices[r][1] = new double[4];
		prMatrices[r][2] = new double[4];
		prMatrices[r][3] = new double[4];

		prMatrices[r][0][0] = 1.0;
		prMatrices[r][0][1] = 0.0;
		prMatrices[r][0][2] = 0.0;
		prMatrices[r][0][3] = 0.0;

		prMatrices[r][1][0] = 0.0;
		prMatrices[r][1][1] = 1.0;
		prMatrices[r][1][2] = 0.0;
		prMatrices[r][1][3] = 0.0;

		prMatrices[r][2][0] = 0.0;
		prMatrices[r][2][1] = 0.0;
		prMatrices[r][2][2] = 1.0;
		prMatrices[r][2][3] = 0.0;

		prMatrices[r][3][0] = 0.0;
		prMatrices[r][3][1] = 0.0;
		prMatrices[r][3][2] = 0.0;
		prMatrices[r][3][3] = 1.0;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Deletes `condLike' array and `prMatrices'.
*/
TreeNodeAttribute::~TreeNodeAttribute()
	{
	assert(condLike != NULL);

	delete [] condLike;

	delete [] prMatrices[0][0];
	delete [] prMatrices[0][1];
	delete [] prMatrices[0][2];
	delete [] prMatrices[0][3];
	delete [] prMatrices[0];
	delete [] prMatrices;
	}

