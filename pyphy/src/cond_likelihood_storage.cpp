#include "pyphy/src/cond_likelihood.hpp"
#include "pyphy/src/cond_likelihood_storage.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Ensures that the `cl_stack' contains at least `capacity' CondLikelihoodShPtr objects.
*/
void CondLikelihoodStorage::fillTo(unsigned capacity)
	{
	PHYCAS_ASSERT(num_patterns > 0);
	PHYCAS_ASSERT(num_rates > 0);
	PHYCAS_ASSERT(num_states > 0);
	unsigned curr_sz = (unsigned)cl_stack.size();
	unsigned num_needed = (capacity > curr_sz ? capacity - curr_sz : 0);
	for (unsigned i = 0; i < num_needed; ++i)
		{
		cl_stack.push(CondLikelihoodShPtr(new CondLikelihood(num_patterns, num_rates, num_states)));
		num_created++;
		}
	}

} // namespace phycas
