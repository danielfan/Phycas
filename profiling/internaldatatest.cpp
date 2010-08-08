// Creates an InternalData object and then immediately deletes it in order to test for memory leaks.

#include <iostream>
#include <cmath>
#include <map>
#include <stack>
#include <vector>
#include <stdint.h>
#include "ncl/nxsallocatematrix.h"          # for ScopedThreeDMatrix declaration
#include "phycas/src/cond_likelihood.hpp"
#include "phycas/src/cond_likelihood_storage.hpp"
#include "phycas/src/internal_data.hpp"
#include "memchk.hpp"

namespace phycas
{

class InternalDataTest
    {
    public:
        unsigned num_patterns;
        unsigned num_rates;
        unsigned num_states;
        unsigned num_reps;
        CondLikelihoodStorageShPtr cla_pool;
        
        InternalDataTest(unsigned npat, unsigned nrat, unsigned nstat, unsigned nreps)
            {
            num_patterns = npat;
            num_rates    = nrat;
            num_states   = nstat;
            num_reps     = nreps;
            cla_pool     = CondLikelihoodStorageShPtr(new CondLikelihoodStorage());
            }
            
        void run()
            {
            unsigned n = 0;
            for (unsigned rep = 0; rep < num_reps; ++rep)
                {
                InternalData * pid = new InternalData(false,    // using unimap
                                    num_patterns,				// number of site patterns
                                    num_rates,					// number of relative rate categories
                                    num_states,					// number of model states
                                    NULL,						// pMat
                                    true,						// managePMatrices
                                    cla_pool);
                                    
                delete pid;
                }
            }
    };

}   // namespace phycas 

int main()
	{
	CREATE_MEMCHK
	unsigned num_patterns = 10;
	unsigned num_rates    = 4;
	unsigned num_states   = 4;
	unsigned num_reps     = 1;
	phycas::InternalDataTest(num_patterns, num_rates, num_states, num_reps).run();
	MEMCHK_REPORT(std::cerr)
	return 0;
	}
