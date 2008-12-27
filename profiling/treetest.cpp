// Repeatedly creates a Tree object and a TreeManip object, builds a tree using TreeManip's equiprobTree function, 
// then deletes the Tree and TreeManip objects. Designed to find memory leaks involved in Tree and TreeNode objects.
// Note that no InternalData or TipData objects are attached to the TreeNode objects in this test, so if the memory
// leak involves these objects this test will not find the leak.
//
// Note: the program is designed to loop indefinitely so that you can use Activity Monitor (MacOSX) or Task Manager
// (Windows) to observe memory usage. You should be able to stop it using Ctrl-C but, failing that, you can force-quit
// it using Activity Monitor or Task Manager.

#include <iostream>
#include <cmath>
#include <map>
#include <stack>
#include <vector>
#include <stdint.h>
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/basic_lot.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/tree_manip.hpp"
#include "memchk.hpp"

namespace phycas
{

class TreeTest
    {
    public:
        unsigned nreps; 
        
        TreeTest(unsigned num_reps)
            {
            nreps = num_reps;
            }
            
        void run()
            {
            unsigned n = 0;
            LotShPtr r = LotShPtr(new Lot());
            ProbDistShPtr prior = ProbDistShPtr(new GammaDistribution(0.5, 2.0));
            for (unsigned rep = 0; rep < nreps; ++rep)
                {
                TreeShPtr t = TreeShPtr(new Tree());
                TreeManip tm(t);
                tm.equiprobTree(10 /*tips*/, r, prior, prior);
                t.reset();
                }
            }
    };
    
}   // namespace phycas 

int main()
	{
	CREATE_MEMCHK
	unsigned num_reps = 1;
	phycas::TreeTest(num_reps).run();
	MEMCHK_REPORT(std::cerr)
	return 0;
	}
