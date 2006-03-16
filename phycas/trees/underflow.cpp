#include "phycas/force_include.h"
#include "phycas/trees/underflow.hpp"
#include "ncl/misc/string_extensions.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	Creates string representation of `vBounces' vector for debugging purposes.
*/
std::string UnderflowManager::DebugCreateBounceVectorRepresentation()
	{
	std::string s;
	std::vector<unsigned>::iterator it;
	for (it = vBounces.begin(); it != vBounces.end(); it++)
		{
		unsigned v = (*it);
		if (v == UINT_MAX)
			StrPrintF(s, ".");
		else
			StrPrintF(s, "[%d]", v);
		}
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns average number of bounces (only using those elements whose value is less than UINT_MAX). Useful for seeing
|	hou much bouncing is going on.
*/
double UnderflowManager::DebugGetAvgNumBounces()
	{
	unsigned n = 0;
	double sum = 0.0;
	std::vector<unsigned>::iterator it;
	for (it = vBounces.begin(); it != vBounces.end(); it++)
		{
		unsigned v = (*it);
		if (v < UINT_MAX)
			{
			n++;
			sum += (double)v;
			}
		}
	if (n > 0)
		return (sum/(double)n);
	else
		return 0.0;
	}

