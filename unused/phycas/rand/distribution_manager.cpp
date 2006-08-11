#include "phycas/force_include.h"
#include "phycas/rand/distribution_manager.hpp"
#include "phycas/command/distribution_settings.hpp"
#include "ncl/output/nxs_output.hpp"
#include "ncl/output/temp_basic_output_operators.hpp"
using ncl::endl;

DistributionManager::DistributionManager()
	:reservedDistributionNames(SplitString("Uniform|Bernoulli|Binomial|Exponential|Gamma|InvGamma|InverseGamma|Dirichlet|Beta", '|'))
	{
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Adds (or replaces) a distribution definition and alerts the user.  Always returns kCmdSucceeded.
*/
CmdResult DistributionManager::HandleNewDistribution(DistributionSettings *s)
	{
	NxsOutputStream * outStream = NxsOutput::GetOutputStreamPtr();
	if (outStream != NULL)
		{
		if (GetDistribution(s->name) == NULL)
			*outStream << "Defining ";
		else
			*outStream << "Redefining ";
		*outStream << "the distribution " << s->name << " as ";
		Emit<kConciseOutStyle>(*outStream, s->distrib) << endl;
		}
	AddDistribution(s->name, s->distrib);
	return kCmdSucceeded;
	}


