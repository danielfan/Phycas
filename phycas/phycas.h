#if !defined(PHYCAS_H)
#define PHYCAS_H

template<typename T>
unsigned GetLongevity(const T * const p);
class DistributionManager;
template<>
inline unsigned GetLongevity<DistributionManager>(const DistributionManager * const)
	{
	return kNxsDistributionManagerLongevity;
	}

class PWD;
template<>
inline unsigned GetLongevity<PWD>(const PWD * const)
	{
	return 1;
	}


#include "ncl/nxs_defs.hpp"
//USING_VARIANCE_FOR_RATEHET as opposed to using the shape parameter of the gamma distribution for rate het.
#define USING_VARIANCE_FOR_RATEHET
//	If USING_NODE_NUMBERS_ONLY is defined, tree Nodes do not have a std::string data member that stores their name
//	an external class must be used to translate from the taxon number to a name.
#define USING_NODE_NUMBERS_ONLY
typedef char NStateInt;	// to be used in code where speed is critical, but we need a integer type big enough to hold the max num states
						// signed because - numbers are used to indicate ambiguity
#define NSTATEINT_MAX CHAR_MAX

//	USING_PROGRAMMER_SPECIFIC_DEFINES added so that we can define 
//	macros in our projects without touching phorest.h (I keep forgetting that I've tweaked
//	this file and committing changes).
// 	the !defined(USING_PROGRAMMER_SPECIFIC_DEFINES) section below provides the 
//	stable defaults.
//
#if !defined(USING_PROGRAMMER_SPECIFIC_DEFINES)
#		define PHYLODIVERSITY_MODULE
#endif
#define SCM_MODULE
#endif
