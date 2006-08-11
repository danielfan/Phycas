/*########## sample_settings.hpp ##########*/
#if !defined (PHYC_SAMPLE_SETTINGS_HPP)
#define PHYC_SAMPLE_SETTINGS_HPP
#include "phycas/rand/distribution_description.hpp"
#include "ncl/output/nxs_output_destination_description.hpp"
class SampleSettings
	{
	public:
		DistributionDescription	distribution;
		unsigned        	nsampled;
		bool            	show;
		unsigned        	seed;
		NxsOutputDestinationDescription	outDestination;
		unsigned        	precision;
		
		SampleSettings()
			:distribution(),
			nsampled(100),
			show(false),
			seed(0),
			outDestination(),
			precision(6)
			{}
	};
#endif // if !defined (PHYC_SAMPLE_SETTINGS_HPP)
/*%%%%%% /sample_settings.hpp %%%%%%*/
