#ifndef PHO_SAMPLER_H
#define PHO_SAMPLER_H

class NxsCommandManager;
class DistributionManager;
class SampleSettings;

/*----------------------------------------------------------------------------------------------------------------------
|	class that handles the Sample command.  Only datamembers are aliases to output and distribution manager
*/
class PhoSampler
	{
	public:
		CmdResult 	HandleSample(SampleSettings *s);
		void 		SetupSampleCommand(NxsCommandManager *cmdMgr);
	private:	
		PhoSampler(){}
		friend PhoSamplerCreator;
	};
#endif
