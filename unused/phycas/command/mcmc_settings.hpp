/*########## mcmc_settings.hpp ##########*/
#if !defined (PHYC_MCMC_SETTINGS_HPP)
#define PHYC_MCMC_SETTINGS_HPP
#include "ncl/output/nxs_output_destination_description.hpp"
class MCMCSettings
	{
	public:
		unsigned        	niterations;
		unsigned        	nchains;
		unsigned        	sampleEvery;
		unsigned        	reportEvery;
		unsigned        	plotEvery;
		NxsOutputDestinationDescription	paramOut;
		NxsOutputDestinationDescription	splitOut;
		unsigned        	rnseed;
		bool            	ignoreData;
		double          	bushMoveWeight;
		unsigned        	gibbsEvery;
		unsigned        	sliceDivisions;
		double          	initUnitWidth;
		double          	scalerMoveWeight;
		double          	localMoveWeight;
		double          	kappaStartingValue;
		double          	kappaPriorMean;
		unsigned        	ncat;
		double          	pTopoPrior;
		bool            	topoPriorFlat;
		bool            	useEdgeLenHyper;
		double          	edgeLenHyperPriorMean;
		double          	edgeLenHyperPriorVar;
		bool            	resClassPrior;
		double          	polytomyLnPriorRatio;
		
		MCMCSettings()
			:niterations(1100),
			nchains(4),
			sampleEvery(10),
			reportEvery(100),
			plotEvery(500),
			paramOut(),
			splitOut(),
			rnseed(0),
			ignoreData(false),
			bushMoveWeight(1),
			gibbsEvery(10),
			sliceDivisions(3),
			initUnitWidth(50.0),
			scalerMoveWeight(0.0),
			localMoveWeight(2.0),
			kappaStartingValue(4.0),
			kappaPriorMean(1.0),
			ncat(1),
			pTopoPrior(0.5),
			topoPriorFlat(true),
			useEdgeLenHyper(true),
			edgeLenHyperPriorMean(1),
			edgeLenHyperPriorVar(10),
			resClassPrior(false),
			polytomyLnPriorRatio(0.0)
			{}
	};
#endif // if !defined (PHYC_MCMC_SETTINGS_HPP)
/*%%%%%% /mcmc_settings.hpp %%%%%%*/
