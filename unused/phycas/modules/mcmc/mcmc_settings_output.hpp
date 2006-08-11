#if !defined(MCMC_SETTINGS_OUTPUT_HPP)
#define MCMC_SETTINGS_OUTPUT_HPP
#include "ncl/output/temp_basic_output_operators.hpp"
#include "phycas/command/mcmc_settings.hpp"
#include "ncl/misc/nxs_file_path.hpp"

template<class OUT_STREAM>
class GenericPrinterClass<kVerboseOutStyle, MCMCSettings, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const MCMCSettings & settings)
			{
			const NxsOutFilePath & fp = settings.paramOut.filePath; //@MTH BUG, should check for other types of output
			outStream << "\n  Number of iterations:               " << settings.niterations;
			outStream << "\n  Number of chains:                   " << settings.nchains;
			if (settings.reportEvery > 0)
				outStream << "\n  Reporting frequency:                " << settings.reportEvery;
			else
				outStream << "\n  Reporting frequency:                " << "never";
			if (settings.sampleEvery > 0)
				outStream << "\n  Sampling frequency:                 " << settings.sampleEvery;
			else
				outStream << "\n  Sampling frequency:                 " << "no samples will be taken";

			outStream << "\n  Parameter file name:                " << fp.GetFileName();
			outStream << " (" << (fp.Exists() ? "exists" : "does not exist") << ')';
			outStream << "\n  Starting random number seed:        " << settings.rnseed;
			outStream << "\n  Local move weight:                  " << settings.localMoveWeight;
			outStream << "\n  Bush move weight:                   " << settings.bushMoveWeight;
			outStream << "\n  Scaler move weight:                 " << settings.scalerMoveWeight;
			if (settings.ignoreData)
				outStream << "\n  Ignoring data (posterior will equal prior)";
			else
				outStream << "\n  Data not ignored";
			}
	};
#endif