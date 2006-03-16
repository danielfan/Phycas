#include "phycas/force_include.h"
#include <fstream>
#include "phycas/rand/probability_distribution.hpp"
#include "phycas/rand/distribution_manager.hpp"
#include "phycas/modules/sampler/sampler.hpp"
#include "phycas/command/sample_settings.hpp"
#include "ncl/output/nxs_output.hpp"
#if defined(USING_SLICESAMPLER)
#	include "phycas/rand/slicesampler.h"
#endif
using std::ofstream;
using std::string;
using std::vector;
using ncl::endl;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::sprintf;
	using std::sqrt;
#endif
#include "ncl/output/temp_basic_output_operators.hpp"
	
class SamplerSummary
	{
	public:
		bool 			showAllSamples;
		unsigned		precision;
		double			expectedMean;
		double			expectedStdDev;
		double			sampleMean;
		double			sampleStdDev;		
		vector<double> 	samples;
	};

class MVSamplerSummary
	{
	public:
		bool 			showAllSamples;
		unsigned		precision;
		VecDbl			expectedMean;
		VecDbl			expectedStdDev;
		VecDbl			sampleMean;
		VecDbl			sampleStdDev;		
		unsigned		nvar;
		vector<VecDbl>	samples;
	};

template<class OUT_STREAM>
class GenericPrinterClass<kVerboseOutStyle, SampleSettings, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const SampleSettings & ss)
			{
			outStream << "  Sampling from: " << ss.distribution.GetValueInStringForm() << '\n';
			outStream << "  Drawing " << ss.nsampled << " samples\n";
			outStream << "  Seed is " << ss.seed << '\n';
			const NxsOutFilePath & fp = ss.outDestination.filePath; //@MTH BUG, should check for other types of output
			if (fp.WasOpened())
				{
				if (fp.IsAppending())
					outStream << "  Appending samples to file " << fp.GetFullName() << '\n';
				else 
					outStream << "  Saving samples to file " << fp.GetFullName() << '\n';
				}
			}
	};
	
template<class OUT_STREAM>
class GenericPrinterClass<kVerboseOutStyle, SamplerSummary, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const SamplerSummary & samples)
			{
			char format[32];
			const UInt precision = samples.precision;
			const unsigned colw = (precision < 9 ? 12 : precision + 3);
			sprintf(format, "%%%dd%%12.%df\n", colw, precision);
			std::string line;
			if (samples.showAllSamples)
				{
				unsigned i = 1;
				for (VecDbl::const_iterator sIt = samples.samples.begin(); sIt != samples.samples.end(); ++sIt, ++i)
					{
					StrPrintF(line, format, i, *sIt);
					outStream << line;
					line.clear();
					}
				}
			StrPrintF(line, "    sample size: %12d\n", (unsigned)samples.samples.size());
			StrPrintF(line, "    sample mean: %*.*f\n", colw, precision, samples.sampleMean);
			StrPrintF(line, "  expected mean: %*.*f\n", colw, precision, samples.expectedMean);
			StrPrintF(line, "    sample s.d.: %*.*f\n", colw, precision, samples.sampleStdDev);
			StrPrintF(line, "  expected s.d.: %*.*f", colw, precision, samples.expectedStdDev);
			outStream << line;
			}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kVerboseOutStyle, MVSamplerSummary, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM &outStream, const MVSamplerSummary &samples)
			{
			VecDbl::const_iterator it;
			char format[32];
			const unsigned precision = samples.precision;
			const unsigned colw = (precision < 9 ? 12 : precision + 3);
			sprintf(format, " %%%d.%df", colw, precision);
			std::string line;
			if (samples.showAllSamples)
				{
				unsigned i = 1;
				for (vector<VecDbl>::const_iterator sIt = samples.samples.begin(); sIt != samples.samples.end(); ++sIt, ++i)
					{
					StrPrintF(line, "%d", i);
					for (it = sIt->begin(); it != sIt->end(); ++it)
						{
						StrPrintF(line, format, *it);
						}
					outStream << line << "\n";
					line.clear();
					}
				}

			StrPrintF(line, "    sample size: %12d\n", (unsigned)samples.samples.size());

			StrPrintF(line, "    sample mean:");
			for (it = samples.sampleMean.begin(); it != samples.sampleMean.end(); ++it)
				{
				StrPrintF(line, format, *it);
				}
			outStream << line << "\n";
			line.clear();

			StrPrintF(line, "  expected mean:");
			for (it = samples.expectedMean.begin(); it != samples.expectedMean.end(); ++it)
				{
				StrPrintF(line, format, *it);
				}
			outStream << line << "\n";
			line.clear();

			StrPrintF(line, "    sample s.d.:");
			for (it = samples.sampleStdDev.begin(); it != samples.sampleStdDev.end(); ++it)
				{
				StrPrintF(line, format, *it);
				}
			outStream << line << "\n";
			line.clear();

			StrPrintF(line, "  expected s.d.:");
			for (it = samples.expectedStdDev.begin(); it != samples.expectedStdDev.end(); ++it)
				{
				StrPrintF(line, format, *it);
				}
			outStream << line << "\n";
			}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kTabSeparatedTabularOutStyle, SamplerSummary, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const SamplerSummary & samples)
			{
			char format[32];
			const UInt precision = samples.precision;
			sprintf(format, "%%d\t%%.%df\n", precision);
			std::string line;
			unsigned i = 1;
			for (VecDbl::const_iterator sIt = samples.samples.begin(); sIt != samples.samples.end(); ++sIt, ++i)
				{
				StrPrintF(line, format, i, *sIt);
				outStream << line;
				line.clear();
				}
			}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kTabSeparatedTabularOutStyle, MVSamplerSummary, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM &outStream, const MVSamplerSummary &samples)
			{
			char format[32];
			const unsigned precision = samples.precision;
			sprintf(format, "\t%%.%df", precision);
			std::string line;
			unsigned i = 1;
			for (vector<VecDbl>::const_iterator sIt = samples.samples.begin(); sIt != samples.samples.end(); ++sIt, ++i)
				{
				StrPrintF(line, "%d", i);

				// sIt is a multivariate sample, so must iterate through all the values
				//
				VecDbl::const_iterator it;
				for (it = sIt->begin(); it != sIt->end(); ++it)
					{
					StrPrintF(line, format, *it);
					}

				outStream << line << "\n";
				line.clear();
				}
			}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kTabSeparatedTabularLabelsOutStyle, SamplerSummary *, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const SamplerSummary * )
			{
			outStream << "#\trep\tvalue\n";
			}
	};

template<class OUT_STREAM>
class GenericPrinterClass<kTabSeparatedTabularLabelsOutStyle, MVSamplerSummary *, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const MVSamplerSummary *s )
			{
			assert(s->nvar > 0);
			outStream << "#\trep";
			for (unsigned i = 0; i < s->nvar; i++)
				{
				outStream << "\tvar_" << (i+1);
				}
			outStream << "\n";
			}
	};

/*----------------------------------------------------------------------------------------------------------------------
|	uses the distribution description (s->distribution) to create a ProbabilityDistribution and draw s->nsampled points
|	from the distribution.
|	reports summary statistics to the user.
*/
CmdResult PhoSampler::HandleSample(SampleSettings *s)
	{	
	NxsOutputStreamWrapperShPtr pdFileShPtr = NxsOutput::GetGenericOutputStreamShPtr(s->outDestination);
	NxsOutputStreamWrapper * outFilePtr = pdFileShPtr.get();
	NxsOutputStream * outStream = NxsOutput::GetOutputStreamPtr();
	Lot r(s->seed);
	if (s->seed == 0)
		s->seed = r.GetInitSeed(); //@@ need to discuss seed == 0 behaviour wrt command parsing
		
	if (outStream != NULL)
		Emit<kVerboseOutStyle>(*outStream, *s) << endl;
	if (outFilePtr)
		Emit<kTabSeparatedTabularLabelsOutStyle, SamplerSummary*>(*outFilePtr, NULL) << endl;

	NxsReturnValStream * returnValStream =  NxsOutput::GetReturnValueStreamPtr();

	if (s->distribution.IsMultivariate())
		{
		// Multivariate distribution 
		//
#		if defined(USING_SLICESAMPLER)
			throw XProbabilityDistribution("multivariate slice sampler not yet implemented");
#		endif

		MVProbDistPtr mvpDistr = s->distribution.CreateMultiVariateProbabilityDistribution();
		mvpDistr->SetLot(&r);

		MVSamplerSummary samples;
		samples.precision		= s->precision;
		samples.showAllSamples	= s->show;
		samples.nvar			= s->distribution.GetNVariates();

		VecDbl sum(samples.nvar, 0.0);
		VecDbl sumsq(samples.nvar, 0.0);
		for (unsigned i = 0; i < s->nsampled; ++i)
			{
			const VecDbl x = mvpDistr->Sample();
			samples.samples.push_back(x);
			for (unsigned k = 0; k < samples.nvar; ++k)
				{
				double x_k = x[k];
				sum[k] += x_k;
				sumsq[k] += x_k*x_k;
				}
			}

		samples.sampleMean.clear();
		samples.sampleStdDev.clear();
		for (unsigned k = 0; k < samples.nvar; ++k)
			{
			double n = (double)(s->nsampled);
			double m = sum[k]/s->nsampled;
			double ss = sumsq[k];
			samples.sampleMean.push_back(m);
			double v = ((n < 2.0) ? 0.0 : (ss - n*m*m)/(n - 1.0));
			samples.sampleStdDev.push_back(sqrt(v));
			}

		samples.expectedMean	= mvpDistr->GetMean();
		samples.expectedStdDev	= mvpDistr->GetStdDev();

		if (returnValStream)
			Emit<kVerboseOutStyle>(*returnValStream, samples) << endl;
			
		if (outFilePtr)
			Emit<kTabSeparatedTabularOutStyle>(*outFilePtr, samples);
		}
	else
		{
		// Univariate distribution 
		//
#		if defined(USING_SLICESAMPLER)
			SliceSampler slicer;
			slicer.SetRandomNumberGenerator(&r);
#		endif

		double sum		= 0.0;
		double sumsq	= 0.0;
		ProbDistPtr pDistr = s->distribution.CreateProbabilityDistribution();

#		if !defined(USING_SLICESAMPLER)
			pDistr->SetLot(&r);	//POL 22Jun2005
#		endif

		SamplerSummary samples;
		samples.precision		= s->precision;
		samples.showAllSamples	= s->show;

		for (unsigned i = 0; i < s->nsampled; ++i)
			{
#			if defined(USING_SLICESAMPLER)
				slicer.AttachDistribution(pDistr);
				const double x = slicer.GetNextSample();
#			else
				const double x = pDistr->Sample();	//POL 22Jun2005
#			endif
			sum		+= x;
			sumsq	+= x*x;
			samples.samples.push_back(x);
			}

		const double mean = sum/s->nsampled;
		const double sd =  sqrt(s->nsampled < 2 ? 0.0 : (sumsq - s->nsampled*mean*mean)/(s->nsampled - 1.0));

		samples.expectedMean	= pDistr->GetMean();
		samples.expectedStdDev	= pDistr->GetStdDev();
		samples.sampleMean		= mean;
		samples.sampleStdDev	= sd;

		if (returnValStream)
			Emit<kVerboseOutStyle>(*returnValStream, samples) << endl;
			
		if (outFilePtr)
			Emit<kTabSeparatedTabularOutStyle>(*outFilePtr, samples);
		}

	s->seed = (unsigned) r.GetSeed();
	return kCmdSucceeded;
	}
