#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle

//#include "phycas/force_include.h"
#if defined(POL_PHYCAS)
#	include "pyphy/src/basic_lot.hpp"
#else
#	include "phycas/rand/lot.hpp"
#endif
#include "pyphy/src/probability_distribution.hpp"
#include "pyphy/src/slice_sampler.hpp"
#include <boost/python.hpp>
#include "pyphy/src/xprobdist.hpp"

using namespace boost::python;

double getEffectiveLnZero()
	{
	return -DBL_MAX;
	}

// The following wrapper struct is needed because we will potentially be deriving Python classes from AdHocDensity
// This struct is thus only necessary if we plan to do something like this in Python:
//		MyFunctor(AdhocDensity):
//			...
// See http://wiki.python.org/moin/boost.python/InternalDataStructures
//
struct AdHocDensityWrapper : AdHocDensity 
	{
	AdHocDensityWrapper(PyObject * p) : self(p) {}
	virtual ~AdHocDensityWrapper()
		{
		//std::cerr << "(ProbDist)AdHocDensityWrapper dying..." << std::endl;
		}
	double operator()(double x) { return boost::python::call_method<double>(self, "__call__", x); }
    PyObject * self;
	};

void translateXProbDist(const XProbDist &e)
	{
    // Use the Python 'C' API to set up an exception object
    PyErr_SetString(PyExc_Exception, e.what());
    }

BOOST_PYTHON_MODULE(_ProbDist)
{
	def("getEffectiveLnZero", getEffectiveLnZero);

	// these lines required by num_util
	import_array();
	numeric::array::set_module_and_type("numarray", "NDArray");

	class_<AdHocDensity, boost::noncopyable, boost::shared_ptr<AdHocDensityWrapper> >("AdHocDensityBase")
		;

	class_<ProbabilityDistribution, bases<AdHocDensity>, boost::shared_ptr<ProbabilityDistribution>, boost::noncopyable>("ProbabilityDistribution", no_init)
		;

	class_<SliceSampler, boost::shared_ptr<SliceSampler> >("SliceSamplerBase")
		.def(init<boost::shared_ptr<phycas::Lot>, boost::shared_ptr<AdHocDensity> >())
		.def("attachFunc", &SliceSampler::AttachFunc)
		.def("sample", &SliceSampler::Sample)
		.def("debugSample", &SliceSampler::DebugSample)
		.def("overrelaxedSample", &SliceSampler::OverrelaxedSample)
		.def("debugOverrelaxedSample", &SliceSampler::DebugOverrelaxedSample)
		.def("attachRandomNumberGenerator", &SliceSampler::AttachRandomNumberGenerator)
		.def("adaptSimple", &SliceSampler::AdaptSimple)
		.def("adaptNeal", &SliceSampler::AdaptNeal)
		.def("getSliceUnitWidth", &SliceSampler::GetSliceUnitWidth)
		.def("getMinX", &SliceSampler::GetMinX)
		.def("getMaxX", &SliceSampler::GetMaxX)
		.def("setMaxUnits", &SliceSampler::SetMaxUnits)
		.def("getOrigLeftEdgeOfSlice", &SliceSampler::GetOrigLeftEdgeOfSlice)
		.def("getOrigRightEdgeOfSlice", &SliceSampler::GetOrigRightEdgeOfSlice)
		.def("getLeftEdgeOfSlice", &SliceSampler::GetLeftEdgeOfSlice)
		.def("getRightEdgeOfSlice", &SliceSampler::GetRightEdgeOfSlice)
		.def("getSliceYValue", &SliceSampler::GetSliceYValue)
		.def("getNumFuncEvals", &SliceSampler::GetNumFuncEvals)
		.def("getNumFailedSamples", &SliceSampler::GetNumFailedSamples)
		.def("getNumUnitsRequired", &SliceSampler::GetNumUnitsRequired)
		.def("getNumSamples", &SliceSampler::GetNumSamples)
		.def("resetDiagnostics", &SliceSampler::ResetDiagnostics)
		.def("adaptYConditional", &SliceSampler::AdaptYConditional)
		.def("calcW", &SliceSampler::CalcW)
		.def("getMode", &SliceSampler::GetMode)
		.def("getLnDensityAtMode", &SliceSampler::GetLnDensityAtMode)
		.def("getLastSampledXValue", &SliceSampler::GetLastSampledXValue)
		.def("getLastSampledYValue", &SliceSampler::GetLastSampledYValue)
		.def("setXValue", &SliceSampler::SetXValue)
		;

//We tell boost::python the smart pointer type we're using, like this:
//class_<DrawableInterface, DrawablePtr, boost::noncopyable>
//("DrawableInterface", no_init)
//    ;
//where DrawablePtr is a typedef for the smart pointer type. 

	class_<phycas::Lot, boost::shared_ptr<phycas::Lot>, boost::noncopyable>("LotBase", init<unsigned>())
		.def("getSeed", &phycas::Lot::GetSeed)
		.def("setSeed", &phycas::Lot::SetSeed)
		.def("getInitSeed", &phycas::Lot::GetInitSeed)
		.def("uniform", &phycas::Lot::Uniform)
		.def("getrandbits", &phycas::Lot::GetRandBits)
		;

	class_<DirichletDistribution>("DirichletDistBase")
		.def(init<const std::vector<double> &>())
		.def("isDiscrete", &DirichletDistribution::IsDiscrete)
		.def("getDistName", &DirichletDistribution::GetDistributionName)
		.def("__str__", &DirichletDistribution::GetDescriptionForPython)
		.def("__repr__", &DirichletDistribution::GetDescriptionForPython)
		.def("setLot", &DirichletDistribution::SetLot)
		.def("setSeed", &DirichletDistribution::SetSeed)
		.def("resetLot", &DirichletDistribution::ResetLot)
		.def("getMean", &DirichletDistribution::GetMean)
		.def("getVar", &DirichletDistribution::GetVar)
		.def("getStdDev", &DirichletDistribution::GetStdDev)
		.def("approxCDF", &DirichletDistribution::ApproxCDF)
		.def("sample", &DirichletDistribution::Sample)
		.def("getLnPDF", &DirichletDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &DirichletDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &DirichletDistribution::SetMeanAndVariance)
		.def("getVarCovarMatrix", &DirichletDistribution::GetVarCovarMatrix)
		.def("getNParams", &DirichletDistribution::GetNParams)
		;
	class_<BetaDistribution, bases<ProbabilityDistribution, AdHocDensity> >("BetaDistBase")
		.def(init<double, double>())
		.def("isDiscrete", &BetaDistribution::IsDiscrete)
		.def("getDistName", &BetaDistribution::GetDistributionName)
		.def("__str__", &BetaDistribution::GetDistributionDescription)
		.def("__repr__", &BetaDistribution::GetDistributionDescription)
		.def("setLot", &BetaDistribution::SetLot)
		.def("setSeed", &BetaDistribution::SetSeed)
		.def("resetLot", &BetaDistribution::ResetLot)
		.def("getMean", &BetaDistribution::GetMean)
		.def("getVar", &BetaDistribution::GetVar)
		.def("getStdDev", &BetaDistribution::GetStdDev)
		.def("getCDF", &BetaDistribution::GetCDF)
		.def("sample", &BetaDistribution::Sample)
		.def("getLnPDF", &BetaDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &BetaDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &BetaDistribution::SetMeanAndVariance)
		;
	class_<BernoulliDistribution, bases<ProbabilityDistribution> >("BernoulliDistBase")
		.def(init<double>())
		.def("isDiscrete", &BernoulliDistribution::IsDiscrete)
		.def("getDistName", &BernoulliDistribution::GetDistributionName)
		.def("__str__", &BernoulliDistribution::GetDistributionDescription)
		.def("__repr__", &BernoulliDistribution::GetDistributionDescription)
		.def("setLot", &BernoulliDistribution::SetLot)
		.def("setSeed", &BernoulliDistribution::SetSeed)
		.def("resetLot", &BernoulliDistribution::ResetLot)
		.def("getMean", &BernoulliDistribution::GetMean)
		.def("getVar", &BernoulliDistribution::GetVar)
		.def("getStdDev", &BernoulliDistribution::GetStdDev)
		.def("getCDF", &BernoulliDistribution::GetCDF)
		.def("sample", &BernoulliDistribution::Sample)
		.def("getLnPDF", &BernoulliDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &BernoulliDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &BernoulliDistribution::SetMeanAndVariance)
		;
	class_<BinomialDistribution, bases<ProbabilityDistribution> >("BinomialDistBase")
		.def(init<double, double>())
		.def("isDiscrete", &BinomialDistribution::IsDiscrete)
		.def("getDistName", &BinomialDistribution::GetDistributionName)
		.def("__str__", &BinomialDistribution::GetDistributionDescription)
		.def("__repr__", &BinomialDistribution::GetDistributionDescription)
		.def("setLot", &BinomialDistribution::SetLot)
		.def("setSeed", &BinomialDistribution::SetSeed)
		.def("resetLot", &BinomialDistribution::ResetLot)
		.def("getMean", &BinomialDistribution::GetMean)
		.def("getVar", &BinomialDistribution::GetVar)
		.def("getStdDev", &BinomialDistribution::GetStdDev)
		.def("getCDF", &BinomialDistribution::GetCDF)
		.def("sample", &BinomialDistribution::Sample)
		.def("getLnPDF", &BinomialDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &BinomialDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &BernoulliDistribution::SetMeanAndVariance)
		;
	class_<UniformDistribution, bases<ProbabilityDistribution, AdHocDensity> >("UniformDistBase")
		.def(init<double, double>())
		.def("isDiscrete", &UniformDistribution::IsDiscrete)
		.def("getDistName", &UniformDistribution::GetDistributionName)
		.def("__str__", &UniformDistribution::GetDistributionDescription)
		.def("__repr__", &UniformDistribution::GetDistributionDescription)
		.def("setLot", &UniformDistribution::SetLot)
		.def("setSeed", &UniformDistribution::SetSeed)
		.def("resetLot", &UniformDistribution::ResetLot)
		.def("getMean", &UniformDistribution::GetMean)
		.def("getVar", &UniformDistribution::GetVar)
		.def("getStdDev", &UniformDistribution::GetStdDev)
		.def("getCDF", &UniformDistribution::GetCDF)
		.def("sample", &UniformDistribution::Sample)
		.def("getLnPDF", &UniformDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &UniformDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &UniformDistribution::SetMeanAndVariance)
		;
	class_<GammaDistribution, bases<ProbabilityDistribution, AdHocDensity> >("GammaDistBase")
		.def(init<double, double>())
		.def("isDiscrete", &GammaDistribution::IsDiscrete)
		.def("getDistName", &GammaDistribution::GetDistributionName)
		.def("__str__", &GammaDistribution::GetDistributionDescription)
		.def("__repr__", &GammaDistribution::GetDistributionDescription)
		.def("setLot", &GammaDistribution::SetLot)
		.def("setSeed", &GammaDistribution::SetSeed)
		.def("resetLot", &GammaDistribution::ResetLot)
		.def("getMean", &GammaDistribution::GetMean)
		.def("getVar", &GammaDistribution::GetVar)
		.def("getStdDev", &GammaDistribution::GetStdDev)
		.def("getCDF", &GammaDistribution::GetCDF)
		.def("sample", &GammaDistribution::Sample)
		.def("getLnPDF", &GammaDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &GammaDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &GammaDistribution::SetMeanAndVariance)
		;
	class_<ExponentialDistribution, bases<GammaDistribution, ProbabilityDistribution, AdHocDensity> >("ExponentialDistBase")
		.def(init<double>())
		.def("isDiscrete", &ExponentialDistribution::IsDiscrete)
		.def("getDistName", &ExponentialDistribution::GetDistributionName)
		.def("__str__", &ExponentialDistribution::GetDistributionDescription)
		.def("__repr__", &ExponentialDistribution::GetDistributionDescription)
		.def("setLot", &ExponentialDistribution::SetLot)
		.def("setSeed", &ExponentialDistribution::SetSeed)
		.def("resetLot", &ExponentialDistribution::ResetLot)
		.def("getMean", &ExponentialDistribution::GetMean)
		.def("getVar", &ExponentialDistribution::GetVar)
		.def("getStdDev", &ExponentialDistribution::GetStdDev)
		.def("getCDF", &ExponentialDistribution::GetCDF)
		.def("sample", &ExponentialDistribution::Sample)
		.def("getLnPDF", &ExponentialDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &ExponentialDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &ExponentialDistribution::SetMeanAndVariance)
		;
	class_<InverseGammaDistribution, bases<ProbabilityDistribution, AdHocDensity> >("InverseGammaDistBase")
		.def(init<double, double>())
		.def("isDiscrete", &InverseGammaDistribution::IsDiscrete)
		.def("getDistName", &InverseGammaDistribution::GetDistributionName)
		.def("__str__", &InverseGammaDistribution::GetDistributionDescription)
		.def("__repr__", &InverseGammaDistribution::GetDistributionDescription)
		.def("setLot", &InverseGammaDistribution::SetLot)
		.def("setSeed", &InverseGammaDistribution::SetSeed)
		.def("resetLot", &InverseGammaDistribution::ResetLot)
		.def("getMean", &InverseGammaDistribution::GetMean)
		.def("getVar", &InverseGammaDistribution::GetVar)
		.def("getStdDev", &InverseGammaDistribution::GetStdDev)
		.def("getCDF", &InverseGammaDistribution::GetCDF)
		.def("sample", &InverseGammaDistribution::Sample)
		.def("getLnPDF", &InverseGammaDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &InverseGammaDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &InverseGammaDistribution::SetMeanAndVariance)
		;
	register_exception_translator<XProbDist>(&translateXProbDist);
}