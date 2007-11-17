/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#if defined(USING_NUMARRAY)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#endif

#include <boost/python.hpp>

#include "basic_lot.hpp"
#include "probability_distribution.hpp"

#if POLPY_NEWWAY
#include "stop_watch.hpp"
#endif

#include "slice_sampler.hpp"
#include "xprobdist.hpp"

using namespace boost::python;
using namespace phycas;

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

#if defined(USING_NUMARRAY)
	// these lines required by num_util
	import_array();
	numeric::array::set_module_and_type("numarray", "NDArray");
#endif

	class_<AdHocDensity, boost::noncopyable, boost::shared_ptr<AdHocDensityWrapper> >("AdHocDensityBase")
		;

	class_<ProbabilityDistribution, bases<AdHocDensity>, boost::shared_ptr<ProbabilityDistribution>, boost::noncopyable>("ProbabilityDistribution", no_init)
#if POLPY_NEWWAY
        .def("lnGamma", &ProbabilityDistribution::LnGamma)
#endif
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
#if POLPY_NEWWAY
		.def("useDoublingMethod", &SliceSampler::UseDoublingMethod)
#endif
		;

#if POLPY_NEWWAY
	class_<phycas::StopWatch, boost::shared_ptr<phycas::StopWatch>, boost::noncopyable>("StopWatchBase")
		.def("start", &phycas::StopWatch::start)
		.def("stop", &phycas::StopWatch::stop)
		.def("reset", &phycas::StopWatch::reset)
		.def("normalize", &phycas::StopWatch::normalize)
		.def("elapsedSeconds", &phycas::StopWatch::elapsedSeconds)
		.def("stopTicks", &phycas::StopWatch::stopTicks)
		.def("doofus", &phycas::StopWatch::doofus)
		;
#endif

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
#if POLPY_NEWWAY
		.def(init<const DirichletDistribution &>())
		.def("clone", &DirichletDistribution::Clone, return_value_policy<manage_new_object>())
#endif
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
#if POLPY_NEWWAY
		.def(init<const BetaDistribution &>())
		.def("clone", &BetaDistribution::Clone, return_value_policy<manage_new_object>())
#endif
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
#if POLPY_NEWWAY
		.def(init<const BernoulliDistribution &>())
		.def("clone", &BernoulliDistribution::Clone, return_value_policy<manage_new_object>())
#endif
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
#if POLPY_NEWWAY
		.def(init<const BinomialDistribution &>())
		.def("clone", &BinomialDistribution::Clone, return_value_policy<manage_new_object>())
#endif
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
	class_<ImproperUniformDistribution, bases<ProbabilityDistribution, AdHocDensity> >("ImproperUniformDistBase")
#if POLPY_NEWWAY
		.def(init<const ImproperUniformDistribution &>())
		.def("clone", &ImproperUniformDistribution::Clone, return_value_policy<manage_new_object>())
#endif
		.def("isDiscrete", &ImproperUniformDistribution::IsDiscrete)
		.def("getDistName", &ImproperUniformDistribution::GetDistributionName)
		.def("__str__", &ImproperUniformDistribution::GetDistributionDescription)
		.def("__repr__", &ImproperUniformDistribution::GetDistributionDescription)
		.def("setLot", &ImproperUniformDistribution::SetLot)
		.def("setSeed", &ImproperUniformDistribution::SetSeed)
		.def("resetLot", &ImproperUniformDistribution::ResetLot)
		.def("getMean", &ImproperUniformDistribution::GetMean)
		.def("getVar", &ImproperUniformDistribution::GetVar)
		.def("getStdDev", &ImproperUniformDistribution::GetStdDev)
		.def("getCDF", &ImproperUniformDistribution::GetCDF)
		.def("sample", &ImproperUniformDistribution::Sample)
		.def("getLnPDF", &ImproperUniformDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &ImproperUniformDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &ImproperUniformDistribution::SetMeanAndVariance)
		;
	class_<UniformDistribution, bases<ProbabilityDistribution, AdHocDensity> >("UniformDistBase")
		.def(init<double, double>())
#if POLPY_NEWWAY
		.def(init<const UniformDistribution &>())
		.def("clone", &UniformDistribution::Clone, return_value_policy<manage_new_object>())
#endif
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
#if POLPY_NEWWAY
		.def(init<const GammaDistribution &>())
		.def("clone", &GammaDistribution::Clone, return_value_policy<manage_new_object>())
#endif
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
#if POLPY_NEWWAY
		.def(init<const ExponentialDistribution &>())
		.def("clone", &ExponentialDistribution::Clone, return_value_policy<manage_new_object>())
#endif
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
#if POLPY_NEWWAY
		.def(init<const InverseGammaDistribution &>())
		.def("clone", &InverseGammaDistribution::Clone, return_value_policy<manage_new_object>())
#endif
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
	class_<NormalDistribution, bases<ProbabilityDistribution, AdHocDensity> >("NormalDistBase")
		.def(init<double, double>())
#if POLPY_NEWWAY
		.def(init<const NormalDistribution &>())
		.def("clone", &NormalDistribution::Clone, return_value_policy<manage_new_object>())
#endif
		.def("isDiscrete", &NormalDistribution::IsDiscrete)
		.def("getDistName", &NormalDistribution::GetDistributionName)
		.def("__str__", &NormalDistribution::GetDistributionDescription)
		.def("__repr__", &NormalDistribution::GetDistributionDescription)
		.def("setLot", &NormalDistribution::SetLot)
		.def("setSeed", &NormalDistribution::SetSeed)
		.def("resetLot", &NormalDistribution::ResetLot)
		.def("getMean", &NormalDistribution::GetMean)
		.def("getVar", &NormalDistribution::GetVar)
		.def("getStdDev", &NormalDistribution::GetStdDev)
		.def("getCDF", &NormalDistribution::GetCDF)
		.def("sample", &NormalDistribution::Sample)
		.def("getLnPDF", &NormalDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &NormalDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &NormalDistribution::SetMeanAndVariance)
		;
	register_exception_translator<XProbDist>(&translateXProbDist);
}
