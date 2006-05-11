#ifndef PHYCAS_SLICE_SAMPLER_H
#define PHYCAS_SLICE_SAMPLER_H

#include <cstdlib>
#include <climits>
#include <cfloat>
#if defined(POL_PYPHY)
#	include "pyphy/prob_dist/basic_lot.hpp"
	typedef boost::shared_ptr<phycas::Lot> LotShPtr;
#else
#	include "phycas/rand/lot.hpp"
	typedef boost::shared_ptr<Lot> LotShPtr;
#endif
#include "phycas/rand/probability_distribution.hpp"
#include "pyphy/prob_dist/xprobdist.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/python/call_method.hpp>
typedef std::pair<double, double> ParamAndLnProb;
typedef std::pair<double, double> SliceInterval;

struct SliceStats
	{
	SliceStats() : nsamples(0), value(0.0), width(0.0), diff(0.0), failed(0.0), evals(0.0) {}

	unsigned nsamples;
	double value;
	double width;
	double diff;
	double failed;
	double evals;
	};

// Note: AdHocDensity is in phycas/rand/probablity_distribution.hpp

typedef boost::shared_ptr<AdHocDensity> FuncToSampleShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Implements the univariate slice sampler described in Neal, Radford M. 2003. Slice sampling. Annals of Statistics 
|	31:705-741.
*/
class SliceSampler
	{
	public:
								SliceSampler();
								SliceSampler(LotShPtr rnd, FuncToSampleShPtr f);
		virtual					~SliceSampler();

		double					Sample();
		VecDbl					DebugSample();

		double					OverrelaxedSample();
		VecDbl					DebugOverrelaxedSample();

		void					AttachFunc(FuncToSampleShPtr f);
		void					AttachRandomNumberGenerator(LotShPtr rnd);

		void					SetXValue(double x);
		void					SetMaxUnits(unsigned umax);

		void					SetSliceUnitWidth(double uwidth);
		double					GetSliceUnitWidth() const;

		double					AdaptSimple(double multiplier);
		double					AdaptNeal(double multiplier);
		void					AdaptYConditional(double from_ends, double multiplier);
		SliceInterval			FindSliceInterval(ParamAndLnProb x0, const double ln_y0, double tol, unsigned max_steps = UINT_MAX) const;
		double					CalcW(double y0) const;

		double					GetMode() const;
		double					GetLnDensityAtMode() const;

		double					GetLastSampledXValue();
		double					GetLastSampledYValue();
		double					GetSliceYValue();

		double					GetLnZero() const {return -DBL_MAX;}

		
		// For diagnosing inefficiency
		//
		double					GetMinX();
		double					GetMaxX();
		double					GetOrigLeftEdgeOfSlice();
		double					GetOrigRightEdgeOfSlice();
		double					GetLeftEdgeOfSlice();
		double					GetRightEdgeOfSlice();
		unsigned				GetNumFuncEvals();
		unsigned				GetNumFailedSamples();
		unsigned				GetNumUnitsRequired();
		unsigned				GetNumSamples();

		SliceStats				SummarizeDiagnostics();
		void					ResetDiagnostics();

	protected:
		
		void					Init();
		ParamAndLnProb			GetNextSample(const ParamAndLnProb);
		ParamAndLnProb			GetNextOverrelaxedSample(const ParamAndLnProb);
		SliceInterval			BisectionSqueeze(double left, double lnf_left, double right, double lnf_right, const double ln_y0, double tol, unsigned max_steps) const;

		FuncToSampleShPtr		func;				/**< is a functor representing the probability distribution to be sampled */
		LotShPtr				r;					/**< is the random number generator */
		ParamAndLnProb			lastSampled;		/**< most recent valid sample and its relative density */
		double					w;					/**< unit size for interval I */
		unsigned				maxUnits;			/**< maximum number of units of size w to use for interval I (set to UINT_MAX for unlimited) */

		// These quantities are used for diagnostic purposes in GetNextSample()
		//
		double					orig_left_edge;		/**< x coordinate of left edge of most recent slice (before cropping) */
		double					orig_right_edge;	/**< x coordinate of right edge of most recent slice (before cropping) */
		double					left_edge;			/**< x coordinate of left edge of most recent slice (after cropping) */
		double					right_edge;			/**< x coordinate of right edge of most recent slice (after cropping) */
		double					ln_y;				/**< log of f(x) representing most recent slice */
		
		// For diagnosing inefficiency, all reset with call to ResetDiagnostics()
		//
		double					min_x;				/**< minimum x value tried */
		double					max_x;				/**< maximum x value tried */
		double					sumValues;			/**< sum of last `num_samples' sampled values */
		double					sumWidths;			/**< sum of last `num_samples' cropped slice widths (where slice width = right_edge - left_edge) */
		double					sumDiffs;			/**< sum of last `num_samples' differences between successively sampled values */
		unsigned				func_evals;			/**< counts number of function evaluations */
		unsigned				failed_samples;		/**< counts number of samples that failed because they were not in the slice */
		unsigned				realized_m;			/**< counts number of units, each of size w, that were required to bracket the slice */
		unsigned				num_samples;		/**< counts number of times GetNextSample called */
		std::vector<ParamAndLnProb>	most_recent;	/**< vector of all (x, lnx) pairs evaluated in most recent sampling effort */

		// These are needed for y-conditional adaptation
		bool					ycond_on;			/**< if true, w chosen anew for each sample based on y-coordinate of slice */
		double					ycond_a;
		double					ycond_b;
		double					ycond_multiplier;
		ParamAndLnProb			mode;

		// These are for overrelaxed sampling
		unsigned				num_overrelaxed_samples;	/**< counts number of times GetNextOverrelaxedSample called */
	};

typedef boost::shared_ptr<SliceSampler> SliceSamplerShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the Init() member function.
*/
inline SliceSampler::SliceSampler()
	{
	Init();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the starting x-value (`lastSampled.first') to 0.1, `left_bound' to -DBL_MAX and `right_bound' to DBL_MAX.
|	In addition, sets data members `r' and `probdist' to the values passed in.
*/
inline SliceSampler::SliceSampler(
  LotShPtr rnd,				/**< is the random number generator object to use */
  FuncToSampleShPtr f)		/**< is the probability distribution to sample */
  : func(f), r(rnd)
	{
	Init();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing to do.
*/
inline SliceSampler::~SliceSampler()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by constructors to initialize the object.
*/
inline void SliceSampler::Init()
	{
	w					= 1.0;
	ln_y				= 0.0;
	maxUnits			= UINT_MAX;

	lastSampled.first	= 0.1;
	lastSampled.second	= 0.0;

	ycond_on			= false;
	ycond_a				= 1.0;
	ycond_b				= 0.0;
	mode.first			= 0.1;
	mode.second			= -DBL_MAX;
	ycond_multiplier	= 1.0;

	ResetDiagnostics();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Attaches a ProbabilityDistribution object representing the distribution to be sampled. Calls Init() to reset the
|	sampler to its just-constructed state.
*/
inline void SliceSampler::AttachFunc(FuncToSampleShPtr f)
	{
	func = f;
	Init();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets diagnostic counters (`func_evals', `failed_samples', and `realized_m') all to 0.
*/
inline void SliceSampler::ResetDiagnostics()
	{
	min_x					= DBL_MAX;
	max_x					= -DBL_MAX;
	num_samples				= 0;
	num_overrelaxed_samples	= 0;
	func_evals				= 0;
	failed_samples			= 0;
	realized_m				= 0;
	sumValues				= 0.0;
	sumWidths				= 0.0;
	sumDiffs				= 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes means of several quantities that have accumulated since the last call to ResetDiagnostics.
*/
inline SliceStats SliceSampler::SummarizeDiagnostics()
	{
	SliceStats ss;
	if (num_samples > 0)
		{
		ss.nsamples	= num_samples;
		ss.value	= sumValues/(double)num_samples;
		ss.width	= sumWidths/(double)num_samples;
		ss.diff		= sumDiffs/(double)num_samples;
		ss.evals	= (double)func_evals/(double)num_samples;
		ss.failed	= (double)failed_samples/(double)num_samples;
		}

	return ss;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of `w'.
*/
inline double SliceSampler::GetSliceUnitWidth() const
	{
	return w;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets `w' to a value equal to the mean width of a cropped slice interval times the multiplier provided.
*/
inline double SliceSampler::AdaptSimple(double multiplier)
	{
	assert(num_samples > 0);
	ycond_on = false;
	w = multiplier*sumWidths/(double)num_samples;
	return w;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	`w' is set to the average distance between sampled values times the multiplier supplied. Suggested by Neal in 
|	section 4.4 (p. 721).
*/
inline double SliceSampler::AdaptNeal(double multiplier)
	{
	assert(num_samples > 0);
	ycond_on = false;
	w = multiplier*sumDiffs/(double)num_samples;
	return w;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Finds edges of the slice at level y0 (the natural log of which is `ln_y0') by first stepping out in units of size 
|	`w' until the edge of the slice is bracketed, then using bisection to locate the slice boundary, starting at `x0'
|	and continuing until the distance between successive points is less than `tol'. If `max_steps' is less than 
|	UINT_MAX, performs exactly `max_steps' bisections and ignores `tol'. Returns outer width of slice, which is always
|	greater than or equal to the exact width of the slice by an amount that depends on `tol' or `max_steps'. 
|	Throws XProbDist exception if it turns out that the mode (which is really an estimate of the mode) is not actually 
|	underneath the density curve at the height y0.
*/
inline SliceInterval SliceSampler::FindSliceInterval(ParamAndLnProb x0, const double ln_y0, double tol, unsigned max_steps) const
	{
	// bracket left edge
	//
	double left_edge = x0.first;
	double curr_lnfx = x0.second;
	double prev_lnfx = x0.second;
	while (curr_lnfx > ln_y0)
		{
		left_edge -= w;
		prev_lnfx = curr_lnfx;
		curr_lnfx = (*func)(left_edge);
		}

	// use bisection to locate left edge
	//
	SliceInterval slice_interval = BisectionSqueeze(left_edge, curr_lnfx, left_edge + w, prev_lnfx, ln_y0, tol, max_steps);
	left_edge = slice_interval.first;

	// bracket right edge
	//
	double right_edge = x0.first;
	curr_lnfx = x0.second;
	prev_lnfx = x0.second;
	while (curr_lnfx > ln_y0)
		{
		right_edge += w;
		prev_lnfx = curr_lnfx;
		curr_lnfx = (*func)(right_edge);
		}

	// use bisection to locate right edge
	//
	slice_interval = BisectionSqueeze(right_edge - w, prev_lnfx, right_edge, curr_lnfx, ln_y0, tol, max_steps);
	right_edge = slice_interval.second;

	return SliceInterval(left_edge, right_edge);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recursive function that uses bisection to locate precisely the point at which the function crosses the value y0
|	(the natural log of which is `1n_y0'). Stops when the distance between `left' and `right' is less than `tol', or
|	when `max_steps' is 0 (each recursive call decrements `max_steps' by 1. Assumes supplied values `left' and `right' 
|	bracket the crossover point.
*/
inline SliceInterval SliceSampler::BisectionSqueeze(double left, double lnf_left, double right, double lnf_right, const double ln_y0, double tol, unsigned max_steps) const
	{
	bool left_below = (lnf_left < ln_y0);
	bool right_below = (lnf_right < ln_y0);
	assert((left_below && !right_below) || (right_below && !left_below));
	double middle = (left + right)/2.0;
	double lnf_middle = (*func)(middle);
	bool middle_below = (lnf_middle < ln_y0);
	bool middle_above = !middle_below;
	double half_width = (right - left)/2.0;
	bool leave_now = (half_width < tol);
	bool go_left = false;
	if (left_below && middle_below)
		go_left = false;
	else if (left_below && middle_above)
		go_left = true;
	else if (right_below && middle_below)
		go_left = true;
	else if (right_below && middle_above)
		go_left = false;
	else 
		assert(0);
	if (leave_now)
		{
		SliceInterval slice_interval;
		if (go_left)
			{
			slice_interval.first = left;
			slice_interval.second = middle;
			}
		else
			{
			slice_interval.first = middle;
			slice_interval.second = right;
			}
		return slice_interval;
		}
	else // half_width not yet smaller than tol
		{
		if (go_left)
			{
			return BisectionSqueeze(left, lnf_left, middle, lnf_middle, ln_y0, tol, max_steps - 1);
			}
		else
			{
			return BisectionSqueeze(middle, lnf_middle, right, lnf_right, ln_y0, tol, max_steps - 1);
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `w' to be used for next sample when using y-conditional adaptation. Uses the simple linear formula
|	w = m*(a + b*y), where a and b are calculated, and m is set, by calling AdaptYConditional, and y is the supplied 
|	value `y0'. Note: `y0' should NOT be on the log scale.
*/
inline double SliceSampler::CalcW(double y0) const
	{
	double w0 = ycond_a + ycond_b*y0;
	return w0*ycond_multiplier;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current estimate of mode based on previous sampling. 
*/
inline double SliceSampler::GetMode() const
	{
	return mode.first;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns density at current estimate of mode based on previous sampling. 
*/
inline double SliceSampler::GetLnDensityAtMode() const
	{
	return mode.second;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	`w' is not set, but instead is recalulated for each sample based on the value of y (height of slice). The purpose
|	of this function is to parameterize the function used to compute `w' given y. The width of the slice is computed 
|	using bisection at two points, one a fraction `from_ends' from the bottom and the other a fraction 
|	`from_ends' from the top of the density at its mode. These two widths define a slope that can be used to 
|	approximate the width of the density at any height. The actual width used for `w' is `multiplier' time the width 
|	returned by the function.
*/
inline void SliceSampler::AdaptYConditional(double from_ends, double multiplier)
	{
	ycond_on = true;
	ycond_multiplier = multiplier;

	double ln_from_ends = std::log(from_ends);

	double ln_y0 = mode.second;
	double ln_y1 = ln_from_ends + ln_y0;
	SliceInterval s1 = FindSliceInterval(mode, ln_y1, 0.01);
	double w1 = s1.second - s1.first;

	double ln_from_ends_complement = std::log(1.0 - from_ends);

	double ln_y2 = ln_from_ends_complement + ln_y0;
	SliceInterval s2 = FindSliceInterval(mode, ln_y2, 0.01);
	double w2 = s2.second - s2.first;

	// Calculate intercept a and slope b using the two equations:
	//   w1 = a + b*y1
	//   w2 = a + b*y2
	// 
	//       (w1 - w2)
	//   b = ---------
	//       (y1 - y2)
	//
	//   a = w1 - b*y1
	//
	double y1 = std::exp(ln_y1);
	double y2 = std::exp(ln_y2);
	ycond_b = (w1 - w2)/(y1 - y2);
	ycond_a = w1 - ycond_b*y1;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of times GetNextSample called. Use the value returned here to compute averages of the other
|	diagnostic counts.
*/
inline unsigned SliceSampler::GetNumSamples()
	{
	return num_samples;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the minimum value of x tried since last call to ResetDiagnostics function, or DBL_MAX if no sampling has
|	been attempted.
*/
inline double SliceSampler::GetMinX()
	{
	return min_x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the maximum value of x tried since last call to ResetDiagnostics function, or -DBL_MAX if no sampling has
|	been attempted.
*/
inline double SliceSampler::GetMaxX()
	{
	return max_x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the log of the y coordinate of the most recent slice.
*/
inline double SliceSampler::GetSliceYValue()
	{
	return ln_y;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the x coordinate of the left edge of the most recent slice, before cropping based on failed samples.
*/
inline double SliceSampler::GetOrigLeftEdgeOfSlice()
	{
	return orig_left_edge;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the x coordinate of the right edge of the most recent slice, before cropping based on failed samples.
*/
inline double SliceSampler::GetOrigRightEdgeOfSlice()
	{
	return orig_right_edge;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the x coordinate of the left edge of the most recent slice, after cropping based on failed samples.
*/
inline double SliceSampler::GetLeftEdgeOfSlice()
	{
	return left_edge;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the x coordinate of the right edge of the most recent slice, after cropping based on failed samples.
*/
inline double SliceSampler::GetRightEdgeOfSlice()
	{
	return right_edge;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of function evaluations required. Divide by GetNumSamples to compute average number of function
|	evaluations required for each call to GetNextSample.
*/
inline unsigned SliceSampler::GetNumFuncEvals()
	{
	return func_evals;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of failed samples. Each time a sample is attempted but fails because the value is outside the slice,
|	the value of the counter `failed_samples' is incremented and the value is used to reduce the size of the sampling 
|	interval. Divide by GetNumSamples to compute average number of failed samples per call to GetNextSample.
*/
inline unsigned SliceSampler::GetNumFailedSamples()
	{
	return failed_samples;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of units, each of size w, required to bracket the slice. Divide by GetNumSamples to compute average
|	number of units required per call to GetNextSample.
*/
inline unsigned SliceSampler::GetNumUnitsRequired()
	{
	return realized_m;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the random number generator pointer `r' to the supplied random number generator pointer `rnd'. Calls Init() to
|	reset the sampler to its just-constructed state.
*/
inline void SliceSampler::AttachRandomNumberGenerator(LotShPtr rnd)
	{
	r = rnd;
	Init();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the starting x-value `lastSampled.first' to the supplied value `x'.
*/
inline void SliceSampler::SetXValue(double x)
	{
#if POLPY_NEWWAY
	if (x != lastSampled.first)
		{
		lastSampled.first = x;
		assert(func);
		lastSampled.second = (*func)(x);
		}
#else
	lastSampled.first = x;
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `lastSampled.first'.
*/
inline double SliceSampler::GetLastSampledXValue()
	{
	return lastSampled.first;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `lastSampled.second'.
*/
inline double SliceSampler::GetLastSampledYValue()
	{
	return lastSampled.second;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the unit size for the interval I used for sampling from the slice S.
*/
inline void SliceSampler::SetSliceUnitWidth(double uwidth)
	{
	w = uwidth;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the maximum number of units of size `w' to use for interval I used for sampling from the slice S. If `umax'
|	equals zero, maxUnits is set to UINT_MAX instead.
*/
inline void SliceSampler::SetMaxUnits(unsigned umax)
	{
	maxUnits = umax;
	if (maxUnits == 0)
		maxUnits = UINT_MAX;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using overrelaxed slice sampling. Current point is `lastSampled.first'
*/
inline double SliceSampler::OverrelaxedSample()
	{
	if (r == NULL)
		{
		throw XProbDist("must attach random number generator to slice sampler before attempting to draw samples");
		}
	if (func == NULL)
		{
		throw XProbDist("must specify density function for slice sampler before attempting to draw samples");
		}

	// We can never be guaranteed that the full conditional density at lastSampled.first has not 
	// changed since the last call during an MCMC run.
	lastSampled.second = (*func)(lastSampled.first);
	++func_evals;

	const ParamAndLnProb currentPoint = lastSampled;
	lastSampled = GetNextOverrelaxedSample(currentPoint);

	// Let the new sampled value be the starting point for the next sample
	//
	return lastSampled.first;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using overrelaxed slice sampling. Same as OverrelaxedSample, but returns
|	a vector the elements of which are:
|		0: sampled x
|		1: x-coord of vertical slice
|		2: y-coord of top of vertical slice (y-coord of bottom of vertical slice always 0.0)
|		3: x-coord of left edge of horizontal slice
|		4: x-coord of right edge of horizontal slice
|		5: y-coord of horizontal slice
*/
inline VecDbl SliceSampler::DebugOverrelaxedSample()
	{
	if (r == NULL)
		{
		throw XProbDist("must attach random number generator to slice sampler before attempting to draw samples");
		}
	if (func == NULL)
		{
		throw XProbDist("must attach a probability distribution to slice sampler before attempting to draw samples");
		}

	double v1 = lastSampled.first; // x-coord of vertical slice

	// We can never be guaranteed that the full conditional density at lastSampled.first has not 
	// changed since the last call during an MCMC run.
	lastSampled.second = (*func)(lastSampled.first);
	++func_evals;

	double v2 = lastSampled.second; // y-coord of top of vertical slice

	const ParamAndLnProb currentPoint = lastSampled;
	lastSampled = GetNextOverrelaxedSample(currentPoint);

	double v0 = lastSampled.first;	// sampled x
	double v3 = left_edge;			// x-coord of left edge of horizontal slice
	double v4 = right_edge;			// x-coord of right edge of horizontal slice
	double v5 = ln_y;				// y-coord of horizontal slice

	VecDbl v;
	v.push_back(v0);
	v.push_back(v1);
	v.push_back(v2);
	v.push_back(v3);
	v.push_back(v4);
	v.push_back(v5);

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using slice sampling. Current point is `lastSampled.first'
*/
inline double SliceSampler::Sample()
	{
	if (r == NULL)
		{
		throw XProbDist("must attach random number generator to slice sampler before attempting to draw samples");
		}
#if 0
	if (probdist == NULL)
		{
		throw XProbDist("must attach a probability distribution to slice sampler before attempting to draw samples");
		}
#endif

	// We can never be guaranteed that the full conditional density at lastSampled.first has not 
	// changed since the last call during an MCMC run.
	lastSampled.second = (*func)(lastSampled.first);
	++func_evals;
	//std::cerr << "~~ should not see this" << std::endl;

	const ParamAndLnProb currentPoint = lastSampled;
	lastSampled = GetNextSample(currentPoint);

	// Let the new sampled value be the starting point for the next sample
	//
	return lastSampled.first;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using slice sampling. Same as Sample, but returns a vector the elements
|	of which are:
|       0: sampled x
|       1: x-coord of vertical slice
|       2: y-coord of top of vertical slice (y-coord of bottom of vertical slice always 0.0)
|       3: x-coord of left edge of horizontal slice
|       4: x-coord of right edge of horizontal slice
|       5: y-coord of horizontal slice
|       6: horizontal slice interval width
|		7+: x-coord of failed sampling attempts (y-coord all equal to element 5)
*/
inline VecDbl SliceSampler::DebugSample()
	{
	if (r == NULL)
		{
		throw XProbDist("must attach random number generator to slice sampler before attempting to draw samples");
		}
	if (func == NULL)
		{
		throw XProbDist("must attach a probability distribution to slice sampler before attempting to draw samples");
		}

	// We can never be guaranteed that the full conditional density at lastSampled.first has not 
	// changed since the last call during an MCMC run.
	lastSampled.second = 0.0;
	lastSampled.second = (*func)(lastSampled.first);
	++func_evals;
	//std::cerr << "~~ re-evaluating density at current point in DebugSample" << std::endl;

	double v1 = lastSampled.first; // x-coord of vertical slice
	double v2 = lastSampled.second; // y-coord of top of vertical slice

	// Let the new sampled value be the starting point for the next sample
	//
	const ParamAndLnProb currentPoint = lastSampled;
	lastSampled = GetNextSample(currentPoint);

	double v0 = lastSampled.first;	// sampled x
	double v3 = left_edge;			// x-coord of left edge of horizontal slice
	double v4 = right_edge;			// x-coord of right edge of horizontal slice
	double v5 = ln_y;				// y-coord of horizontal slice
	double v6 = w;					// horizontal slice interval width

	VecDbl v;
	v.push_back(v0);
	v.push_back(v1);
	v.push_back(v2);
	v.push_back(v3);
	v.push_back(v4);
	v.push_back(v5);
	v.push_back(v6);

	// let the last items in the vector be the x-coordinates of the sampling attempts (very last element
	// is the x-coord of the successful sample, which should be the same as v0)
	for (std::vector<ParamAndLnProb>::const_iterator it = most_recent.begin(); it != most_recent.end(); ++it)
		{
		v.push_back(it->first);
		}

	return v;
	}

#endif
