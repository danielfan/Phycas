#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

//#include "phycas/force_include.h"
//#include <cassert>
#include <cmath>
#include "slice_sampler.hpp"
using std::ofstream;
using std::ios;
using std::vector;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::sprintf;
	using std::fabs;
	using std::log;
#endif

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using overrelaxed slice sampling. Throws XPyPhy.
*/
ParamAndLnProb SliceSampler::GetNextOverrelaxedSample(const ParamAndLnProb initialPair)
	{
	assert(r != NULL);
	assert(func != NULL);

	++num_overrelaxed_samples;

	// Step 1: choose a y value for this slice uniformly from 0.0 to f(initialX)
	// where f(initialX) is the density at the original value initialX. To avoid underflow,
	// instead choose ln_y by subtracting an exponential deviate from the 
	// value log(f(initialX)). If a new point x is in the slice, the density curve will 
	// lie above y (and, equivalently, the log-density will be above ln_y)
	//
	const double initialLnfX = initialPair.second;
	const double initialX = initialPair.first;
	double exponential_deviate = -log(r->Uniform(FILE_AND_LINE));
	ln_y = initialLnfX - exponential_deviate;

	// Step 2: Find slice interval at height y
	//
	SliceInterval si = FindSliceInterval(initialPair, ln_y, 1.e-6);
	left_edge = si.first;
	right_edge = si.second;

	// Step 3: Find and return new sampled point
	//
	double x = si.first + si.second - initialX;

	double ln_fx = (*func)(x);
	++func_evals;

	ParamAndLnProb p;
	p.first  = x;
	p.second = ln_fx;
	if (ln_fx < ln_y)
		{
		// Outside slice, not a valid sample
		// Only two reasons for this: 1) bimodal distribution and 2) inadequate precision locating slice boundaries
		// XPyPhy exception thrown in both of these cases
		//
		throw XProbDist("overrelaxed sample failed: possibly not unimodal distribution, or poor precision locating slice boundaries");
		}
	else
		{
		// About to return a sample
		//
		if (x < min_x)
			min_x = x;
		if (x > max_x)
			max_x = x;
		sumValues += x;
		sumWidths += (si.second - si.first);
		sumDiffs += std::fabs(x - initialX);

		// keep track of the sample that is closest to the mode
		if (ln_fx > mode.second)
			{
			mode.first = x;
			mode.second = ln_fx;
			}

		return p;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using slice sampling. Takes <tt>pair<double>(x, ln(f(x)))</tt> and
|	returns <tt>pair<double>(newPoint,ln(f(newPoint)))</tt>
*/
ParamAndLnProb SliceSampler::GetNextSample(const ParamAndLnProb initialPair)
	{
	assert(r != NULL);
	assert(func != NULL);

	++num_samples;

	// Start new vector of x values tried
	//@ most_recent vector is useful for debugging, but may slow things down too much
	most_recent.clear();
	//most_recent.push_back(initialPair);
	
	// Step 1: choose a y value for this slice uniformly from 0.0 to f(initialX)
	// where f(initialX) is the density at the original value initialX. To avoid underflow,
	// instead choose ln_y by subtracting an exponential deviate from the 
	// value log(f(initialX)). If a new point x is in the slice, the density curve will 
	// lie above y (and, equivalently, the log-density will be above ln_y)
	//
	const double initialLnfX = initialPair.second;
	const double initialX = initialPair.first;
	double exponential_deviate = -log(r->Uniform(FILE_AND_LINE));
	ln_y = initialLnfX - exponential_deviate;
	//std::cerr << "~~ ln_y = " << ln_y << std::endl;

	// Step 1.5: if using y-conditional adaptation, choose a new value for w
	//
	if (ycond_on)
		{
		w = CalcW(std::exp(ln_y));
		}

	// Step 2: randomly place interval of width w around initialX (see Fig. 3, p. 715)
	//
	double U	= r->Uniform(FILE_AND_LINE);
	left_edge	= initialX - w*U;
	right_edge	= left_edge + w;
	++realized_m;

	// Step 3: choose maximum number of units on left (J) and right (K)
	// by dividing up the maxUnits total units randomly (see Fig. 3, p. 715)
	//
	double V	= r->Uniform(FILE_AND_LINE);
	unsigned J	= (unsigned)(maxUnits*V);
	unsigned K	= (maxUnits - 1) - J;

	//std::cerr << "maxUnits = " << maxUnits << "\n";
	//std::cerr << "V        = " << V << "\n";
	//std::cerr << "J        = " << J << "\n";
	//std::cerr << "K        = " << K << "\n";
	//std::cerr << "-DBL_MAX = " << (-DBL_MAX) << "\n";
	//std::cerr << std::endl;
	//int ch = scanf("%c");

	// Step 4: Grow interval to the left until left edge is not in the slice
	// or J units have been added, whichever comes first
	//
	for (;;)
		{
		double left_edge_ln_y = (*func)(left_edge);
		//std::cerr << "~~ left edge = " << left_edge << ", left_edge_ln_y = " << left_edge_ln_y << std::endl;
		++func_evals;
		if (left_edge_ln_y < ln_y || J == 0)
			break;

		left_edge -= w;
		--J;
		++realized_m;
		}

	// Step 5: Grow interval to the right until right edge is not in the slice 
	// or K units have been added, whichever comes first
	//
	for (;;)
		{
		double right_edge_ln_y = (*func)(right_edge);
		++func_evals;
		//std::cerr << "~~ right edge = " << right_edge << ", right_edge_ln_y = " << right_edge_ln_y << std::endl;
		if (right_edge_ln_y < ln_y || K == 0)
			break;

		right_edge += w;
		--K;
		++realized_m;
		}

#if 0
	if (K == 0)
		{
		std::cerr << "\n*** K = 0 in SliceSampler::GetNextSample ***\n" << std::endl;
		}
	else if (J == 0)
		{
		std::cerr << "\n*** J = 0 in SliceSampler::GetNextSample ***\n" << std::endl;
		}
#endif

	orig_left_edge = left_edge;
	orig_right_edge = right_edge;
		
	// Step 6: Choose new x value uniformly from interval. If chosen x value 
	// is in the slice, return this as the sampled value. If outside the slice, 
	// adjust the interval accordingly and repeat until a sampled value is found.
	//
	ParamAndLnProb p;
	for(;;)
		{
		double x = left_edge + ((right_edge - left_edge) * r->Uniform(FILE_AND_LINE));

		double ln_fx = (*func)(x);
		//std::cerr << "~~ x = " << x << ", ln_fx = " << ln_fx << std::endl;
		++func_evals;

		p.first  = x;
		p.second = ln_fx;
		if (ln_fx < ln_y)
			{
			// Outside slice, not a valid sample
			// Adjust edge so that it is closer to slice boundary and try again
			//
			++failed_samples;
			most_recent.push_back(p);
			if (x > initialX)
				{
				right_edge = x;
				//std::cerr << "~~ failed: right_edge cropped to " << right_edge << std::endl;
				}
			else
				{
				left_edge = x;
				//std::cerr << "~~ failed: left_edge cropped to " << left_edge << std::endl;
				}
			}
		else
			{
			// About to return a sample
			//
			if (x < min_x)
				min_x = x;
			if (x > max_x)
				max_x = x;
			sumValues += x;
			sumWidths += (right_edge - left_edge);
			sumDiffs += std::fabs(x - initialX);

			// keep track of the sample that is closest to the mode
			if (ln_fx > mode.second)
				{
				mode.first = x;
				mode.second = ln_fx;
				}

			//std::cerr << "*** success: ln_y = " << ln_y << ", left = " << left_edge << ", right = " << right_edge << std::endl;
			return p;
			}
		}
	}
