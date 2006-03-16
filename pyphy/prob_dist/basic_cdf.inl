#if ! defined(BASIC_CDF_INL)
#define BASIC_CDF_INL

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the gamma function for the value `a'.
*/
inline double CDF::LnGamma(double a) const
	{
	double alpha = a;
	return gamln(&alpha);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative beta distribution function for value `x' and parameters `alpha' and `beta'. This is the integral
|	of the beta probability density function from 0 to x.
*/
inline double CDF::CumBeta(
  double x,		/**< is the integral of the probability density function from 0 to x will be returned */
  double alpha,	/**< is the first parameter of the Beta distribution */
  double beta)	const /**< is the second parameter of the Beta distribution */
	{
	double p	= x;
	double q	= 1.0 - x;
	double a	= alpha;
	double b	= beta;
	double cum	= 0.0;
	double ccum	= 0.0;

	cumbet(&p, &q, &a, &b, &cum, &ccum);
	return cum;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative gamma distribution function for the value `x', shape parameter `alpha' and scale parameter 
|	`beta'. This is the integral of the gamma probability density function from 0 up to `x'. The gamma density as 
|	defined here is
|>
|	       x^{alpha - 1} e^{-x/beta}
|	f(x) = -------------------------
|	       beta^{alpha} Gamma(alpha)
|>
|	Note: the CDF library for which this CDF class is a wrapper defines the scale parameter as the inverse of beta!
*/
inline double CDF::CumGamma(
  double x,					/**< is the upper limit of the interval starting at 0.0 for which the integral will be computed */
  double alpha,				/**< is the shape parameter of the gamma distribution */
  double beta) const		/**< is the scale parameter of the gamma distribution */
	{
	double X		= x;
	int status		= 0;
	int which		= 1;	// compute p given x, shape and scale
	double P		= 0.0;
	double Q		= 0.0;
	double shape	= alpha;
	double scale	= 1.0/beta; // CDF library uses inverse of normal definition of scale parameter
	double bound	= 0.0;
	cdfgam(&which, &P, &Q, &X, &shape, &scale, &status, &bound);
	assert(status == 0);
	return P;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `x' such that the integral of the gamma probability density function from 0 up to `x' is equal to
|	the supplied parameter `p' (with shape parameter `alpha' and scale parameter `beta' provided). Useful for simulating
|	draws from a gamma distribution. The gamma density as defined here is
|>
|	       x^{alpha - 1} e^{-x/beta}
|	f(x) = -------------------------
|	       beta^{alpha} Gamma(alpha)
|>
|	Note: the CDF library for which this CDF class is a wrapper defines the scale parameter as the inverse of beta!
*/
inline double CDF::SampleGamma(
  double p,				/**< is the integral of the gamma probability density function from 0 up to `x' */
  double alpha,			/**< is the shape parameter of the gamma distribution */
  double beta) const	/**< is the scale parameter of the gamma distribution */
	{
	double X = 0.0;
	int status = 0;
	int which = 2;	// compute x given p, shape and scale
	double P = p;
	double Q = 1.0 - P;
	double shape = alpha;
	double scale = 1.0 / beta; // CDF library uses inverse of normal definition of scale parameter
	double bound = 0.0;
	cdfgam(&which, &P, &Q, &X, &shape, &scale, &status, &bound);
	if (status == 10)
		{
		// If alpha is tiny (alpha < 0.001) and beta is huge (beta > 1000.0), then cdfgam may 
		// fail with status 10, in which case we return 0.0
		assert(alpha < 0.001 && beta > 1000.0);
		//@POL not sure about the validity of this
		return 0.0;
		}
	return X;
	}

} // namespace phycas

#endif
