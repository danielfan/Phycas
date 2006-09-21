#if !defined(PROBABILITY_DISTRIBUTION_INL)
#define PROBABILITY_DISTRIBUTION_INL

#include "probability_distribution.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	The () operator calls the member function GetRelativeLnPDF.
*/
inline double ProbabilityDistribution::operator()(double x)
	{
	return GetRelativeLnPDF(x);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the random number generator used with the ProbabilityDistribution::Sample member function. The original 
|	random number generator	(data member `myLot') can be replaced by calling the ProbabilityDistribution::ResetLot 
|	function. Note that this object does not take ownership of the Lot object whose pointer is specified as `other'. 
|	It is assumed that `other' is non-NULL.
*/
inline void ProbabilityDistribution::SetLot(
	Lot * other) /**< is a pointer to the random number generator object to be used subsequently by Sample */
	{
	if (other == NULL)
		throw XProbDist("attempt made to install a non-existent pseudorandom number generator");
	lot = other;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the random number seed for the `myLot' data member. Note that if ProbabilityDistribution::SetLot has been 
|	called, calling ProbabilityDistribution::SetSeed is pointless because you will not be setting the seed for the 
|	correct random number generator!
*/
inline void	ProbabilityDistribution::SetSeed(
  unsigned rnseed)	/**< is the new seed value */
	{
	myLot.SetSeed(rnseed);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Makes the data member `lot' (which is used as the random number generator by the member function 
|	ProbabilityDistribution::Sample) point to the local data member myLot. This function only needs to be called if 
|	ProbabilityDistribution::SetLot has been called previously to replace the random number generator used by 
|	ProbabilityDistribution::Sample.
*/
inline void ProbabilityDistribution::ResetLot()
	{
	lot = &myLot;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This version of the function returns the sum of the relative log probability densities for all of the points in an
|	array. 
*/
inline double ProbabilityDistribution::GetRelativeLnPDFArray(
  double *x,	/* the array of values for which the density function is to be evaluated */
  int arrLen) const	/* the number of elements in the array x */
	{
	double retval = 0.0;
	for (int i = 0; i < arrLen; ++i)
		{
		double tmp = GetRelativeLnPDF(x[i]);
		if (tmp == -DBL_MAX)
			return -DBL_MAX;
		retval += tmp;
		}
	return retval;
	}

//############################################################################################
//###### NORMAL DISTRIBUTION INLINED FUNCTIONS ###############################################
//############################################################################################

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `mean' to 0.0 and `sd' to 1.0. Also initializes `pi_const' and `sqrt2_const'. The following is relevant
|	to the calculation of `pi_const':
|>
|	      /|    
|	     / |     sin(theta) = y/r
|	  r /  |     sin(90 deg)  = sin(pi/2) = 1.0
|	   /   | y   asin(1.0) = pi/2
|	  /    |     2.0*asin(1.0) = pi
|	 /_____|    
|	    x
|>
*/
inline NormalDistribution::NormalDistribution()
  	{
	mean = 0.0;
	sd = 1.0;
	ComputeLnConst();
	pi_const = 2.0*std::asin(1.0);
	sqrt2_const = std::sqrt(2.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the `mean' and `sd' data members to `m' and `s', respectively. Also initializes `pi_const' and 
|	`sqrt2_const'. The following is relevant to the calculation of `pi_const':
|>
|	      /|    
|	     / |     sin(theta) = y/r
|	  r /  |     sin(90 deg)  = sin(pi/2) = 1.0
|	   /   | y   asin(1.0) = pi/2
|	  /    |     2.0*asin(1.0) = pi
|	 /_____|    
|	    x
|>
*/
inline NormalDistribution::NormalDistribution(
  double m,		/* the mean parameter */
  double s)		/* the standard deviation parameter */
  	{
	PHYCAS_ASSERT(s > 0.0);
	mean = m;
	sd = s;
	pi_const = 2.0*std::asin(1.0);
	sqrt2_const = std::sqrt(2.0);
	ComputeLnConst();
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
inline NormalDistribution::~NormalDistribution()
	{
	//std::cerr << "Deleting a NormalDistribution object" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false because the univariate normal is a continuous distribution.
*/
inline bool NormalDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Normal", which is the name of this distribution.
*/
inline std::string NormalDistribution::GetDistributionName() const
	{
	return "Normal";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Normal(<mean>, <stddev>)", with <mean> and <stddev> replaced with the current values of the 
|	`mean' and `sd' data members.
*/
inline std::string NormalDistribution::GetDistributionDescription() const
	{
	return str(boost::format("NormalDist(%#.5f, %#.5f)") % mean % sd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected mean of the normal distribution as currently specified, which is simply the value of the `mean'
|	data member.
*/
inline double NormalDistribution::GetMean() const
	{
	return mean;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected variance of the normal distribution as currently specified, which is simply the square of the 
|	value of the data member `sd'.
*/
inline double NormalDistribution::GetVar() const
	{
	return std::pow(sd, 2.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected standard deviation of the normal distribution as currently specified, which is simply the value
|	of the `sd' data member.
*/
inline double NormalDistribution::GetStdDev() const
	{
	return sd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative normal distribution for x (integral of normal density function from negative infinity to x). 
*/
inline double NormalDistribution::GetCDF(
  double x)	 const	/**< is the value for which the cumulative distribution function is to be evaluated */
	{
	return cdf.CumNorm(x, mean, sd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a sampled value from a normal distribution having mean `mean' and standard deviation `sd'.
*/
inline double NormalDistribution::Sample() const
	{
	return cdf.SampleNorm(lot->Uniform(FILE_AND_LINE), mean, sd); 
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The probability density function of the normal distribution is 
|>
|	              1            /  (x - mean)^2 \
|	f(x) = --------------- exp| --------------- |
|	       sd (2 pi)^{0.5}     \     2 sd^2    /
|>
|	This function returns the natural log of the density function at `x':
|>
|	ln[f(x)] = (x - mu)^2/(2 sd^2) - ln[sd] - ln[2*pi]/2.0
|<
|	The sum of the last two terms are precalculated and available as the variable `ln_const', but this function only
|	returns the first term because all we need is something proportional to the density at `x'.
*/
inline double NormalDistribution::GetRelativeLnPDF(
  double x)   const	/* the value for which the density function is to be evaluated */
	{
	double term1 = x - mean;
	double term2 = sqrt2_const*sd;
	double term3 = term1/term2;
	double lnpdf = pow(term3, 2.0);
	return lnpdf;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns NormalDistribution::GetRelativeLnPDF plus `ln_const' (log of the part of the density function that depends
|	only on the standard deviation parameter, and which is recalculated when `sd' changes).
*/
inline double NormalDistribution::GetLnPDF(
  double x)   const /* the value for which the density function is to be evaluated */
	{
	return (GetRelativeLnPDF(x) + ln_const);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows distribution to be specified in terms of the mean and variance rather than the shape and standard deviation.
|	This function sets `mean' to the supplied `mu' and `sd' to the square root of the supplied `var'. Assumes `var' is 
|	greater than zero.
*/
inline void NormalDistribution::SetMeanAndVariance(
  double mu, 	/* the mean of the gamma distribution */
  double var)	/* the variance of the gamma distribution */
  	{
	PHYCAS_ASSERT(var > 0.0);
	mean = mu;
  	sd = sqrt(var);
	ComputeLnConst();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the value of the data member `ln_const' (the log of that part of the normal density function that depends
|	only on the standard deviation parameter, and which is recalculated whenever `sd' changes).
*/
inline void NormalDistribution::ComputeLnConst()
  	{
	ln_const = -log(sd) - log(2.0*pi_const)/2.0;
	}

//############################################################################################
//###### BERNOULLI DISTRIBUTION INLINED FUNCTIONS ############################################
//############################################################################################

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes p to 0.5.
*/
inline BernoulliDistribution::BernoulliDistribution()
  	{
	p = 0.5;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
inline BernoulliDistribution::~BernoulliDistribution()
	{
	//std::cerr << "Deleting a BernoulliDistribution object" << std::endl;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes parameter p to prob_success. Assumes p is greater than or equa to zero and less than or equal to 1.0.
*/
inline BernoulliDistribution::BernoulliDistribution(
  double prob_success)	/* probability of success on any given trial */
  	{
	if (prob_success < 0.0 || prob_success > 1.0)
		throw XProbDist("success probability supplied to BernoulliDist constructor out of bounds (should be between 0.0 and 1.0)");
	p = prob_success;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns true.
*/
inline bool BernoulliDistribution::IsDiscrete() const
	{
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Bernoulli", which is the name of this distribution.
*/
inline std::string BernoulliDistribution::GetDistributionName() const
	{
	return "Bernoulli";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Bernoulli(<psuccess>)", with <psuccess> replaced with actual probability of success.
*/
inline std::string BernoulliDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("BernoulliDist(%.5f)", p);
	return str(boost::format("BernoulliDist(%#.5f)") % p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected mean, which is simply the p parameter.
*/
inline double BernoulliDistribution::GetMean() const
  	{
	return p;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected variance, which is p*(1-p).
*/
inline double BernoulliDistribution::GetVar() const
  	{
	return p*(1.0 - p);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected standard deviation, which is sqrt(p*(1.0 - p)).
*/
inline double BernoulliDistribution::GetStdDev() const
  	{
	return std::sqrt(p*(1.0 - p));
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the cumulative distribution function evaluated at the specified value x. Assumes x is either 0.0 or 1.0.
*/
inline double BernoulliDistribution::GetCDF(
  double x) const	/* the value for which the cumulative distribution function is to be evaluated; specify either 0.0 or 1.0 */
  	{
	if (x != 0.0 && x != 1.0)
		throw XProbDist("only 0 or 1 are acceptable values for the BernoulliDist getCDF function");
	return (x == 0.0 ? 1.0 - p : 1.0);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns a sampled value from a bernoulli distribution having probabilty of success parameter p.
*/
inline double BernoulliDistribution::Sample() const
  	{
	return (lot->Uniform(FILE_AND_LINE) <= p ? 1.0 : 0.0);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	This function returns the natural log of the probability of x, which is p if x equals 1 and 1-p if x equals 0. If p
|	is exactly 1.0 and x is zero, or if p is exactly 0.0 and x is 1, the probability is zero and the value returned 
|	-DBL_MAX (which is as close as possible to negative infinity).
*/
inline double BernoulliDistribution::GetLnPDF(
  double x) const	/**< is the value for which the density function is to be evaluated; should be either 0.0 or 1.0 */
	{
	if (x != 0.0 && x != 1.0)
		throw XProbDist("only 0 or 1 are acceptable values for the BernoulliDist getLnPDF function");
	double lnp = 0.0;
	if (x == 0.0 && p < 1.0)
		lnp = std::log(1.0 - p);
	else if (x == 1.0 && p > 0.0)
		lnp = std::log(p);
	else
		lnp = -DBL_MAX;
	return lnp;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls GetLnPDF to do the work.
*/
inline double BernoulliDistribution::GetRelativeLnPDF(
  double x) const	/**< is the value for which the density function is to be evaluated; should be either 0.0 or 1.0 */
	{
	return GetLnPDF(x);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows distribution to be specified in terms of the mean rather than the probability of success on any given trial.
|	Because there is only one parameter, only the mean is used. The parameter p is set to the specified mean. The mean 
|	is assumed to lie in [0.0, 1.0]. The specified variance is ignored.
*/
inline void BernoulliDistribution::SetMeanAndVariance(
  double mean,	/* the mean of the bernoulli distribution */
  double var)	/* ignored */
  	{
 #if defined (HAVE_PRAGMA_UNUSED)
 #	pragma unused (var)
 #endif
	if (mean < 0.0 || mean > 1.0)
		throw XProbDist("the mean of BernoulliDist must be between 0.0 and 1.0");
  	p = mean;
	}

//############################################################################################
//###### BINOMIAL DISTRIBUTION INLINED FUNCTIONS #############################################
//############################################################################################

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes p to 0.5.
*/
inline BinomialDistribution::BinomialDistribution()
  	{
	n = 1.0;
	p = 0.5;
	if (p > 0.0)
		lnp = std::log(p);
	else
		lnp = -DBL_MAX;
	q = 1.0 - p;
	if (q > 0.0)
		lnq = std::log(q);
	else
		lnq = -DBL_MAX;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes n to sample_size and parameter p to prob_success. Assumes p is greater than or equal to zero and less 
|	than or equal to 1.0, and that n is greater than or equal to zero.
*/
inline BinomialDistribution::BinomialDistribution(
  double sample_size,	/* sample size (number of trials) */
  double prob_success)	/* probability of success on any given trial */
  	{
	if (sample_size < 0.0)
		throw XProbDist("the supplied sample size must be positive when creating BinomialDist objects");
	if (prob_success < 0.0 || prob_success > 1.0)
		throw XProbDist("the supplied probability of success must be between 0.0 and 1.0 when creating BinomialDist objects");
	n = sample_size;
	p = prob_success;
	if (p > 0.0)
		lnp = std::log(p);
	else
		lnp = -DBL_MAX;
	q = 1.0 - p;
	if (q > 0.0)
		lnq = std::log(q);
	else
		lnq = -DBL_MAX;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
inline BinomialDistribution::~BinomialDistribution()
	{
	//std::cerr << "Deleting a BinomialDistribution object" << std::endl;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns true.
*/
inline bool BinomialDistribution::IsDiscrete() const
	{
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Binomial", which is the name of this distribution.
*/
inline std::string BinomialDistribution::GetDistributionName() const
	{
	return "Binomial";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Binomial(n, p)", where n is the number of trials and p the probability of success on each trial.
*/
inline std::string BinomialDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("BinomialDist(%d, %.5f)", (unsigned)n, p);
	unsigned nn = (unsigned)n;
	return str(boost::format("BinomialDist(%d, %#.5f)") % nn % p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected mean, which is n*p.
*/
inline double BinomialDistribution::GetMean() const
  	{
	return n*p;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected variance, which is n*p*(1-p).
*/
inline double BinomialDistribution::GetVar() const
  	{
	return n*p*q;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected standard deviation, which is sqrt(n*p*(1.0 - p)).
*/
inline double BinomialDistribution::GetStdDev() const
  	{
	return std::sqrt(n*p*q);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the cumulative distribution function evaluated at the specified value x. Assumes x is either 0.0 or a 
|	positive integer less than or equal to n.
*/
inline double BinomialDistribution::GetCDF(
  double x) const	/* the value for which the cumulative distribution function is to be evaluated; specify either 0.0 or a positive integer */
  	{
	if (x < 0.0)
		throw XProbDist("the value supplied to the BinomialDist getCDF function must be greater than or equal to 0.0");
	if (x > n)
		throw XProbDist("the value supplied to the BinomialDist getCDF function must be less than or equal to the sample size");
	if (x - (unsigned)x != 0.0)
		throw XProbDist("fractional values are not allowed as arguments to the BinomialDist getCDF function");

	double cum_prob = 0.0;
	double n_only = n + 1.0;

	double term1 = cdf.LnGamma(n_only);

	for (unsigned i = 0; i <= x; ++i)
		{
		double x = (double)i;
		double term2 = x*lnp + (n - x)*lnq;
		double x_only = x + 1.0;
		double term3 = cdf.LnGamma(x_only);
		double n_minus_x = n - x + 1.0;
		double term4 = cdf.LnGamma(n_minus_x);
		double lnprob = term1 + term2 - term3 - term4;
		cum_prob += std::exp(lnprob);
		}

	return cum_prob;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns a sampled value from a binomial distribution having probabilty of success parameter p. 
*/
inline double BinomialDistribution::Sample() const
  	{
	//@POL if p is small, a faster method may be to add together x geometric random numbers
	// until the sum is greater than n-x. The number of such geometric random numbers is
	// a sample from a binomial distribution. From Evans, Hastings, Peacock Statistical Distributions
	// book. Generate a geometric random number using [ln(u)/ln(1-p)] - 1, rounded up to the next
	// larger integer.
	//
	unsigned cum = 0;
	for (unsigned i = 0; i < n; ++i)
		{
		if (lot->Uniform(FILE_AND_LINE) <= p)
			++cum;
		}
	return cum;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	This function returns the natural log of the probability of x, which is n! p^x (1-p)^(n-p)/(x! (n-x)!). If p is 
|	exactly 1.0 and x is zero, or if p is exactly 0.0 and x is 1, the probability is zero and the value returned 
|	-DBL_MAX (which is as close as possible to negative infinity).
*/
inline double BinomialDistribution::GetLnPDF(
  double x) const	/**< is the value for which the density function is to be evaluated; should be either 0.0 or a positive integer */
	{
	// The GetRelativeLnPDF function checks for validity of x
	//
	double lnprob = GetRelativeLnPDF(x);

	// The portion below is only necessary if normalized probabilities are desired
	// because these terms involve only the data x and the sample size n, both of 
	// which are constant during likelihood or bayesian analyses
	//

	double n_only		= n + 1.0;
	double x_only		= x + 1.0;
	double n_minus_x	= n - x + 1.0;

	lnprob += cdf.LnGamma(n_only);
	lnprob -= cdf.LnGamma(x_only);
	lnprob -= cdf.LnGamma(n_minus_x);

	return lnprob;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function returns the natural log of the probability of x, which is n! p^x (1-p)^(n-p)/(x! (n-x)!). If p is 
|	exactly 1.0 and x is zero, or if p is exactly 0.0 and x is 1, the probability is zero and the value returned 
|	-DBL_MAX (which is as close as possible to negative infinity).
*/
inline double BinomialDistribution::GetRelativeLnPDF(
  double x) const	/* the value for which the density function is to be evaluated; should be either 0.0 or a positive integer */
	{
	//@POL these XProbDist throws should be eliminated in favor of returning -DBL_MAX
	if (x < 0.0)
		throw XProbDist("the value supplied to the BinomialDist getRelativeLnPDF function must be greater than or equal to 0.0");
	if (x > n)
		throw XProbDist("the value supplied to the BinomialDist getRelativeLnPDF function must be less than or equal to the sample size");
	if (x - (unsigned)x != 0.0)
		throw XProbDist("fractional values are not allowed as arguments to the BinomialDist getRelativeLnPDF function");

	double lnprob = -DBL_MAX;
	if ((p == 0.0 && x > 0.0) || (p == 1.0 && x == 0.0))
		return lnprob;

	// Note: do not change this function without also verifying GetLnPDF (which calls this function)
	//
	lnprob = (x*lnp);
	lnprob += ((n-x)*lnq);

	return lnprob;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows distribution to be specified in terms of the mean rather than the probability of success on any given trial.
|	The argument mean is assumed to be greater than 0.0, and var is ignored because it is determined entirely by the 
|	sample size and mean.
*/
inline void BinomialDistribution::SetMeanAndVariance(
  double mean,	/**< the mean of the binomial distribution */
  double var)	/**< ignored for the Binomial distribution */
  	{
	if (mean <= 0.0 || mean > n)
		throw XProbDist("mean must be greater than 0.0 and less than or equal to sample size for BinomialDist");

	p = mean/n;
	}

//############################################################################################
//###### UNIFORM DISTRIBUTION INLINED FUNCTIONS ##############################################
//############################################################################################

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the left bound to 0.0 and right bound to 1.0.
*/
inline UniformDistribution::UniformDistribution()
  	{
	a = 0.0;
	b = 1.0;
	log_density = 0.0;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the left bound (a) and right bound (b) parameters to the specified values. Assumes right_bound is 
|	greater than left_bound.
*/
inline UniformDistribution::UniformDistribution(
  double left_bound,	/* left bound */
  double right_bound)	/* right bound */
  	{
	if (right_bound <= left_bound)
		throw XProbDist("right bound must exceed left bound for UniformDist");
	a = left_bound;
	b = right_bound;
	log_density = -1.0*std::log(b-a);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
inline UniformDistribution::~UniformDistribution()
	{
	//std::cerr << "Deleting a UniformDistribution object" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false.
*/
inline bool UniformDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Uniform", which is the name of this distribution.
*/
inline std::string UniformDistribution::GetDistributionName() const
	{
	return "Uniform";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Uniform(<lower>, <upper>)", with <lower> and <upper> replaced with actual lower and upper bounds.
*/
inline std::string UniformDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("UniformDist(%.5f, %.5f)", a, b);
	return str(boost::format("UniformDist(%#.5f, %#.5f)") % a % b);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected mean, which is simply the arithmetic average of the a and b parameters.
*/
inline double UniformDistribution::GetMean() const
  	{
	return ((a + b)/2.0);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected variance, which is (b - a)^2/12.
*/
inline double UniformDistribution::GetVar() const
  	{
	return ((b - a)*(b - a)/12.0);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected standard deviation, which is fabs(b - a)/sqrt(12).
*/
inline double UniformDistribution::GetStdDev() const
  	{
	return (std::fabs(b - a)/std::sqrt(12.0));
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the cumulative distribution function evaluated at the specified value x. If GetCDF(x) returns 0.4, this 
|	means that 40% of sampled values of X, where X is a uniform random variable, will be less than or equal to the 
|	value x. Assumes x is in the interval [a, b].
*/
inline double UniformDistribution::GetCDF(
  double x)	 const/* the value for which the cumulative distribution function is to be evaluated */
  	{
	if (x < a)
		return 0.0;
	else if (x > b)
		return 1.0;
	else
		return ((x - a)/(b - a));
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns a sampled value from a uniform distribution having left and right bounds a and b, respectively. The value 
|	returned is a + (b - a)*r.Uniform(FILE_AND_LINE).
*/
inline double UniformDistribution::Sample() const
  	{
	double u = lot->Uniform(FILE_AND_LINE);
	return (a + (b - a)*u);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	The probability density function of the uniform distribution is 
|>
|		f(y) = 1 / (b - a)
|<
|	This function returns the natural log of the normalized density function at x, i.e. -log(b-a)
|	If x is outside the interval [a, b], the value -DBL_MAX is returned.
*/
inline double UniformDistribution::GetLnPDF(
  double x) const	/* the value for which the density function is to be evaluated */
	{
#	if defined (NDEBUG) && defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(x)
#	endif 

	if (x < a || x > b)
		return -DBL_MAX;
	else
		return log_density;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The probability density function of the uniform distribution is 
|>
|		f(y) = 1 / (b - a)
|<
|	This function returns the natural log of the non-normalized density function at x, which is just 0.0 because the 
|	denominator is the normalizing factor and is ignored. If x is outside the interval [a, b], the value -DBL_MAX is 
|	returned.
*/
inline double UniformDistribution::GetRelativeLnPDF(
  double x) const	/* the value for which the density function is to be evaluated */
	{
#	if defined (NDEBUG) && defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(x)
#	endif 

	if (x < a || x > b)
		return -DBL_MAX;
	else
		return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows distribution to be specified in terms of the mean and variance rather than the left and right bounds. 
|	Letting c equal the quantity sqrt(12*var), this function sets the left bound parameter a = mean - c/2 and the right 
|	bound parameter b = a + c. Assumes var > 0.0.
*/
inline void UniformDistribution::SetMeanAndVariance(
  double mean,	/* the mean of the uniform distribution */
  double var)	/* the variance of the uniform distribution */
  	{
	if (var <= 0.0)
		throw XProbDist("specified variance must be greater than 0.0 for UniformDist");
	double c = std::sqrt(12.0*var);
	a = mean - (c/2.0);
  	b = a + c;
	log_density = -1.0*std::log(c);
	}

//############################################################################################
//###### GAMMA DISTRIBUTION INLINED FUNCTIONS ################################################
//############################################################################################

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the shape and scale parameters to 1.0 and 1.0, respectively.
*/
inline GammaDistribution::GammaDistribution()
  	{
	alpha = 1.0;
	beta = 1.0;
	ln_const = 0.0;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
inline GammaDistribution::~GammaDistribution()
	{
	//std::cerr << "Deleting a GammaDistribution object" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the shape and scale parameters to the specified values.
*/
inline GammaDistribution::GammaDistribution(
  double shape,		/* the shape parameter */
  double scale)		/* the scale parameter */
  	{
	PHYCAS_ASSERT(shape > 0.0);
	PHYCAS_ASSERT(scale > 0.0);
	alpha = shape;
	beta = scale;
	ComputeLnConst();
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns false.
*/
inline bool GammaDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Gamma", which is the name of this distribution.
*/
inline std::string GammaDistribution::GetDistributionName() const
	{
	return "Gamma";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Gamma(<shape>, <scale>)", with <shape> and <scale> replaced with actual shape and scale 
|	parameters.
*/
inline std::string GammaDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("GammaDist(%.5f, %.5f)", alpha, beta);
	return str(boost::format("GammaDist(%#.5f, %#.5f)") % alpha % beta);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected mean of the gamma distribution as currently specified, which is the product of the shape and 
|	scale parameters.
*/
inline double GammaDistribution::GetMean() const
	{
	return alpha*beta;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected variance of the gamma distribution as currently specified, which is the product of the shape
|	parameters and the square of the scale parameter.
*/
inline double GammaDistribution::GetVar() const
	{
	return alpha*beta*beta;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected standard deviation of the gamma distribution as currently specified. Returns beta times the 
|	square root of alpha.
*/
inline double GammaDistribution::GetStdDev() const
	{
	return std::fabs(beta)*std::sqrt(alpha);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative gamma distribution for x (integral of gamma density function from 0.0 to x). Assumes x is 
|	greater than or equal to zero.
*/
inline double GammaDistribution::GetCDF(
  double x)	 const	/**< is the value for which the cumulative distribution function is to be evaluated */
	{
	if (x <= 0.0)
		return 0.0;
	else
		return cdf.CumGamma(x, alpha, beta);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a sampled value from a gamma distribution having shape alpha and scale beta.
*/
inline double GammaDistribution::Sample() const
	{
	return cdf.SampleGamma(lot->Uniform(FILE_AND_LINE), alpha, beta); 
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The probability density function of the gamma is 
|>
|	       x^(alpha-1) exp(-x/beta)
|	f(x) = ------------------------
|	       beta^alpha Gamma(alpha)
|<
|	This function returns the natural log of the density function at `x':
|>
|	ln[f(x)] = (alpha - 1)*ln(x) - x/beta - alpha*ln(beta) - lnGamma(alpha)
|<
|	The sum of the last two terms are precalculated and available as the variable `ln_const'.
|	Returns -DBL_MAX if PDF is zero, which happens if x < 0.0, or if (x = 0.0 and alpha > 1.0).
*/
inline double GammaDistribution::GetRelativeLnPDF(
  double x)   const/* the value for which the density function is to be evaluated */
	{
	if (x < 0.0 || (x == 0.0 && alpha > 1.0))
		return -DBL_MAX;
	else
		{
		double term1 = (alpha - 1.0)*std::log(x);
		double term2 = -x/beta;
		double lnpdf = term1 + term2;

		return lnpdf;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns GammaDistribution::GetRelativeLnPDF plus `ln_const' (log of the part of the density function that depends
|	only on the shape and scale parameters, and which is recalculated when either shape or scale is changed).
*/
inline double GammaDistribution::GetLnPDF(
  double x)   const /* the value for which the density function is to be evaluated */
	{
	return (GetRelativeLnPDF(x) + ln_const);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows distribution to be specified in terms of the mean and variance rather than the shape and scale. This function 
|	sets the shape parameter to the mean squared divided by the variance, and sets the scale parameter to the variance 
|	divided by the mean. Assumes var is greater than zero.
*/
inline void GammaDistribution::SetMeanAndVariance(
  double mean, 	/* the mean of the gamma distribution */
  double var)	/* the variance of the gamma distribution */
  	{
	PHYCAS_ASSERT(mean > 0.0);
	PHYCAS_ASSERT(var > 0.0);
	alpha = mean*mean/var;
  	beta = var/mean;
	ComputeLnConst();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the value of the data member `ln_const' (the log of that part of the gamma density function that depends
|	only on the shape and scale parameter, and which is recalculated whenever either the shape or scale is changed).
*/
inline void GammaDistribution::ComputeLnConst()
  	{
	double a = alpha;
	double gammalnalpha = cdf.LnGamma(a);
	ln_const = -1.0*(alpha*std::log(beta) + gammalnalpha);
	}

//############################################################################################
//###### EXPONENTIAL DISTRIBUTION INLINED FUNCTIONS ##########################################
//############################################################################################

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the shape and scale parameters of the underlying gamma distribution to 1.0 and 1.0, respectively.
*/
inline ExponentialDistribution::ExponentialDistribution()
  : GammaDistribution()
  	{
  	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the shape and scale parameters of the underlying gamma distribution to 1.0 and 1.0/lambda, respectively.
*/
inline ExponentialDistribution::ExponentialDistribution(
  double lambda)	/* the single hazard rate parameter of the exponential distribution */ 
  : GammaDistribution(1.0, 1.0/lambda)
  	{
	if (lambda <= 0.0)
		throw XProbDist("specified hazard parameter must be greater than 0.0 for ExponentialDist");
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
inline ExponentialDistribution::~ExponentialDistribution()
	{
	//std::cerr << "ExponentialDistribution dying..." << std::endl;	
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false.
*/
inline bool ExponentialDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Exponential", which is the name of this distribution.
*/
inline std::string ExponentialDistribution::GetDistributionName() const
	{
	return "Exponential";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Exponential(<hazard>)", with <hazard> replaced with actual hazard rate parameter.
*/
inline std::string ExponentialDistribution::GetDistributionDescription() const
	{
	double hazard_rate = 1.0/beta;
	//return MakeStrPrintF("ExponentialDist(%.5f)", hazard_rate);
	return str(boost::format("ExponentialDist(%#.5f)") % hazard_rate);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows distribution to be specified in terms of the mean and variance rather than the hazard rate parameter. This 
|	function sets shape = 1.0 and scale = mean. Throws XProbDist if mean is not greater than 0.0. This
|	distribution is entirely determined by the mean, so the variance argument is ignored.
*/
inline void ExponentialDistribution::SetMeanAndVariance(
  double mean,	/* the mean of the exponential distribution */
  double var)	/* ignored */
  	{
	if (mean <= 0.0)
		throw XProbDist("specified mean must be greater than 0.0 for ExponentialDist");
#	if defined (NDEBUG) && defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(var)
#	endif 
	alpha = 1.0;
	beta = mean;
	}

//############################################################################################
//###### INVERSE GAMMA DISTRIBUTION INLINED FUNCTIONS ########################################
//############################################################################################

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the shape and scale parameters to 3.0 and 0.5, respectively, which yields a distribution with mean and 
|	variance 1.0.
*/
inline InverseGammaDistribution::InverseGammaDistribution()
  	{
	alpha = 3.0;
	beta = 0.5;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the shape and scale parameters to the specified values. Assumes shape is greater than 2.0 (variance is 
|	not defined unless alpha is greater than 2.0)
*/
inline InverseGammaDistribution::InverseGammaDistribution(
  double shape,		/* the shape parameter */
  double scale)		/* the scale parameter */
  	{
	if (shape <= 2.0)
		throw XProbDist("variance undefined for shape less than or equal to 2");
	alpha = shape;
	beta = scale;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
inline InverseGammaDistribution::~InverseGammaDistribution()
	{
	//std::cerr << "InverseGammaDistribution dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false.
*/
inline bool InverseGammaDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "InverseGamma", which is the name of this distribution.
*/
inline std::string InverseGammaDistribution::GetDistributionName() const
	{
	return "InverseGamma";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "InverseGamma(<shape>, <scale>)", with <shape> and <scale> replaced with actual shape and scale 
|	parameters.
*/
inline std::string InverseGammaDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("InverseGammaDist(%.5f, %.5f)", alpha, beta);
	return str(boost::format("InverseGammaDist(%#.5f, %#.5f)") % alpha % beta);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected mean of the inverse gamma distribution as currently specified, which is 1/(beta * (alpha - 1)).
*/
inline double InverseGammaDistribution::GetMean() const
	{
	double denom = beta * (alpha - 1.0);
	return (1.0/denom);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected variance of the inverse gamma distribution as currently specified, which is 
|	1/(beta^2 * (alpha - 1)^2 * (alpha - 2)).
*/
inline double InverseGammaDistribution::GetVar() const
	{
	double term1 = beta*beta;
	double term2 = (alpha - 1.0)*(alpha - 1.0);
	double denom = term1*term2*(alpha - 2.0);
	return (1.0/denom);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected standard deviation of the gamma distribution as currently specified. Returns beta times the 
|	square root of alpha.
*/
inline double InverseGammaDistribution::GetStdDev() const
	{
	double denom = beta*(alpha - 1.0)*std::sqrt(alpha - 2.0);
	return 1.0/denom;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative inverse gamma distribution for x (integral of inverse gamma density function from 0.0 to x). 
|	Assumes x is greater than or equal to zero.
*/
inline double InverseGammaDistribution::GetCDF(
  double x)	 const/* the value for which the cumulative distribution function is to be evaluated */
	{
	if (x <= 0.0)
		return 0.0;
	else
		{
		double y = 1.0/x;

		double Fy = cdf.CumGamma(y, alpha, beta);

		return 1.0 - Fy;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a sampled value from an inverse gamma distribution having shape parameter alpha and scale parameter beta. 
|	If the next uniform random deviate generated using the supplied random number generator just happens to be such that
|	the gamma deviate is exactly zero, the value DBL_MAX is returned (closest we can get to positive infinity.
*/
inline double InverseGammaDistribution::Sample() const
	{
	double gamma_deviate = cdf.SampleGamma(lot->Uniform(FILE_AND_LINE), alpha, beta);

	if (gamma_deviate == 0.0)
		return DBL_MAX;
	else
		return 1.0/gamma_deviate;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The probability density function of the inverse gamma distribution is 
|>
|	                    exp(-1/(x*beta))
|		f(x) = -----------------------------------
|	           x^(alpha+1) beta^alpha Gamma(alpha)
|<
|	This function returns the natural log of the density function at x. If x is less than or equal to 0.0, returns 
|	-DBL_MAX.
*/
inline double InverseGammaDistribution::GetLnPDF(
  double x)   const/* the value for which the density function is to be evaluated */
	{
	if (x <= 0.0)
		return -DBL_MAX;
	else
		{
		double lnpdf = -1.0/(x*beta);
		lnpdf -= ((alpha + 1.0)*std::log(x));
		lnpdf -= alpha*std::log(beta);
		double a = alpha;
		lnpdf -= cdf.LnGamma(a);
		return lnpdf;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The probability density function of the inverse gamma distribution is 
|>
|	                    exp(-1/(x*beta))
|		f(x) = -----------------------------------
|	           x^(alpha+1) beta^alpha Gamma(alpha)
|<
|	This function returns the natural log of the density function at x. Since this will be used in the context of MCMC,
|	we can omit terms not containing x because they would cancel anyway when computing the acceptance ratio, making the
|	log of the relative density at x simply
|>
|		ln[f(x)] = -1/(x*beta) - (alpha + 1)*ln(x)
|<
|	If x is less than or equal to 0.0, returns -DBL_MAX.
*/
inline double InverseGammaDistribution::GetRelativeLnPDF(
  double x)   const/* the value for which the density function is to be evaluated */
	{
	if (x <= 0.0)
		return -DBL_MAX;
	else
		{
		double term1 = (alpha + 1.0) * std::log(x);
		double term2 = 1.0/(x*beta);
		return -(term1 + term2);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows inverse gamma distribution to be specified in terms of the mean and variance rather than the shape and scale.
|	The mean of an inverse gamma distribution equals 1/(beta*(alpha - 1)) and the variance equals mean^2/(alpha - 2).
|	Inverting these, alpha = 2 + mean^2/var and beta = var/(mean*(mean^2 + var)). Throws XProbDist if
|	var is greater than zero. Letting epsilon = mean^2/var, alpha = 2 + epsilon and beta = 1.0/(mean + mean*epsilon).
*/
inline void InverseGammaDistribution::SetMeanAndVariance(
  double mean,	/* the mean of the inverse gamma distribution */
  double var)	/* the variance of the inverse gamma distribution */
  	{
	if (var <= 0.0)
		throw XProbDist("specified variance must be greater than 0.0 for InverseGammaDist");
	double epsilon = mean*mean/var;
	alpha	= 2.0 + epsilon;
  	beta	= 1.0/(mean + mean*epsilon);
	}

inline const GammaDistribution &DirichletDistribution::GetDistributionOnParameter(
  unsigned i) const
  	{
  	PHYCAS_ASSERT( i < paramDistributions.size());
	return paramDistributions[i];
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	Construct the equivalent of a flat Beta distribution by default.
*/
inline DirichletDistribution::DirichletDistribution() 
	{
	dirParams.clear();
	scratchSpace.clear();
	paramDistributions.clear();
	dirParams.push_back(1.0);
	dirParams.push_back(1.0);
	scratchSpace.push_back(0.0);
	scratchSpace.push_back(0.0);
	paramDistributions.push_back(GammaDistribution(1.0, 1.0));
	paramDistributions.push_back(GammaDistribution(1.0, 1.0));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the random number generator used with the MultivariateProbabilityDistribution::Sample member function. 
|	The original random number generator (data member `myLot') can be replaced by calling the 
|	MultivariateProbabilityDistribution::ResetLot function. Note that this object does not take ownership of the Lot 
|	object whose pointer is specified as `other'. It is assumed that `other' is non-NULL.
|	This overridden member function differs from the base class version ProbabilityDistribution::SetLot in that it
|	calls SetLot for all its component distributions so that the same random number generator is used throughout.
*/
inline void DirichletDistribution::SetLot(
	Lot * other) /**< is a pointer to the random number generator object to be used subsequently by Sample */
	{
	if (other == NULL)
		throw XProbDist("attempt made to install a non-existent pseudorandom number generator");
	lot = other;
	for (unsigned i = 0; i < paramDistributions.size(); ++i)
		{
		paramDistributions[i].SetLot(lot);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the random number seed for the `myLot' data member. Note that if DirichletDistribution::SetLot has been 
|	called, calling DirichletDistribution::SetSeed is pointless because you will not be setting the seed for the 
|	correct random number generator!
*/
inline void	DirichletDistribution::SetSeed(
  unsigned rnseed)	/**< is the new seed value */
	{
	myLot.SetSeed(rnseed);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Makes the data member `lot' (which is used as the random number generator by the member function 
|	MultivariateProbabilityDistribution::Sample) point to the local data member `myLot'. This function only needs to be
|	called if MultivariateProbabilityDistribution::SetLot has been called previously to replace the random number 
|	generator used by MultivariateProbabilityDistribution::Sample. This overridden member function differs from the base
|	class version ProbabilityDistribution::ResetLot in that it calls SetLot for all its component distributions so that
|	the same random number generator is used throughout.
*/
inline void DirichletDistribution::ResetLot()
	{
	lot = &myLot;
	for (unsigned i = 0; i < paramDistributions.size(); ++i)
		{
		paramDistributions[i].SetLot(&myLot);
		}
	}

inline bool DirichletDistribution::IsDiscrete() const 
	{
	return false;
	}

inline std::string DirichletDistribution::GetDistributionName() const 
	{
	if (GetNParams() == 2)
		return "Beta";
	return "Dirichlet";
	}

inline std::string DirichletDistribution::GetDistributionDescription() const 
	{
	PHYCAS_ASSERT(GetNParams() > 1);
	std::string s;
	s << "DirichletDist(";
	//s << MakeStrPrintF("%.5f", dirParams[0]);
	s << str(boost::format("%#.5f") % dirParams[0]);
	for (unsigned i = 1; i < GetNParams(); ++i) 
		{
		//s << MakeStrPrintF(",%.5f", dirParams[i]);
		s << str(boost::format(",%#.5f") % dirParams[i]);
		}
	s << ')';
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function is needed because in Python the parameters of the distribution are supplied to the constructor via
|	a tuple, which requires an extra set of parentheses. For example, in Phycas, one could sample from a flat Dirichlet
|	as follows:
|	
|	begin phycas;
|		sample dist=dirichlet(1.0,1.0,1.0,1.0);
|	end;
|
|	Note only one set of parentheses surround the vector of Dirichlet parameters. In Python, we have:
|
|	>>> DirichletDist((1.0,1.0,1.0,1.0)).sample()
|
|	Note the extra set of parentheses needed in the Python representation.
*/
inline std::string DirichletDistribution::GetDescriptionForPython() const 
	{
	PHYCAS_ASSERT(GetNParams() > 1);
	std::string s;
	s << "DirichletDist((";
	//s << MakeStrPrintF("%.5f", dirParams[0]);
	s << str(boost::format("%#.5f") % dirParams[0]);
	for (unsigned i = 1; i < GetNParams(); ++i) 
		{
		//s << MakeStrPrintF(",%.5f", dirParams[i]);
		s << str(boost::format(", %#.5f") % dirParams[i]);
		}
	s << "))";
	return s;
	}

inline VecDbl DirichletDistribution::GetMean() const 
	{
	VecDbl retvect;

	double sum = 0.0;
	for (std::vector<double>::const_iterator dIt = dirParams.begin(); dIt != dirParams.end(); ++dIt) 
		sum += *dIt;

	for (unsigned i = 0; i < GetNParams(); ++i)
		{
		double x = dirParams[i]/sum;
		retvect.push_back(x);
		}

	return retvect;
	}

inline VecDbl DirichletDistribution::GetVar() const 
	{
	VecDbl retvect;

	double c = 0.0;
	for (std::vector<double>::const_iterator dIt = dirParams.begin(); dIt != dirParams.end(); ++dIt) 
		c += *dIt;

	double denom = c*c*(c + 1.0);

	for (unsigned i = 0; i < GetNParams(); ++i)
		{
		double c_i = dirParams[i];
		double v = c_i*(c - c_i)/denom;
		retvect.push_back(v);
		}

	return retvect;
	}

inline VecDbl DirichletDistribution::GetStdDev() const 
	{
	VecDbl retvect = GetVar();

	for (unsigned i = 0; i < GetNParams(); ++i)
		{
		retvect[i] = sqrt(retvect[i]);
		}

	return retvect;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Approximates the cumulative distribution function evaluated at the supplied point `x'. The precision of the 
|	approximation is controlled by `nsamples'. The approximation is done using a brute force approach: `nsamples'
|	samples are drawn from this Dirichlet distribution, and the proportion of those samples that are inside the region
|	defined by `x' is returned as the approximated CDF. The supplied point should be a vector of length k, where k is
|	one fewer than the number of parameters of the Dirichlet distribution. If `x' has length greater than k, all 
|	elements of `x' after the (k-1)th. will be ignored.
*/
inline double DirichletDistribution::ApproxCDF(
  const VecDbl &x,		/**< is the point to be evaluated */
  unsigned nsamples)	/**< is the number of samples to use in approximating the CDF */
  const 
	{
	PHYCAS_ASSERT(nsamples > 0);
	unsigned nparams = GetNParams();
	PHYCAS_ASSERT(x.size() >= nparams - 1);
	unsigned num_inside = 0;

	for (unsigned k = 0; k < nsamples; ++k)
		{
		double sum = 0.0;
		for (unsigned i = 0; i < nparams; ++i)
			{
			const double y = paramDistributions[i].Sample();
			scratchSpace[i] = y;
			sum += y;
			}

		bool inside = true;
		for (unsigned j = 0; j < nparams - 1; ++j)
			{
			const double p = scratchSpace[j]/sum;
			if (p > x[j])
				{
				inside = false;
				break;
				}
			}

		if (inside)
			num_inside++;
		}

	return (double)num_inside/(double)nsamples;
	}

inline VecDbl DirichletDistribution::Sample() const 
	{
	VecDbl x;

	double sum = 0.0;
	for (unsigned i = 0; i < GetNParams(); ++i)
		{
		const double y = paramDistributions[i].Sample();
		scratchSpace[i] = y;
		sum += y;
		}

	for (unsigned j = 0; j < GetNParams(); ++j)
		x.push_back(scratchSpace[j]/sum);

	return x;
	}

inline double DirichletDistribution::GetLnPDF(
  const VecDbl &x)
  const 
	{
	PHYCAS_ASSERT(GetNParams() == (unsigned)x.size());
	double retVal = 0.0;

	double c = 0.0;
	for (unsigned i = 0; i < GetNParams(); ++i)
		{
		double c_i = dirParams[i];
		c += c_i;
		retVal += (c_i - 1.0) * std::log(x[i]);
		retVal -= cdf.LnGamma(c_i);
		}

	retVal += cdf.LnGamma(c);

	return retVal;
	}

inline double DirichletDistribution::GetRelativeLnPDF(
  const VecDbl &x)
  const 
	{
	PHYCAS_ASSERT(GetNParams() == (unsigned)x.size());
	double retVal = 0.0;

	for (unsigned i = 0; i < GetNParams(); ++i)
		retVal += (dirParams[i] - 1.0) * std::log(x[i]);

	return retVal;
	}

#if defined(PYTHON_ONLY)
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the variance-covariance matrix in the form of a numarray object.
*/
inline boost::python::numeric::array DirichletDistribution::GetVarCovarMatrix()
	{
	unsigned i, j;
	unsigned dim = (unsigned)dirParams.size();

	std::vector<int> dims;
	dims.push_back((int)dim);
	dims.push_back((int)dim);

	double c = 0.0;
	for (i = 0; i < dim; ++i)
		{
		c += dirParams[i];
		}
	double denom = c*c*(c + 1.0);

	VecDbl V;
	for (i = 0; i < dim; ++i)
		{
		for (j = 0; j < dim; ++j)
			{
			if (i == j)
				{
				double var_i = dirParams[i]*(c - dirParams[i])/denom;
				V.push_back(var_i);
				}
			else
				{
				double cov_ij = -dirParams[i]*dirParams[j]/denom;
				V.push_back(cov_ij);
				}
			}
		}

	return num_util::makeNum(&V[0], dims);
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Sets parameters of the distribution from the vector of means and vector of variances (note: v is not the variance
|	covariance matrix, but just a single-dimensional array of variances - the covariances are not needed).
*/
inline void DirichletDistribution::SetMeanAndVariance(const VecDbl &m, const VecDbl &v)
	{
	unsigned m_length = (unsigned)m.size();

	// check to make sure mean and variance vectors are the same length
	PHYCAS_ASSERT(m_length == (unsigned)v.size());

	// check to make sure user isn't trying to change the dimension of the distribution
	PHYCAS_ASSERT(m_length == GetNParams());

	// create iterators
	// get pointers to the underlying C arrays
	double *mean     = (double *)(&m[0]);
	double *variance = (double *)(&v[0]);

	// calculate the sum of variances and the sum of squared means
	unsigned i;
	double sum_of_variances = 0.0;
	double sum_of_squared_means = 0.0;
	for (i = 0; i < m_length; ++i)
		{
		sum_of_squared_means += (mean[i]*mean[i]);
		sum_of_variances += variance[i];
		}

	// calculate c, the sum of parameters
	if (sum_of_variances <= 0.0)
		throw XProbDist("sum of variances supplied to DirichletDist.setMeanAndVariance must be positive");
	double c = (1.0 - sum_of_squared_means - sum_of_variances)/sum_of_variances;

	// set parameters of the distribution
	dirParams.clear();
	paramDistributions.clear();
	for (unsigned i = 0; i < m_length; ++i)
		{
		double c_i = c*mean[i];
		dirParams.push_back(c_i);
		paramDistributions.push_back(GammaDistribution(c_i, 1.0));
		}
	}

#if defined(PYTHON_ONLY)
/*----------------------------------------------------------------------------------------------------------------------
|	Identical to SetMeanAndVariance, but uses numarray objects rather than std::vector. This is temporary, only serving
|	as a template for how to import a numarray and access its data within a member function.
*/
inline void DirichletDistribution::AltSetMeanAndVariance(boost::python::numeric::array m, boost::python::numeric::array v)
	{
	unsigned m_length = (unsigned)m.nelements();

	// check to make sure mean and variance vectors are the same length
	PHYCAS_ASSERT(m_length == (unsigned)v.nelements());

	// check to make sure user isn't trying to change the dimension of the distribution
	PHYCAS_ASSERT(m_length == GetNParams());

	// get pointers to the underlying C arrays
	double *mean     = (double *)num_util::data(m);
	double *variance = (double *)num_util::data(v);

	// calculate the sum of variances and the sum of squared means
	unsigned i;
	double sum_of_variances = 0.0;
	double sum_of_squared_means = 0.0;
	for (i = 0; i < m_length; ++i)
		{
		sum_of_squared_means += (mean[i]*mean[i]);
		sum_of_variances += variance[i];
		}

	// calculate c, the sum of parameters
	if (sum_of_variances <= 0.0)
		throw XProbDist("sum of variances supplied to DirichletDist.setMeanAndVariance must be positive");
	double c = (1.0 - sum_of_squared_means - sum_of_variances)/sum_of_variances;

	// set parameters of the distribution
	dirParams.clear();
	paramDistributions.clear();
	for (i = 0; i < m_length; ++i)
		{
		double c_i = c*mean[i];
		dirParams.push_back(c_i);
		paramDistributions.push_back(GammaDistribution(c_i, 1.0));
		}
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of parameters.
*/
inline unsigned DirichletDistribution::GetNParams() const 
	{
	return (unsigned)dirParams.size();
	}

inline BetaDistribution::~BetaDistribution()
	{
	//std::cerr << "Deleting a BetaDistribution object" << std::endl;
	}

inline BetaDistribution::BetaDistribution(double a, double b)
	: alphaParam(a), betaParam(b)
	{}
	
inline bool BetaDistribution::IsDiscrete() const	
	{
	return false;
	}
	
inline std::string BetaDistribution::GetDistributionName() const	
	{
	return "Beta";
	}
	
inline std::string BetaDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("BetaDist(%.5f, %.5f)", alphaParam, betaParam);
	return str(boost::format("BetaDist(%#.5f, %#.5f)") % alphaParam % betaParam);
	}
	
inline double BetaDistribution::GetMean() const
	{
	return (alphaParam / (alphaParam+betaParam));
	}
	
inline double BetaDistribution::GetVar() const
	{
	const double ab = alphaParam + betaParam;
	return (alphaParam * betaParam)/ (ab * ab * (ab + 1.0));
	}
	
inline double BetaDistribution::GetStdDev() const
	{
	return std::sqrt(GetVar());
	}
	
inline double BetaDistribution::GetCDF(double x) const
	{
	if (x == 0.0)
		return 0.0;

	return cdf.CumBeta(x, alphaParam, betaParam);
	}
	
inline double BetaDistribution::Sample() const
	{
	double x_1 = cdf.SampleGamma(lot->Uniform(FILE_AND_LINE), alphaParam, 1.0); 
	double x_2 = cdf.SampleGamma(lot->Uniform(FILE_AND_LINE), betaParam, 1.0); 

	return x_1/(x_1 + x_2);
	}
	
inline double BetaDistribution::GetLnPDF(double x) const
	{
	if (x <= 0.0 || x >= 1.0)
		return -DBL_MAX;
	else
		{
		double lnpdf = ((alphaParam - 1.0)*std::log(x)) + ((betaParam - 1.0)*std::log(1.0 - x));
		double a_only = alphaParam;
		double b_only = betaParam;
		double a_plus_b = alphaParam + betaParam;

		lnpdf += cdf.LnGamma(a_plus_b);
		lnpdf -= cdf.LnGamma(a_only);
		lnpdf -= cdf.LnGamma(b_only);

		return lnpdf;
		}
	}
	
inline double BetaDistribution::GetRelativeLnPDF(double x) const
	{
	if (x <= 0.0 || x >= 1.0)
		return -DBL_MAX;
	else
		return ((alphaParam - 1.0)*std::log(x)) + ((betaParam - 1.0)*std::log(1.0 - x));
	}
	
inline void BetaDistribution::SetMeanAndVariance(double m, double v)
	{
	if (m < 0.0 || m > 1.0)
		throw XProbDist("specified mean is out of bounds, should be between 0.0 and 1.0");
	if (v < 0.0 || v > 1.0)
		throw XProbDist("specified variances is out of bounds, should be between 0.0 and 1.0");
	double temp = m - m*m - v;
	if (temp < 0.0)
		throw XProbDist("the mean/variance combination specified are incompatible with a Beta distribution");
	alphaParam = m*temp/v;
	betaParam = (1.0 - m)*temp/v;
	}

//############################################################################################
//###### EXPONENTIAL DISTRIBUTION INLINED FUNCTIONS ##########################################
//############################################################################################

/*----------------------------------------------------------------------------------------------------------------------
|	The probability density function of the exponential distrbution is 
|>
|	f(x) = (1/mu) exp(-x/mu)
|<
|	where mu is the mean of the exponential. Since the exponential is really a gamma distribution, with shape alpha
|	equal to 1 and scale beta equal to mu, beta can be substituted in for mu in the above equation. This function 
|	returns the natural log of the density function at `x', which is 
|>
|	ln[f(x)] = -x/mu - ln(mu)
|<
|	Once again, beta is really the mean, and can be substituted in for mu. If x < 0.0, returns -DBL_MAX.
*/
inline double ExponentialDistribution::GetLnPDF(
  double x)   const	/**< is the value for which the density function is to be evaluated */
	{
	if (x < 0.0)
		return -DBL_MAX;
	else
		{
		return -x/beta - std::log(beta);
		}
	}

} //namespace phycas

#endif

