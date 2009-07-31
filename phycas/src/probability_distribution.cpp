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

#if defined(USING_NUMARRAY)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#	define NO_IMPORT_ARRAY
#endif

#include "probability_distribution.hpp"
#include "ncl/nxsexception.h"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|   Virtual destructor.
*/
ProbabilityDistribution::~ProbabilityDistribution()
	{
	//std::cerr << "\n>>>>> ProbabilityDistribution dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the gamma function evaluated at the supplied value `x'. This function is a wrapper
|   around CDF::LnGamma.
*/
double ProbabilityDistribution::LnGamma(
  double x)
    {
    return cdf.LnGamma(x);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the shape and scale parameters of the underlying gamma distribution to those of the supplied 
|   distribution `other'.
*/
ExponentialDistribution::ExponentialDistribution(
  const ExponentialDistribution & other)	/* the exponential distribution to clone */ 
  : GammaDistribution(other.alpha, other.beta)
  	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object. 
*/
ExponentialDistribution * ExponentialDistribution::Clone() const
	{
    double mean = alpha*beta;
    PHYCAS_ASSERT(mean > 0.0);
    double mean_inv = 1.0/mean;
    return new ExponentialDistribution(mean_inv);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes parameter `p' to `other.p'.
*/
BernoulliDistribution::BernoulliDistribution(
  const BernoulliDistribution & other)	/* the bernoulli distribution to clone */
  : p(other.p)
  	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object. 
*/
BernoulliDistribution * BernoulliDistribution::Clone() const
	{
    return new BernoulliDistribution(p);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes parameter `p' to `other.p' and `n' to `other.n'. Also initializes `lnp' and `lnq' accordingly.
*/
BinomialDistribution::BinomialDistribution(
  const BinomialDistribution & other)	/* the binomial distribution to clone */
  : BernoulliDistribution(other.p), n(other.n)
  	{
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
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object. 
*/
BinomialDistribution * BinomialDistribution::Clone() const
	{
    return new BinomialDistribution(n, p);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes parameter `alphaParam' to `other.alphaParam' and `betaParam' to `other.betaParam'.
*/
BetaDistribution::BetaDistribution(
  const BetaDistribution & other)	/* the beta distribution to clone */
  : alphaParam(other.alphaParam), betaParam(other.betaParam)
	{}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes parameter `alphaParam' to `other.alphaParam' and `betaParam' to `other.betaParam'.
*/
BetaDistribution * BetaDistribution::Clone() const
	{
    return new BetaDistribution(alphaParam, betaParam);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes parameter `alphaParam' to `other.alphaParam' and `betaParam' to `other.betaParam'.
*/
BetaPrimeDistribution::BetaPrimeDistribution(
  const BetaPrimeDistribution & other)	/* the beta prime distribution to clone */
  : alphaParam(other.alphaParam), betaParam(other.betaParam)
	{}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes parameter `alphaParam' to `other.alphaParam' and `betaParam' to `other.betaParam'.
*/
BetaPrimeDistribution * BetaPrimeDistribution::Clone() const
	{
    return new BetaPrimeDistribution(alphaParam, betaParam);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing to be done in constructor.
*/
ImproperUniformDistribution::ImproperUniformDistribution(
  const ImproperUniformDistribution & other)	/* the improper uniform distribution to clone */
  	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object. 
*/
ImproperUniformDistribution * ImproperUniformDistribution::Clone() const
	{
    return new ImproperUniformDistribution();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the left bound `a' to `other.a' and the right bound `b' to `other.b'.
*/
UniformDistribution::UniformDistribution(
  const UniformDistribution & other)	/* the uniform distribution to clone */
  : a(other.a), b(other.b)
  	{
	log_density = -1.0*std::log(b-a);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object. 
*/
UniformDistribution * UniformDistribution::Clone() const
	{
    return new UniformDistribution(a, b);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `alpha' to `other.alpha' and `beta' to `other.beta'.
*/
GammaDistribution::GammaDistribution(
  const GammaDistribution & other)	/* the gamma distribution to clone */
  : alpha(other.alpha), beta(other.beta)
  	{
	ComputeLnConst();
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object. 
*/
GammaDistribution * GammaDistribution::Clone() const
	{
    return new GammaDistribution(alpha, beta);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `alpha' to `other.alpha' and `beta' to `other.beta'.
*/
InverseGammaDistribution::InverseGammaDistribution(
  const InverseGammaDistribution & other)	/* the inverse gamma distribution to clone */
  : GammaDistribution(other.alpha, other.beta)
  	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object. 
*/
InverseGammaDistribution * InverseGammaDistribution::Clone() const
	{
    return new InverseGammaDistribution(alpha, beta);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the `mean' and `sd' data members to `other.mean' and `other.sd', respectively. Also initializes 
|   `pi_const' and `sqrt2_const'. The following is relevant to the calculation of `pi_const':
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
NormalDistribution::NormalDistribution(
  const NormalDistribution & other)	/* the normal distribution to clone */
  : mean(other.mean), sd(other.sd)
  	{
	pi_const = 2.0*std::asin(1.0);
	sqrt2_const = std::sqrt(2.0);
	ComputeLnConst();
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object. 
*/
NormalDistribution * NormalDistribution::Clone() const
	{
    return new NormalDistribution(mean, sd);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `dirParams', `scratchSpace' and `paramDistributions' vectors to the values contained in their 
|   counterparts in `other'.
*/
DirichletDistribution::DirichletDistribution(
  const DirichletDistribution & other)	/* the dirichlet distribution to clone */
	{
	dirParams.resize(other.dirParams.size());
	std::copy(other.dirParams.begin(), other.dirParams.end(), dirParams.begin());

	scratchSpace.resize(other.scratchSpace.size());
	std::copy(other.scratchSpace.begin(), other.scratchSpace.end(), scratchSpace.begin());

    paramDistributions.clear();
    for (std::vector<GammaDistribution>::const_iterator it = other.paramDistributions.begin(); it != other.paramDistributions.end(); ++it)
	    paramDistributions.push_back(GammaDistribution(*it));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object. 
*/
DirichletDistribution * DirichletDistribution::Clone() const
	{
    return new DirichletDistribution(dirParams);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	The () operator calls the member function GetRelativeLnPDF.
*/
double ProbabilityDistribution::operator()(double x)
	{
	return GetRelativeLnPDF(x);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the random number generator used with the ProbabilityDistribution::Sample member function. The original 
|	random number generator	(data member `myLot') can be replaced by calling the ProbabilityDistribution::ResetLot 
|	function. Note that this object does not take ownership of the Lot object whose pointer is specified as `other'. 
|	It is assumed that `other' is non-NULL.
*/
void ProbabilityDistribution::SetLot(
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
void	ProbabilityDistribution::SetSeed(
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
void ProbabilityDistribution::ResetLot()
	{
	lot = &myLot;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This version of the function returns the sum of the relative log probability densities for all of the points in an
|	array. 
*/
double ProbabilityDistribution::GetRelativeLnPDFArray(
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
NormalDistribution::NormalDistribution()
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
NormalDistribution::NormalDistribution(
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
NormalDistribution::~NormalDistribution()
	{
	//std::cerr << "Deleting a NormalDistribution object" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false because the univariate normal is a continuous distribution.
*/
bool NormalDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Normal", which is the name of this distribution.
*/
std::string NormalDistribution::GetDistributionName() const
	{
	return "Normal";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Normal(<mean>, <stddev>)", with <mean> and <stddev> replaced with the current values of the 
|	`mean' and `sd' data members.
*/
std::string NormalDistribution::GetDistributionDescription() const
	{
	return str(boost::format("Normal(%#.5f, %#.5f)") % mean % sd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected mean of the normal distribution as currently specified, which is simply the value of the `mean'
|	data member.
*/
double NormalDistribution::GetMean() const
	{
	return mean;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected variance of the normal distribution as currently specified, which is simply the square of the 
|	value of the data member `sd'.
*/
double NormalDistribution::GetVar() const
	{
	return std::pow(sd, 2.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected standard deviation of the normal distribution as currently specified, which is simply the value
|	of the `sd' data member.
*/
double NormalDistribution::GetStdDev() const
	{
	return sd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative normal distribution for x (integral of normal density function from negative infinity to x). 
*/
double NormalDistribution::GetCDF(
  double x)	 const	/**< is the value for which the cumulative distribution function is to be evaluated */
	{
	return cdf.CumNorm(x, mean, sd);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a sampled value from a normal distribution having mean `mean' and standard deviation `sd'.
*/
double NormalDistribution::Sample() const
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
double NormalDistribution::GetRelativeLnPDF(
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
double NormalDistribution::GetLnPDF(
  double x)   const /* the value for which the density function is to be evaluated */
	{
	return (GetRelativeLnPDF(x) + ln_const);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows distribution to be specified in terms of the mean and variance rather than the shape and standard deviation.
|	This function sets `mean' to the supplied `mu' and `sd' to the square root of the supplied `var'. Assumes `var' is 
|	greater than zero.
*/
void NormalDistribution::SetMeanAndVariance(
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
void NormalDistribution::ComputeLnConst()
  	{
	ln_const = -log(sd) - log(2.0*pi_const)/2.0;
	}

//############################################################################################
//###### BERNOULLI DISTRIBUTION INLINED FUNCTIONS ############################################
//############################################################################################

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes p to 0.5.
*/
BernoulliDistribution::BernoulliDistribution()
  	{
	p = 0.5;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
BernoulliDistribution::~BernoulliDistribution()
	{
	//std::cerr << "Deleting a BernoulliDistribution object" << std::endl;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes parameter p to prob_success. Assumes p is greater than or equa to zero and less than or equal to 1.0.
*/
BernoulliDistribution::BernoulliDistribution(
  double prob_success)	/* probability of success on any given trial */
  	{
	if (prob_success < 0.0 || prob_success > 1.0)
		throw XProbDist("success probability supplied to BernoulliDist constructor out of bounds (should be between 0.0 and 1.0)");
	p = prob_success;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns true.
*/
bool BernoulliDistribution::IsDiscrete() const
	{
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Bernoulli", which is the name of this distribution.
*/
std::string BernoulliDistribution::GetDistributionName() const
	{
	return "Bernoulli";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Bernoulli(<psuccess>)", with <psuccess> replaced with actual probability of success.
*/
std::string BernoulliDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("Bernoulli(%.5f)", p);
	return str(boost::format("Bernoulli(%#.5f)") % p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected mean, which is simply the p parameter.
*/
double BernoulliDistribution::GetMean() const
  	{
	return p;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected variance, which is p*(1-p).
*/
double BernoulliDistribution::GetVar() const
  	{
	return p*(1.0 - p);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected standard deviation, which is sqrt(p*(1.0 - p)).
*/
double BernoulliDistribution::GetStdDev() const
  	{
	return std::sqrt(p*(1.0 - p));
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the cumulative distribution function evaluated at the specified value x. Assumes x is either 0.0 or 1.0.
*/
double BernoulliDistribution::GetCDF(
  double x) const	/* the value for which the cumulative distribution function is to be evaluated; specify either 0.0 or 1.0 */
  	{
	if (x != 0.0 && x != 1.0)
		throw XProbDist("only 0 or 1 are acceptable values for the BernoulliDist getCDF function");
	return (x == 0.0 ? 1.0 - p : 1.0);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns a sampled value from a bernoulli distribution having probabilty of success parameter p.
*/
double BernoulliDistribution::Sample() const
  	{
	return (lot->Uniform(FILE_AND_LINE) <= p ? 1.0 : 0.0);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	This function returns the natural log of the probability of x, which is p if x equals 1 and 1-p if x equals 0. If p
|	is exactly 1.0 and x is zero, or if p is exactly 0.0 and x is 1, the probability is zero and the value returned 
|	-DBL_MAX (which is as close as possible to negative infinity).
*/
double BernoulliDistribution::GetLnPDF(
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
double BernoulliDistribution::GetRelativeLnPDF(
  double x) const	/**< is the value for which the density function is to be evaluated; should be either 0.0 or 1.0 */
	{
	return GetLnPDF(x);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows distribution to be specified in terms of the mean rather than the probability of success on any given trial.
|	Because there is only one parameter, only the mean is used. The parameter p is set to the specified mean. The mean 
|	is assumed to lie in [0.0, 1.0]. The specified variance is ignored.
*/
void BernoulliDistribution::SetMeanAndVariance(
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
BinomialDistribution::BinomialDistribution()
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
BinomialDistribution::BinomialDistribution(
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
BinomialDistribution::~BinomialDistribution()
	{
	//std::cerr << "Deleting a BinomialDistribution object" << std::endl;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns true.
*/
bool BinomialDistribution::IsDiscrete() const
	{
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Binomial", which is the name of this distribution.
*/
std::string BinomialDistribution::GetDistributionName() const
	{
	return "Binomial";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Binomial(n, p)", where n is the number of trials and p the probability of success on each trial.
*/
std::string BinomialDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("Binomial(%d, %.5f)", (unsigned)n, p);
	unsigned nn = (unsigned)n;
	return str(boost::format("Binomial(%d, %#.5f)") % nn % p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected mean, which is n*p.
*/
double BinomialDistribution::GetMean() const
  	{
	return n*p;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected variance, which is n*p*(1-p).
*/
double BinomialDistribution::GetVar() const
  	{
	return n*p*q;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected standard deviation, which is sqrt(n*p*(1.0 - p)).
*/
double BinomialDistribution::GetStdDev() const
  	{
	return std::sqrt(n*p*q);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the cumulative distribution function evaluated at the specified value x. Assumes x is either 0.0 or a 
|	positive integer less than or equal to n.
*/
double BinomialDistribution::GetCDF(
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
double BinomialDistribution::Sample() const
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
double BinomialDistribution::GetLnPDF(
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
double BinomialDistribution::GetRelativeLnPDF(
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
void BinomialDistribution::SetMeanAndVariance(
  double mean,	/**< the mean of the binomial distribution */
  double var)	/**< ignored for the Binomial distribution */
  	{
	if (mean <= 0.0 || mean > n)
		throw XProbDist("mean must be greater than 0.0 and less than or equal to sample size for BinomialDist");

	p = mean/n;
	}

//############################################################################################
//###### IMPROPER UNIFORM DISTRIBUTION INLINED FUNCTIONS #####################################
//############################################################################################

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing to be done.
*/
ImproperUniformDistribution::ImproperUniformDistribution()
  	{
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Nothing to be done.
*/
ImproperUniformDistribution::~ImproperUniformDistribution()
	{
	//std::cerr << "Deleting a ImproperUniformDistribution object" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false.
*/
bool ImproperUniformDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Improper Uniform", which is the name of this distribution.
*/
std::string ImproperUniformDistribution::GetDistributionName() const
	{
	return "Improper Uniform";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "ImproperUniformDist()", which is the string that would be necessary on the Python side to 
|	create an object of this type.
*/
std::string ImproperUniformDistribution::GetDistributionDescription() const
	{
	return "ImproperUniform()";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The mean is not defined in this distribution, so an exception is thrown.
*/
double ImproperUniformDistribution::GetMean() const
  	{
	throw XProbDist("the mean is undefined for an improper uniform distribution");
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	The variance is not defined in this distribution, so an exception is thrown.
*/
double ImproperUniformDistribution::GetVar() const
  	{
	throw XProbDist("the variance is undefined for an improper uniform distribution");
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	The standard deviation is not defined in this distribution, so an exception is thrown.
*/
double ImproperUniformDistribution::GetStdDev() const
  	{
	throw XProbDist("the standard deviation is undefined for an improper uniform distribution");
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	The cumulative distribution function is not defined in this distribution, so an exception is thrown.
*/
double ImproperUniformDistribution::GetCDF(
  double x)	 const /**> is the value for which the cumulative distribution function is to be evaluated */
  	{
	throw XProbDist("the cumulative distribution function is undefined for an improper uniform distribution");
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	It is not possible to sample from an improper uniform distribution, so an exception is thrown.
*/
double ImproperUniformDistribution::Sample() const
  	{
	throw XProbDist("it is not possible to sample from an improper uniform distribution");
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	The probability density function of the improper uniform distribution is infinitesimal, so the log of this value is 
|	undefined. Therefore, an exception is thrown.
*/
double ImproperUniformDistribution::GetLnPDF(
  double x) const	/* the value for which the density function is to be evaluated */
	{
	throw XProbDist("the probability density function of an improper uniform distribution is not defined");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The relative probability density of any point in an improper uniform distribution can be any constant. The value 
|	1.0 is arbitrarily chosen here, making the log 0.0.
*/
double ImproperUniformDistribution::GetRelativeLnPDF(
  double x) const	/* the value for which the density function is to be evaluated */
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The mean and variance are both undefined for an improper uniform distribution, so an exception is thrown.
*/
void ImproperUniformDistribution::SetMeanAndVariance(
  double mean,	/* the mean of the uniform distribution */
  double var)	/* the variance of the uniform distribution */
  	{
	throw XProbDist("neither the mean nor the variance of an improper uniform distribution is defined");
	}

//############################################################################################
//###### UNIFORM DISTRIBUTION INLINED FUNCTIONS ##############################################
//############################################################################################

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the left bound to 0.0 and right bound to 1.0.
*/
UniformDistribution::UniformDistribution()
  	{
	a = 0.0;
	b = 1.0;
	log_density = 0.0;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the left bound (a) and right bound (b) parameters to the specified values. Assumes right_bound is 
|	greater than left_bound.
*/
UniformDistribution::UniformDistribution(
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
UniformDistribution::~UniformDistribution()
	{
	//std::cerr << "Deleting a UniformDistribution object" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false.
*/
bool UniformDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Uniform", which is the name of this distribution.
*/
std::string UniformDistribution::GetDistributionName() const
	{
	return "Uniform";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Uniform(<lower>, <upper>)", with <lower> and <upper> replaced with actual lower and upper bounds.
*/
std::string UniformDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("Uniform(%.5f, %.5f)", a, b);
	return str(boost::format("Uniform(%#.5f, %#.5f)") % a % b);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected mean, which is simply the arithmetic average of the a and b parameters.
*/
double UniformDistribution::GetMean() const
  	{
	return ((a + b)/2.0);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected variance, which is (b - a)^2/12.
*/
double UniformDistribution::GetVar() const
  	{
	return ((b - a)*(b - a)/12.0);
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns expected standard deviation, which is fabs(b - a)/sqrt(12).
*/
double UniformDistribution::GetStdDev() const
  	{
	return (std::fabs(b - a)/std::sqrt(12.0));
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the cumulative distribution function evaluated at the specified value x. If GetCDF(x) returns 0.4, this 
|	means that 40% of sampled values of X, where X is a uniform random variable, will be less than or equal to the 
|	value x. Assumes x is in the interval [a, b].
*/
double UniformDistribution::GetCDF(
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
double UniformDistribution::Sample() const
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
double UniformDistribution::GetLnPDF(
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
double UniformDistribution::GetRelativeLnPDF(
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
void UniformDistribution::SetMeanAndVariance(
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
GammaDistribution::GammaDistribution()
  	{
	alpha = 1.0;
	beta = 1.0;
	ln_const = 0.0;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
GammaDistribution::~GammaDistribution()
	{
	//std::cerr << "Deleting a GammaDistribution object" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the shape and scale parameters to the specified values.
*/
GammaDistribution::GammaDistribution(
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
bool GammaDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Gamma", which is the name of this distribution.
*/
std::string GammaDistribution::GetDistributionName() const
	{
	return "Gamma";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Gamma(<shape>, <scale>)", with <shape> and <scale> replaced with actual shape and scale 
|	parameters.
*/
std::string GammaDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("Gamma(%.5f, %.5f)", alpha, beta);
	return str(boost::format("Gamma(%#.5f, %#.5f)") % alpha % beta);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected mean of the gamma distribution as currently specified, which is the product of the shape and 
|	scale parameters.
*/
double GammaDistribution::GetMean() const
	{
	return alpha*beta;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected variance of the gamma distribution as currently specified, which is the product of the shape
|	parameters and the square of the scale parameter.
*/
double GammaDistribution::GetVar() const
	{
	return alpha*beta*beta;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected standard deviation of the gamma distribution as currently specified. Returns beta times the 
|	square root of alpha.
*/
double GammaDistribution::GetStdDev() const
	{
	return std::fabs(beta)*std::sqrt(alpha);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative gamma distribution for x (integral of gamma density function from 0.0 to x). Assumes x is 
|	greater than or equal to zero.
*/
double GammaDistribution::GetCDF(
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
double GammaDistribution::Sample() const
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
double GammaDistribution::GetRelativeLnPDF(
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
double GammaDistribution::GetLnPDF(
  double x)   const /* the value for which the density function is to be evaluated */
	{
	return (GetRelativeLnPDF(x) + ln_const);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows distribution to be specified in terms of the mean and variance rather than the shape and scale. This function 
|	sets the shape parameter to the mean squared divided by the variance, and sets the scale parameter to the variance 
|	divided by the mean. Assumes var is greater than zero.
*/
void GammaDistribution::SetMeanAndVariance(
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
void GammaDistribution::ComputeLnConst()
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
ExponentialDistribution::ExponentialDistribution()
  : GammaDistribution()
  	{
  	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the shape and scale parameters of the underlying gamma distribution to 1.0 and 1.0/lambda, respectively.
*/
ExponentialDistribution::ExponentialDistribution(
  double lambda)	/* the single hazard rate parameter of the exponential distribution */ 
  : GammaDistribution(1.0, 1.0/lambda)
  	{
	if (lambda <= 0.0)
		throw XProbDist("specified hazard parameter must be greater than 0.0 for ExponentialDist");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
ExponentialDistribution::~ExponentialDistribution()
	{
	//std::cerr << "ExponentialDistribution dying..." << std::endl;	
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false.
*/
bool ExponentialDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Exponential", which is the name of this distribution.
*/
std::string ExponentialDistribution::GetDistributionName() const
	{
	return "Exponential";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "Exponential(<hazard>)", with <hazard> replaced with actual hazard rate parameter.
*/
std::string ExponentialDistribution::GetDistributionDescription() const
	{
	double hazard_rate = 1.0/beta;
	//return MakeStrPrintF("Exponential(%.5f)", hazard_rate);
	return str(boost::format("Exponential(%#.5f)") % hazard_rate);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows distribution to be specified in terms of the mean and variance rather than the hazard rate parameter. This 
|	function sets shape = 1.0 and scale = mean. Throws XProbDist if mean is not greater than 0.0. This
|	distribution is entirely determined by the mean, so the variance argument is ignored.
*/
void ExponentialDistribution::SetMeanAndVariance(
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
InverseGammaDistribution::InverseGammaDistribution()
  	{
	alpha = 3.0;
	beta = 0.5;
	}
 
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes the shape and scale parameters to the specified values. Assumes shape is greater than 2.0 (variance is 
|	not defined unless alpha is greater than 2.0)
*/
InverseGammaDistribution::InverseGammaDistribution(
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
InverseGammaDistribution::~InverseGammaDistribution()
	{
	//std::cerr << "InverseGammaDistribution dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns false.
*/
bool InverseGammaDistribution::IsDiscrete() const
	{
	return false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "InverseGamma", which is the name of this distribution.
*/
std::string InverseGammaDistribution::GetDistributionName() const
	{
	return "InverseGamma";
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the string "InverseGamma(<shape>, <scale>)", with <shape> and <scale> replaced with actual shape and scale 
|	parameters.
*/
std::string InverseGammaDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("InverseGamma(%.5f, %.5f)", alpha, beta);
	return str(boost::format("InverseGamma(%#.5f, %#.5f)") % alpha % beta);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected mean of the inverse gamma distribution as currently specified, which is 1/(beta * (alpha - 1)).
*/
double InverseGammaDistribution::GetMean() const
	{
	double denom = beta * (alpha - 1.0);
	return (1.0/denom);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the expected variance of the inverse gamma distribution as currently specified, which is 
|	1/(beta^2 * (alpha - 1)^2 * (alpha - 2)).
*/
double InverseGammaDistribution::GetVar() const
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
double InverseGammaDistribution::GetStdDev() const
	{
	double denom = beta*(alpha - 1.0)*std::sqrt(alpha - 2.0);
	return 1.0/denom;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative inverse gamma distribution for x (integral of inverse gamma density function from 0.0 to x). 
|	Assumes x is greater than or equal to zero.
*/
double InverseGammaDistribution::GetCDF(
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
double InverseGammaDistribution::Sample() const
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
double InverseGammaDistribution::GetLnPDF(
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
double InverseGammaDistribution::GetRelativeLnPDF(
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
void InverseGammaDistribution::SetMeanAndVariance(
  double mean,	/* the mean of the inverse gamma distribution */
  double var)	/* the variance of the inverse gamma distribution */
  	{
	if (var <= 0.0)
		throw XProbDist("specified variance must be greater than 0.0 for InverseGammaDist");
	double epsilon = mean*mean/var;
	alpha	= 2.0 + epsilon;
  	beta	= 1.0/(mean + mean*epsilon);
	}

//#######################################################################################
//###### DIRICHLET DISTRIBUTION MEMBER FUNCTIONS ########################################
//#######################################################################################

const GammaDistribution &DirichletDistribution::GetDistributionOnParameter(
  unsigned i) const
  	{
  	PHYCAS_ASSERT( i < paramDistributions.size());
	return paramDistributions[i];
  	}

/*----------------------------------------------------------------------------------------------------------------------
|	Construct the equivalent of a flat Beta distribution by default.
*/
DirichletDistribution::DirichletDistribution() 
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
	for (unsigned i = 0; i < paramDistributions.size(); ++i)
		{
		paramDistributions[i].SetLot(&myLot);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a DirichletDistribution based on the parameter values supplied in the vector `params'.
*/
DirichletDistribution::DirichletDistribution(const std::vector<double> &params) 
	{
	unsigned params_size = (unsigned)params.size();
	//std::cerr << "params_size = " << params_size << std::endl;
	PHYCAS_ASSERT(params_size > 1);

	// Reserve necessary space to save copying
	//
	dirParams.reserve(params_size);
	scratchSpace.reserve(params_size);
	paramDistributions.reserve(params_size);

#if 0   //@POL Mark, this causes an "out of keys" compiler error on VC
#	if defined NCL_NXS_THROW_UNDEFINED
		if (params.size() < 2)
			throw NxsX_UndefinedException("Illegal Dirichlet", __FILE__, __LINE__);
#	endif
#endif

	std::vector<double>::const_iterator iter;
	for (iter = params.begin(); iter != params.end(); ++iter)
		{
		double c_i = *iter;
		dirParams.push_back(c_i);
		scratchSpace.push_back(0.0);
    	paramDistributions.push_back(GammaDistribution(c_i, 1.0));
		}

	for (unsigned i = 0; i < paramDistributions.size(); ++i)
		{
		paramDistributions[i].SetLot(&myLot);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the random number generator used with the MultivariateProbabilityDistribution::Sample member function. 
|	The original random number generator (data member `myLot') can be replaced by calling the 
|	MultivariateProbabilityDistribution::ResetLot function. Note that this object does not take ownership of the Lot 
|	object whose pointer is specified as `other'. It is assumed that `other' is non-NULL.
|	This overridden member function differs from the base class version ProbabilityDistribution::SetLot in that it
|	calls SetLot for all its component distributions so that the same random number generator is used throughout.
*/
void DirichletDistribution::SetLot(
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
void	DirichletDistribution::SetSeed(
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
void DirichletDistribution::ResetLot()
	{
	lot = &myLot;
	for (unsigned i = 0; i < paramDistributions.size(); ++i)
		{
		paramDistributions[i].SetLot(&myLot);
		}
	}

bool DirichletDistribution::IsDiscrete() const 
	{
	return false;
	}

std::string DirichletDistribution::GetDistributionName() const 
	{
	if (GetNParams() == 2)
		return "Beta";
	return "Dirichlet";
	}

std::string DirichletDistribution::GetDistributionDescription() const 
	{
	PHYCAS_ASSERT(GetNParams() > 1);
	std::string s;
	s << "Dirichlet(";
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
std::string DirichletDistribution::GetDescriptionForPython() const 
	{
	PHYCAS_ASSERT(GetNParams() > 1);
	std::string s;
	s << "Dirichlet((";
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

std::vector<double> DirichletDistribution::GetMean() const 
	{
	std::vector<double> retvect;

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

std::vector<double> DirichletDistribution::GetVar() const 
	{
	std::vector<double> retvect;

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

std::vector<double> DirichletDistribution::GetStdDev() const 
	{
	std::vector<double> retvect = GetVar();

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
double DirichletDistribution::ApproxCDF(
  const std::vector<double> &x,		/**< is the point to be evaluated */
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

std::vector<double> DirichletDistribution::Sample() const 
	{
	std::vector<double> x;

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

double DirichletDistribution::GetLnPDF(
  const std::vector<double> &x)
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

double DirichletDistribution::GetRelativeLnPDF(
  const std::vector<double> &x)
  const 
	{
	PHYCAS_ASSERT(GetNParams() == (unsigned)x.size());
	double retVal = 0.0;

	for (unsigned i = 0; i < GetNParams(); ++i)
		retVal += (dirParams[i] - 1.0) * std::log(x[i]);

	return retVal;
	}

#if defined(PYTHON_ONLY)
#if defined(USING_NUMARRAY)
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the variance-covariance matrix in the form of a numarray object.
*/
boost::python::numeric::array DirichletDistribution::GetVarCovarMatrix()
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

	std::vector<double> V;
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
#else	// USING_NUMARRAY not #defined
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the variance-covariance matrix in the form of a 1-dimensional vector (first row, then second row, etc.)
*/
std::vector<double> DirichletDistribution::GetVarCovarMatrix()
	{
	unsigned i, j;
	unsigned dim = (unsigned)dirParams.size();

	double c = 0.0;
	for (i = 0; i < dim; ++i)
		{
		c += dirParams[i];
		}
	double denom = c*c*(c + 1.0);

	std::vector<double> V;
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

	return V;
	}
#endif
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Sets parameters of the distribution from the vector of means and vector of variances (note: v is not the variance
|	covariance matrix, but just a single-dimensional array of variances - the covariances are not needed).
*/
void DirichletDistribution::SetMeanAndVariance(const std::vector<double> &m, const std::vector<double> &v)
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

#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
/*----------------------------------------------------------------------------------------------------------------------
|	Identical to SetMeanAndVariance, but uses numarray objects rather than std::vector. This is temporary, only serving
|	as a template for how to import a numarray and access its data within a member function.
*/
void DirichletDistribution::AltSetMeanAndVariance(boost::python::numeric::array m, boost::python::numeric::array v)
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
unsigned DirichletDistribution::GetNParams() const 
	{
	return (unsigned)dirParams.size();
	}

BetaDistribution::~BetaDistribution()
	{
	//std::cerr << "Deleting a BetaDistribution object" << std::endl;
	}

BetaDistribution::BetaDistribution(double a, double b)
	: alphaParam(a), betaParam(b)
	{}
	
bool BetaDistribution::IsDiscrete() const	
	{
	return false;
	}
	
std::string BetaDistribution::GetDistributionName() const	
	{
	return "Beta";
	}
	
std::string BetaDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("Beta(%.5f, %.5f)", alphaParam, betaParam);
	return str(boost::format("Beta(%#.5f, %#.5f)") % alphaParam % betaParam);
	}
	
double BetaDistribution::GetMean() const
	{
	return (alphaParam / (alphaParam+betaParam));
	}
	
double BetaDistribution::GetVar() const
	{
	const double ab = alphaParam + betaParam;
	return (alphaParam * betaParam)/ (ab * ab * (ab + 1.0));
	}
	
double BetaDistribution::GetStdDev() const
	{
	return std::sqrt(GetVar());
	}
	
double BetaDistribution::GetCDF(double x) const
	{
	if (x == 0.0)
		return 0.0;

	return cdf.CumBeta(x, alphaParam, betaParam);
	}
	
double BetaDistribution::GetQuantile(
  double p) const   /**< is the integral from 0.0 up to x (x is the value returned) */
	{
	if (p <= 0.0)
		return 0.0;
    else if (p >= 1.0)
        return 1.0;

	return cdf.BetaQuantile(p, alphaParam, betaParam);
	}
	
double BetaDistribution::Sample() const
	{
	double x_1 = cdf.SampleGamma(lot->Uniform(FILE_AND_LINE), alphaParam, 1.0); 
	double x_2 = cdf.SampleGamma(lot->Uniform(FILE_AND_LINE), betaParam, 1.0); 

	return x_1/(x_1 + x_2);
	}
	
double BetaDistribution::GetLnPDF(double x) const
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
	
double BetaDistribution::GetRelativeLnPDF(double x) const
	{
	if (x <= 0.0 || x >= 1.0)
		return -DBL_MAX;
	else
		return ((alphaParam - 1.0)*std::log(x)) + ((betaParam - 1.0)*std::log(1.0 - x));
	}
	
void BetaDistribution::SetMeanAndVariance(double m, double v)
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

////////////////// beta prime start
BetaPrimeDistribution::~BetaPrimeDistribution()
	{
	//std::cerr << "Deleting a BetaPrimeDistribution object" << std::endl;
	}

BetaPrimeDistribution::BetaPrimeDistribution(double a, double b)
	: alphaParam(a), betaParam(b)
	{}
	
bool BetaPrimeDistribution::IsDiscrete() const	
	{
	return false;
	}
	
std::string BetaPrimeDistribution::GetDistributionName() const	
	{
	return "BetaPrime";
	}
	
std::string BetaPrimeDistribution::GetDistributionDescription() const
	{
	//return MakeStrPrintF("BetaPrime(%.5f, %.5f)", alphaParam, betaParam);
	return str(boost::format("BetaPrime(%#.5f, %#.5f)") % alphaParam % betaParam);
	}
	
double BetaPrimeDistribution::GetMean() const
	{
	if (betaParam <= 1.0)
		throw XProbDist("mean of BetaPrime distribution undefined");

	return (alphaParam/(betaParam - 1.0));
	}
	
double BetaPrimeDistribution::GetVar() const
	{
	if (betaParam <= 2.0)
		throw XProbDist("variance of BetaPrime distribution undefined");
	const double a = alphaParam;
	const double b = betaParam;
	double var = a*(a + b - 1.0)/((b - 2.0)*(b - 1.0)*(b - 1.0));
	return var;
	}
	
double BetaPrimeDistribution::GetStdDev() const
	{
	return std::sqrt(GetVar());
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|   Assume X is a BetaPrime(a,b) random variable and Y is a Beta(a,b) random variable. Since F_X(x) = F_Y(y), the 
|	cumulative distribution of a BetaPrime(a,b) random variable X can be obtained as the cumulative distribution of the 
|	corresponding Beta(a,b) random variable Y = X/(X+1).
*/
double BetaPrimeDistribution::GetCDF(double x) const
	{
	if (x <= 0.0)
		return 0.0;

	double y = x/(1.0 + x);
	return cdf.CumBeta(y, alphaParam, betaParam);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|   Assume X is a BetaPrime(a,b) random variable and Y is a Beta(a,b) random variable. Since F_X(x) = F_Y(y), the 
|	quantiles of a BetaPrime(a,b) distribution can be obtained by transforming the corresponding quantile of a Beta(a,b)
|	distribution. If x is the value such that F_X(x) = p, then y = x/(1 - x) is the value such that F_Y(y) = p.
*/
double BetaPrimeDistribution::GetQuantile(
  double p) const   /**< is the integral from 0.0 up to x (x is the value returned) */
	{
	if (p <= 0.0)
		return 0.0;
		
	if (p >= 1.0)
		throw XProbDist("argument supplied to BetaPrimeDistribution::GetQuantile out of range: must be less than 1.0");

	double x = cdf.BetaQuantile(p, alphaParam, betaParam);
	double y = x/(1.0 - x);
	return y;
	}
	
double BetaPrimeDistribution::Sample() const
	{
	double x_1 = cdf.SampleGamma(lot->Uniform(FILE_AND_LINE), alphaParam, 1.0);
	double x_2 = 0.0;
	while (x_2 == 0.0)
		x_2 = cdf.SampleGamma(lot->Uniform(FILE_AND_LINE), betaParam, 1.0); 

	return x_1/x_2;
	}
	
double BetaPrimeDistribution::GetLnPDF(double x) const
	{
	if (x <= 0.0)
		return -DBL_MAX;
	else
		{
		double a_only = alphaParam;
		double b_only = betaParam;
		double a_plus_b = alphaParam + betaParam;
		double lnpdf = (a_only - 1.0)*std::log(x);
		lnpdf -= a_plus_b*std::log(1.0 + x);
		lnpdf += cdf.LnGamma(a_plus_b);
		lnpdf -= cdf.LnGamma(a_only);
		lnpdf -= cdf.LnGamma(b_only);

		return lnpdf;
		}
	}
	
double BetaPrimeDistribution::GetRelativeLnPDF(double x) const
	{
	if (x <= 0.0)
		return -DBL_MAX;
	else
		return ((alphaParam - 1.0)*std::log(x)) - ((alphaParam + betaParam - 1.0)*std::log(1.0 + x));
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Mean of a BetaPrime(a,b) random variable is a/(b - 1) and variance is a*(a + b - 1)/[(b - 2)*(b - 1)^2]. Thus, the 
|	mean is undefined if b <= 1 and the variance is undefined if b <= 2. Given mean m and variance v, one can find a and
|	b as follows. Letting c = m*(m+1)/v, a = m*(c+1) and b = c+2.
*/
void BetaPrimeDistribution::SetMeanAndVariance(double m, double v)
	{
	if (m <= 0.0)
		throw XProbDist("specified mean is out of bounds for a BetaPrime distribution; should be greater than 0.0");
	if (v <= 0.0)
		throw XProbDist("specified variance is out of bounds, should be greater than 0.0");
	double c   = m*(m + 1.0)/v;
	alphaParam = m*(c + 1.0);
	betaParam  = c + 2.0;
	}
////////////////// beta prime end

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
double ExponentialDistribution::GetLnPDF(
  double x)   const	/**< is the value for which the density function is to be evaluated */
	{
	if (x < 0.0)
		return -DBL_MAX;
	else
		{
		return -x/beta - std::log(beta);
		}
	}

} // namespace phycas
