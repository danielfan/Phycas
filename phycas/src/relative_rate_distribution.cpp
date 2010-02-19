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

#include <numeric>
#include <functional>
#include <cmath>
#include <boost/lambda/lambda.hpp>

#include "phycas/src/relative_rate_distribution.hpp"
//#include "ncl/nxsexception.h"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Construct the equivalent of a flat Beta distribution by default.
*/
RelativeRateDistribution::RelativeRateDistribution() 
  : DirichletDistribution(), dim(2)
	{
	sum_params = std::accumulate(dirParams.begin(), dirParams.end(), 0.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a RelativeRateDistribution based on the parameter values supplied in the vector `params'.
*/
RelativeRateDistribution::RelativeRateDistribution(
  const double_vect_t & params) /**< is the vector of parameters, which are assumed to have mean 1.0 (i.e. the sum of the elements of the `params' vector is equal to the length of the `params' vector) */
  : DirichletDistribution(params), dim(2)
	{
	PHYCAS_ASSERT(params.size() > 1);
	dim = (unsigned)params.size();
	sum_params = std::accumulate(dirParams.begin(), dirParams.end(), 0.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `dirParams', `scratchSpace' and `paramDistributions' vectors to the values contained in their 
|   counterparts in `other'.
*/
RelativeRateDistribution::RelativeRateDistribution(
  const RelativeRateDistribution & other)	/* the relative rate distribution to clone */
	{
	dirParams.resize(other.dirParams.size());
	std::copy(other.dirParams.begin(), other.dirParams.end(), dirParams.begin());

	scratchSpace.resize(other.scratchSpace.size());
	std::copy(other.scratchSpace.begin(), other.scratchSpace.end(), scratchSpace.begin());

    paramDistributions.clear();
    for (std::vector<GammaDistribution>::const_iterator it = other.paramDistributions.begin(); it != other.paramDistributions.end(); ++it)
	    paramDistributions.push_back(GammaDistribution(*it));
		
	dim = other.dim;
	sum_params = other.sum_params;
	}

std::string RelativeRateDistribution::GetDistributionName() const 
	{
	return "RelativeRateDistribution";
	}

std::string RelativeRateDistribution::GetDistributionDescription() const 
	{
	PHYCAS_ASSERT(dim > 1);
	std::string s;
	s << "RelativeRateDistribution(";
	for (double_vect_t::const_iterator it = dirParams.begin(); it != dirParams.end(); ++it)
		{
		if (it == dirParams.begin())
			s << boost::str(boost::format("%#.5f") % (*it));
		else 
			s << boost::str(boost::format(",%#.5f") % (*it));
		}
	s << ')';
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function is needed because in Python the parameters of the distribution are supplied to the constructor via
|	a tuple, which requires an extra set of parentheses. For example, in Phycas, one could sample from a 5-parameter
|	relative rate distribution as follows:
|	
|	>>> RelativeRateDistribution((0.8,0.9,1.0,1.1,1.2)).sample()
|
|	Note the extra set of parentheses needed in the Python representation.
*/
std::string RelativeRateDistribution::GetDescriptionForPython() const 
	{
	PHYCAS_ASSERT(dim > 1);
	std::string s;
	s << "RelativeRateDistribution((";
	for (double_vect_t::const_iterator it = dirParams.begin(); it != dirParams.end(); ++it)
		{
		if (it == dirParams.begin())
			s << boost::str(boost::format("%#.5f") % (*it));
		else 
			s << boost::str(boost::format(",%#.5f") % (*it));
		}
	s << "))";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of doubles representing the mean of this RelativeRateDistribution. The mean is `dim' times the
|	mean vector for a Dirichlet random variable with the same parameters.
*/
double_vect_t RelativeRateDistribution::GetMean() const 
	{
	double_vect_t retvect(dim, 0.0);
	double x = (double)dim;
	
	//std::cerr << "@@@@@@@@@@ RelativeRateDistribution::GetMean(): dirParams @@@@@@@@@@" << std::endl;
	//std::copy(dirParams.begin(), dirParams.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << "@@@@@@@@@@ sum_params = " << sum_params << std::endl;
	
	std::transform(dirParams.begin(), dirParams.end(), retvect.begin(), x*boost::lambda::_1/sum_params);

	//std::cerr << "@@@@@@@@@@ RelativeRateDistribution::GetMean(): retvect @@@@@@@@@@" << std::endl;
	//std::copy(retvect.begin(), retvect.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << "@@@@@@@@@@ sum_params = " << sum_params << std::endl;
	
	return retvect;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of doubles representing the variance of this RelativeRateDistribution. The variance is `dim' 
|	squared times the variance vector for a Dirichlet random variable with the same parameters.
*/
double_vect_t RelativeRateDistribution::GetVar() const 
	{
	double_vect_t retvect(dim, 0.0);

	double x = (double)dim;
	double c = sum_params;
	double denom = c*c*(c + 1.0);
	std::transform(dirParams.begin(), dirParams.end(), retvect.begin(), x*x*(boost::lambda::_1)*(c - boost::lambda::_1)/denom);

	return retvect;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of doubles representing the square root of this RelativeRateDistribution. The square root is `dim' 
|	times the square root vector for a Dirichlet random variable with the same parameters.
*/
double_vect_t RelativeRateDistribution::GetStdDev() const 
	{
	double_vect_t retvect(dim, 0.0);

	double x = (double)dim;
	double c = sum_params;
	double denom = c*c*(c + 1.0);
	double_vect_t::iterator dest = retvect.begin();
	for (double_vect_t::const_iterator it = dirParams.begin(); it != dirParams.end(); ++it)
		{
		double v = (*it);
		double dirvar = v*(c - v)/denom;
		*dest++ = x*std::sqrt(dirvar);
		}
	return retvect;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Approximates the cumulative distribution function evaluated at the supplied point `x'. The precision of the 
|	approximation is controlled by `nsamples'. The approximation is done using a brute force approach: `nsamples'
|	samples are drawn from this RelativeRateDistribution, and the proportion of those samples that are inside the region
|	defined by `x' is returned as the approximated CDF. The supplied point should be a vector of length k, where k is
|	one fewer than the number of parameters of the RelativeRateDistribution. If `x' has length greater than k, all 
|	elements of `x' after the (k-1)th. will be ignored.
*/
double RelativeRateDistribution::ApproxCDF(
  const double_vect_t & x,		/**< is the point to be evaluated */
  unsigned nsamples)			/**< is the number of samples to use in approximating the CDF */
  const 
	{
	PHYCAS_ASSERT(nsamples > 0);
	PHYCAS_ASSERT(x.size() >= dim - 1);
	double_vect_t tmp(x.size(), 0.0);
	double z = (double)dim;
	std::transform(x.begin(), x.end(), tmp.begin(), boost::lambda::_1/z);
	return DirichletDistribution::ApproxCDF(tmp, nsamples);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Samples and returns a `dim'-dimensional point from this distribution.
*/
double_vect_t RelativeRateDistribution::Sample() const 
	{
	double_vect_t x(dim, 0.0);

	double sum = 0.0;
	for (unsigned i = 0; i < dim; ++i)
		{
		const double y = paramDistributions[i].Sample();
		scratchSpace[i] = y;
		sum += y;
		}

	double z = (double)dim;
	std::transform(scratchSpace.begin(), scratchSpace.end(), x.begin(), z*boost::lambda::_1/sum);
	return x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the probability density function evaluated at the supplied point `x'. To derive
|	the relative rate distribution density, imagine a random variable X ~ Beta(a,b) and let Y=2X. Y now has a relative
|	rate distribution with density
|>
|	f_Y(y) = f_X(y/2) dX/dY
|	dX/dY = 1/2
|>
|	Thus, the density of Y will be one half the Beta density evaluated at y/2. This generalizes to the Dirichlet: the
|	density of a value y of a random variable Y having an n-parameter relative rate distribution will be 1/n times the
|	corresponding Dirichlet pdf evaluated at the point y/n.
*/
double RelativeRateDistribution::GetLnPDF(
  const double_vect_t & x) const 	/**< is the point at which to evaluate the density */
	{
	//std::cerr << "@@@@@@@@@@ RelativeRateDistribution::GetLnPDF(): x @@@@@@@@@@" << std::endl;
	//std::copy(x.begin(), x.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << "@@@@@@@@@@ sum_params = " << sum_params << std::endl;

	double_vect_t tmp(dim, 0.0);
	std::transform(x.begin(), x.end(), tmp.begin(), boost::lambda::_1/(double)dim);

	//std::cerr << "@@@@@@@@@@ RelativeRateDistribution::GetLnPDF(): tmp @@@@@@@@@@" << std::endl;
	//std::copy(tmp.begin(), tmp.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << "@@@@@@@@@@ sum_params = " << sum_params << std::endl;

	double ln_pdf = DirichletDistribution::GetLnPDF(tmp);
	ln_pdf -= log((double)dim);
	return ln_pdf;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural logarithm of the probability density function evaluated at the supplied point `x'.
*/
double RelativeRateDistribution::GetRelativeLnPDF(
  const double_vect_t & x) const 
	{
	return GetLnPDF(x);
	}

#if defined(PYTHON_ONLY)
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the variance-covariance matrix in the form of a 1-dimensional vector (first row, then second row, etc.)
*/
double_vect_t RelativeRateDistribution::GetVarCovarMatrix()
	{
	unsigned i, j;

	double x = (double)dim;
	double c = sum_params;
	double denom = c*c*(c + 1.0);

	double_vect_t V;
	for (i = 0; i < dim; ++i)
		{
		for (j = 0; j < dim; ++j)
			{
			if (i == j)
				{
				double var_i = x*x*dirParams[i]*(c - dirParams[i])/denom;
				V.push_back(var_i);
				}
			else
				{
				double cov_ij = -x*x*dirParams[i]*dirParams[j]/denom;
				V.push_back(cov_ij);
				}
			}
		}

	return V;
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Sets parameters of the distribution from the vector of means and vector of variances (note: v is not the variance
|	covariance matrix, but just a single-dimensional array of variances - the covariances are not needed).
*/
void RelativeRateDistribution::SetMeanAndVariance(
  const double_vect_t & m, 	/**< is the vector of means */
  const double_vect_t & v) 	/**< is the vector of variances */
	{
	if (dim != (unsigned)m.size())
		throw XProbDist(boost::str(boost::format("Mean vector supplied to SetMeanAndVariance has dimension %d but RelativeRateDistribution has dimension %d") % (unsigned)m.size() % dim));
	if (dim != (unsigned)v.size())
		throw XProbDist(boost::str(boost::format("Variance vector supplied to SetMeanAndVariance has dimension %d but RelativeRateDistribution has dimension %d") % (unsigned)v.size() % dim));
		
	double x = (double)dim;

	//std::cerr << "@@@@@@@@@@ RelativeRateDistribution::SetMeanAndVariance(): m @@@@@@@@@@" << std::endl;
	//std::copy(m.begin(), m.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << "@@@@@@@@@@ sum_params = " << sum_params << std::endl;

	double_vect_t transformed_means(dim);
	std::transform(m.begin(), m.end(), transformed_means.begin(), boost::lambda::_1/x);
	
	//std::cerr << "@@@@@@@@@@ RelativeRateDistribution::SetMeanAndVariance(): transformed_means @@@@@@@@@@" << std::endl;
	//std::copy(transformed_means.begin(), transformed_means.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << "@@@@@@@@@@ sum_params = " << sum_params << std::endl;
	
	//std::cerr << "@@@@@@@@@@ RelativeRateDistribution::SetMeanAndVariance(): v @@@@@@@@@@" << std::endl;
	//std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << "@@@@@@@@@@ sum_params = " << sum_params << std::endl;
	
	double xsquared = x*x;
	double_vect_t transformed_variances(dim);
	std::transform(v.begin(), v.end(), transformed_variances.begin(), boost::lambda::_1/xsquared);
	
	//std::cerr << "@@@@@@@@@@ RelativeRateDistribution::SetMeanAndVariance(): transformed_variances @@@@@@@@@@" << std::endl;
	//std::copy(transformed_variances.begin(), transformed_variances.end(), std::ostream_iterator<double>(std::cerr, " "));
	//std::cerr << "@@@@@@@@@@ sum_params = " << sum_params << std::endl;
	
	DirichletDistribution::SetMeanAndVariance(transformed_means, transformed_variances);
	sum_params = std::accumulate(dirParams.begin(), dirParams.end(), 0.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of parameters.
*/
unsigned RelativeRateDistribution::GetNParams() const 
	{
	return dim;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object, calls the new objects SetLot member function (passing the 
|	supplied Lot object `other'), and returns a pointer to it. The caller is expected to manage the new object. 
*/
RelativeRateDistribution * RelativeRateDistribution::cloneAndSetLot(Lot * other) const
	{
    RelativeRateDistribution * clone = new RelativeRateDistribution(*this);
	clone->SetLot(other);
	return clone;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new object that is a clone of this object and returns a pointer to it. Caller is expected to manage the
|   new object. 
*/
RelativeRateDistribution * RelativeRateDistribution::Clone() const
	{
    return new RelativeRateDistribution(*this);
    }

} // namespace phycas


