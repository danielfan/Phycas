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

#if !defined(PROBABILITY_DISTRIBUTION_HPP)
#define PROBABILITY_DISTRIBUTION_HPP

#if defined(_MSC_VER)
#	pragma warning(disable: 4267)	// warning about loss of data when converting size_t to int
#endif

#include <cmath>
#include "ncl/nxs_defs.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>
#include "basic_cdf.hpp"
#include "basic_lot.hpp"
#include "phycas_string.hpp"
#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
#	include <boost/python/tuple.hpp>
#	include <boost/python/numeric.hpp>
#	include "thirdparty/num_util/num_util.h"
#endif
#include "xprobdist.hpp"
class XUnderflow{};

namespace phycas
{

struct AdHocDensity 
	{
	virtual ~AdHocDensity() 
		{
		//std::cerr << "AdHocDensity dying..." << std::endl;
		}
    virtual double operator()(double) = 0;
	};

class ProbabilityDistribution : public AdHocDensity
	{
	public:
							ProbabilityDistribution() {lot = &myLot;}
							virtual ~ProbabilityDistribution() {}

		virtual void		SetLot(Lot * other); //@POL seems like this should be a shared pointer
		virtual Lot *		GetLot() {return lot;} // temporary!
		virtual void		ResetLot();
		virtual void		SetSeed(unsigned rnseed);

		virtual	bool		IsDiscrete() const						= 0;
		virtual std::string GetDistributionName() const				= 0;
		virtual std::string GetDistributionDescription() const		= 0;
		virtual double 		GetMean() const							= 0;
		virtual double 		GetVar() const							= 0;
		virtual double 		GetStdDev() const						= 0;
		virtual double		GetCDF(double x) const					= 0; 
		virtual double		Sample() const							= 0; 
		virtual double		GetLnPDF(double x) const				= 0; 
		virtual double		GetRelativeLnPDF(double x) const		= 0; 
		double 				GetRelativeLnPDFArray(double *x, int arrLen) const;
		virtual void 		SetMeanAndVariance(double m, double v)	= 0;

		virtual double		operator()(double x);

		CDF			cdf;
		Lot			myLot;
		Lot *		lot;
	};

typedef boost::shared_ptr<ProbabilityDistribution> ProbDistShPtr;

class MultivariateProbabilityDistribution
	{
	public:
									MultivariateProbabilityDistribution() {lot = &myLot;}
							virtual ~MultivariateProbabilityDistribution() {}

		virtual void				SetLot(Lot *other)									= 0;
		virtual void				ResetLot()													= 0;
		virtual void				SetSeed(unsigned rnseed)									= 0;

		virtual	bool				IsDiscrete() const											= 0;
		virtual std::string 		GetDistributionName() const									= 0;
		virtual std::string 		GetDistributionDescription() const							= 0;
		virtual VecDbl				GetMean() const												= 0;
		virtual VecDbl 				GetVar() const												= 0;
		virtual VecDbl 				GetStdDev() const											= 0;
		virtual VecDbl				Sample() const												= 0; 
		virtual double				ApproxCDF(const VecDbl &x, unsigned nsamples = 10000) const	= 0;
		virtual double				GetLnPDF(const VecDbl &) const								= 0; 
		virtual double				GetRelativeLnPDF(const VecDbl &) const						= 0; 
		virtual void 				SetMeanAndVariance(const VecDbl &m, const VecDbl &v)		= 0;
#		if defined(PYTHON_ONLY)
#			if defined(USING_NUMARRAY)
				virtual void 			AltSetMeanAndVariance(boost::python::numeric::array m, boost::python::numeric::array v) = 0;
				virtual boost::python::numeric::array	GetVarCovarMatrix()						= 0;
#			else
				//virtual void 			AltSetMeanAndVariance(VecDbl m, VecDbl v)				= 0;
				virtual VecDbl			GetVarCovarMatrix()										= 0;
#			endif
#		endif
		virtual unsigned			GetNParams()	const										= 0;

		CDF					cdf;
		Lot					myLot;
		Lot *				lot;
	};

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Encapsulates the discrete bernoulli probability distribution with parameter p, the probability of success.
*/
class BernoulliDistribution : public ProbabilityDistribution
	{
	protected:
		double p;	/* the probability of success on any given trial */

	public:
							BernoulliDistribution();
							BernoulliDistribution(double prob_success);
							~BernoulliDistribution();

		virtual	bool		IsDiscrete() const;
		virtual std::string	GetDistributionName() const;
		virtual std::string	GetDistributionDescription() const;
		virtual double		GetMean() const;
		virtual double		GetVar() const;
		virtual double		GetStdDev() const;
		virtual double		GetCDF(double x) const;
		virtual double		Sample() const; 
		virtual double		GetLnPDF(double x) const; 
		virtual double		GetRelativeLnPDF(double x) const; 
		virtual void		SetMeanAndVariance(double mean, double var);
	};

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Encapsulates the discrete binomial probability distribution with parameter p, the probability of success.
*/
class BinomialDistribution : public BernoulliDistribution
	{
	private:
		double q;
		double lnp;
		double lnq;

	protected:
		double n;

	public:
							BinomialDistribution();
							BinomialDistribution(double sample_size, double prob_success);
							~BinomialDistribution();

		virtual	bool		IsDiscrete() const;
		virtual std::string	GetDistributionName() const;
		virtual std::string	GetDistributionDescription() const;
		virtual double		GetMean() const;
		virtual double		GetVar() const;
		virtual double		GetStdDev() const;
		virtual double		GetCDF(double x) const;
		virtual double		Sample() const; 
		virtual double		GetLnPDF(double x) const; 
		virtual double		GetRelativeLnPDF(double x) const; 
		virtual void		SetMeanAndVariance(double mean, double var);
	};
	
/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Encapsulates the continuous Beta probability distribution with parameters alpha and beta.
*/
class BetaDistribution  : public ProbabilityDistribution
	{
		double alphaParam;
		double betaParam;

	public:
						BetaDistribution() : alphaParam(1.0), betaParam(1.0) {}
						BetaDistribution(double a, double b);
						~BetaDistribution();

		bool			IsDiscrete() const;
		std::string 	GetDistributionName() const;
		std::string 	GetDistributionDescription() const;
		double 			GetMean() const;
		double 			GetVar() const;
		double 			GetStdDev() const;
		double			GetCDF(double x) const;
		double			Sample() const;
		double			GetLnPDF(double x) const;
		double			GetRelativeLnPDF(double x) const;
		void 			SetMeanAndVariance(double m, double v);
	};

/*----------------------------------------------------------------------------------------------------------------------
|	A uniform probability distribution with left bound 0.0 and right bound infinity. This is an improper distribution 
|	(the area under its density curve is infinite).
*/
class ImproperUniformDistribution : public ProbabilityDistribution
	{
	public:
					ImproperUniformDistribution();
					~ImproperUniformDistribution();

		bool		IsDiscrete() const;
		std::string	GetDistributionName() const;
		std::string	GetDistributionDescription() const;
		double		GetMean() const;
		double		GetVar() const;
		double		GetStdDev() const;
		double		GetCDF(double x) const;
		double		Sample() const; 
		double		GetLnPDF(double x) const; 
		double		GetRelativeLnPDF(double x) const; 
		void		SetMeanAndVariance(double mean, double var);
	};
	
/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Encapsulates the continuous uniform probability distribution with left bound a and right bound b.
*/
class UniformDistribution : public ProbabilityDistribution
	{
	public:
					UniformDistribution();
					UniformDistribution(double left_bound, double right_bound);
					~UniformDistribution();

		bool		IsDiscrete() const;
		std::string	GetDistributionName() const;
		std::string	GetDistributionDescription() const;
		double		GetMean() const;
		double		GetVar() const;
		double		GetStdDev() const;
		double		GetCDF(double x) const;
		double		Sample() const; 
		double		GetLnPDF(double x) const; 
		double		GetRelativeLnPDF(double x) const; 
		void		SetMeanAndVariance(double mean, double var);

	protected:
		double a;				/**< the left bound */
		double b;				/**< the right bound */
		double log_density;		/**< the precalculated log of the density function */

	};
	
/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	Encapsulates the gamma probability distribution with shape parameter (alpha) and scale parameter (beta).
*/
class GammaDistribution : public ProbabilityDistribution
	{
	protected:
		double alpha;		/* the shape parameter */
		double beta;		/* the scale parameter */
		double ln_const;	/* the natural logarithm of the constant part of the density function */

	public:
						GammaDistribution();
						GammaDistribution(double shape, double scale);
						~GammaDistribution();

		bool			IsDiscrete() const;
		std::string		GetDistributionName() const;
		std::string		GetDistributionDescription() const;
		double			GetMean() const;
		double			GetVar() const;
		double			GetStdDev() const;
		double			GetCDF(double x) const;
		double			Sample() const; 
		virtual double	GetLnPDF(double x) const; 
		double			GetRelativeLnPDF(double x) const; 
		void			SetMeanAndVariance(double mean, double var);

	protected:
		void			ComputeLnConst();
	};
	
/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	The inverse (or inverted) gamma distribution with parameters alpha and beta is the distribution of 1/X where X is a gamma distributed random variable with shape
|	parameter alpha and scale parameter beta. 
*/
class InverseGammaDistribution : public GammaDistribution
	{
	public:
					InverseGammaDistribution();
					InverseGammaDistribution(double shape, double scale);
					~InverseGammaDistribution();

		bool		IsDiscrete() const;
		std::string	GetDistributionName() const;
		std::string	GetDistributionDescription() const;
		double		GetMean() const;
		double		GetVar() const;
		double		GetStdDev() const;
		double		GetCDF(double x) const;
		double		Sample() const; 
		double		GetLnPDF(double x) const; 
		double		GetRelativeLnPDF(double x) const; 
		void		SetMeanAndVariance(double mean, double var);
	};
	
/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	This is a special case of the gamma distribution in which the shape parameter equals the mean and the scale parameter is 1.0. 
*/
class ExponentialDistribution : public GammaDistribution
	{
	public :
					ExponentialDistribution();
					ExponentialDistribution(double lambda);
					~ExponentialDistribution();

		bool		IsDiscrete() const;
		std::string	GetDistributionName() const;
		std::string	GetDistributionDescription() const;
		void		SetMeanAndVariance(double mean, double var = 0.0);
		double		GetLnPDF(double x) const; 
	};

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------
|	The Normal distribution with two parameters, the mean and standard deviation. 
*/
class NormalDistribution : public ProbabilityDistribution
	{
	protected:
		double		mean;			/**< the mean parameter of the normal distribution */
		double		sd;				/**< the standard deviation parameter of the normal distribution */
		double		ln_const;		/**< the natural logarithm of the constant part of the density function */
		double		pi_const;		/**< precalculated (in constructor) value of pi */
		double		sqrt2_const;	/**< precalculated (in constructor) value of sqrt(2.0) */

	public:
					NormalDistribution();
					NormalDistribution(double mean, double stddev);
					~NormalDistribution();

		bool		IsDiscrete() const;
		std::string	GetDistributionName() const;
		std::string	GetDistributionDescription() const;
		double		GetMean() const;
		double		GetVar() const;
		double		GetStdDev() const;
		double		GetCDF(double x) const;
		double		Sample() const; 
		double		GetLnPDF(double x) const; 
		double		GetRelativeLnPDF(double x) const; 
		void		SetMeanAndVariance(double mean, double var);

	protected:
		void			ComputeLnConst();
	};
	
class DirichletDistribution : public MultivariateProbabilityDistribution
	{
	private:

		std::vector<double>					dirParams;
		std::vector<GammaDistribution>		paramDistributions;
		mutable std::vector<double>			scratchSpace;
	
	public:
											DirichletDistribution();
											DirichletDistribution(const std::vector<double> &params);
											~DirichletDistribution() {}
	
		void								SetLot(Lot * other);
		void								ResetLot();
		void								SetSeed(unsigned rnseed);

		bool								IsDiscrete() const;
		std::string 						GetDistributionName() const;
		std::string 						GetDistributionDescription() const;
		std::string 						GetDescriptionForPython() const;
		VecDbl								GetMean() const;
		VecDbl 								GetVar() const;
		VecDbl 								GetStdDev() const;
		VecDbl								Sample() const;
		double								ApproxCDF(const VecDbl &x, unsigned nsamples = 10000) const;
		double								GetLnPDF(const VecDbl &x) const;
		double								GetRelativeLnPDF(const VecDbl &x) const;
		void 								SetMeanAndVariance(const VecDbl &m, const VecDbl &v);
#		if defined(PYTHON_ONLY)
#			if defined(USING_NUMARRAY)
				void							AltSetMeanAndVariance(boost::python::numeric::array m, boost::python::numeric::array v);
				boost::python::numeric::array	GetVarCovarMatrix();
#			else
				void							AltSetMeanAndVariance(VecDbl m, VecDbl v);
				VecDbl							GetVarCovarMatrix();
#			endif
#		endif

		unsigned							GetNParams() const;
		const GammaDistribution 			&GetDistributionOnParameter(unsigned i) const;
	};

} // namespace phycas

#include "probability_distribution.inl"

#endif
