#ifndef PHO_CDF_H
#define PHO_CDF_H

#include <cassert>
#include <cstdio>
#include <cmath>
#include "ncl/nxs_exception.hpp"

#if defined(__MWERKS__)
#	if __ide_target("Console-Mac")
#		define NEW_CDF_WAY 1
#	else
#		define NEW_CDF_WAY 0
#	endif
#else
#	define NEW_CDF_WAY 0
#endif
#	define OLD_CDF_WAY !NEW_CDF_WAY

#if (NEW_CDF_WAY)
namespace CDF
	{
	class X_CDF_IllegalParam : public NxsException
		{
		public:
			X_CDF_IllegalParam(const std::string &s): NxsException(s){}
		};
	class X_CDF_CouldNotCompute : public NxsException
		{
		public:
			X_CDF_CouldNotCompute(const std::string &s): NxsException(s){}
		};
		
	DblPair		CumBeta(DblPair xyPair, double alpha, double beta);
	DblPair		CumChiSquare(double x, double df);
	DblPair 	CumF(const double f, const double dfn, const double dfd);
	DblPair		CumGamma(double x, double alpha, double beta);
	DblPair 	CumNormal(const double x, const double mean, const double sd); //throw X_CDF
	DblPair 	CumNonCentralChiSq(const double x, const double df, const double pnonc);
	DblPair 	CumNonCentralF(const double f, const double dfn, const double dfd, const double pnonc);
	DblPair 	CumNonCentralT(const double t, const double df, const double pnonc);
	DblPair 	CumStdNormal(const double x);
	DblPair 	CumT(const double t,const double df);

	DblPair 	InverseCumBeta(const DblPair prob, const double a, const double b);
	double 		InverseCumNormal(const DblPair prob, const double mean, const double sd);
	//double 		InverseCumNonCentralChiSq(
	double 		InverseCumStdNormal(const DblPair prob);

	
	double		FindAlphaOfBeta(const DblPair prob, const DblPair xy, const double b);
	double		FindBetaOfBeta(const DblPair prob, const DblPair xy, const double a);
	double 		FindMeanOfNormal(const DblPair prob, const double x, const double sd);
	double 		FindStdDevOfNormal(const DblPair prob, const double x, const double mean);

	DblPair		IncompleteBeta(const double a, const double b, const DblPair xy);
	DblPair		IncompleteGamma(const double x, const double aShape);
	
	double 	alnrel(const double);
	double 	algdiv(const double,const double);
	double 	BetaLn(const double, const double);
	double 	CDFALnGamma(const double);
	double 	erf1(const double x);
	double 	erfc1(const int, const double);
	double 	esum(const int ,const double);
	double 	EvalPoly(const double *a,const unsigned len, const double x);
	double 	CDFLnGamma(const double);
	double 	gamln1(double);
	double 	gam1(const double);
	double 	Gam1SumHelper(const double a, const double b);
	DblPair 	grat1(const double a, const double z,const double r,const double eps);
	double 	InitInvCumT(const DblPair pq, const double df);
	double	InvNormNewtonStartVal(const double p);
	double 	psi(const double xx);
	double 	rlog1(double);
	double 	Xgamm(const double);


	/*
	void 	cdfbet(int*,double*,double*,double*,double*,double*,double*, int*,double*);
	void 	cdfbin(int*,double*,double*,double*,double*,double*,double*, int*,double*);
	void 	cdfchi(int*,double*,double*,double*,double*,int*,double*);
	void 	cdfchn(int*,double*,double*,double*,double*,double*,int*,double*);
	void 	cdff(int*,double*,double*,double*,double*,double*,int*,double*);
	void 	cdffnc(int*,double*,double*,double*,double*,double*,double*, int*s,double*);
	void 	cdfgam(int*,double*,double*,double*,double*,double*,int*,double*);
	void 	cdfnbn(int*,double*,double*,double*,double*,double*,double*, int*,double*);
	void 	cdfnor(int*,double*,double*,double*,double*,double*,int*,double*);
	void 	cdfpoi(int*,double*,double*,double*,double*,int*,double*);
	void 	cdft(int*,double*,double*,double*,double*,int*,double*);
	void 	cdftnc(int*,double*,double*,double*,double*,double*,int*,double*);
	*/
	void 	cumbet(double*,double*,double*,double*,double*,double*);
	void 	cumbin(double*,double*,double*,double*,double*,double*);
	void 	cumchi(double*,double*,double*,double*);
	void 	cumgam(double*,double*,double*,double*); //IncompleteGamma
	void 	cumnbn(double*,double*,double*,double*,double*,double*);
	void 	cumpoi(double*,double*,double*,double*);
	void 	dinvr(int*,double*,double*,unsigned long*,unsigned long*);
	void 	dstinv(double*,double*,double*,double*,double*,double*, double*);
	void 	dzror(int*,double*,double*,double*,double *, unsigned long*,unsigned long*);
	void 	dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl);
	double 	exparg(int*);
	void 	gaminv(double*,double*,double*,double*,double*,int*);
	void 	gratio(double*,double*,double*,double*,int*);
	double 	gsumln(double a,double b);
	double 	rcomp(double*,double*);
	double 	rexp(double*);
	double 	rlog(double*);
	double 	spmpar(int*);
	double 	fifdint(double);
	double 	fifdsign(double,double);
	long 	fifmod(long,long);
	void 	ftnstop(char*);
	double			SampleGamma(double p, double alpha, double beta);
	double			SampleChiSquare(double p, double df);
	double			SampleNormal(double p, double mean, double var);
	double 			SampleBinomial(double cumProb, double p, double nTrials);

	// These were declared in the original distribution
	//
	void E0000(int,int*,double*,double*,unsigned long*,
					  unsigned long*,double*,double*,double*,
					  double*,double*,double*,double*);
	void E0001(int,int*,double*,double*,double*,double*,
					  unsigned long*,unsigned long*,double*,double*,
					  double*,double*);

	/*--------------------------------------------------------------------------------------------------------------------------
	|	returns DBL_EPSILON
	|	Converted from: dcdflib::spmpar(1) 
	*/
	inline double MachinePrec() 
		{
		return DBL_EPSILON;
	    }
	    
	/*--------------------------------------------------------------------------------------------------------------------------
	|	returns DBL_MIN
	|	Converted from: dcdflib::spmpar(2) 
	*/
	inline double SmallestMagnitude() //spmpar(2)
		{
		return DBL_MIN;
		}
		
	/*--------------------------------------------------------------------------------------------------------------------------
	|	returns DBL_MAX
	|	Converted from: dcdflib::spmpar(3) 
	*/
	inline double LargestMagnitude() //spmpar(3)
		{
		return  DBL_MAX;
		}
	} // namespace CDF

	/*--------------------------------------------------------------------------------------------------------------------------
	|	Calculates the cdf to xyPair of the incomplete beta distribution with parameters a and b.  
	|	LATEX 	\int_0^x \frac{t^{a-1}(1-t)^{b-1}}{B(a, b)}
	|	Redirected To: 	IncompleteBeta
	|	Converted from:	dcdflib::cumbet()
	|	Equivalent to:	dcdflib::cdfbet(which = 1)
	*/
	inline DblPair CDF::CumBeta(
	  DblPair xyPair, 	/*first element is the upper limit of integration */
	  double  alpha, 	/*First parameter of the beta distribution.*/
	  double  beta)		/*Second parameter of the beta distribution.*/
		{
		return IncompleteBeta(alpha, beta, xyPair);
		}

	/*--------------------------------------------------------------------------------------------------------------------------
	|	Calculates of the Chi-Square distribution with df degrees of freedom, up to x 
	|	Redirected To: 	IncompleteGamma
	|	Converted from:	dcdflib::cumchi() 
	|	Equivalent to:	dcdflib::cdfchi(which = 1) 
	*/
	inline DblPair CDF::CumChiSquare(
	  double x, /*Upper limit of integration*/
	  double df) /*Degrees of freedom */
		{
		return IncompleteGamma(x/2.0, df/2.0);
		}


#else

/*----------------------------------------------------------------------------------------------------------------------
|	CDF is a wrapper around version 1.1 of DCDFLIB. All functions from the original distribution are now given either
|	protected or private status (with private status reserved for functions declared STATIC_DATA_FUNC (i.e., file scope). These
|	oriignal functions still bear the hallmarks of the original Fortran code from which they were translated (e.g. 
|	many goto statements and labels), and no attempt has been made to either document these (beyond supplying the 
|	original comments) or clean them up to be more C++-like. The functionality needed in Phorest has been added using
|	public member functions that call the private or protected DCDFLIB member functions. The DCDFLIB functions all take
|	pointers to variables, and a set of commonly used variables has been added to the class and declared private because
|	they are only used internally to provide a path for communication between the public member functions to the DCDFLIB
|	member functions. The CDF class also contains the function ipmpar, which occupied a separate file in the DCDFLIB
|	distribution. The original documentation for the DCDFLIB library has been appended at the bottom of the cdf.cpp 
|	file. The DCDFLIB library was written/compiled by Barry W. Brown, James Lovato and Kathy Russell, Department of 
|	Biomathematics, Box 237, The University of Texas, M. D. Anderson Cancer Center, Houston, Texas 77030. The 
|	source code used here was downloaded Nov. 17, 2002, from ftp://odin.mdacc.tmc.edu/pub/source/dcdflib.c-1.1.tar.gz
|	but is also available on StatLib at http://lib.stat.cmu.edu/general/Utexas/
|
|	Note: MTH added mutable to all of the data members and const to all of the functions.  This (at first glance, silly)
|	modification was made because the internal variables all appear to be reset with each public function call (so mutable
|	is fine), but we ProbabilityDistribution is derived from CDF so we need the CDF::functions that are called by const
|	ProbabilityDistribution::functions  o be declared const.
*/
class CDF
	{
	public:

				CDF();
				~CDF(){}
			double			CumBeta(double x, double alpha, double beta) const;
			double			CumChiSquare(double x, double df);
			double			CumGamma(double x, double alpha, double beta) const;
			
			STATIC_DATA_FUNC double	InverseCumStdNormal(const DblPair prob);
			STATIC_DATA_FUNC DblPair 	CumStdNormal(double );
			
			STATIC_DATA_FUNC double 	alnrel(double*);
			STATIC_DATA_FUNC double 	algdiv(double*,double*);
			STATIC_DATA_FUNC double 	alngam(double*);
			STATIC_DATA_FUNC double 	apser(double*,double*,double*,double*);
			STATIC_DATA_FUNC double 	bcorr(double*,double*);
			STATIC_DATA_FUNC double 	betaln(double*,double*);
			STATIC_DATA_FUNC double 	basym(double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	bgrat(double*,double*,double*,double*,double*,double*,int*i);
			STATIC_DATA_FUNC double 	bfrac(double*,double*,double*,double*,double*,double*);
			STATIC_DATA_FUNC double 	bpser(double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	bratio(double*,double*,double*,double*,double*,double*,int*);
			STATIC_DATA_FUNC double	brcmp1(int*,double*,double*,double*,double*);
			STATIC_DATA_FUNC double 	brcomp(double*,double*,double*,double*);
			STATIC_DATA_FUNC double 	bup(double*,double*,double*,double*,int*,double*);
			STATIC_DATA_FUNC void 	cumchn(double*,double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	cumf(double*,double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	cumfnc(double*,double*,double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	cumnor(double*,double*,double*);
			STATIC_DATA_FUNC void 	cumt(double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	cumtnc(double*,double*,double*,double*,double*);
			STATIC_DATA_FUNC double 	dinvnr(double *p,double *q);
			STATIC_DATA_FUNC double 	devlpl(double [],int*,double*);
			STATIC_DATA_FUNC double 	dt1(double*,double*,double*);
			STATIC_DATA_FUNC double 	erf1(const double);
			STATIC_DATA_FUNC double 	esum(int*,double*);
			STATIC_DATA_FUNC double 	fpser(double*,double*,double*,double*);
			STATIC_DATA_FUNC double 	gam1(double*);
			STATIC_DATA_FUNC double 	CDFLnGamma(double*);
			STATIC_DATA_FUNC double 	gamln1(double*);
			STATIC_DATA_FUNC void 	grat1(double*,double*,double*,double*,double*,double*);
			STATIC_DATA_FUNC double 	fifdmax1(double,double);
			STATIC_DATA_FUNC double 	fifdmin1(double,double);
			STATIC_DATA_FUNC long 	fifidint(double);
			STATIC_DATA_FUNC double 	psi(double*);
			STATIC_DATA_FUNC double 	rlog1(double*);
			STATIC_DATA_FUNC double 	stvaln(double*);
			STATIC_DATA_FUNC double 	Xgamm(double*);

			STATIC_DATA_FUNC void 	cdfbet(int*,double*,double*,double*,double*,double*,double*, int*,double*);
			STATIC_DATA_FUNC void 	cdfbin(int*,double*,double*,double*,double*,double*,double*, int*,double*);
			STATIC_DATA_FUNC void 	cdfchi(int*,double*,double*,double*,double*,int*,double*);
			STATIC_DATA_FUNC void 	cdfchn(int*,double*,double*,double*,double*,double*,int*,double*);
			STATIC_DATA_FUNC void 	cdff(int*,double*,double*,double*,double*,double*,int*,double*);
			STATIC_DATA_FUNC void 	cdffnc(int*,double*,double*,double*,double*,double*,double*, int*s,double*);
			STATIC_DATA_FUNC void 	cdfgam(int*,double*,double*,double*,double*,double*,int*,double*);
			STATIC_DATA_FUNC void 	cdfnbn(int*,double*,double*,double*,double*,double*,double*, int*,double*);
			STATIC_DATA_FUNC void 	cdfnor(int*,double*,double*,double*,double*,double*,int*,double*);
			STATIC_DATA_FUNC void 	cdfpoi(int*,double*,double*,double*,double*,int*,double*);
			STATIC_DATA_FUNC void 	cdft(int*,double*,double*,double*,double*,int*,double*);
			STATIC_DATA_FUNC void 	cdftnc(int*,double*,double*,double*,double*,double*,int*,double*);
			STATIC_DATA_FUNC void 	cumbet(double*,double*,double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	cumbin(double*,double*,double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	cumchi(double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	cumgam(double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	cumnbn(double*,double*,double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	cumpoi(double*,double*,double*,double*);
			STATIC_DATA_FUNC void 	dinvr(int*,double*,double*,unsigned long*,unsigned long*);
			STATIC_DATA_FUNC void 	dstinv(double*,double*,double*,double*,double*,double*, double*);
			STATIC_DATA_FUNC void 	dzror(int*,double*,double*,double*,double *, unsigned long*,unsigned long*);
			STATIC_DATA_FUNC void 	dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl);
			STATIC_DATA_FUNC double 	erfc1(int*,double*);
			STATIC_DATA_FUNC double 	exparg(int*);
			STATIC_DATA_FUNC void 	gaminv(double*,double*,double*,double*,double*,int*);
			STATIC_DATA_FUNC void 	gratio(double*,double*,double*,double*,int*);
			STATIC_DATA_FUNC double 	gsumln(double*,double*);
			STATIC_DATA_FUNC double 	rcomp(double*,double*);
			STATIC_DATA_FUNC double 	rexp(double*);
			STATIC_DATA_FUNC double 	rlog(double*);
			STATIC_DATA_FUNC double 	spmpar(int*);
			STATIC_DATA_FUNC double 	fifdint(double);
			STATIC_DATA_FUNC double 	fifdsign(double,double);
			STATIC_DATA_FUNC long 	fifmod(long,long);
			STATIC_DATA_FUNC void 	ftnstop(char*);

			// The following was declared extern in the original distribution
			//
			STATIC_DATA_FUNC int		ipmpar(int *i);
		
		
		double			SampleGamma(double p, double alpha, double beta) const;
		double			SampleChiSquare(double p, double df) const;
		double			SampleNormal(double p, double mean, double var) const;
		double 			SampleBinomial(double cumProb, double p, double nTrials);

	protected:

	private:
		// These were declared STATIC_DATA_FUNC in the original distribution
		//
		STATIC_DATA_FUNC void E0000(int,int*,double*,double*,unsigned long*,
						  unsigned long*,double*,double*,double*,
						  double*,double*,double*,double*);
		STATIC_DATA_FUNC void E0001(int,int*,double*,double*,double*,double*,
						  unsigned long*,unsigned long*,double*,double*,
						  double*,double*);

	private:
			struct CDFDataStruct
				{
				int		which;
				double	p;
				double	q;
				double	x;
				double	y;
				double	alpha;
				double	beta;
				double	df;
				int		status;
				double	bound;
				double	mean;
				double	var;
				};
	};
/*----------------------------------------------------------------------------------------------------------------------
|
*/
inline CDF::CDF()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative gamma distribution function for value x, shape parameter alpha and scale parameter beta. This is 
|	the integral of the gamma probability density function from 0 up to x. The gamma density as defined here is
|>
|	       x^{alpha - 1} e^{-x/beta}
|	f(x) = -------------------------
|	       beta^{alpha} Gamma(alpha)
|>
|	Note: the CDF library for which this CDF class is a wrapper defines the scale parameter as the inverse of beta!
|	Only the private member functions use this definition however, the CDF interface uses the density as defined above.
*/
inline double CDF::CumGamma(
  double x,		/* the integral of the probability density function from 0 to x will be returned */
  double alpha,	/* the shape parameter of the gamma distribution */
  double beta) const	/* the scale parameter of the gamma distribution */
	{
	double X = x;
	int status = 0;
	int which = 1;	// compute p given x, shape and scale
	double P = 0.0;
	double Q = 0.0;
	double shape = alpha;
	double scale = 1.0 / beta; // CDF library uses inverse of normal definition of scale parameter
	double bound = 0.0;
	cdfgam(&which, &P, &Q, &X, &shape, &scale, &status, &bound);
	assert(status == 0);
	return P;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of x such that the integral of the gamma probability density function from 0 up to x is equal to
|	the supplied parameter p (with shape parameter alpha and scale parameter beta provided). Useful for simulating draws
|	from a gamma distribution. The gamma density as defined here is
|>
|	       x^{alpha - 1} e^{-x/beta}
|	f(x) = -------------------------
|	       beta^{alpha} Gamma(alpha)
|>
|	Note: the CDF library for which this CDF class is a wrapper defines the scale parameter as the inverse of beta!
|	Only the private member functions use this definition however, the CDF interface uses the density as defined above.
*/
inline double CDF::SampleGamma(
  double p,		/* the integral of the gamma probability density function from 0 up to x, where x represents the value computed and returned */
  double alpha,	/* the shape parameter */
  double beta) const	/* the scale parameter */
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
	assert(status == 0); //@POL should do slice sampling here
	return X;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of x such that the integral of the chi-square probability density function from 0 up to x is equal to
|	the supplied parameter p (with df degrees of freedom). Useful for simulating draws from a chi-square distribution.
*/
inline double CDF::SampleChiSquare(
  double p,		/* the integral of the chi-square probability density function from 0 up to x, where x represents the value computed and returned */
  double df) const	/* the degrees of freedom */
	{
	CDFDataStruct dcdflib;
	dcdflib.which	= 2;	// calculate x given p, q and df
	dcdflib.p		= p;
	dcdflib.q		= 1.0 - p;
	dcdflib.x		= 0.0;
	dcdflib.df		= df;
	dcdflib.status	= 0;
	dcdflib.bound	= 0.0;

	cdfchi(&dcdflib.which, &dcdflib.p, &dcdflib.q, &dcdflib.x, &dcdflib.df, &dcdflib.status, &dcdflib.bound);
	assert(dcdflib.status == 0);
	return dcdflib.x;
	}

inline double	CDF::SampleNormal(
  double p, 
  double mean, 
  double var) const
  	{
  	CDFDataStruct dcdflib;
	dcdflib.which	= 2;	// calculate x given p, q and df
	dcdflib.p		= p;	// p is the cumulative prob
	dcdflib.q		= 1.0 - p;
	dcdflib.x		= 0.0;	// result
	dcdflib.mean	= mean;
	dcdflib.var		= std::sqrt(var);
	dcdflib.status	= 0;	//result code
	dcdflib.bound	= 0.0;	// second part of result code
	if (var == 0.0)
		return mean;
	cdfnor(&dcdflib.which, &dcdflib.p, &dcdflib.q, &dcdflib.x, &dcdflib.mean, &dcdflib.var, &dcdflib.status, &dcdflib.bound);
	assert(dcdflib.status == 0);
	return dcdflib.x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns cumulative beta distribution function for value x and parameters alpha and beta. This is the integral of
|	the beta probability density function from 0 up to x.
*/
inline double CDF::CumBeta(
  double x,		/* the integral of the probability density function from 0 to x will be returned */
  double alpha,	/* the first parameter */
  double beta)	const /* the second parameter */
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
|	Returns cumulative chi-square distribution function for value x and degrees of freedom df. This is the integral of
|	the chi-square probability density function from 0 up to x.
*/
inline double CDF::CumChiSquare(
  double x,		/* the integral of the probability density function from 0 to x will be returned */
  double df)	/* the degrees of freedom */
	{
  	CDFDataStruct dcdflib;
	dcdflib.which	= 1;	// calculate p and q given x and df
	dcdflib.p		= 0.0;
	dcdflib.q		= 0.0;
	dcdflib.x		= x;
	dcdflib.df		= df;
	dcdflib.status	= 0;
	dcdflib.bound	= 0.0;

	cdfchi(&dcdflib.which, &dcdflib.p, &dcdflib.q, &dcdflib.x, &dcdflib.df, &dcdflib.status, &dcdflib.bound);
	assert(dcdflib.status == 0);
	return dcdflib.p;
	}


#endif //NEW_CDF_WAY

#if (NEW_CDF_WAY)
#endif
#endif
