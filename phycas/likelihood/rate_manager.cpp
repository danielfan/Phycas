#include "phycas/force_include.h"
#include "phycas/likelihood/rate_manager.hpp"
#include "ncl/misc/algorithm_extensions.hpp"


const double kMaxGammaShape = 300.0;
const double kMinGammaShape = 1.0e-6;
const double kMinGammaVar = 1.0/kMaxGammaShape;
const double kMaxGammaVar = 1.0/kMinGammaShape;

const unsigned kGammFuncItmax = 100;
const double kGammFuncEpsilon = 3.0e-7;

#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::log;
	using std::exp;
	using std::fabs;
	using std::sqrt;
	using std::pow;
#endif
		
/*----------------------------------------------------------------------------------------------------------------------
|	Returns z so that Prob{X <= z} = pr where X ~ N(0,1) and (1e-12) < pr < 1-(1e-12). Returns (-9999) if in error.
|	From Odeh, R. E. and J. O. Evans. 1974. The percentage points of the normal distribution. Applied Statistics, 22:
|	96-97 (AS70). Newer methods include:
|	- Wichura, M. J. 1988. Algorithm AS 241: The percentage points of the normal distribution. 37: 477-484.
|	- Beasley, JD and S. G. Springer. 1977. Algorithm AS 111: The percentage points of the normal distribution. 26: 
|		118-121.
|  From Ziheng Yang via John Huelsenbeck via Dave Swofford.
*/
double PhoRateManager::point_normal(
  double pr)	/**< is the cumulative probability */
	{
	double 		a0 = -0.322232431088, a1 = -1.0, a2 = -0.342242088547, a3 = -0.0204231210245,
			a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495,
			b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634,
			y, z, p = pr, p1;

	p1 = (p< 0.5 ? p : 1-p);
	if (p1<1e-20)
		{
		return (-9999);
		}
	y = sqrt (log(1/(p1*p1)));
	z = y + ((((y*a4 + a3)*y + a2)*y + a1)*y + a0) / ((((y*b4 + b3)*y + b2)*y + b1)*y + b0);
	return (p< 0.5 ? -z : z);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	From the Numerical Recipes (Press et al.) function of the same name. 
|	\todo Should switch to using an equivalent not from the Numerical Recipes book.
*/
double PhoRateManager::gammln(const double xx)
	{
	const double cof[6] = {76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	double x = xx - 1.0;
	const double tmp = -x - 5.5 + (x + 6.0)*log(x + 5.5);
	double ser = 1.0;
	for (unsigned int j = 0; j <= 5; ++j) 
		{
		x += 1.0;
		ser += cof[j]/x;
		}
	return tmp + log(2.50662827465*ser);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns z so that Prob{X <= z} = pr where X is Chi2 distributed with df = v. Returns -1 if in error. 
|	0.000002 < pr < 0.999998.
|
|	From RATNEST FORTRAN by Best, D. J. and D. E. Roberts. 1975. The percentage points of the Chi2 distribution. Applied 
|	Statistics 24:385-388. (AS91)
|
|	Converted into C by Ziheng Yang, Oct. 1993.
|
|	From Ziheng Yang via John Huelsenbeck via Dave Swofford.
*/
double PhoRateManager::point_chi2 (
  double pr,	/**< is the cumulative probability */
  double v)		/**< is the degrees of freedom */
	{
	double e = 0.5e-6, aa = 0.6931471805, p = pr, g,
		xx, c, ch, a, q, p1, p2, t,
		x, b, s1, s2, s3, s4, s5, s6;

	if (p < 0.000002 || p > 0.999998 || v <= 0.0)
		return (-1.0);
	g = gammln(v / 2.0);
	xx=v/2.0;
	c=xx-1.0;
	if (v >= -1.24*log(p))
		{
		goto L1;
		}
	ch = pow((p*xx*exp(g + xx*aa)), 1.0/xx);
	if (ch - e< 0)
		return (ch);
	goto L4;
	L1:
		if (v > 0.32)
			goto L3;
		ch = 0.4;
		a = log(1.0-p);
	L2:
		q = ch;
		p1 = 1.0 + ch*(4.67 + ch);
		p2 = ch*(6.73 + ch*(6.66 + ch));
		t = -0.5 + (4.67 + 2.0*ch)/p1 - (6.73 + ch*(13.32 + 3.0*ch))/p2;
		ch -= (1.0 - exp(a + g + 0.5*ch + c*aa)*p2/p1)/t;
		if (fabs(q/ch - 1.0) - 0.01 <= 0.0)
			goto L4;
		else
			goto L2;
	L3:
		x = point_normal (p);
		p1 = 0.222222/v;
		ch = v*pow((x*sqrt(p1) + 1.0 - p1), 3.0);
		if (ch > 2.2*v + 6.0)
			ch = -2.0*(log(1.0 - p) - c*log(0.5*ch) + g);
	L4:
		q = ch;
		p1 = 0.5*ch;
		t = 1.0 - gammq(xx, p1);
		// TO DO: pass in lnGamma since we've already calculated it
		assert(t >= 0.0);
		p2 = p - t;
		t = p2*exp(xx*aa + g + p1 - c*log(ch));
		b = t/ch;
		a = 0.5*t-b*c;
		s1 = (210.0 + a*(140.0 + a*(105.0 + a*(84.0 + a*(70.0 + 60.0*a))))) / 420.0;
		s2 = (420.0 + a*(735.0 + a*(966.0 + a*(1141.0 + 1278.0*a))))/2520.0;
		s3 = (210.0 + a*(462.0 + a*(707.0 + 932.0*a)))/2520.0;
		s4 = (252.0 + a*(672.0 + 1182.0*a) + c*(294.0 + a*(889.0 + 1740.0*a)))/5040.0;
		s5 = (84.0 + 264.0*a + c*(175.0 + 606.0*a))/2520.0;
		s6 = (120.0 + c*(346.0 + 127.0*c))/5040.0;
		ch += t*(1 + 0.5*t*s1 - b*c*(s1 - b*(s2 - b*(s3 - b*(s4 - b*(s5 - b*s6))))));
		if (fabs(q/ch - 1.0) > e)
			goto L4;
		return (ch);
	}


/*----------------------------------------------------------------------------------------------------------------------
|	From the Numerical Recipes (Press et al.) function of the same name. \todo Should switch to non-Numerical-Recipes
|	code.
*/
inline double PhoRateManager::gammq(
  double a,	/**< the a parameter */ 
  double x)	/**< the x parameter */
	{
	
	//if (x < 0.0 || a <= 0.0) {
	//	ofstream tmpf("tmp.txt");
	//	tmpf << " x < 0.0 or a <= 0.0 in PhoRateManager::gammq" << endl;
	//	tmpf.close();
	//}
	double gln;
	assert(x >= 0.0 && a > 0.0);
	if (x < (a + 1.0)) 
		{
		double gamser;
		gser(gamser, a, x, gln);
		return 1.0 - gamser;
		}
	double gammcf;
	gcf(gammcf, a, x, gln);
	return gammcf;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	From the Numerical Recipes (Press et al.) function of the same name. \todo Should switch to non-Numerical-Recipes 
|	code.
*/
void PhoRateManager::gcf(
  double &gammcf, 
  double a, 
  double x, 
  double &gln)
	{
	double gold = 0.0,g,fac = 1.0,b1 = 1.0;
	double b0 = 0.0,anf,ana,an,a1,a0=1.0;

	gln=gammln(a);
	a1=x;
	for (unsigned n = 1; n <= kGammFuncItmax; ++n) 
		{
		an = (double) n;
		ana = an - a;
		a0  = (a1 + a0*ana)*fac;
		b0 = (b1 + b0*ana)*fac;
		anf = an*fac;
		a1 = x*a0 + anf*a1;
		b1 = x*b0 + anf*b1;
		if (a1) 
			{
			fac = 1.0/a1;
			g = b1*fac;
			if (fabs((g-gold)/g) < kGammFuncEpsilon) 
				{
				gammcf = exp(-x + a * log(x) - gln) * g;
				return;
				}
			gold=  g;
			}
		}
	throw XIterMax();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	From the Numerical Recipes (Press et al.) function of the same name. \todo Should switch to non-Numerical-Recipes 
|	code.
*/
void PhoRateManager::gser(
  double &gamser, 
  double a, 
  double x, 
  double &gln)
	{
	gln=gammln(a);
	if (x <= 0.0) 
		{
		assert(x == 0.0);
		gamser = 0.0;
		return;
		} 
	else 
		{
		double ap = a;
		double sum,del;
		del = sum = 1.0/a;
		for (unsigned n = 1; n <= kGammFuncItmax; ++n) 
			{
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*kGammFuncEpsilon) 
				{
				gamser = sum * exp(-x + a*log(x) - gln);
				return;
				}
			}
		throw XIterMax();
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates the rate for each category of variable sites based on `pInvar' and the gamma distribution. Using mean or
|	median rates (depending on the value of `useMeanRate'). The rates are calculated from the gamma and then scaled by
|	`pinvBrCorrection' if `pInvar' > 0. Each of the rate categories is of equal size (in terms of the probability mass).
|	If the gamma shape parameter is above `kMaxGammaShape' (or the variance is below `kMinGammaVar'), all of the rates 
|	are set to `pinvBrCorrection'.
*/
void PhoRateManager::CalculateRates()
	{
	assert(UsingPInvar() || pinvBrCorrection == 1.0);

	if (UsingPInvar())
		pinvBrCorrection = 1.0 / (1.0 - pInvar->GetValue()) ;

#	if defined (USING_VARIANCE_FOR_RATEHET)
		const bool noGammaRateHet = (nVarRateCats == 1 || gammaParam->GetValue() < kMinGammaVar );
#	else
		const bool noGammaRateHet = (nVarRateCats == 1 || gammaParam->GetValue() > kMaxGammaShape);
#	endif

	if (noGammaRateHet)
		{
		const double theRate = pinvBrCorrection * (meanRate ? meanRate->GetValue() : 1.0);
		replace_all(rates, rates + nVarRateCats, theRate);
		return;
		}
		
	const double gamVariable = gammaParam->GetValue();

#	if defined (USING_VARIANCE_FOR_RATEHET)
		const double gamShape = 1.0 / (gamVariable > kMaxGammaVar ? kMaxGammaVar : gamVariable);
#	else
		const double gamShape = (gamVariable < kMinGammaShape ? kMinGammaShape : gammaParam->val);
#	endif
	
	const double twoShape = 2.0 * gamShape;

	double a, b, x;
	const double N = (double) nVarRateCats;
	

//	if (!useMeanRate) //@ is this a bug? we were trying to use the meanRate (as opposed to the median) did we comment out
						// the wrong section of code?  I think that useMeanRate is not the correct variable name and we are OK
//		{
		const double invK = 1.0 / N;
		b = 0.0;
		x = 1.0;  /* (in case user sets nRateCategs=1) */
		for (unsigned i = 0; i < nVarRateCats - 1; ++i)
			{
			a = b;
			b = point_chi2(invK * (i+1), twoShape) / twoShape;
			x = gammq(gamShape + 1.0, b * gamShape);
			double y = (a > 0.0 ? gammq(gamShape + 1.0, a * gamShape) : 1.0);
			rates[i] = (y - x);	// mult by N that was done here now done below in MultiplyRatesBy
			}
		rates[nVarRateCats - 1] = x; // mult by N that was done here now done below in MultiplyRatesBy
		const double rateMult =  N * pinvBrCorrection * (meanRate ?  meanRate->GetValue() : 1.0 );
		MultiplyRatesBy(rateMult);
	
//		}
//	else CODE FOR RATES AS THE MEDIAN OF THE RATE CATEGORY
//		{
//		invK = 0.5 / N;
//		for (i = 0; i < nVarRateCats; ++i)
//			{
//			rates[i] = point_chi2(invK * (i*2+1), twoShape) / twoShape;
//			}
//		ScaleRatesToMeanOfOne();
//		const double rateMult =  (meanRate != NULL ? pinvBrCorrection * meanRate->GetValue() : pinvBrCorrection);
//		MultiplyRatesBy(rateMult);
//		}
//		

	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of parameters that are owned and being used by the rate manager.
*/
unsigned PhoRateManager::GetNumParameters() const
	{
	unsigned np =0;
	if (gammaParam && nVarRateCats > 1)
		++np;
	if (UsingPInvar())
		++np;
	if (meanRate)
		++np;
	return np;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Deletes all allocated fields (right now only rates)
*/
void PhoRateManager::FlushRates()
	{
	delete [] rates;
	rates = NULL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates all pointer fields (right now only rates)
*/
void PhoRateManager::ReallocateRates()
	{
	FlushRates();
	rates = new double[nVarRateCats];	//deleted in Flush
	}

#if 0
/*----------------------------------------------------------------------------------------------------------------------
|	Initializes all fields to default values. Need to call InitRateManager() if this constructor is used.
*/
PhoRateManager::PhoRateManager()
  : rates(NULL)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes all fields (allocating memory for rates)
*/
void PhoRateManager::InitRateManager(
  ConstParamShPtr pinv, 		/**< is the initial value of `pInvar' parameter (with the correct initial value; this is ignored if the next field is false) */
  ConstParamShPtr gamParam, 	/**< is the shape parameter of the gamma distribution */
  ConstParamShPtr mr,			/**< is the mean rate */
  unsigned nGamCats)			/**< is the number of VARIABLE rate categories */
	{
	meanRate		= mr;
	gammaParam		= gamParam;
	pInvar			= pinv;
	nVarRateCats	= nGamCats;
	assert(nVarRateCats != 0 && nVarRateCats != UINT_MAX);
	assert(meanRate == ConstParamShPtr() || meanRate->GetValue() > 0.0);

	if (pinv)
		{
		assert(pinv->GetValue() < 1.0 && pinv->GetValue() >= 0.0);
		pinvBrCorrection = 1.0 / (1.0 - pInvar->GetValue()) ; // number that must be multiplied by the branch length so that the exp. # changes = br len
		}
	else
		pinvBrCorrection = 1.0;

	ReallocateRates();

	if (nVarRateCats > 1)
		{
		assert(gamParam);
#		if defined (USING_VARIANCE_FOR_RATEHET)
			assert(gamParam->GetValue() < kMaxGammaVar);
#		else
			assert(gamParam->GetValue() < kMaxGammaShape);
#		endif
		CalculateRates();
		}
	else if (meanRate)
		rates[0] = meanRate->GetValue() * pinvBrCorrection;
	else
		rates[0] = pinvBrCorrection;
	}
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes all fields (allocating memory for rates)
*/
PhoRateManager::PhoRateManager(
  ConstParamShPtr pinv, 		/**< is the initial value of `pInvar' parameter (with the correct initial value; this is ignored if the next field is false) */
  ConstParamShPtr gamParam, 	/**< is the shape parameter of the gamma distribution */
  ConstParamShPtr mr,			/**< is the mean rate */
  unsigned nGamCats)			/**< is the number of VARIABLE rate categories */
	:pInvar(pinv), 
	gammaParam(gamParam),
	meanRate(mr), 
	rates(NULL), 
	nVarRateCats(nGamCats)
	{
	assert(nVarRateCats != 0 && nVarRateCats != UINT_MAX);
	assert(meanRate == ConstParamShPtr() || meanRate->GetValue() > 0.0);
	if (pinv)
		{
		assert(pinv->GetValue() < 1.0 && pinv->GetValue() >= 0.0);
		pinvBrCorrection = 1.0 / (1.0 - pInvar->GetValue()) ; // number that must be multiplied by the branch length so that the exp. # changes = br len
		}
	else
		pinvBrCorrection = 1.0;
	ReallocateRates();
	if (nVarRateCats > 1)
		{
		assert(gamParam);
#		if defined (USING_VARIANCE_FOR_RATEHET)
			assert(gamParam->GetValue() < kMaxGammaVar);
#		else
			assert(gamParam->GetValue() < kMaxGammaShape);
#		endif
		CalculateRates();
		}
	else if (meanRate)
		rates[0] = meanRate->GetValue() * pinvBrCorrection;
	else
		rates[0] = pinvBrCorrection;
	}

#if 0
#	include "lot.h"

	int PhoRateManager::itmax = 100;
	double PhoRateManager::epsilon = 3.0e-7;

	/*------------------------------------------------------------------------------------------------------------------
	|	Calls Flush to delete all allocated fields (right now only rates)
	*/
	PhoRateManager::~PhoRateManager()
		{
		Flush();
		}

	/*------------------------------------------------------------------------------------------------------------------
	| 	this function is used in simulations and will fills in an array with rates for len random sites
	*/
	void PhoRateManager::FillWithRandomRate(
	  double *dest,
	  charIndex len,
	  Lot *rnd)
	  	{
	  	double multiplier = 1.0;
	  	if (usingPInvar)
	  		multiplier = pinvBrCorrection;
	  	if (meanRate)
	  		multiplier *= meanRate->val;
	  	if (usingPInvar)
	  		{
	  		double pV = pInvar->val;
	  		if (nVarRateCats > 1)
	  			{
#				if defined (USING_VARIANCE_FOR_RATEHET)
	  				double invGS = gammaParam->val;
	  				double gs = 1.0/invGS;
#				else
		  			double gs = gammaParam->val;
		  			double invGS = 1.0/gs;
#				endif
	  			for (charIndex i = 0;i < len; ++i) 
	  				{
	  				if ( rnd->Uniform() < pV)
	  					dest[i] = 0.0;
	  				else
	  					dest[i] = multiplier * rnd->Gamma(gs, invGS);
	  				}
	  			}
	  		else
	  			{
	  			for (charIndex i = 0;i < len; ++i) 
	  				{
	  				if ( rnd->Uniform() < pV)
	  					dest[i] = 0.0;
	  				else
	  					dest[i] = multiplier;
	  				}
	  			}
	  		}
	  	else
	  		{
	  		if (nVarRateCats > 1)
	  			{
#				if defined (USING_VARIANCE_FOR_RATEHET)
		  			double invGS = gammaParam->val;
		  			double gs = 1.0/invGS;
#				else
		  			double gs = gammaParam->val;
		  			double invGS = 1.0/gs;
#				endif
	  			for (charIndex i = 0;i < len; ++i) 
	  				dest[i] = multiplier * rnd->Gamma(gs, invGS);
	  			}
	  		else
	  			{
	  			for (charIndex i = 0;i < len; ++i) 
	  				dest[i] = multiplier;
	  			}
	  		}
	  }
	  
	/*------------------------------------------------------------------------------------------------------------------
	| 	this function is used in simulations and will return a rate for a random site, based on the current settings
	*/
	double PhoRateManager::DrawRandomRate(
	 Lot *rnd)
	  	{
	  	double multiplier = 1.0;
	  	if (usingPInvar)
	  		{
	  		if ( rnd->Uniform() < pInvar->val)
	  			return 0.0;
	  		multiplier = pinvBrCorrection;
	  		}
	  	if (meanRate)
	  		multiplier *= meanRate->val;
	  	if (nVarRateCats > 1)
#			if defined (USING_VARIANCE_FOR_RATEHET)
  				return multiplier*(rnd->Gamma(1.0/gammaParam->val, gammaParam->val));
#			else
  				return multiplier*(rnd->Gamma(gammaParam->val, 1.0/gammaParam->val));
#			endif
	  	return multiplier;
	  	}

#	if	defined (DEBUGGING_LIKELIHOOD_CALC)
		void PhoRateManager::AppendToPAUPLSetCmd(
		  string &s)
			{
			if (nVarRateCats > 1)
				{
				s << " rates = gamma ncat = ";
				s << nVarRateCats;
				s << " shape = ";
#				if defined (USING_VARIANCE_FOR_RATEHET)
					s << (1.0 / gammaParam->val);
#				else
					s << gammaParam->val;
#				endif
				}
			else
				s << " rates = equal";
			if (usingPInvar)
				{
				s << " pinv = ";
				s << pInvar->val;
				}
			else
				s << " pinv = 0.0";
			if (meanRate)
				{
				s << " [!The mean rate is ";
				s << meanRate->val;
				s << " so the branch lengths should be scaled appropriately]";
				}
			}
		void PhoRateManager::AppendToPhorestModelCmd(
		  string &s)
		  	{
			if (nVarRateCats > 1)
				{
				if (nVarRateCats != 4)
					{
					assert(0);	//ncat option not written yet
					s << " rates = gamma ncat = ";
					s << nVarRateCats;
					s << " shape = ";
					}
				else
#					if defined (USING_VARIANCE_FOR_RATEHET)
						s << " rates = gamma RateVariance = ";
#					else
						s << " rates = gamma shape = ";
#					endif
				s << gammaParam->val;
				}
			else
				s << " rates = equal";
			if (usingPInvar)
				{
				s << " pinv = ";
				s << pInvar->val;
				}
			else
				s << " pinv = 0.0";
			if (meanRate)
				{
				assert(0);	// mean rate option not written yet
				s << " [!The mean rate is ";
				s << meanRate->val;
				s << " so the branch lengths should be scaled appropriately]";
				}
		  	}
#	endif // defined (DEBUGGING_LIKELIHOOD_CALC)
#endif //0

