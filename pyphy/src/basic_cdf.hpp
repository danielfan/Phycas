#if ! defined(BASIC_CDF_HPP)
#define BASIC_CDF_HPP

#include <cstdio>
#include <cmath>

extern "C"
{
#include "thirdparty/dcdflib/src/cdflib.h"
}

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	CDF is a wrapper around version 1.1 of DCDFLIB. The DCDFLIB library was written/compiled by Barry W. Brown, James 
|	Lovato and Kathy Russell, Department of Biomathematics, Box 237, The University of Texas, M. D. Anderson Cancer 
|	Center, Houston, Texas 77030. The source code used here was downloaded Nov. 17, 2002, from 
|	ftp://odin.mdacc.tmc.edu/pub/source/dcdflib.c-1.1.tar.gz but is also available on StatLib at 
|	http://lib.stat.cmu.edu/general/Utexas/.
|
|	This file contains the minimum functionality needed for Phycas and does not attempt to wrap all the functions in the
|	original DCDFLIB library. The original source code and documentation for DCDFLIB is in the directory 
|	phypy/prob_dist/dcdflib.
*/
class CDF
	{
	public:
		double 						LnGamma(double x) const;
		double						CumBeta(double x, double alpha, double beta) const;
		double						CumGamma(double x, double alpha, double beta) const;
		double						SampleGamma(double p, double alpha, double beta) const;
		double						CumNorm(double x, double mu, double sigma) const;
		double						SampleNorm(double p, double alpha, double beta) const;
	};

} // namespace phycas

#include "basic_cdf.inl"

#endif
