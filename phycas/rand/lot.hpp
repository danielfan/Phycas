#ifndef PHO_LOT_H
#define PHO_LOT_H

#include <cassert>
#include <boost/shared_ptr.hpp>

#define SLIGHTLY_BETTER_BINOMIAL 0
#define NEW_SEED_WAY 1
#define OLD_SEED_WAY !NEW_SEED_WAY
class CDF;
unsigned FindIndexFromFreq(double prob, const double *freqs, unsigned nCateg);
const unsigned kMaxLotUInt = 2147483647; // 2^31 - 1 

/*----------------------------------------------------------------------------------------------------------------------
|	This class was called Lot because the noun lot is defined as "an object used in deciding something by chance" 
|	according to The New Merriam-Webster Dictionary.
*/
class Lot
	{
	public:
		STATIC_DATA_FUNC double LotLnGamma(double x);
		
					Lot();
					Lot(unsigned);
					Lot(const Lot &);
					~Lot();

		// State Accessors
		unsigned	 GetSeed() const
							{ 
							return currSeed;
							}	/**< Returns current seed (which changes after each draw) */
		unsigned 		GetInitSeed() const 
							{
							return lastSeedSetting;
							}	/**< Returns seed used to initialize generator */
		// Manipulations
		void 			UseClockToSeed();
		void 			UseClockToSeedThenSpin(unsigned spin = 100);
		void 			SetSeed(unsigned s)
							{
							assert(s > 0 && s <= kMaxLotUInt);
							currSeed = lastSeedSetting = s;
							}
		void 			Spin(unsigned spin = 100);

		unsigned long	GetRandBits(unsigned nbits);

		// Discrete Samplers
		unsigned 		Binomial(double p, const unsigned nSamples);
		unsigned char 	Bernoulli(double p)
							{
							return (unsigned char)((Uniform() < p) ? 1 : 0);
							}
#		if defined(USE_UNTESTED_RNG_CODE) && (USE_UNTESTED_RNG_CODE)
		unsigned 		RandUniRaw(); // from dls
#		endif
		unsigned 		SampleUInt(unsigned);
		long 			SampleLong(long);
		
							/*--------------------------------------------------------------------------------------------------------------------------
							|	This function takes an array of double (freqs) of length nCateg, and selects an index [0 - nCateg) such that 
							|	Pr(x) = freq[x].
							|	It assumes that the sum of the freqs array is 1.0, and all rounding error is "given" to the last category (nCateg-1) 
							*/
		unsigned		SingleMulitinomial(const double *freqs, unsigned nCateg)
							{
							return FindIndexFromFreq(Uniform(), freqs, nCateg);
							}
							/*--------------------------------------------------------------------------------------------------------------------------
							|	Temp. slow sampling from a multinomial by repeatinf nSamples calls to SingleMulitinomial.
							|	the seleceted elements of arr have their values incremented (The elements of arr are NOT initialized to 0)
							*/
		void			SampleMultinomial(const double *freqs, unsigned nCateg, unsigned nSamples, unsigned *arr)
							{
							for (unsigned i = 0; i < nSamples; ++i)
								++arr[SingleMulitinomial(freqs, nCateg)];
							}
		//				
		// Continuous Samplers
		double 			Beta( double a, double b )
							{ 
							double x = Gamma(a, 1.0); 
							double y = Gamma(b, 1.0); 
							return ( x / (x+y) ); 
							}
		double 			ChiSquare(double df);
		double 			Exponential( double lambda )
							{
							return ( lambda * Expon() ); 
							}
		double 			Gamma( double a, double b);
		double 			Normal(double mean, double var);
		unsigned		Poisson(double mean);
		double 			Uniform();
#		if defined(USE_UNTESTED_RNG_CODE) && (USE_UNTESTED_RNG_CODE)
		double 			RandUni(); // from dls
#		endif					
					
		unsigned	 	BruteForceBinomial(const double p, const unsigned nSamples);
					
	private:    	
		double 			Expon();
		double 			Unbounded(double alpha);
		double 			Bounded(double alpha);
		
		unsigned 	lastSeedSetting;
		unsigned	currSeed;
		CDF		*cdfConverter;
	};

typedef boost::shared_ptr<Lot> LotShPtr;

#endif
