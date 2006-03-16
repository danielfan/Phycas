#define POLPY_SPIN

#include "phycas/force_include.h"
#include "phycas/rand/lot.hpp"
#include "phycas/rand/cdf.hpp"
#include <ctime>
using std::max;
using std::min;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::exp;
	using std::log;
	using std::pow;
	using std::sqrt;
#endif
#define USE_FATTENED_BINOMIAL 1
const double kSingleOnlyBinomialProbCutoff = 1e-8; 	/// cutoff for when a probability of a "hit" in binomial sampling is so low that 
												/// we can ignore double hits
const double kSmallBinomialProbCutoff = kSingleOnlyBinomialProbCutoff/2.0;	///cutoff for using the SmallProbBinomial instead of the ModProbBinomial routine
const unsigned MASKSIGNBIT = 0x80000000;


const double kExpOne = exp(1.0);

#if !(defined(HAVE_64BIT_INTS) || defined(THINK_C) || (defined(__MWERKS__) && macintosh) || defined(OS_MAC))
#	define NEED_LONGMUL
#endif

double LnChoose(const unsigned n, const unsigned k);


/// returns the lowest n from [0,nCateg) such that 
///		$prob \leq \sum_{i=0}^n freq[i]$
///	Used in sampling from a multinomial (if prob is a U(0,1) variate and freqs is
///	an array of the probabilities of each category.
unsigned FindIndexFromFreq(double prob, const double * freqs, unsigned nCateg)
	{
	assert(nCateg > 0);
	--nCateg;
	for(unsigned index = 0; index < nCateg; ++index)
		{
		if (prob <= *freqs)
			return index;
		prob -= *freqs++;
		}
	return nCateg;
	}
	
#if defined (USE_UNTESTED_RNG_CODE) &&  (USE_UNTESTED_RNG_CODE)	
#	if defined(MAC_PHOREST)
#		define TARGET_CPU_PPC 1
#	endif
#	if defined(NEED_LONGMUL)
		template<typename T>
		class Int64Bit
			{
			public:
				T	hiLong;
				T	loLong;
			};
		
		template <typename T>
		void	LongMul(T, T, Int64Bit<T> *);
		/*----------------------------------------------------------------------------------------------------------------------
		|
		|	LongMul
		|
		|	Multiplies the two 32-bit integers 'ab' and 'cd' and stores the result in 'mn' (high-order 32 bits) and 'op' (low
		|	order 32 bits); needed when it is not possible to store a 64-bit result directly.
		|
		|	Note:  This routine is a substitute for the Macintosh toolbox utility of the same name, however it has been
		|	written specifically for random number generation and assumes that its input arguments are nonnegative, even
		|	though they are stored as signed quantities.  Thus, the scope of valid inputs is [0, 2^31 - 1]); no checking is
		|	performed.
		|
		|	The algorithm used is to break each 32-bit integer into 16-bit components and perform the multiplication as:
		|
		|
		|		  ab		(a = high 16 bits, b = low 16 bits of 1st arg)
		|		x cd		(c = high 16 bits, d = low 16 bits of 2nd arg)
		|		----
		|		  ef		(d*b)
		|		 gh			(d*a)
		|		 ij			(c*b)
		|		kl			(c*a)
		|		----
		|		mnop		('mn' returned as result->hiLong; 'op' as result->loLong)
		*/
		template<typename T>
		void LongMul(T ab, T cd, Int64Bit<T> *result)
			{
			const unsigned a = (unsigned)(ab >> 16);
			const unsigned b = (unsigned)(ab & 0x0000FFFF);
			const unsigned c = (unsigned)(cd >> 16);
			const unsigned d = (unsigned)(cd & 0x0000FFFF);

			const unsigned ef = d * b;
			const unsigned gh = d * a;
			const unsigned ij = c * b;
			const unsigned kl = c * a;

			const unsigned o = (ef >> 16) + (gh & 0x0000FFFF) + (ij & 0x0000FFFF);

			result->loLong = (T)((o << 16) | (ef & 0x0000FFFF));
			result->hiLong = (T)((o >> 16 /* carry */) + (gh >> 16) + (ij >> 16) + kl);
			}

#	endif //defined(NEED_LONGMUL)

/*----------------------------------------------------------------------------------------------------------------------
|
|	RandUni
|
|	Returns a pseudorandom number on U[0,1] via the multiplicative congruential method:
|
|		x(i) = A * x(i-1) mod (2^31 - 1).
|
|	It is recommended that the multiplier A be set to 397204094.  This is the smallest multiplier that performed well
|	in Fishman and Moore's (1982, JASA 77:129-136) comprehensive evaluation.  (Smaller multipliers are much faster,
|	and may be "good enough" depending on the purpose.)  Do not change the multiplier unless you know what you're
|	doing!
|
|	Pass a pointer to an integer variable containing the previous number in the sequence (or the starting seed, if this
|	is the first call).  The seed must be greater than 0 and less than 2^31 - 1.  It is up to the caller to make sure
|	that the seed falls within these bounds.
|
|	The algorithm used here is based on that of Payne, Rabung, and Bogyo (1969; Comm. ACM 12:85-86).
|
|	This routine is portable, so long as the default integer size is at least 32 bits.  (The code can easily be
|	modified, and simplified as well, for machines with long's of 64 bits or longer.)  Note that the "overflow"
|	in the comments below simply means that the value exceeded 2^31 - 1 (but it will never exceed 2^32 - 1).
|
|	Dave Swofford 9/21/87
|	Modified 1/8/93
*/
const unsigned kRandUniM = 2147483647;				// modulus = 2^31 - 1  used in RandUni and RandUniRaw
double Lot::RandUni()
	{
	/* Get next random number in sequence */
	unsigned u = RandUniRaw();
	const double x = (1.0 / (kRandUniM-2)) * (u - 1);
	
	/* Force value into [0,1] due to possibility of roundoff error */
	if (x < 0.0)
		return 0.0;
	else if (x > 1.0)
		return 1.0;
	return x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The guts of RandUni.  If speed is critical, floating-point operations can be avoided by calling this
|	function directly.
|
|	Returns a pseudorandom integer from U[1,2^31 - 2].  The caller should use this value as the next seed.
*/
unsigned Lot::RandUniRaw()
	{
	const unsigned kRandUniA = 397204094;					/* multiplier */
	const unsigned MASK31BITS = 0x7FFFFFFF;

#	if !defined(HAVE_64BIT_INTS)
#		if defined (TARGET_CPU_PPC) && TARGET_CPU_PPC
			wide Ax;
			(void)WideMultiply(currSeed, kRandUniA, &Ax);	/* kRandUniA * x */
			currSeed = Ax.lo;								/* move to register variables */
			unsigned y = (unsigned)Ax.hi;					/*	for speed				  */
#		else //TARGET_CPU_PPC
			Int64Bit<unsigned> Ax;
			LongMul(currSeed, kRandUniA, &Ax);		/* kRandUniA * x */
			currSeed = (unsigned)Ax.loLong;			/* move to register variables */
			unsigned y = (unsigned)Ax.hiLong;		/*	for speed				  */
#		endif	 //TARGET_CPU_PPC
		y = (y << 1) | (currSeed >> 31);			/* isolate high-order 31 bits */
		currSeed &= MASK31BITS;						/* isolate low-order 31 bits */
		currSeed += y;								/* x'(i + 1) unless overflows */
		if (currSeed & MASKSIGNBIT) 				/* overflow check */
			currSeed -= kRandUniM;					/* deal with overflow */

#	else //defined(HAVE_64BIT_INTS)

		/* IMPORTANT: Use this 64-bit version only if multiplication of two 32-bit operands yields
		              a result with 64-bit precision, which is not necessarily the case!
		*/
#		if defined(__ALPHA) || defined(__alpha) || defined(__alpha__)

			/* This is faster on Dec Alpha, which seems to have a fast modulo operation */
			currSeed = (unsigned)(((uint64_t)kRandUniA * currSeed) % (uint64_t)kRandUniM);

#		else	/* This version is probably faster on most processors */

			const uint64_t long kMASK32BITS 0x00000000FFFFFFFFL
			uint64_t w = (uint64_t)kRandUniA * currSeed;
			currSeed = (unsigned)(w & kMASK32BITS);
			y = (unsigned)(w >> 32);
			y = (y << 1) | (currSeed >> 31);		/* isolate high-order 31 bits */
			currSeed &= MASK31BITS;					/* isolate low-order 31 bits */
			currSeed += y;							/* x'(i + 1) unless overflows */
			if (currSeed & MASKSIGNBIT) 			/* overflow check */
				currSeed -= kRandUniM;				/* deal with overflow */
#		endif
#	endif //defined(HAVE_64BIT_INTS)

	return currSeed;
	}

#endif //(USE_UNTESTED_RNG_CODE)	

double Lot::Gamma( double a, double b) 
	{
	// a and b are the parameters of the gamma distribution
	//	 mean = a * b
	//	variance = mean * b
		assert( a > 0.0 );
		if( a == 1.0 )
			return ( b * Expon() );
		else if( a < 1.0 )
			return ( b * Unbounded(a) );
		else
			return( b * Bounded(a) );
	}
		

unsigned  Lot::BruteForceBinomial(const double p, const unsigned nSamples)
	{
	assert(p > DBL_EPSILON && p <=.5);
	unsigned currCount = 0;
	for (unsigned i = 0; i < nSamples; ++i)
  		{
  		if (Uniform() < p)
  			++currCount;
  		}
  	return currCount;
	} 
	
unsigned HelperSampleFromCategs(const double *nHitProb, const unsigned nReps, const double probLowerPart, const unsigned nCatsInLowerPart, const unsigned nCatsInUpperPart, Lot &rnd);
template<int>
unsigned FattenedBinomial(const double p, const unsigned nSamples, Lot &rnd);

inline unsigned HelperSampleFromCategs(const double *nHitProb, const unsigned nReps, const double probLowerPart, const unsigned  nCatsInLowerPart, const unsigned nCatsInUpperPart, Lot &rnd)
	{
	const double *probUpperHalf = nHitProb + nCatsInLowerPart;
	const unsigned validSecondHalf = nCatsInUpperPart - 1;
	unsigned currCount = 0;
	for (unsigned k = 1; k <= nReps; ++k)	//starting at 1 so we can decrement if we need to skip a rep.
		{
		double prob = rnd.Uniform();
		if (prob < probLowerPart)
			{
			const unsigned nHits = FindIndexFromFreq(prob, nHitProb, nCatsInLowerPart);
			currCount += nHits;
			}
		else
			{
			const unsigned nHits = FindIndexFromFreq(prob - probLowerPart, probUpperHalf, nCatsInUpperPart);
			if (nHits <= validSecondHalf)// last categ is for rounding error (sign that we need to resimulate)
				currCount += nCatsInLowerPart + nHits;
			else
				--k; // round off error , don't count this rep this should be very rare
			}
		}
	return currCount;
	}

template<>
inline unsigned FattenedBinomial<20>(const double p, const unsigned nSamples, Lot &rnd)
        {
        assert(p <=.5 && p > 0.0);
        double nHitProb[22];
        const double kChooseN[] = {1, 20, 190, 1140, 4845, 15504, 38760, 77520, 125970, 167960, 184756, 167960, 125970, 77520, 38760, 15504, 4845, 1140, 190, 20, 1};
        double currPMult = p;
        nHitProb[0] = 1.0;
        for (unsigned i = 1; i < 21; ++i)
                {
                nHitProb[i] = currPMult*kChooseN[i];
                currPMult *= p; 
                }
        const double omp = 1.0 - p;
        double currOMPMult = omp;
        for (unsigned j = 19;; --j)
                {
                nHitProb[j] *= currOMPMult;
                if (j == 0)
                        break;
                currOMPMult *= omp;
                }
        nHitProb[21]= 0.0;      //last category is for roundoff error
        unsigned maxCateg = 0;  //find max using the fact that the distribution is unimodal (max will usually be a small index because p <.5)
        double sumBelowMax = 0.0;
        for (; (nHitProb[maxCateg] < nHitProb[maxCateg+1] && maxCateg < 21); ++maxCateg)
                {
                sumBelowMax += nHitProb[maxCateg];
                }
        const unsigned nReps = nSamples/20;
        unsigned currCount = 0;
        if (maxCateg == 0 || maxCateg == 1)
                {
                register const double zeroCatProb = (maxCateg == 0) ? nHitProb[0] : sumBelowMax;
                for (unsigned k = 1; k <= nReps; ++k)   //starting at 1 so we can decrement if we need to skip a rep.
                        {
                        double prob = rnd.Uniform();
                        if (prob > zeroCatProb)
                                {
                                const unsigned nHits = 1 + FindIndexFromFreq(prob - zeroCatProb, (nHitProb +1), 11);
                                if (nHits <= 20)// last categ is for rounding error (sign that we need to resimulate)
                                        currCount += nHits;
                                else
                                        --k; // round off error , don't count this rep this should be very rare
                                }
                        }
                }
        else
                currCount =  HelperSampleFromCategs(nHitProb, nReps, sumBelowMax, maxCateg, 22 - maxCateg, rnd);
        currCount += rnd.BruteForceBinomial(p, nSamples % 20);
        return currCount;
        }

template<>
inline unsigned FattenedBinomial<10>(const double p, const unsigned nSamples, Lot &rnd)
	{
	assert(p <=.5 && p > 0.0);
	double nHitProb[12];
	const double kChooseN[] = {1, 10,  45, 120, 210, 252, 210, 120, 45, 10, 1};
	double currPMult = p;
	nHitProb[0] = 1.0;
	for (unsigned i = 1; i < 11; ++i)
		{
		nHitProb[i] = currPMult*kChooseN[i];
		currPMult *= p; 
		}
	const double omp = 1.0 - p;
	double currOMPMult = omp;
	for (unsigned j = 9;; --j)
		{
		nHitProb[j] *= currOMPMult;
		if (j == 0)
			break;
		currOMPMult *= omp;
		}
	nHitProb[11]= 0.0;	//last category is for roundoff error	
	unsigned maxCateg = 0;	//find max using the fact that the distribution is unimodal (max will usually be a small index because p <.5)
	double sumBelowMax = 0.0;
	for (; (nHitProb[maxCateg] < nHitProb[maxCateg+1]); ++maxCateg)
		{
		assert(maxCateg < 11);
		sumBelowMax += nHitProb[maxCateg];
		}
		
	const unsigned nReps = nSamples/10;
	unsigned currCount = 0;
	if (maxCateg == 0 || maxCateg == 1)
		{
		register const double zeroCatProb = (maxCateg == 0) ? nHitProb[0] : sumBelowMax;
		for (unsigned k = 1; k <= nReps; ++k)	//starting at 1 so we can decrement if we need to skip a rep.
			{
			double prob = rnd.Uniform();
			if (prob > zeroCatProb)
				{
				const unsigned nHits = 1 + FindIndexFromFreq(prob - zeroCatProb, (nHitProb +1), 11);
				if (nHits <= 10)// last categ is for rounding error (sign that we need to resimulate)
					currCount += nHits;
				else
					--k; // round off error , don't count this rep this should be very rare
				}
			}
		}
	else
		currCount =  HelperSampleFromCategs(nHitProb, nReps, sumBelowMax, maxCateg, 12 - maxCateg, rnd);
	currCount += rnd.BruteForceBinomial(p, nSamples % 10);
	return currCount;
	}

	
const unsigned kNSimulBinomialCategs = 10; // don't change kNSimulBinomialCategs without changing kChooseRatio
#if ! defined(USE_FATTENED_BINOMIAL)
const double kChooseRatio[] = 
	{
		0.0,  			// unused 
		10.0, 			// (kNSimulBinomialCategs choose 1) / (kNSimulBinomialCategs choose 0)
		4.5,			// (kNSimulBinomialCategs choose 2) / (kNSimulBinomialCategs choose 1)
		120.0/45.0,		// (kNSimulBinomialCategs choose 3) / (kNSimulBinomialCategs choose 2)
		210.0/120.0,	// (kNSimulBinomialCategs choose 4) / (kNSimulBinomialCategs choose 3)
		252.0/210.0,	// (kNSimulBinomialCategs choose 5) / (kNSimulBinomialCategs choose 4)
		210.0/252.0,	// (kNSimulBinomialCategs choose 6) / (kNSimulBinomialCategs choose 5)
		120.0/210.0,	// (kNSimulBinomialCategs choose 7) / (kNSimulBinomialCategs choose 6)
		45.0/120.0,		// (kNSimulBinomialCategs choose 8) / (kNSimulBinomialCategs choose 7)
		10.0/45.0,		// (kNSimulBinomialCategs choose 9) / (kNSimulBinomialCategs choose 8)
		0.1				// (kNSimulBinomialCategs choose 10) / (kNSimulBinomialCategs choose 9)
	};
#endif
unsigned Lot::Binomial(double p, const unsigned nSamples)
	{
	if (p > .5)
  		return nSamples - Binomial(1.0-p, nSamples);
	if (p <= 0.0)
  		return 0;
	if (nSamples <= 3*kNSimulBinomialCategs)  //@@@temp bench
		return BruteForceBinomial(p, nSamples);
	else
		{
# 		if	USE_FATTENED_BINOMIAL
		return FattenedBinomial<10>(p, nSamples, *this);
#		else
		unsigned currCount = 0;
		const double omp = 1.0 - p;
		const double ratio = p/omp;

		double probNHits[kNSimulBinomialCategs + 2]; // last categ is for rounding error (sign that we need to resimulate)
		probNHits[0] = pow(omp, kNSimulBinomialCategs);
		for (unsigned i = 1; i <= kNSimulBinomialCategs; ++i)
			{
			//const double chooseRatio = exp(LnChoose(kNSimulBinomialCategs, i) - LnChoose(kNSimulBinomialCategs, i - 1));
			probNHits[i] = probNHits[i-1]*ratio*kChooseRatio[i];
			}
		probNHits[kNSimulBinomialCategs] = 0.0;
		unsigned nReps = nSamples/kNSimulBinomialCategs;
		for (unsigned j = 0; j < nReps;)
			{
			const unsigned nHits = SingleMulitinomial(probNHits, kNSimulBinomialCategs + 2);
			if (nHits <= kNSimulBinomialCategs)// last categ is for rounding error (sign that we need to resimulate)
				{
				currCount += nHits;
				++j;
				}
			}
		currCount += BinomialImpl(p, nSamples % kNSimulBinomialCategs);
		return currCount;
#		endif
		}

	}

	
#if 0
unsigned Lot::SmallProbBinomial(
  const double p, 
  const unsigned nSamples)
	{
	assert(p < kSmallBinomialProbCutoff);
	if (p <= 0.0)
		return 0;
	if (p >= kSingleOnlyBinomialProbCutoff)
		return BinomialImpl(p, nSamples);
	unsigned count = 0;
	unsigned remaining = nSamples;
	
	if (p < kSingleOnlyBinomialProbCutoff/2.0)
		{
		if (p*(double)nSamples < kSingleOnlyBinomialProbCutoff)
			return (Uniform() < p*(double)nSamples ? 1U : 0U);
		unsigned multToUse = (unsigned) (kSingleOnlyBinomialProbCutoff/p);
		if (multToUse < 2)
			multToUse = 2;
		remaining = nSamples % multToUse;
		const unsigned newNSamp = nSamples/multToUse;
		const double newProb = p * (double)multToUse;
		count = BinomialImpl(newProb, newNSamp);
		}
	count += BinomialImpl(p, remaining);
	return count;
	}
	
#endif		

Lot::Lot() : lastSeedSetting(1U), currSeed(1U)
	{
	//std::cerr << "In Lot default constructor (" << this << ")" << std::endl;
	cdfConverter = new CDF();
#if defined(POLPY_SPIN)
	UseClockToSeedThenSpin();
#else
	UseClockToSeed();
#endif
	}

Lot::Lot(unsigned rndSeed) : lastSeedSetting(1U), currSeed(1U)
	{
	//std::cerr << "In Lot(unsigned) constructor (" << this << ")" << std::endl;
	cdfConverter = new CDF();
	if (rndSeed == 0 || rndSeed == LONG_MAX)
#if defined(POLPY_SPIN)
		UseClockToSeedThenSpin();
#else
		UseClockToSeed();
#endif
	else
		SetSeed(rndSeed);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copy constructor needed if Lot objects are to be used in std::vector: for example, without copy constructor, this
|	fails:
|	
|	std::vector<Lot> lots;
|	lots.push_back(Lot());
|	
|	Fails (I think) because cdfConverter of temporary object created by Lot() is deleted before push_back call is 
|	finished.
*/
Lot::Lot(const Lot &other) : lastSeedSetting(other.lastSeedSetting), currSeed(other.currSeed)
	{
	//std::cerr << "In Lot copy constructor (" << this << ")" << std::endl;
	cdfConverter = new CDF();
	}

Lot::~Lot() 
	{
	//std::cerr << "In Lot destructor (" << this << ")" << std::endl;
	delete cdfConverter;
	}

double Lot::Normal(
  double mean, 
  double var)
  	{
  	return cdfConverter->SampleNormal(Uniform(), mean, var);
  	}
  
double Lot::ChiSquare( 
  double df)
	{
	return cdfConverter->SampleChiSquare(Uniform(), df);
	}

double LnChoose(const unsigned n, const unsigned k)
	{
	assert(n >= k);
	if (k == 0 || k == n)
		return 0.0;
	const unsigned largerK = max(k, n-k);
	const unsigned smallerK = min(k, n-k);
	double sum = 0.0;
	for (unsigned i = n; i > largerK; --i)
		sum += log((double)i);
	for (unsigned j = smallerK; j > 1; --j)
		sum -= log((double)j);
	return sum;
	}
	
#if 0
double ChooseRatio(unsigned kOne, unsigned kTwo, unsigned n)
	{
	if (kOne == kTwo)
		return 1.0;
	if (kOne < kTwo)
		return 1/ChooseRatio(kOne, kTwo, n);
	double ratio = 1.0;
	unsigned lowerN = kTwo;
	unsigned higherN = n - kTwo;
	for (; lowerN <= kOne; ++lowerN, --higherN)
		{
		
		ratio *= 
		}
	}
#endif
	
void Lot::Spin(unsigned spin)
	{
	for(unsigned k = 0; k < spin; ++k )
		Uniform();
	}

void Lot::UseClockToSeed()
	{
#	if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::time;
	using std::time_t;
#	endif
	time_t timer;
	currSeed = (unsigned) time(&timer);
	lastSeedSetting = currSeed;
	}

void Lot::UseClockToSeedThenSpin(unsigned spin)
	{
#	if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::time;
	using std::time_t;
#	endif
	time_t timer;
	currSeed = (unsigned) time(&timer);
#	if (OLD_SEED_WAY)
		lastSeedSetting = currSeed; // why are we resetting  lastSeedSetting here ?
        Spin( spin );
#	else
		Spin( spin );
		lastSeedSetting = currSeed;
#	endif
	}

long Lot::SampleLong( long max )
{
	long return_val = max;

	while( return_val == max )
		return_val = (long)( (double)max * Uniform() );

	return return_val;
}

//	returns an index from [0-max) in which all indices are equiprobable
//	INFINITE LOOP if max is zero!!!
unsigned Lot::SampleUInt( unsigned max )
	{
	assert(max > 0);
	unsigned return_val = max;

	while( return_val == max ) 
		{
		double r = Uniform();
		return_val = (unsigned)( (double)max * r );
		}

	return return_val;
	}

unsigned long Lot::GetRandBits(unsigned nbits)
	{
	assert(nbits > 0);
	assert(nbits < 32);
	double term1 = log(2.0)*(double)nbits;
	double u = Uniform();
	assert(u > 0.0);
	double term2 = log(u);
	double term3 = exp(term1 + term2);
	return (unsigned long)floor(term3);
	}

// Uniform pseudorandom number generator
// Provided by J. Monahan, Statistics Dept., N.C. State University
//   From Schrage, ACM Trans. Math. Software 5:132-138 (1979)
// Translated to C by Paul O. Lewis, Dec. 10, 1992
//
double Lot::Uniform()
	{
	//std::cerr << "Lot::Uniform()..." << GetInitSeed() << std::endl;

	const unsigned a = 16807U;
	const unsigned b15 = 32768U;
	const unsigned b16 = 65536U;
	const unsigned p = 2147483647U;
	const unsigned xhi = currSeed / b16;
	const unsigned xalo = (currSeed - xhi * b16) * a;
	const unsigned leftlo = xalo / b16;
	const unsigned fhi = xhi * a + leftlo;
	const unsigned k = fhi / b15;
	currSeed = (((xalo - leftlo * b16) - p) + (fhi - k * b15) * b16) + k;
	if (currSeed & MASKSIGNBIT) 
		currSeed = unsigned (((int)currSeed) + p);
	return currSeed * 4.6566128575e-10;
	}


#if 0
	// No longer used - now use the expon function from Ripley's book
	// which doesn't involve taking logarithms
	double rng::Exponential(double lambda)
	{
		double x = 0.0;

		while( x <= 0.0 || x > 1.0 )
			x = 1.0 - Uniform();
		x = -log(x) / lambda;

		return x;
	}
#endif

// Translated from Fortran to C++ by Paul O. Lewis
// September 16, 1999.  From Appendix B.6 (algorithm 3.7)
// p. 230 in: Ripley, Brian D. 1987. Stochastic Simulation.
// John Wiley & Sons, Inc.
//
double Lot::Expon()
{
	double u, u0, ustar;
	double exprv = 0.0;

	bool done = false;
	double A = 0.0;
	
	while( !done )
	{
		u = u0 = Uniform();
	
		for(;;)
		{
			ustar = Uniform();
			if( u < ustar ) {
				exprv = A + u0;
				done = true;
				break;		
			}
			
			u = Uniform();
			if( u >= ustar ) 
				break;
		}
		
		A += 1.0;
	}
	
	return exprv;
}

// Translated from Fortran to C++ by Paul O. Lewis
// September 16, 1999.  From Appendix B.7 (algorithm 3.19)
// pp. 230-231 in: Ripley, Brian D. 1987. Stochastic Simulation.
// John Wiley & Sons, Inc.
//
double Lot::Unbounded( double alpha )
	{
	double p, x;
	double b = (alpha + kExpOne)/kExpOne;
	
	for(;;)
	{
		p = b * Uniform();
		if( p > 1.0 ) {
			double tmp = ( b - p ) / alpha;
			x = -log(tmp);
			tmp = pow( x, alpha - 1.0 );
			double u = Uniform();
			if( tmp >= u )
				break;
		}
		else {
			x = pow( p, 1.0 / alpha );
			double u = Uniform();
			if( x <= -log(u) )
				break;
		}
	}	
	
	return x;
	}

// Translated from Fortran to C++ by Paul O. Lewis
// September 16, 1999.  From Appendix B.7 (algorithm 3.20)
// p. 231 in: Ripley, Brian D. 1987. Stochastic Simulation.
// John Wiley & Sons, Inc.
//
double Lot::Bounded(double alpha)
	{
	STATIC_DATA double c1 = 0.0;
	STATIC_DATA double prev_alpha = 0.0;
	STATIC_DATA double c2 = 0.0;
	STATIC_DATA double c3 = 0.0;
	STATIC_DATA double c4 = 0.0;
	STATIC_DATA double c5 = 0.0;
	double aa, u1, u2, w;
	if (alpha != prev_alpha ) 
		{
		c1 = alpha - 1.0;
		aa = 1.0 / c1;
		c2 = aa * (alpha - 1.0/(6.0 * alpha));
		c3 = 2.0 * aa;
		c4 = c3 + 2.0;
		if (alpha > 2.5)
			c5 = 1.0 / sqrt(alpha);
		}

	for(;;)
		{
		u1 = Uniform();
		u2 = Uniform();

		if( alpha > 2.5 ) {
		u1 = u2 + c5 * ( 1.0 - 1.86 * u1 );
		if( u1 <= 0.0 || u1 >= 1.0 )
		continue;
		}

		w = c2 * u2 / u1;

		if( c3*u1 + w + 1.0/w < c4 )
			break;		

		if( c3*log(u1) - log(w) + w < 1.0 )
			break;
		}

	prev_alpha = alpha;
	return c1 * w;
	}


double Lot::LotLnGamma(double x)
	{
    // ====================================================================== 
    // NIST Guide to Available Math Software. 
    // Source for module GAMLN from package CMLIB. 
    // Retrieved from TIBER on Wed Apr 29 17:30:20 1998. 
    // ====================================================================== 
    //     WRITTEN BY D. E. AMOS, SEPTEMBER, 1977. 
    //
    //     REFERENCES 
    //         SAND-77-1518 
    //
    //         COMPUTER APPROXIMATIONS BY J.F.HART, ET.AL., SIAM SERIES IN 
    //         APPLIED MATHEMATICS, WILEY, 1968, P.135-136. 
    //
    //         NBS HANDBOOK OF MATHEMATICAL FUNCTIONS, AMS 55, BY 
    //         M. ABRAMOWITZ AND I.A. STEGUN, DECEMBER. 1955, P.257. 
    //
    //     ABSTRACT 
    //         GAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR 
    //         X.GT.0. A RATIONAL CHEBYSHEV APPROXIMATION IS USED ON 
    //         8.LT.X.LT.1000., THE ASYMPTOTIC EXPANSION FOR X.GE.1000. AND 
    //         A RATIONAL CHEBYSHEV APPROXIMATION ON 2.LT.X.LT.3. FOR 
    //         0.LT.X.LT.8. AND X NON-INTEGRAL, FORWARD OR BACKWARD 
    //         RECURSION FILLS IN THE INTERVALS  0.LT.X.LT.2 AND 
    //         3.LT.X.LT.8. FOR X=1.,2.,...,100., GAMLN IS SET TO 
    //         NATURAL LOGS OF FACTORIALS. 
    //
    //     DESCRIPTION OF ARGUMENTS 
    //
    //         INPUT 
    //           X      - X.GT.0 
    //
    //         OUTPUT 
    //           GAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT X 
    //
    //     ERROR CONDITIONS 
    //         IMPROPER INPUT ARGUMENT - A FATAL ERROR 

    STATIC_DATA double xlim1 = (double)8.;
    STATIC_DATA double xlim2 = (double)1e3;
    STATIC_DATA double rtwpil = (double).918938533204673;
    STATIC_DATA double p[5] = { (double)7.66345188e-4,(double)-5.9409561052e-4,(double)
	    7.936431104845e-4,(double)-.00277777775657725,(double)
	    .0833333333333169 };
    STATIC_DATA double q[2] = { (double)-.00277777777777778,(double).0833333333333333 }
	    ;
    STATIC_DATA double pcoe[9] = { (double).00297378664481017,(double)
	    .0092381945590276,(double).109311595671044,(double).398067131020357,
	    (double)2.15994312846059,(double)6.33806799938727,(double)
	    20.7824725317921,(double)36.0367725300248,(double)62.0038380071273 }
	    ;
    STATIC_DATA double qcoe[4] = { (double)1.,(double)-8.90601665949746,(double)
	    9.82252110471399,(double)62.003838007127 };
    STATIC_DATA double gln[100] = { (double)0.,(double)0.,(double).693147180559945,(
	    double)1.79175946922806,(double)3.17805383034795,(double)
	    4.78749174278205,(double)6.5792512120101,(double)8.52516136106541,(
	    double)10.6046029027453,(double)12.8018274800815,(double)
	    15.1044125730755,(double)17.5023078458739,(double)19.9872144956619,(
	    double)22.5521638531234,(double)25.1912211827387,(double)
	    27.8992713838409,(double)30.6718601060807,(double)33.5050734501369,(
	    double)36.3954452080331,(double)39.3398841871995,(double)
	    42.3356164607535,(double)45.3801388984769,(double)48.4711813518352,(
	    double)51.6066755677644,(double)54.7847293981123,(double)
	    58.0036052229805,(double)61.261701761002,(double)64.5575386270063,(
	    double)67.8897431371815,(double)71.257038967168,(double)
	    74.6582363488302,(double)78.0922235533153,(double)81.557959456115,(
	    double)85.0544670175815,(double)88.5808275421977,(double)
	    92.1361756036871,(double)95.7196945421432,(double)99.3306124547874,(
	    double)102.968198614514,(double)106.631760260643,(double)
	    110.320639714757,(double)114.034211781462,(double)117.771881399745,(
	    double)121.533081515439,(double)125.317271149357,(double)
	    129.123933639127,(double)132.952575035616,(double)136.802722637326,(
	    double)140.673923648234,(double)144.565743946345,(double)
	    148.477766951773,(double)152.409592584497,(double)156.360836303079,(
	    double)160.331128216631,(double)164.320112263195,(double)
	    168.327445448428,(double)172.352797139163,(double)176.395848406997,(
	    double)180.456291417544,(double)184.533828861449,(double)
	    188.628173423672,(double)192.739047287845,(double)196.86618167289,(
	    double)201.009316399282,(double)205.168199482641,(double)
	    209.342586752537,(double)213.532241494563,(double)217.736934113954,(
	    double)221.95644181913,(double)226.190548323728,(double)
	    230.439043565777,(double)234.701723442818,(double)238.978389561834,(
	    double)243.268849002983,(double)247.572914096187,(double)
	    251.890402209723,(double)256.22113555001,(double)260.564940971863,(
	    double)264.921649798553,(double)269.29109765102,(double)
	    273.673124285694,(double)278.067573440366,(double)282.47429268763,(
	    double)286.893133295427,(double)291.32395009427,(double)
	    295.766601350761,(double)300.220948647014,(double)304.686856765669,(
	    double)309.164193580147,(double)313.652829949879,(double)
	    318.152639620209,(double)322.663499126726,(double)327.185287703775,(
	    double)331.717887196928,(double)336.261181979198,(double)
	    340.815058870799,(double)345.379407062267,(double)349.95411804077,(
	    double)354.539085519441,(double)359.134205369575 };

    /* System generated locals */
    long int i__1;
    double ret_val = 0.0;

    /* Local variables */
    STATIC_DATA double dgam;
    STATIC_DATA long int i__;
    STATIC_DATA double t, dx, px, qx, rx, xx;
    STATIC_DATA long int ndx, nxm;
    STATIC_DATA double sum, rxx;

    if ( x <= (double)0.) {
		goto L90;
    } else {
		goto L5;
    }
L5:
    ndx = (long int)x; //POL added cast
    t = x - (double) ndx;
    if (t == (double)0.) {
		goto L51;
    }
    dx = xlim1 - x;
    if (dx < (double)0.) {
		goto L40;
    }

/*     RATIONAL CHEBYSHEV APPROXIMATION ON 2.LT.X.LT.3 FOR GAMMA(X) */

    nxm = ndx - 2;
    px = pcoe[0];
    for (i__ = 2; i__ <= 9; ++i__) {
/* L10: */
	px = t * px + pcoe[i__ - 1];
    }
    qx = qcoe[0];
    for (i__ = 2; i__ <= 4; ++i__) {
/* L15: */
	qx = t * qx + qcoe[i__ - 1];
    }
    dgam = px / qx;
    if (nxm > 0) {
	goto L22;
    }
    if (nxm == 0) {
	goto L25;
    }

/*     BACKWARD RECURSION FOR 0.LT.X.LT.2 */

    dgam /= t + (double)1.;
    if (nxm == -1) {
	goto L25;
    }
    dgam /= t;
    ret_val = log(dgam);
    return ret_val;

/*     FORWARD RECURSION FOR 3.LT.X.LT.8 */

L22:
    xx = t + (double)2.;
    i__1 = nxm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dgam *= xx;
/* L24: */
	xx += (double)1.;
    }
L25:
    ret_val = log(dgam);
    return ret_val;

/*     X.GT.XLIM1 */

L40:
    rx = (double)1. / x;
    rxx = rx * rx;
    if (x - xlim2 < (double)0.) {
	goto L41;
    }
    px = q[0] * rxx + q[1];
    ret_val = px * rx + (x - (double).5) * log(x) - x + rtwpil;
    return ret_val;

/*     X.LT.XLIM2 */

L41:
    px = p[0];
    sum = (x - (double).5) * log(x) - x;
    for (i__ = 2; i__ <= 5; ++i__) {
	px = px * rxx + p[i__ - 1];
/* L42: */
    }
    ret_val = px * rx + sum + rtwpil;
    return ret_val;

/*     TABLE LOOK UP FOR INTEGER ARGUMENTS LESS THAN OR EQUAL 100. */

L51:
    if (ndx > 100) {
	goto L40;
    }
    ret_val = gln[ndx - 1];
    return ret_val;
L90:
    //cerr << "GAMLN  ARGUMENT IS LESS THAN OR EQUAL TO ZERO " << endl;
    return ret_val;
}

//@@ Yang's factorial (added because we are currently using Yang's YangsLnGamma and rndpoisson
long factorial(int n);
long factorial(int n)
	{
   long f, i;
   assert(n <10); //@should throw
   for (i=2,f=1; i<=(long)n; ++i)
   	f*=i;
   return (f);
	}
//@@ Yang's YangsLnGamma (added because we are currently using Yang's rndpoisson
double YangsLnGamma (double x);
double YangsLnGamma (double x)
	{
/* returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.

   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double f=0, fneg=0, z, lng;
   int nx=(int)x-1;

   if((double)nx==x && nx>0 && nx<10)
      lng=log((double)factorial(nx));
   else {
      	assert(x>0); //@should throw
      if (x<7) {
         f=1;  z=x-1;
         while (++z<7)  f*=z;
         x=z;   f=-log(f);
      }
      z = 1/(x*x);
      lng = fneg+ f + (x-0.5)*log(x) - x + .918938533204673 
             + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
                  +.083333333333333)/x;
   }
   return  lng;
}

//@@using Yang's rndpoisson
unsigned Lot::Poisson (double m)
	{
#	if defined(C_FUNCS_IN_STD_NAMESPACE)
		using std::tan;
		using std::floor;
#	endif
	/* m is the rate parameter of the poisson
	   Numerical Recipes in C, 2nd ed. pp. 293-295
	*/
	   STATIC_DATA double sq, alm, g, oldm=-1;
	   double em, t, y;

	/* search from the origin
	   if (m<5) { 
	      if (m!=oldm) { oldm=m; g=exp(-m); }
	      y=rndu();  sq=alm=g;
	      for (em=0; ; ) {
	         if (y<sq) break;
	         sq+= (alm*=m/ ++em);
	      }
	   }
	*/
	   if (m<12) { 
	      if (m!=oldm) { oldm=m; g=exp(-m); }
	      em=-1; t=1;
	      for (; ;) {
	         ++em;
	         t*=Uniform();
	         if (t<=g)
	         	break;
	      }
	   }
	   else {
	     if (m!=oldm) {
	        oldm=m;  sq=sqrt(2*m);  alm=log(m);
	        g=m*alm-YangsLnGamma(m+1);
	     }
	     do {
	        do {
	           y=tan(3.141592654*Uniform());
	           em=sq*y+m;
	        } while (em<0);
	        em=floor(em);
	        t=0.9*(1+y*y)*exp(em*alm-YangsLnGamma(em+1)-g);
	     } while (Uniform()>t);
	   }
	   return ((unsigned) em);
	}

#if defined(PHYCAS_EXTENDING_PYTHON)
	#include <boost/python/class.hpp>
	#include <boost/python/module.hpp>
	#include <boost/python/def.hpp>
	
	BOOST_PYTHON_MODULE(temp_jam)
	{
		using namespace boost::python;
		class_<Lot>("Lot", init<unsigned int>())
			// Add a regular member function.
			.def("binomial", &Lot::Binomial)
			;
		
	}
#endif
