#include "pyphy/prob_dist/basic_lot.hpp"

const unsigned MASKSIGNBIT = 0x80000000;

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Default constructor. Sets `last_seed_setting' and `curr_seed' both to 1U, creates a new CDF object and stored the
|	pointer in `cdf_converter', then calls the UseClockToSeed function.
*/
Lot::Lot() : last_seed_setting(1U), curr_seed(1U)
#if POLPY_NEWWAY
	, num_seeds_generated(0)
#endif
	{
	//cdf_converter = new CDF();
	UseClockToSeed();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor taking a starting seed. Initializes `last_seed_setting' and `curr_seed' both to 1U, creates a new CDF 
|	object and stores the pointer in `cdf_converter', then calls the UseClockToSeed function if `rnd_seed' is zero or
|	UINT_MAX, or the SetSeed function if `rnd_seed' is any other value.
*/
Lot::Lot(unsigned rnd_seed) : last_seed_setting(1U), curr_seed(1U)
#if POLPY_NEWWAY
	, num_seeds_generated(0)
#endif
	{
	//cdf_converter = new CDF();
	if (rnd_seed == 0 || rnd_seed == UINT_MAX)
		UseClockToSeed();
	else
		SetSeed(rnd_seed);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor deletes `cdf_converter'.
*/
Lot::~Lot() 
	{
	//delete cdf_converter;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Member function returning uniform deviate. Provided by J. Monahan, Statistics Dept., North Carolina State 
|	University. Originally from Schrage, ACM Trans. Math. Software 5:132-138 (1979). Translated to C by Paul O. Lewis, 
|	Dec. 10, 1992.
*/
#if POLPY_NEWWAY
#	if defined(NDEBUG)
		double Lot::Uniform()
#	else
		double Lot::Uniform(const char * file, int line)
#	endif
#else
	double Lot::Uniform()
#endif
	{
	const unsigned a = 16807U;
	const unsigned b15 = 32768U;
	const unsigned b16 = 65536U;
	const unsigned p = 2147483647U;
	const unsigned xhi = curr_seed / b16;
	const unsigned xalo = (curr_seed - xhi * b16) * a;
	const unsigned leftlo = xalo / b16;
	const unsigned fhi = xhi * a + leftlo;
	const unsigned k = fhi / b15;
	curr_seed = (((xalo - leftlo * b16) - p) + (fhi - k * b15) * b16) + k;
	if (curr_seed & MASKSIGNBIT) 
		curr_seed = unsigned (((int)curr_seed) + p);
#if 0 & POLPY_NEWWAY
#	if defined(NDEBUG)
		return curr_seed * 4.6566128575e-10;
#	else
		double retval = curr_seed * 4.6566128575e-10;
		std::cerr << num_seeds_generated << " -> " << retval << "(" << file << ":" << line << ")" << std::endl;
		num_seeds_generated++;
		return retval;
#	endif
#else
	return curr_seed * 4.6566128575e-10;
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns an unsigned integer in [0-`max') in which all values are equiprobable.
*/
unsigned Lot::SampleUInt(unsigned max)
	{
	assert(max > 0);

	unsigned samples_uint = max;
	while(samples_uint == max) 
		{
		double r = Uniform(FILE_AND_LINE);
		samples_uint = (unsigned)((double)max*r);
		}

	return samples_uint;
	}

