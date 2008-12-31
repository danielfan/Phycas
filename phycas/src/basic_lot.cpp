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

#include "phycas/src/basic_lot.hpp"

const unsigned MASKSIGNBIT = 0x80000000;

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Default constructor. Sets `last_seed_setting' and `curr_seed' both to 1U, creates a new CDF object and stored the
|	pointer in `cdf_converter', then calls the UseClockToSeed function.
*/
Lot::Lot() : last_seed_setting(1U), curr_seed(1U), num_seeds_generated(0)
	{
	UseClockToSeed();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor taking a starting seed. Initializes `last_seed_setting' and `curr_seed' both to 1U, creates a new CDF 
|	object and stores the pointer in `cdf_converter', then calls the UseClockToSeed function if `rnd_seed' is zero or
|	UINT_MAX, or the SetSeed function if `rnd_seed' is any other value.
*/
Lot::Lot(unsigned rnd_seed) : last_seed_setting(1U), curr_seed(1U), num_seeds_generated(0)
	{
	//cdf_converter = new CDF();
	if (rnd_seed == 0 || rnd_seed == UINT_MAX)
		UseClockToSeed();
	else
		SetSeed(rnd_seed);
	}

#include <iostream>
/*----------------------------------------------------------------------------------------------------------------------
|	Destructor deletes `cdf_converter'.
*/
Lot::~Lot() 
	{
	std::cerr << "\n>>>>> Lot object dying..." << std::endl;
	//delete cdf_converter;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Member function returning uniform deviate. Provided by J. Monahan, Statistics Dept., North Carolina State 
|	University. Originally from Schrage, ACM Trans. Math. Software 5:132-138 (1979). Translated to C by Paul O. Lewis, 
|	Dec. 10, 1992.
*/
#if defined(LOG_LOT_UNIFORM_CALLS)
double Lot::Uniform(const char * file, const int line)
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
#if defined(LOG_LOT_UNIFORM_CALLS)
#   pragma TODO("undef LOG_LOT_UNIFORM_CALLS in std_force_include.hpp for release version")
	double retval = curr_seed * 4.6566128575e-10;
    std::ofstream tmpf("uniform.txt", std::ios::out | std::ios::app);
    tmpf << num_seeds_generated << " -> " << retval << " (" << file << ":" << line << ")" << std::endl;
    tmpf.close();
	num_seeds_generated++;
	return retval;
#else
	return curr_seed * 4.6566128575e-10;
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns an unsigned integer in [0-`max') in which all values are equiprobable.
*/
unsigned Lot::SampleUInt(unsigned max)
	{
	PHYCAS_ASSERT(max > 0);

	unsigned samples_uint = max;
	while(samples_uint == max) 
		{
		double r = Uniform(FILE_AND_LINE);
		samples_uint = (unsigned)((double)max*r);
		}

	return samples_uint;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a multinomial deviate in [0,n) given the bin probabilities in probs. The sum of the first n elements of 
|   probs is assumed to equal 1.0.
*/
unsigned Lot::MultinomialDraw(const double * probs, unsigned n, double totalProb)
	{
    PHYCAS_ASSERT(probs != NULL);
    PHYCAS_ASSERT(n > 0);
    double u = totalProb*Uniform(FILE_AND_LINE);
    for (unsigned i = 0; i < n; ++i)
        {
        u -= probs[i];
        if (u < 0.0)
            return i;
        }
    return n-1;
	}

// below here formerly in basic_lot.inl

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `curr_seed', which stores the current seed (which changes after each draw).
*/
unsigned Lot::GetSeed() const
	{
	return curr_seed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `last_seed_setting', which stores the seed used to initialize generator.
*/
unsigned Lot::GetInitSeed() const 
	{
	return last_seed_setting;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `last_seed_setting', which stores the seed used to initialize generator.
*/
void Lot::SetSeed(unsigned s)
	{
	PHYCAS_ASSERT(s > 0 && s < UINT_MAX);
	curr_seed = last_seed_setting = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `last_seed_setting', which stores the seed used to initialize generator.
*/
void Lot::UseClockToSeed()
	{
	time_t timer;
	curr_seed = (unsigned)time(&timer);
	last_seed_setting = curr_seed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Designed to emulate Python function random.getrandbits(k). Returns an unsigned value in which the low `nbits' bits
|	are randomly chosen. Assumes `nbits' is greater than 0 and less than 32.
|>
|	nbits   decimal             binary
|	-----------------------------------------------------------------------
|	1       0, 1                0 or 1
|	2		0, 1, 2, 3	        00, 01, 10, 11
|	3       0, 1, 2, ..., 7     000, 001, 010, 011, 100, 101, 110, 111
|	-----------------------------------------------------------------------
|>
|	In general, specifying `nbits' to n results in a random choice of unsigned values in the range [0, 1, ..., (2^n)-1].  
*/
unsigned Lot::GetRandBits(unsigned nbits)
	{
	PHYCAS_ASSERT(nbits > 0);
	PHYCAS_ASSERT(nbits < 32);

	double u = Uniform(FILE_AND_LINE);

	if (u == 0.0)
		{
		return 0;
		}
	else
		{
		// e.g. for nbits = 3, u = 0.99
		//   term1 = log(8)
		//   term2 = log(0.99)
		//   term3 = 0.99*8 = 7.92
		//   return floor(7.92) = 7
		double term1 = log(2.0)*(double)nbits;
		double term2 = log(u);
		double term3 = exp(term1 + term2);
		return (unsigned)floor(term3);
		}
	}

