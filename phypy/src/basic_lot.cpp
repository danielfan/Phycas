/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

#include "phypy/src/basic_lot.hpp"

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
#if defined(NDEBUG)
double Lot::Uniform()
#else
double Lot::Uniform(const char * file, int line)
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
#if 1 || defined(NDEBUG)
	return curr_seed * 4.6566128575e-10;
#else
	double retval = curr_seed * 4.6566128575e-10;
	std::cerr << num_seeds_generated << " -> " << retval << "(" << file << ":" << line << ")" << std::endl;
	num_seeds_generated++;
	return retval;
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

