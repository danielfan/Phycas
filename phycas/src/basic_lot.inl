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

#if ! defined(BASIC_LOT_INL)
#define BASIC_LOT_INL

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `curr_seed', which stores the current seed (which changes after each draw).
*/
inline unsigned Lot::GetSeed() const
	{
	return curr_seed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `last_seed_setting', which stores the seed used to initialize generator.
*/
inline unsigned Lot::GetInitSeed() const 
	{
	return last_seed_setting;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `last_seed_setting', which stores the seed used to initialize generator.
*/
inline void Lot::SetSeed(unsigned s)
	{
	PHYCAS_ASSERT(s > 0 && s < UINT_MAX);
	curr_seed = last_seed_setting = s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `last_seed_setting', which stores the seed used to initialize generator.
*/
inline void Lot::UseClockToSeed()
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
inline unsigned Lot::GetRandBits(unsigned nbits)
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

}  // namespace phycas

#endif
