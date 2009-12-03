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


#include "phycas/src/cond_likelihood.hpp"

namespace phycas
{

#if POLPY_NEWWAY
unsigned CondLikelihood::calcCLALength(
  const uint_vect_t & npatterns, 	/**< is a vector containing the number of data patterns for each partition subset */
  const uint_vect_t & nrates, 		/**< is a vector containing the number of among-site relative rate categories for each partition subset */
  const uint_vect_t & nstates) 		/**< is a vector containing the number of states for each partition subset */
	{
	unsigned sz = (unsigned)npatterns.size();
	PHYCAS_ASSERT(nrates.size() == sz);
	PHYCAS_ASSERT(nstates.size() == sz);
	
	unsigned total = 0;
	for (unsigned i = 0; i < sz; ++i)
		{
		total += (npatterns[i]*nrates[i]*nstates[i]);
		}
	return total;
	}
#endif

} //namespace phycas

