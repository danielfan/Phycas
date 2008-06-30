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
#if POLPY_NEWWAY
#include "phycas/src/univents.hpp"

namespace phycas {




void Univents::setEndStates(const int8_t *p)
{
	PHYCAS_ASSERT(p);
	is_valid = false;
	const unsigned n = end_states_vec.size();
	PHYCAS_ASSERT(n > 0);
	for (unsigned i = 0; i < n ; ++i)
		end_states_vec[i] = p[i];
}


void Univents::resize(unsigned n)
	{
	univents.resize(n);
	times.resize(n);
	end_states_vec.resize(n);
	mdot = UINT_MAX;
	is_valid = false;
	times_valid = false;
	}	

void Univents::swap(Univents & other)
	{
	univents.swap(other.univents);
	times.swap(other.times);
	end_states_vec.swap(other.end_states_vec);
	std::swap(mdot, other.mdot);
	std::swap(is_valid, other.is_valid);
	std::swap(times_valid, other.times_valid);
	}	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of univents for site `site'. Assumes `unimap' is true and `size' is less than the length of the
|	`univents' vector.
*/
unsigned Univents::getNumEvents(
  unsigned site) const		/**< is the site of interest */
	{
	PHYCAS_ASSERT(this->is_valid);
	PHYCAS_ASSERT(site < univents.size());
	return univents.at(site).size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of univent states for site `site'. Assumes `unimap' is true and `size' is less than the length of 
|	the `state_time' vector. This function is not particularly efficient, and it intended primarily for transferring
|	univent states to Python code for debugging purposes.
*/
std::vector<unsigned> Univents::getEventsVec(
  unsigned site) const		/**< is the site of interest */
	{
	std::vector<unsigned> v(getNumEvents(site));
	const StateMapping & u = univents.at(site);
	unsigned i = 0;
	for (StateMapping::const_iterator it = u.begin(); it != u.end(); ++it, ++i)
		v[i] = (unsigned)*it;
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of univent times for site `site'. Assumes `unimap' is true and `size' is less than the length of 
|	the `state_time' vector. This function is not particularly efficient, and it intended primarily for transferring
|	univent times to Python code for debugging purposes.
*/
std::vector<double> Univents::getTimes(
  unsigned site) const		/**< is the site of interest */
	{
	PHYCAS_ASSERT(times_valid);
	std::vector<double> v(times.at(site));
	return v;
	}


void Univents::fillStateCodeArray(int8_t * tipSpecificStateCode) const
	{
	const unsigned n = end_states_vec.size();
	PHYCAS_ASSERT(n > 0);
	for (unsigned i = 0 ; i < n; ++i)
		tipSpecificStateCode[i] = end_states_vec[i];
	}

} // namespace phycas {

#endif //POLPY_NEWWAY

