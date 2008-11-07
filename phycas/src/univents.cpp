/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2008 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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
#include "phycas/src/univents.hpp"

namespace phycas {

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor simply sets initial values for `mdot' (UINT_MAX) and `is_valid' (false).
*/
Univents::Univents()
  : mdot(UINT_MAX), /**< is the total number of univents */
  is_valid(false)   /**< is true if univent mapping is in a valid state, false if it needs to be refreshed */
    {
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Copies the end states in the supplied vector `p' to the data member `end_states_vec'. This function sets 'is_valid'
|   to false because changing the end states invalidates any univents currently mapped.
*/
void Univents::setEndStates(
  const int8_t * p)  /**< is the vector of end states to be copied to `end_states_vec' */
    {
	PHYCAS_ASSERT(p);
	is_valid = false;   //@POL shouldn't we also set times_valid to false here?
	const unsigned n = (const unsigned)end_states_vec.size();
	PHYCAS_ASSERT(n > 0);
	for (unsigned i = 0; i < n; ++i)
		end_states_vec[i] = p[i];   //@POL should assert that p is long enough for this
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Resizes the `univents', `times' and `end_states_vec' vectors, and resets `mdot' to UINT_MAX. This function sets 
|   'is_valid' and `times_valid' to false because changing the sizes of these vectors invalidates any currently mapped
|   univents.
*/
void Univents::resize(unsigned n)
	{
	univents.resize(n);
	times.resize(n);
	end_states_vec.resize(n);
	mdot = UINT_MAX;
	is_valid = false;
	times_valid = false;
	}	

/*----------------------------------------------------------------------------------------------------------------------
|   Swaps all data members with `other'.
*/
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
|	Returns the number of univents for site `site'. Assumes `is_valid' is true and `site' is less than the length of the
|	`univents' vector.
*/
unsigned Univents::getNumEvents(
  unsigned site) const		/**< is the site of interest */
	{
	PHYCAS_ASSERT(this->is_valid);
	PHYCAS_ASSERT(site < univents.size());
	return (unsigned)univents.at(site).size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of univent states for site `site'. Assumes `is_valid' is true and `site' is less than the length of 
|	the `univents' vector. This function is not particularly efficient, and it intended primarily for transferring
|	univent states to Python code for debugging purposes.
*/
std::vector<unsigned> Univents::getEventsVec(
  unsigned site) const		/**< is the site of interest */
	{
	PHYCAS_ASSERT(is_valid);
	std::vector<unsigned> v(getNumEvents(site));
	const StateMapping & u = univents.at(site);
	unsigned i = 0;
	for (StateMapping::const_iterator it = u.begin(); it != u.end(); ++it, ++i)
		v[i] = (unsigned)*it;
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of univent times for site `site'. Assumes `times_valid' is true and `site' is less than the length 
|	the `times' vector. This function is not particularly efficient, and it intended primarily for transferring
|	univent times to Python code for debugging purposes.
*/
std::vector<double> Univents::getTimes(
  unsigned site) const		/**< is the site of interest */
	{
	PHYCAS_ASSERT(times_valid);
	std::vector<double> v(times.at(site));
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Fills the supplied `tipSpecificStateCode' array with the elements of the `end_states_vec' vector. Before calling 
|   this function, ensure that `tipSpecificStateCode' is long enough (this assumption is not checked in this function).
*/
void Univents::fillStateCodeArray(int8_t * tipSpecificStateCode) const
	{
	const unsigned n = (const unsigned)end_states_vec.size();
	PHYCAS_ASSERT(n > 0);
	for (unsigned i = 0 ; i < n; ++i)
		tipSpecificStateCode[i] = end_states_vec[i];
	}

} // namespace phycas

