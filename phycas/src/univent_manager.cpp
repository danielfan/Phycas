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

#include "phycas/src/univent_manager.hpp"

namespace phycas
{
#if 0
/*----------------------------------------------------------------------------------------------------------------------
|	The constructor dymanically allocates the `state_time' data member and sets `map' to 0.
*/
UniventManager::UniventManager()
    {
    mdot = 0;
    state_time = new StateTimeListVect();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	The destructor deletes the dynamically-allocated `state_time' vector.
*/
UniventManager::~UniventManager()
    {
    PHYCAS_ASSERT(state_time != NULL);
    delete state_time;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the `state_time' vector.
*/
StateTimeListVect & UniventManager::getVect()
    {
    PHYCAS_ASSERT(state_time != NULL);
    return *state_time;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the length of the `state_time' vector.
*/
unsigned UniventManager::size() const
    {
    PHYCAS_ASSERT(state_time != NULL);
    return state_time->size();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns an iterator to (one beyond) the end of the `state_time' vector.
*/
void UniventManager::resize(
  unsigned sz)  /**< is the new size of the `state_time' vector */
    {
    PHYCAS_ASSERT(state_time != NULL);
    state_time->resize(sz);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the `i'th. element of the `state_time' vector.
*/
StateTimeList & UniventManager::operator[](
  unsigned i)   /**< is the index of the element to return */
    {
    PHYCAS_ASSERT(state_time != NULL);
    return (*state_time)[i];
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a reference to the (const) `i'th. element of the `state_time' vector.
*/
const StateTimeList & UniventManager::operator[](
  unsigned i) const   /**< is the index of the element to return */
    {
    PHYCAS_ASSERT(state_time != NULL);
    return (*state_time)[i];
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Swaps `state_time' and `map' data members with `other'.
*/
void UniventManager::swap(
  UniventManager & other) /**< is the other UniventManager involved in the swap */
    {
    unsigned tmpmdot = mdot;
    mdot = other.mdot;
    other.mdot = tmpmdot;

    StateTimeListVect * tmp = state_time;
    state_time = other.state_time;
    other.state_time = tmp;
    }
#endif
}   // namespace phycas
