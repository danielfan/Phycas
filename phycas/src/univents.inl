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

#if ! defined(UNIVENTS_INL)
#define UNIVENTS_INL

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the length of the `univents' vector.
*/
inline unsigned Univents::size() const 
    {
    return univents.size();
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the current value of `is_valid'.
*/
inline bool Univents::isValid() const 
    {
    return is_valid;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns true if both `is_valid' and `times_valid' are true, and returns false if either or both of these are false.
*/
inline bool Univents::timesValid() const 
    {
    return is_valid && times_valid;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the current value of `mdot'. Assumes `is_valid' is true.
*/
inline unsigned Univents::getMDot() const
    {
    PHYCAS_ASSERT(is_valid);
    return mdot;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a const reference to the `end_states_vec' vector.
*/
inline const std::vector<int8_t> & Univents::getEndStatesVecConstRef() const 
    {
    return end_states_vec;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a non-const reference to the `end_states_vec' vector.
*/
inline std::vector<int8_t> & Univents::getEndStatesVecRef() 
    {
    return end_states_vec;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a const reference to the `univents' vector.
*/
inline const std::vector< StateMapping > & Univents::getVecEventsVecConstRef() const 
    {
    return univents;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a non-const reference to the `univents' vector.
*/
inline std::vector< StateMapping > & Univents::getVecEventsVecRef() 
    {
    return univents;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a const reference to the StateMapping object for site `i'.
*/
inline const StateMapping & Univents::getEventsVecConstRef(
  unsigned i) const  
    {
    return univents.at(i);
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns a non-const reference to the StateMapping object for site `i'.
*/
inline StateMapping & Univents::getEventsVecRef(
  unsigned i) 
    {
    return univents.at(i);
    }

} // namespace phycas
#endif
