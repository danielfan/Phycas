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

#if ! defined(UNIVENT_MANAGER_HPP)
#define UNIVENT_MANAGER_HPP

#include "phycas/src/states_patterns.hpp"

namespace phycas
{
#if 1
	typedef StateTimeListVect UniventManager;
#else
class UniventManager
    {
    public:

                                        UniventManager();
                                        ~UniventManager();

        unsigned                        size() const;
        void                            resize(unsigned sz);
        StateTimeList &                 operator[](unsigned i);
        const StateTimeList &           operator[](unsigned i) const;
        StateTimeListVect &             getVect();

        void                            swap(UniventManager & other);

    private:

        unsigned            mdot;           /**< the total number of univents over all sites on the edge owned by this node */
        StateTimeListVect   state_time;     /**< state_time[i][j].first holds a state for univent j at site i, whereas state_time[i][j].second holds the fraction of the edgelen representing the time at which the univent occurred */
    };
#endif
}
#endif
