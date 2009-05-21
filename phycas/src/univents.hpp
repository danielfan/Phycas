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
#if ! defined(UNIVENTS_HPP)
#define UNIVENTS_HPP

#include <vector>

#include "phycas/src/states_patterns.hpp"

namespace phycas
{
class UniventProbMgr;
class TreeNode;
class Univents;

class Univents
    {
	public:

		                                    Univents();
		unsigned                            size() const;
		void                                resize(unsigned n);
		void                                swap(Univents & other);

		void                                setEndStates(const int8_t * p);

		bool                                isValid() const;
		void								setValid(bool v)
		{
		is_valid = v;
		}
		bool                                timesValid() const;

		const std::vector<int8_t> &         getEndStatesVecConstRef() const;
		std::vector<int8_t> &               getEndStatesVecRef();

		const std::vector< StateMapping > & getVecEventsVecConstRef() const;
		std::vector< StateMapping > &       getVecEventsVecRef();
		
		const StateMapping &                getEventsVecConstRef(unsigned i) const;
		StateMapping &                      getEventsVecRef(unsigned i);
		
		unsigned                            getNumEvents(unsigned) const;
		std::vector<unsigned>               getEventsVec(unsigned) const;
		std::vector<double>                 getTimes(unsigned) const;

		void                                fillStateCodeArray(int8_t * tipSpecificStateCode) const;
		unsigned                            getMDot() const;

	private:
		/************ remember to add any new data members to swap(); *************/
		std::vector< StateMapping > 	    univents;		/**< univents[i][j] holds a state for univent j at site i */
		std::vector< std::vector<double> >  times;		    /**< times[i][j] holds the fraction of the edgelen representing the time at which the univent occurred */
		std::vector<int8_t>				    end_states_vec; /**<* end_states_vec[i] holds the end state for site i */
		unsigned						    mdot;			/**< the total number of univents over all sites on the edge owned by this node */
		bool							    is_valid;       /**< */
		bool							    times_valid;    /**< */
		/************ remember to add any new data members to swap(); *************/
	friend class UniventProbMgr;
    };

} // phycas

#include "phycas/src/univents.inl"

#endif
