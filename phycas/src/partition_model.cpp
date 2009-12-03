/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2009 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#include "phycas/src/partition_model.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	PartitionModel constructor.
*/
PartitionModel::PartitionModel()
    {
	}

/*----------------------------------------------------------------------------------------------------------------------
|	PartitionModel destructor.
*/
PartitionModel::~PartitionModel()
	{
	//std::cerr << "\n>>>>> PartitionModel dying..." << std::endl;
	} 
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns current number of subsets by returning the length of the `subset_model' vector.
*/
unsigned PartitionModel::getNumSubsets()
    {
    // In debug mode check to make sure everyone agrees about the number of subsets
    PHYCAS_ASSERT(subset_num_patterns.size() == subset_model.size());
    PHYCAS_ASSERT(subset_num_states.size() == subset_model.size());
    PHYCAS_ASSERT(subset_num_rates.size() == subset_model.size());
    
    return (unsigned)subset_model.size();
    }
    
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the total number of patterns (sum over all partition subsets).
*/
unsigned PartitionModel::getTotalNumPatterns()
    {
    //@POL need to make total_num_patterns data member if this function is called often!
    return (unsigned)std::accumulate(subset_num_patterns.begin(), subset_num_patterns.end(), 0);
    }
    
} // namespace phycas