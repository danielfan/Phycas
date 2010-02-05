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
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_model' vector.
*/
const ModelVect & PartitionModel::getModelsVect() const
	{
	return subset_model;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_num_patterns' vector.
*/
const std::vector<unsigned> & PartitionModel::getNumPatternsVect() const
	{
	return subset_num_patterns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_num_states' vector.
*/
const std::vector<unsigned> & PartitionModel::getNumStatesVect() const
	{
	return subset_num_states;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_num_rates' vector.
*/
const std::vector<unsigned> & PartitionModel::getNumRatesVect() const
	{
	return subset_num_rates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds model `m' to the `subset_model' vector.
*/
void PartitionModel::addModel(ModelShPtr m)
	{
	PHYCAS_ASSERT(m);
	subset_model.push_back(m);
	
	// Update subset_num_rates and subset_num_states vectors
	subset_num_states.push_back(m->getNStates());
	subset_num_rates.push_back(m->getNRatesTotal());
	subset_num_patterns.push_back(0);	// can't get this from model
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `subset_model' to a copy of the supplied vector `models'.
*/
void PartitionModel::setModelsVect(const ModelVect & models)
	{
	unsigned new_size = (unsigned)models.size();
	PHYCAS_ASSERT(new_size > 0);
	subset_model.resize(new_size);
	std::copy(models.begin(), models.end(), subset_model.begin());
	
	// Go ahead and setup subset_num_rates and subset_num_states vectors
	// according to the models
	subset_num_rates.resize(new_size);
	subset_num_states.resize(new_size);
	subset_num_patterns.resize(new_size);
	for (unsigned i = 0; i < new_size; ++i)
		{
		subset_num_states[i] = subset_model[i]->getNStates();
		subset_num_rates[i] = subset_model[i]->getNRatesTotal();
		subset_num_patterns[i] = 0;	// can't get this from model
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `subset_num_patterns' to a copy of the supplied vector `npatterns'.
*/
void PartitionModel::setNumPatternsVect(const std::vector<unsigned> & npatterns)
	{
	unsigned new_size = (unsigned)npatterns.size();
	PHYCAS_ASSERT(new_size > 0);
	subset_num_patterns.resize(new_size);
	std::copy(npatterns.begin(), npatterns.end(), subset_num_patterns.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `subset_num_states' to a copy of the supplied vector `nstates'.
*/
void PartitionModel::setNumStatesVect(const std::vector<unsigned> & nstates)
	{
	unsigned new_size = (unsigned)nstates.size();
	PHYCAS_ASSERT(new_size > 0);
	subset_num_states.resize(new_size);
	std::copy(nstates.begin(), nstates.end(), subset_num_states.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `subset_num_rates' to a copy of the supplied vector `nrates'.
*/
void PartitionModel::setNumRatesVect(const std::vector<unsigned> & nrates)
	{
	unsigned new_size = (unsigned)nrates.size();
	PHYCAS_ASSERT(new_size > 0);
	subset_num_rates.resize(new_size);
	std::copy(nrates.begin(), nrates.end(), subset_num_rates.begin());
	}
    
} // namespace phycas