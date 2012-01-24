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
unsigned PartitionModel::getNumSubsets() const
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
unsigned PartitionModel::getTotalNumPatterns() const
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
|	Returns const reference to `subset_relrates' vector.
*/
double PartitionModel::getSubsetRelRate(unsigned i) const
	{
	PHYCAS_ASSERT(subset_relrates.size() > i);
	return subset_relrates[i];
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_relrates' vector.
*/
const std::vector<double> & PartitionModel::getSubsetRelRatesVect() const
	{
	return subset_relrates;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_relrates' vector.
*/
MultivarProbDistShPtr PartitionModel::getSubsetRelRatePrior() const
	{
	return subset_relrate_prior;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the values in the data member vector `subset_relrates' with those in the supplied `rrates' vector. Note 
|   that we cannot assume that the values supplied in `rrates' are normalized. To normalize, we need to solve for the
|   value `x' in the following example involving 2 subsets, one of which contains 4/5 of the sites and an unnormalized
|   relative rate of 1, with the other subset containing 1/5 of the sites and having an unnormalized relative rate of
|   10:
|>
|   (0.8)*(1/x) + (0.2)*(10/x) = 1.0  ==> x = (0.8)*(1) + (0.2)*(10) = 2.8
|>
*/
void PartitionModel::setSubsetRelRatesVect(
  const std::vector<double> & rrates)	/**< is the vector of relative rates */
	{
	PHYCAS_ASSERT(rrates.size() == getNumSubsets());
	PHYCAS_ASSERT(subset_relrates.size() == getNumSubsets());
    
    double num_sites_total = (double)std::accumulate(subset_num_sites.begin(), subset_num_sites.end(), 0.0);
    //std::cerr << "num_sites_total:" << num_sites_total << std::endl;
    
    std::vector<double> relrate_proportion;
    relrate_proportion.resize(rrates.size());
    std::transform(rrates.begin(), rrates.end(), subset_num_sites.begin(), relrate_proportion.begin(), (boost::lambda::_1)*(boost::lambda::_2)/num_sites_total);
    double x = (double)std::accumulate(relrate_proportion.begin(), relrate_proportion.end(), 0.0);
    //std::cerr << "x:" << x << std::endl;

    std::transform(rrates.begin(), rrates.end(), subset_relrates.begin(), boost::lambda::_1/x);

    //std::cerr << ">>>>>>>>>> PartitionModel::setSubsetRelRatesVect <<<<<<<<<<<" << std::endl;
    //std::cerr << "subset_relrates:";
	//std::copy(subset_relrates.begin(), subset_relrates.end(), std::ostream_iterator<double>(std::cerr, " "));
    //std::cerr << std::endl;
    
    //std::cerr << "subset_num_sites:";
	//std::copy(subset_num_sites.begin(), subset_num_sites.end(), std::ostream_iterator<unsigned>(std::cerr, " "));
    //std::cerr << std::endl;
    
    //std::cerr << "relrate_proportion:";
	//std::copy(relrate_proportion.begin(), relrate_proportion.end(), std::ostream_iterator<double>(std::cerr, " "));
    //std::cerr << std::endl;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the multivariate probability distribution (`subset_relrate_prior') used for subset relative rates.
*/
void PartitionModel::setSubsetRelRatePrior(
  MultivarProbDistShPtr rrate_prior)	/**< is the new relative rate prior */
	{
	subset_relrate_prior = rrate_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_num_patterns' vector.
*/
const std::vector<unsigned> & PartitionModel::getNumPatternsVect() const
	{
	return subset_num_patterns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of sites in partition subset `i' by summing pattern counts over patterns in that subset.
*/
unsigned PartitionModel::getNumSites(unsigned i) const
	{
	return std::accumulate(subset_num_sites.begin(), subset_num_sites.end(), 0);
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
	subset_relrates.push_back(1.0);
	subset_num_states.push_back(m->getNumStates());
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
	
	subset_relrates.assign(new_size, 1.0);
	
	// Go ahead and setup subset_num_rates and subset_num_states vectors
	// according to the models
	subset_num_rates.resize(new_size);
	subset_num_states.resize(new_size);
	subset_num_patterns.resize(new_size);
	for (unsigned i = 0; i < new_size; ++i)
		{
		subset_num_states[i] = subset_model[i]->getNumStates();
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
|	Sets `subset_num_sites' to a copy of the supplied vector `nsites'.
*/
void PartitionModel::setNumSitesVect(const std::vector<unsigned> & nsites)
	{
	unsigned new_size = (unsigned)nsites.size();
	PHYCAS_ASSERT(new_size > 0);
	subset_num_sites.resize(new_size);
	std::copy(nsites.begin(), nsites.end(), subset_num_sites.begin());
	
	//tmp
	//std::cerr << "\nDebug output: here are the number of sites in each subset:\n**********\n";
	//std::copy(subset_num_sites.begin(), subset_num_sites.end(), std::ostream_iterator<unsigned>(std::cerr, " "));
	//std::cerr << "\n**********\n" << std::endl;
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
    
/*----------------------------------------------------------------------------------------------------------------------
|	Sets `site_assignments' to a copy of the supplied vector `v'.
*/
void PartitionModel::setSiteAssignments(const std::vector<unsigned> & v)
	{
	unsigned new_size = (unsigned)v.size();
	PHYCAS_ASSERT(new_size > 0);
	site_assignments.resize(new_size);
	std::copy(v.begin(), v.end(), site_assignments.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `site_assignments' to a copy of the supplied vector `v'.
*/
const std::vector<unsigned> & PartitionModel::getSiteAssignments() const
	{
    return site_assignments;
	}
        
} // namespace phycas
