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

#if ! defined(PARTITION_MODEL_HPP)
#define PARTITION_MODEL_HPP

#include "phycas/src/likelihood_models.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Container for information about a data partition: stores number of subsets, a vector of the number of rate 
|   categories for each subset, the number of states for each subset, and a vector of shared pointers to the model 
|   applied to each subset.
*/
class PartitionModel
{
	friend class TreeLikelihood;
	friend class TipData;
	friend class InternalData;
	
public:

        								PartitionModel();
    virtual								~PartitionModel();
    	
	// Utilities
    unsigned                        	getNumSubsets();
    unsigned                        	getTotalNumPatterns();
	
	// Accessors
	ModelShPtr							getModel(unsigned i) {return subset_model[i];}
	const ModelVect					&	getModelsVect() const;
	const std::vector<unsigned>		&	getNumPatternsVect() const;
	const std::vector<unsigned>		&	getNumStatesVect() const;
	const std::vector<unsigned>		&	getNumRatesVect() const;
	
	// Setters
	void								addModel(ModelShPtr m);
	void								setModelsVect(const std::vector<ModelShPtr> & models);
	void								setNumPatternsVect(const std::vector<unsigned> & npatterns);
	void								setNumStatesVect(const std::vector<unsigned> & nstates);
	void								setNumRatesVect(const std::vector<unsigned> & nrates);

private:
    
    std::vector<unsigned>           	subset_num_patterns;	/**< The number of site patterns for each partition subset */
    std::vector<unsigned>           	subset_num_states;		/**< The number of states for each partition subset */
    std::vector<unsigned>           	subset_num_rates;		/**< The number of relative rate categories for each partition subset */
    ModelVect							subset_model;			/**< The substitution model for each partition subset */
};  // class PartitionModel

typedef boost::shared_ptr<PartitionModel> PartitionModelShPtr;

} // namespace phycas

#endif
