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

//#include "phycas/src/cipres/CipresDataMatrixHelper.h"
//#include "phycas/src/probability_distribution.hpp"
//#include "phycas/src/likelihood_models.hpp"
//#include "phycas/src/basic_tree_node.hpp"
#include "phycas/src/tree_likelihood.hpp"
//#include "phycas/src/xlikelihood.hpp"
//#include "phycas/src/mcmc_chain_manager.hpp"
//#include "phycas/src/basic_tree.hpp"
//#include "boost/format.hpp"
#include "phycas/src/mapping_move.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|   MappingMove constructor. Nothing to be done because this class only manipulates the TreeLikelihood object
|   and has no data members of its own.
*/
MappingMove::MappingMove()
	{
    }

/*----------------------------------------------------------------------------------------------------------------------
|   MappingMove destructor.
*/
MappingMove::~MappingMove()
	{
    //std::cerr << "MappingMove dying..." << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new univent mapping for all characters on the tree. This move regenerates only latent variables, and is 
|   always accepted.
*/
bool MappingMove::update()
	{
    PHYCAS_ASSERT(likelihood->isUsingUnimap());
    likelihood->fullRemapping(tree, rng, true);
    return true;
    }

}   // phycas namespace
