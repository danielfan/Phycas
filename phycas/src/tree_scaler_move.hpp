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

#if ! defined(TREE_SCALER_MOVE_HPP)
#define TREE_SCALER_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
//#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "phycas/src/mcmc_updater.hpp"		// for base class MCMCUpdater

namespace phycas
{

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates a whole-tree-scaling move. Each update scales all edge lengths by a common factor. Useful for moving
|   quickly to an appropriate scale for the tree. Depending on the LargetSimonMove alone sometimes requires a long
|   burn-in period due to the fact that each LargetSimonMove update only affects three edges.
|>
*/
class TreeScalerMove : public MCMCUpdater
	{
	public:
						TreeScalerMove();
						virtual ~TreeScalerMove() 
							{
							//std::cerr << "TreeScalerMove dying..." << std::endl;
							}

        void            proposeNewState();

		// These are virtual functions in the MCMCUpdater base class
        virtual void    setPosteriorTuningParam(double x);
        virtual void    setPriorTuningParam(double x);
		virtual void	setBoldness(double x);
		virtual bool	update();
		virtual double	recalcPrior();			// override virtual from MCMCUpdater base class
		virtual void	revert();
		virtual void	accept();
		virtual double	getLnHastingsRatio() const;
		
		double			recalcWorkingPriorForMove() const;	// recomputes working prior for purpose of this move only (not used to recompute joint working prior)

    private:

		double			boldness;		/**< Ranges from 0 to 100 and determines the boldness of the move */
        double          lambda;         /**< the tuning parameter used for this move */
        double          min_lambda;     /**< the tuning parameter used for exploring the posterior distribution */
        double          max_lambda;     /**< the tuning parameter used for exploring the prior distribution */
        unsigned        n;              /**< the number of edges in the tree */
        double          m;              /**< the preproposal tree length */
        double          mstar;          /**< the proposed new tree length */
        double          forward_scaler; /**< the factor by which the original edge lengths should be multiplied to create proposed new edge lengths */
        double          reverse_scaler; /**< the factor by which the new edge lengths should be multiplied to return them to their original lengths*/
	};

} // namespace phycas

#include "phycas/src/tree_scaler_move.inl"

#endif
