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

#if ! defined(UNIMAP_EDGE_MOVE_HPP)
#define UNIMAP_EDGE_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
#include "phycas/src/mcmc_updater.hpp"		// for base class MCMCUpdater

namespace phycas
{

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	An UnimapEdgeMove changes the length of just one randomly-chosen edge in the tree. 
|	Needs work
*/
class UnimapEdgeMove : public MCMCUpdater
	{
	public:
									UnimapEdgeMove();
									virtual ~UnimapEdgeMove() 
										{
										//std::cerr << "UnimapEdgeMove dying..." << std::endl;
										}


		// Accessors
		double						getLambda() const;

		// Modifiers
		void						setLambda(double x);

		// Utilities
		void						reset();

		// These are virtual functions in the MCMCUpdater base class
        virtual void                setPosteriorTuningParam(double x);
        virtual void                setPriorTuningParam(double x);
		virtual void	            setBoldness(double x);
		virtual bool				update();
		virtual double				getLnHastingsRatio() const;
		virtual double				getLnJacobian() const;
		virtual void				proposeNewState();
		virtual void				revert();
		virtual void				accept();

	private:

		UnimapEdgeMove &			operator=(const UnimapEdgeMove &);	// never use - don't define

	private:

		double			            boldness;		/**< Ranges from 0 to 100 and determines the boldness of the move */
        double                      lambda;         /**< the tuning parameter used for this move */
        double                      min_lambda;     /**< the tuning parameter used for exploring the posterior distribution */
        double                      max_lambda;     /**< the tuning parameter used for exploring the prior distribution */
        
		double						origEdgelen;	/**< Length of modified edge saved (in case revert is necessary) */
		TreeNode *					origNode;		/**< Node owning the modified edge (in case revert is necessary) */
//		std::vector<unsigned>		mdot;			/**< Number of univents over all sites for this edge */
		unsigned					mdot;			/**< Number of univents over all sites for this edge */
		double						r;				/**< ratio of new edgelen to original edgelen */
	};

} // namespace phycas

#include "phycas/src/unimap_edge_move.inl"

#endif
