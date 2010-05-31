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

#if ! defined(UNIMAP_NODE_SLIDE_MOVE_HPP)
#define UNIMAP_NODE_SLIDE_MOVE_HPP

#include <vector>							// for std::vector
#include <boost/shared_ptr.hpp>				// for boost::shared_ptr
#include <boost/weak_ptr.hpp>				// for boost::weak_ptr
#include "phycas/src/mcmc_updater.hpp"		// for base class MCMCUpdater
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/unimap_nni_move.hpp"

namespace phycas
{

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager> ChainManagerWkPtr;

/*----------------------------------------------------------------------------------------------------------------------
|   
*/
class UnimapNodeSlideMove : public UnimapTopoMove
	{
	public:
						        UnimapNodeSlideMove(TreeLikeShPtr treeLike);
						        virtual ~UnimapNodeSlideMove() {}

		virtual double			getLnHastingsRatio() const;
		virtual double			getLnJacobian() const;
		virtual void			accept();
		double getWindowWidth() const {return windowSize;} 
		void setWindowWidth(double v) {windowSize = v;} 
	protected:
		enum NodeSlideMoveEnum {
			MOVE_X_TOWARD_W,
			MOVE_X_TOWARD_Z,
			MOVE_W_TOWARD_Y,
			MOVE_Z_TOWARD_Y
			};

		virtual void			ProposeStateWithTemporaries(ChainManagerShPtr &);
		
		NodeSlideMoveEnum moveType;

		// finalDist and startingDist are the distance to the sliding node from either:
		//		node C (if we are moving along the AC path), or
		//		node D (if we are moving along the AD path)

		double startingDist; 
		double finalDist; 
		double nniPartOfPath; // the portion of the path that delimits an NNI move from a move that just changes branch lengths
		double windowSize;

	};


} // namespace phycas

#endif
