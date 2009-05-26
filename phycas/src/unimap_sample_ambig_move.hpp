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

#if ! defined(UNIMAP_SAMPLE_AMBIG_MOVE_HPP)
#define UNIMAP_SAMPLE_AMBIG_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
#include "phycas/src/mcmc_updater.hpp"		// for base class MCMCUpdater

namespace phycas
{

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	UnimapSampleAmbigMove
*/
class UnimapSampleAmbigMove : public MCMCUpdater
	{
	public:
		UnimapSampleAmbigMove(TreeLikeShPtr treeLike, TreeShPtr t, ModelShPtr model, unsigned weight);
		virtual ~UnimapSampleAmbigMove() 
			{
			}

		// Utilities
		void						reset();

		// These are virtual functions in the MCMCUpdater base class
		virtual bool				update();
		virtual double				getLnHastingsRatio() const;
		virtual double				getLnJacobian() const;
		virtual void				proposeNewState();
		virtual void				revert();
		virtual void				accept();

		// resamples all of the tips without considering the edge length or neighbor nodes (called before there is a mapping at internal nodes).
		void sampleTipsAsDisconnected();
		// returns the number of nodes that have at least one ambiguous site
		unsigned getNumAmbigNodes() const 
			{
			return ambigTipToAmbigCol.size();
			}
		
	private:

		UnimapSampleAmbigMove &			operator=(const UnimapSampleAmbigMove &);	// never use - don't define

	protected:
	
	
		typedef std::map<TreeNode *, std::vector<unsigned> > AmbigTipMap;
		AmbigTipMap ambigTipToAmbigCol;
		
		void proposeNewStateArrayForNode(TreeNode * nd, const std::vector<unsigned> & ambigInds);
		void sampleNewStateArrayForNodeAsDisconnected(TreeNode * nd, const std::vector<unsigned> & ambigInds);


	};

} // namespace phycas

#endif
