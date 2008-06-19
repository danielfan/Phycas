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

#if ! defined(UNIMAP_NNI_MOVE_HPP)
#define UNIMAP_NNI_MOVE_HPP

#include <vector>							// for std::vector
#include <boost/shared_ptr.hpp>				// for boost::shared_ptr
#include <boost/weak_ptr.hpp>				// for boost::weak_ptr
#include "phycas/src/mcmc_updater.hpp"		// for base class MCMCUpdater

namespace phycas
{

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager> ChainManagerWkPtr;
UniventManager * GetUniventManager(TreeNode * nd);

/*----------------------------------------------------------------------------------------------------------------------
|   
*/
class UnimapNNIMove : public MCMCUpdater
	{
	public:
						UnimapNNIMove();
						virtual ~UnimapNNIMove() 
							{
							//std::cerr << "UnimapNNIMove dying..." << std::endl;
							}

		// These are virtual functions in the MCMCUpdater base class
		//
		virtual bool	update();
		virtual double	getLnHastingsRatio() const;
		virtual double	getLnJacobian() const;
		virtual void	proposeNewState();
		virtual void	revert();
		virtual void	accept();
	protected:
		TreeNode * randomInternalAboveSubroot();
		TreeNode * x;
		TreeNode * y;
		TreeNode * z;

		TipData * xTipData;
		TipData * yTipData;
		TipData * zTipData;
		TipData * wTipData;

		bool x_is_left;
		
		double prev_ln_prior;	/**< The log prior of the starting state */
		double prev_ln_like;	/**< The log likelihood of the starting state */

		TipData * allocTipDataFromUnivents(TreeNode * nd , bool use_last);
		double FourTaxonLnL(TreeNode * nd);
		double FourTaxonLnLFromCorrectTipDataMembers(TreeNode * nd);
		double HarvestLnLikeFromCondLikePar(ConstCondLikelihoodShPtr focalCondLike, ConstCondLikelihoodShPtr neighborCondLike, const double * const * childPMatrix);

        void DebugSaveNexusFile(TipData * xtd, TipData * ytd, TipData * ztd, TipData * wtd, double lnlike);
	};

} // namespace phycas

#include "phycas/src/unimap_nni_move.inl"

#endif
