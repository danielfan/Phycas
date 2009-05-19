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

/*----------------------------------------------------------------------------------------------------------------------
|   
*/
class UnimapNNIMove : public MCMCUpdater
	{
	public:
						        UnimapNNIMove();
						        virtual ~UnimapNNIMove();

		// These are virtual functions in the MCMCUpdater base class
		//
		virtual bool	        update();
		virtual double	        getLnHastingsRatio() const;
		virtual double	        getLnJacobian() const;
		virtual void	        proposeNewState();
		virtual void	        revert();
		virtual void	        accept();
		virtual void	        setLot(LotShPtr p);

		bool                    getDoSampleUnivents() const;
		void                    setDoSampleUnivents(bool v);

	protected:

		void                    calculatePairwiseDistances();
		void                    calculateProposalDist(bool);

		double                  calcProposalLnDensity(double mean, double x);
		double                  proposeEdgeLen(double mean);

		TipData *               createTipDataFromUnivents(const Univents &, TipData *);
		double                  FourTaxonLnLBeforeMove(TreeNode * nd);
		double                  FourTaxonLnLFromCorrectTipDataMembers(TreeNode * nd);
		double                  HarvestLnLikeFromCondLikePar(CondLikelihoodShPtr focalCondLike, ConstCondLikelihoodShPtr neighborCondLike, const double * const * childPMatrix);
		void                    storePMatTransposed(double **& cached, const double *** p_mat_array);
        void                    DebugSaveNexusFile(TipData * xtd, TipData * ytd, TipData * ztd, TipData * wtd, double lnlike, TreeNode * nd);
	    TreeNode *              randomInternalAboveSubroot();

	protected:

        TreeNode *              x;                          /**< xxxx */
	    TreeNode *              y;                          /**< xxxx */
	    TreeNode *              z;                          /**< xxxx */
	    TreeNode *              wSis;                       /**< xxxx */
	    TreeNode *              ySis;                       /**< xxxx */
		bool                    x_is_left;                  /**< xxxx */

		TipData *               ySisTipData;                /**< xxxx */
		TipData *               yTipData;                   /**< xxxx */
		TipData *               wSisTipData;                /**< xxxx */
		TipData *               wTipData;                   /**< xxxx */
		
		double  * *             pre_x_pmat_transposed;      /**< xxxx */
		double  * *             pre_y_pmat_transposed;      /**< xxxx */
		double  * *             pre_w_pmat_transposed;      /**< xxxx */
		double  * *             pre_z_pmat_transposed;      /**< xxxx */

		bool                    doSampleUnivents;           /**< xxxx */
		GammaDistribution       gammaDist;                  /**< xxxx */
		double                  min_edge_len_mean;          /**< xxxx */
        double                  edge_len_prop_cv;           /**< xxxx */
		double                  ln_density_reverse_move;    /**< xxxx */
        double                  ln_density_forward_move;    /**< xxxx */
		double                  prev_x_len;                 /**< xxxx */
        double                  prev_y_len;                 /**< xxxx */
        double                  prev_z_len;                 /**< xxxx */
        double                  prev_nd_len;                /**< xxxx */
        double                  prev_ndP_len;               /**< xxxx */

		double                  dXY;                        /**< xxxx */
		double                  dWX;                        /**< xxxx */
		double                  dXZ;                        /**< xxxx */
		double                  dWY;                        /**< xxxx */
		double                  dYZ;                        /**< xxxx */
		double                  dWZ;                        /**< xxxx */

		double                  propMeanX;                  /**< xxxx */
		double                  propMeanY;                  /**< xxxx */
		double                  propMeanZ;                  /**< xxxx */
		double                  propMeanW;                  /**< xxxx */
		double                  propMeanInternal;           /**< xxxx */
		double                  prev_ln_prior;	            /**< The log prior of the starting state */
		double                  prev_ln_like;	            /**< The log likelihood of the starting state */
		CondLikelihoodShPtr     pre_root_posterior;         /**< xxxx */
		CondLikelihoodShPtr     pre_cla;                    /**< xxxx */
		double * * *            pre_p_mat;                  /**< xxxx */
		CondLikelihoodShPtr     post_root_posterior;        /**< xxxx */
		CondLikelihoodShPtr     post_cla;                   /**< xxxx */
		double * * *            post_p_mat;                 /**< xxxx */
		bool                    scoringBeforeMove;          /**< xxxx */
	};

} // namespace phycas

#include "phycas/src/unimap_nni_move.inl"

#endif
