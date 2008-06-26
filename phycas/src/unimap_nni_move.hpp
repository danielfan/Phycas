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
StateTimeListVect * GetStateTimeListVect(TreeNode * nd);

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
		virtual bool	update();
		virtual double	getLnHastingsRatio() const;
		virtual double	getLnJacobian() const;
		virtual void	proposeNewState();
		virtual void	revert();
		virtual void	accept();
		virtual void	setLot(LotShPtr p);
	protected:
		TreeNode * randomInternalAboveSubroot();
		TreeNode * x;
		TreeNode * y;
		TreeNode * wSis;
		TreeNode * ySis;
		TreeNode * z;

		TipData * ySisTipData;
		TipData * yTipData;
		TipData * wSisTipData;
		TipData * wTipData;
		
		double  * * pre_x_pmat_transposed;
		double  * * pre_y_pmat_transposed;
		double  * * pre_w_pmat_transposed;
		double  * * pre_z_pmat_transposed;

		GammaDistribution gammaDist; 
		double min_edge_len_mean, edge_len_prop_cv;
		double ln_density_reverse_move, ln_density_forward_move;
		double prev_x_len, prev_y_len, prev_z_len, prev_nd_len, prev_ndP_len;
		bool x_is_left;
		double dXY, dWX, dXZ, dWY, dYZ, dWZ;
		double propMeanX, propMeanY, propMeanZ, propMeanW, propMeanInternal; 
		double prev_ln_prior;	/**< The log prior of the starting state */
		double prev_ln_like;	/**< The log likelihood of the starting state */
		CondLikelihoodShPtr pre_root_posterior;
		CondLikelihoodShPtr pre_cla;
		double * * * pre_p_mat;
		CondLikelihoodShPtr post_root_posterior;
		CondLikelihoodShPtr post_cla;
		double * * * post_p_mat;
		bool scoringBeforeMove;
		void calculatePairwiseDistances();
		void calculateProposalDist(bool);

		double calcProposalLnDensity(double mean, double x);
		double proposeEdgeLen(double mean);

		void FillStateCodeArray(const StateTimeListVect * um, int8_t * tipSpecificStateCode, bool);
		TipData * createTipDataFromUnivents(TreeNode * nd , bool use_last, TipData *);
		double FourTaxonLnLBeforeMove(TreeNode * nd);
		double FourTaxonLnLFromCorrectTipDataMembers(TreeNode * nd);
		double HarvestLnLikeFromCondLikePar(CondLikelihoodShPtr focalCondLike, ConstCondLikelihoodShPtr neighborCondLike, const double * const * childPMatrix);
		void storePMatTransposed(double **& cached, const double *** p_mat_array);
        void DebugSaveNexusFile(TipData * xtd, TipData * ytd, TipData * ztd, TipData * wtd, double lnlike);

		void sampleDescendantStates(const unsigned num_patterns, int8_t * nd_states, const double ** p_mat, const LikeFltType * des_cla, const int8_t * parent_states);
		void sampleRootStates(const unsigned num_patterns, int8_t * nd_states, LikeFltType * rootStatePosterior);
		void sampleUniventsKeepEndStates(TreeNode * nd, const int8_t * par_states, const double * * p_mat_transposed);
		void sampleUnivents(TreeNode * nd, const int8_t * par_states, const int8_t * des_states, const double * * p_mat);
		void sampleUniventsKeepBegStates(TreeNode * nd, const int8_t * des_states, const double * * p_mat_transposed);
	};

} // namespace phycas

#include "phycas/src/unimap_nni_move.inl"

#endif
