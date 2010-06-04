/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
|																			  |
|  This program is free software; you can redistribute it and/or modify		  |
|  it under the terms of the GNU General Public License as published by		  |
|  the Free Software Foundation; either version 2 of the License, or		  |
|  (at your option) any later version.										  |
|																			  |
|  This program is distributed in the hope that it will be useful,			  |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of			  |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the			  |
|  GNU General Public License for more details.								  |
|																			  |
|  You should have received a copy of the GNU General Public License along	  |
|  with this program; if not, write to the Free Software Foundation, Inc.,	  |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.				  |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#if ! defined(UNIMAP_NNI_MOVE_HPP)
#define UNIMAP_NNI_MOVE_HPP
#include <queue>							

#include <vector>							// for std::vector
#include <boost/shared_ptr.hpp>				// for boost::shared_ptr
#include <boost/weak_ptr.hpp>				// for boost::weak_ptr
#include "phycas/src/mcmc_updater.hpp"		// for base class MCMCUpdater
#include "phycas/src/mcmc_chain_manager.hpp"

namespace phycas
{

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager> ChainManagerWkPtr;



/*----------------------------------------------------------------------------------------------------------------------
|	
*/
class UnimapTopoMove : public MCMCUpdater
	{
	public:
								UnimapTopoMove(TreeLikeShPtr treeLike);
								virtual ~UnimapTopoMove();

		// These are the virtual functions in the MCMCUpdater base class that we do not overload in UnimapTopoMove
		//
		virtual double			getLnHastingsRatio() const = 0;
		virtual double			getLnJacobian() const = 0;
		virtual void			accept() = 0;

		// These are virtual functions in the MCMCUpdater base class
		//
		virtual bool			update();
		virtual void			proposeNewState();
		virtual void			revert();
		virtual void			setLot(LotShPtr p);
		virtual void			setTreeLikelihood(TreeLikeShPtr treeLike);


		void queueNode(TreeNode *nd) {ndQ.push(nd);}
	protected:

		void DebugSaveNexusFile(std::ostream & nxsf, double lnlike, unsigned subsetIndex);
		bool CheckWithPaup(double lnlike, unsigned);

		// ProposeStateWithTemporaries is called by proposeNewState and must be defined by subclasses
		virtual void			ProposeStateWithTemporaries(ChainManagerShPtr &) = 0;

		TipData *				createTipDataFromUnivents(const std::vector<Univents> &, TipData *);
		double					FourTaxonLnLBeforeMove();
		double					FourTaxonLnLFromCorrectTipDataMembers(unsigned subsetIndex);
		double					HarvestLnLikeFromCondLikePar(CondLikelihoodShPtr focalCondLike, ConstCondLikelihoodShPtr neighborCondLike, const double * const * childPMatrix, unsigned subsetIndex);
		void					storePMatTransposed(double **& cached, const double *** p_mat_array, unsigned subsetIndex);
		void					resampleInternalNodeStates(const LikeFltType * root_state_posterior, const LikeFltType * des_cla, unsigned subsetIndex);


		std::queue<TreeNode *>  ndQ;
		TreeNode *				origNode;
		TreeNode *				origNodePar; /**< */


		TreeNode *				x;							/**< xxxx */
		TreeNode *				y;							/**< xxxx */
		TreeNode *				z;							/**< xxxx */
		TreeNode *				w;							/**< xxxx */

		TreeNode *				a;							/**< xxxx */
		TreeNode *				b;							/**< xxxx */
		TreeNode *				c;							/**< xxxx */
		TreeNode *				d;							/**< xxxx */
		TreeNode *				aLenNd;						/**< xxxx */
		TreeNode *				bLenNd;						/**< xxxx */
		TreeNode *				cLenNd;						/**< xxxx */
		TreeNode *				dLenNd;						/**< xxxx */
		TipData *				aTipData;					/**< xxxx */
		TipData *				bTipData;					/**< xxxx */
		TipData *				cTipData;					/**< xxxx */
		TipData *				dTipData;					/**< xxxx */
		
		std::vector<double	* *	>			pre_x_pmat_transposed;		/**< xxxx */
		std::vector<double	* *	>			pre_y_pmat_transposed;		/**< xxxx */
		std::vector<double	* *	>			pre_w_pmat_transposed;		/**< xxxx */
		std::vector<double	* *	>			pre_z_pmat_transposed;		/**< xxxx */

		bool					doSampleInternalStates;

		bool					x_is_left;					/**< xxxx */
		double					prev_x_len;					/**< xxxx */
		double					prev_y_len;					/**< xxxx */
		double					prev_z_len;					/**< xxxx */
		double					prev_nd_len;				/**< xxxx */
		double					prev_ndP_len;				/**< xxxx */

		double					prev_ln_prior;				/**< The log prior of the starting state */
		double					prev_ln_like;				/**< The log likelihood of the starting state */
		CondLikelihoodShPtr	pre_root_posterior;			/**< xxxx */
		CondLikelihoodShPtr	pre_cla;					/**< xxxx */
		CondLikelihoodShPtr	post_root_posterior;		/**< xxxx */
		CondLikelihoodShPtr	post_cla;					/**< xxxx */
		std::vector<double * * *>			pre_p_mat;					/**< xxxx */
		std::vector<double * * *>			post_p_mat;					/**< xxxx */
		bool					scoringBeforeMove;			/**< xxxx */
		
		bool isFirstTime;
		double cached_posterior;

	};

class UnimapNNIMove : public UnimapTopoMove
	{
	public:
								UnimapNNIMove(TreeLikeShPtr treeLike);
								virtual ~UnimapNNIMove() {}

		// These are virtual functions in the UnimapTopoMove base class
		//
		virtual double			getLnHastingsRatio() const;
		virtual double			getLnJacobian() const;
		virtual void			accept();

		virtual void			setLot(LotShPtr p);

	protected:
		virtual void			ProposeStateWithTemporaries(ChainManagerShPtr &);

		void					calculatePairwiseDistances();
		void					calculatePairwiseDistancesForSubset(unsigned subsetIndex);
		void					calculateProposalDist(bool);

		virtual double			calcProposalLnDensity(double mean, double x);
		double					proposeEdgeLen(double mean);

		GammaDistribution		gammaDist;					/**< xxxx */
		double					min_edge_len_mean;			/**< xxxx */
		double					edge_len_prop_cv;			/**< xxxx */
		double					ln_density_reverse_move;	/**< xxxx */
		double					ln_density_forward_move;	/**< xxxx */

		double					dXY;						/**< xxxx */
		double					dWX;						/**< xxxx */
		double					dXZ;						/**< xxxx */
		double					dWY;						/**< xxxx */
		double					dYZ;						/**< xxxx */
		double					dWZ;						/**< xxxx */

		double					propMeanX;					/**< xxxx */
		double					propMeanY;					/**< xxxx */
		double					propMeanZ;					/**< xxxx */
		double					propMeanW;					/**< xxxx */
		double					propMeanInternal;			/**< xxxx */
	};


class UnimapLargetSimonMove : public UnimapTopoMove
	{
	public:
								UnimapLargetSimonMove(TreeLikeShPtr treeLike);
								virtual ~UnimapLargetSimonMove() {}

		// These are virtual functions in the UnimapTopoMove base class
		//
		virtual double			getLnHastingsRatio() const;
		virtual double			getLnJacobian() const;
		virtual void			accept();

		virtual void			setLot(LotShPtr p);
		void setLambda(double v) {lambda = v;}
		double getLambda() const {return lambda;}
	protected:
		virtual void			ProposeStateWithTemporaries(ChainManagerShPtr &);

		double lambda; /* lambda in LS-terminology */
		bool wOnThreeEdgePath;
		double expandContractFactor;
		bool detachUpperNode;
		bool swapYwithZ;
		bool topoChanging;
	};
	

class UnimapTopoMoveSpreader: public MCMCUpdater
	{
	public:
		UnimapTopoMoveSpreader() {}
		void addTopoMoveSpreader(UnimapTopoMove &m) {topoMoves.insert(&m);}
		void debugShowSelectedSubtrees();
		
		bool update();
	private:
		std::set<UnimapTopoMove *> topoMoves;
		std::set<TreeNode*> selectedNodes;
		std::vector<TreeNode *> queued;
		bool conflictsWithPrevious(TreeNode *) const ;
	};

} // namespace phycas

#include "phycas/src/unimap_nni_move.inl"

#endif
