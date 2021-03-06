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

#if ! defined(UNIMAP_FAST_NNI_MOVE_HPP)
#define UNIMAP_FAST_NNI_MOVE_HPP

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
class UnimapFastNNIMove : public MCMCUpdater
	{
	public:
								UnimapFastNNIMove();
								virtual ~UnimapFastNNIMove();

		// These are the virtual functions in the MCMCUpdater base class that we do not overload here
		//
		virtual double			getLnHastingsRatio() const;
		virtual double			getLnJacobian() const;
		virtual void			accept();
		virtual bool            computesTopologyPrior() const {return true;}

		// These are virtual functions in the MCMCUpdater base class that should be overloaded
		//
		virtual bool			update();
		virtual void			proposeNewState();
		virtual void			revert();
		virtual void			setLot(LotShPtr p);

	protected:
		double 					calcFourTaxonLogLikelihood();
		void 					storeOrigEdgeInfo();
		double 					calcEdgeLenLnPrior(const TreeNode & x, double edge_len, ChainManagerShPtr & chain_mgr);
		void 					addUniventsOneEdge(SquareMatrix & smat, TreeNode * nd);
		TreeNode * 				randomInternalAboveSubroot();
		TipData *				createTipDataFromUnivents(const Univents &, TipData *);
		void					DebugSaveNexusFile(TipData * xtd, TipData * ytd, TipData * ztd, TipData * wtd, double lnlike);
		void 					debugCheckUnivents();
		
	protected:
		//      a  b
		//       \/
		//   c   x <- randomly-chosen non-subroot internal node
		//    \ /
		//     y
		//    /
		//   d
		TreeNode *				x;	/**< randomly-chosen internal node (parent is `y', children are `a' and `b') */
		TreeNode *				y;	/**< internal node whose two children are `c' and `x' */
		TreeNode *				a;	/**< one of the two children of `x' (other one is `b') */
		TreeNode *				b;	/**< one of the two children of `x' (other one is `a') */
		TreeNode *				c;	/**< one child of `y' (the other is `x') */
		TreeNode *				d;	/**< the parent of `y' (may be the tip node that roots the entire tree) */
		
		unsigned				which_case;					/**< keeps track of which of the 8 cases was used for the proposal (needed for reverting) */
		double					tuning_factor;				/**< determines boldness of move (larger values mean bolder moves); if m is proposed 3-edge segment before move, m' = m*exp{-tuning_factor*(u-0.5)} will be proposed new length, where u is a uniform(0,1) deviate */
		double					three_edge_length_ratio;	/**< ratio of proposed new length of 3-edge segment to current length of 3-edge segment */
		
		bool					a_is_top;					/**< if true, `x'-`a' is top segment of the 3-edge path; if false, `x'-`b' is the top segment */
		bool					c_is_bottom;				/**< if true, `y'-`c' is bottom segment of the 3-edge path; if false, `y'-`d' is the bottom segment */
		bool					sliding_x;					/**< if true, `x' slides to a random location along the 3-edge path; if false, `y' does the sliding */
		
		unsigned				num_states;
		unsigned				num_sites;
		unsigned				num_moves_attempted;
		unsigned				num_moves_accepted;
		
		double 					lnL_before;
		double					lnL_after;
		
		double 					ln_density_reverse_move;
		double					ln_density_forward_move;

		double					tlen_before;
		double					alen_before;
		double					blen_before;
		double					clen_before;
		double					xlen_before;
		double					ylen_before;
		
		unsigned				mdota_before;
		unsigned				mdotb_before;
		unsigned				mdotc_before;
		unsigned				mdotx_before;
		unsigned				mdoty_before;
		
		SquareMatrix			smat_before;
		SquareMatrix			smat_after;
		SquareMatrix			log_umat_caretaker;
		double * * 				log_umat;
		
	};

} // namespace phycas

#include "phycas/src/unimap_fast_nni_move.inl"

#endif
