/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2008 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
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
#if ! defined(UNIVENT_PROB_MGR_HPP)
#define UNIVENT_PROB_MGR_HPP

#include <vector>

#include "phycas/src/cond_likelihood.hpp"
#include "phycas/src/square_matrix.hpp"
#include "phycas/src/univents.hpp"
#include "phycas/src/likelihood_models.hpp"

namespace phycas
{
class Lot;
class Tree; 

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
class UniventProbMgr
    {
	public:
		                                    UniventProbMgr(ModelShPtr);
		
		void                                sampleDescendantStates(Univents & u, const double * const * p_mat, const LikeFltType * des_cla, const int8_t * parent_states, Lot & rng) const;
		void                                sampleRootStates(Univents & u, const LikeFltType * rootStatePosterior, Lot & rng, bool posteriors_normalized, unsigned * obs_state_counts = NULL) const;

		void                                sampleDescendantStatesImpl(const unsigned num_patterns, int8_t * nd_states, const double * const * p_mat, const LikeFltType * des_cla, const int8_t * parent_states, Lot & rng) const;
		void                                sampleRootStatesImpl(const unsigned num_patterns, int8_t * nd_states, const LikeFltType * rootStatePosterior, Lot & rng, bool posteriors_normalized,  unsigned * obs_state_counts = NULL) const;

		void                                sampleUniventsKeepEndStates(Univents & u, const double edgelen, const int8_t * par_states, const double * * p_mat_transposed, Lot & rng) const;
		void                                sampleUnivents(Univents & u,  const double edgelen, const int8_t * par_states, const double * const * p_mat, Lot & rng, unsigned ** s_mat) const;
		void                                sampleUniventsImpl(Univents & u, const double edgelen, const int8_t * par_states, const int8_t * des_states, const double * const * p_mat, Lot & rng, unsigned ** s_mat) const;
		
		const double * const *			    getUMatConst(unsigned m) const;
		//double * *					    getUMat(unsigned m);
		void							    recalcUMat();
		double							    calcUnimapLnL(const Tree & t, const unsigned num_patterns, const unsigned * obs_state_counts, const unsigned * const  * sMat);
		
		bool								isMappingValid() const {return isMappingValidVar;}
		void								setIsMappingValid(bool v) {isMappingValidVar = v;}
	private:

		void                                unimapEdgeOneSite(Univents &u, unsigned index, int8_t start_state, int8_t end_state, double transition_prob, double edgelen, bool sampleTimes, Lot & rng) const;
		unsigned                            sampleM(int8_t start_state, int8_t end_state, double transition_prob, double edgelen, Lot & rng) const;

		void							    recalcUMatVect() const;
		void							    expandUMatVect(unsigned) const;

		SquareMatrix					    lnUMatMemMgt;   /**< deletes lnUMat when it goes out of scope */
		mutable std::vector<SquareMatrix>	uMatVect;       /**< uMat[m][i][j] is the (i,j)th element of the matrix obtained by raising the uniformized transition matrix to the power m */				  
		double * *						    lnUMat;         /**< lnUMat[i][j] is the log of the corresponding element of the uniformized transition matrix element */
		mutable std::vector<double>		    logmfact;       /**< contains the log of the factorial of m for m = 0..maxm */
		mutable double					    lambda;	        /**< rate at which uniformization events, or univents, occur */
		mutable unsigned 				    maxm;           /**< */
		unsigned						    numStates;      /**< */
		ModelShPtr						    model;          /**< */
		bool							    sampleTimes;    /**< */
		mutable SquareMatrix			    scratchMatOne;  /**< */
		mutable SquareMatrix			    scratchMatTwo;  /**< */
		bool							    storeUnivents;  /**< */
		bool								isMappingValidVar;
    };

} // phycas namespace

#include "phycas/src/univent_prob_mgr.inl"

#endif
