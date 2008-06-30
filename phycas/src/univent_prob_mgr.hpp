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
#if ! defined(UNIVENT_PROB_MGR_HPP)
#define UNIVENT_PROB_MGR_HPP

#include <vector>

#if POLPY_NEWWAY
#include "phycas/src/cond_likelihood.hpp"
#include "phycas/src/square_matrix.hpp"
#include "phycas/src/univents.hpp"
#include "phycas/src/likelihood_models.hpp"


namespace phycas
{
class Lot;
class Tree; 
class UniventProbMgr
{
	public:
		UniventProbMgr(ModelShPtr);
		
		void sampleDescendantStates(Univents & u, const double * const * p_mat, const LikeFltType * des_cla, const int8_t * parent_states, Lot & rng) const
			{
			PHYCAS_ASSERT(u.size() == u.end_states_vec.size());
			u.is_valid = false;
			sampleDescendantStates(u.size(), &u.end_states_vec[0], p_mat, des_cla, parent_states, rng);
			}

		void sampleRootStates(Univents & u, const LikeFltType * rootStatePosterior, Lot & rng, bool posteriors_normalized, unsigned * obs_state_counts = NULL) const
			{
			PHYCAS_ASSERT(u.size() == u.end_states_vec.size());
			u.is_valid = false;
			sampleRootStates(u.size(), &u.end_states_vec[0], rootStatePosterior, rng, posteriors_normalized, obs_state_counts);
			}

		void sampleDescendantStates(const unsigned num_patterns, int8_t * nd_states, const double * const * p_mat, const LikeFltType * des_cla, const int8_t * parent_states, Lot & rng) const;
		void sampleRootStates(const unsigned num_patterns, int8_t * nd_states, const LikeFltType * rootStatePosterior, Lot & rng, bool posteriors_normalized,  unsigned * obs_state_counts = NULL) const;

		void sampleUniventsKeepEndStates(Univents & u, const double edgelen, const int8_t * par_states, const double * * p_mat_transposed, Lot & rng) const;
		void sampleUnivents(Univents & u,  const double edgelen, const int8_t * par_states, const double * const * p_mat, Lot & rng, unsigned ** s_mat) const
			{
			const std::vector<int8_t> & u_end_states_vec = u.getEndStatesVecRef();
			const int8_t * u_states_ptr = &u_end_states_vec[0];
			sampleUnivents(u, edgelen, par_states, u_states_ptr,  p_mat, rng, s_mat);
			}
		void sampleUnivents(Univents & u, const double edgelen, const int8_t * par_states, const int8_t * des_states, const double * const * p_mat, Lot & rng, unsigned ** s_mat) const;
		
		
		const double * const *			getUMatConst(unsigned m) const;
		//double * *					getUMat(unsigned m);
		void							recalcUMat();
		double							calcUnimapLnL(const Tree & t, const unsigned num_patterns, const unsigned * obs_state_counts, const unsigned * const  * sMat);
		
		
	private:
		void unimapEdgeOneSite(Univents &u, unsigned index, int8_t start_state, int8_t end_state, double transition_prob, double edgelen, bool sampleTimes, Lot & rng) const;
		unsigned sampleM(int8_t start_state, int8_t end_state, double transition_prob, double edgelen, Lot & rng) const;

		void							recalcUMatVect() const;
		void							expandUMatVect(unsigned) const;

		SquareMatrix					lnUMatMemMgt; /* deletes lnUMat when it goes out of scope */
		mutable std::vector<SquareMatrix>	uMatVect;				/**< uMat[m][i][j] is the (i,j)th element of the matrix obtained by raising the uniformized transition matrix to the power m */				  
		double * *						lnUMat;				/**< lnUMat[i][j] is the log of the corresponding element of the uniformized transition matrix element */
		mutable std::vector<double>		logmfact;				/**< contains the log of the factorial of m for m = 0..maxm */
		mutable double					lambda;					/**< rate at which uniformization events, or univents, occur */
		mutable unsigned 				maxm;
		unsigned						numStates;
		ModelShPtr						model;
		bool							sampleTimes;
		mutable SquareMatrix			scratchMatOne;
		mutable SquareMatrix			scratchMatTwo;
		bool							storeUnivents;
};

/*----------------------------------------------------------------------------------------------------------------------
|   Returns requested element of uMatVect.
*/
inline const double * const *  UniventProbMgr::getUMatConst(unsigned m) const
    {
    if (m >= uMatVect.size())
    	expandUMatVect(m);
    return const_cast<const double * const *>( uMatVect[m].GetMatrix());
    }

} // phycas


#endif //POLPY_NEWWAY
#endif
