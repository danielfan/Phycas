/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2008 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#if ! defined(UNIVENT_PROB_MGR_INL)
#define UNIVENT_PROB_MGR_INL

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline void UniventProbMgr::sampleDescendantStates(Univents & u, const double * const * p_mat, const LikeFltType * des_cla, const int8_t * parent_states, Lot & rng) const
	{
	PHYCAS_ASSERT(u.size() == u.end_states_vec.size());
	u.setValid(false);
	sampleDescendantStatesImpl(u.size(), &u.end_states_vec[0], p_mat, des_cla, parent_states, rng);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline void UniventProbMgr::sampleRootStates(Univents & u, const LikeFltType * rootStatePosterior, Lot & rng, bool posteriors_normalized, unsigned * obs_state_counts) const
	{
	PHYCAS_ASSERT(u.size() == u.end_states_vec.size());
	u.setValid(false);
	sampleRootStatesImpl(u.size(), &u.end_states_vec[0], rootStatePosterior, rng, posteriors_normalized, obs_state_counts);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline void UniventProbMgr::sampleUnivents(Univents & u,  const double edgelen, const int8_t * par_states, const double * const * p_mat, Lot & rng, unsigned ** s_mat) const
	{
	const std::vector<int8_t> & u_end_states_vec = u.getEndStatesVecRef();
	const int8_t * u_states_ptr = &u_end_states_vec[0];
	sampleUniventsImpl(u, edgelen, par_states, u_states_ptr,  p_mat, rng, s_mat);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns requested element of uMatVect.
*/
inline const double * const *  UniventProbMgr::getUMatConst(unsigned m) const
    {
    if (m >= uMatVect.size())
    	expandUMatVect(m);
    return const_cast<const double * const *>( uMatVect[m].GetMatrixAsRawPointer());
    }

} // namespace phycas
#endif
