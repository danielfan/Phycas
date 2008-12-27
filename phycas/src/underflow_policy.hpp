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

#if !defined(UNDERFLOW_POLICY_HPP)
#define UNDERFLOW_POLICY_HPP

#include <cmath>
#include "phycas/src/states_patterns.hpp"
#include "phycas/src/cond_likelihood.hpp"
#include "phycas/src/cond_likelihood_storage.hpp"

//POL The original intention was to create several underflow policy classes in this file, then 
// make TreeLikelihood a template taking one of these policy classes. It was more complicated than
// it first appeared, however (more than TreeLikelihood would need to be templatized, and then there
// was the issue of TreeLikelihoodShPtr, which can be defined to only point to one instantiation!), 
// so the templatization step has been postponed. For now, TreeLikelihood has a SimpleUnderflowPolicy
// data member, which could be replaced by one of the others, but that would require recompliation.

//@POL I was surprised to discover that NaiveUnderflowPolicy, which has only empty member functions,
// was not nearly as fast as a version of the program that deletes all calls to the underflow policy 
// object. I thought that the inlined member function bodies would be optimized out entirely, but 
// I guess my intuition was wrong again.

namespace phycas
{


/*----------------------------------------------------------------------------------------------------------------------
|	Underflow policy that provides no underflow protection. Useful for small problems where speed is critical.
*/
class NaiveUnderflowPolicy
	{
	public:
									NaiveUnderflowPolicy() {}

		void						setTriggerSensitivity(unsigned nedges) {}
		void						setCorrectToValue(double maxval) {}
		void						setDimensions(unsigned np, unsigned nr, unsigned ns) {}

		void						twoTips(CondLikelihood & cond_like) const {}
		void						oneTip(CondLikelihood & cond_like, const CondLikelihood & other_cond_like, const CountVectorType & counts) const {}
		void						noTips(CondLikelihood & cond_like, const CondLikelihood & left_cond_like, const CondLikelihood & right_cond_like, const CountVectorType & counts) const {}

		void						correctSiteLike(double & site_like, unsigned pat, ConstCondLikelihoodShPtr condlike_shptr) const {}
		void						correctLnLike(double & ln_like, ConstCondLikelihoodShPtr condlike_shptr) const {}
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Default underflow policy for use with TreeLikelihood class.
*/
class SimpleUnderflowPolicy
	{
	public:
									SimpleUnderflowPolicy() {}

		void						setTriggerSensitivity(unsigned nedges);
		void						setCorrectToValue(double maxval);
		void						setDimensions(unsigned np, unsigned nr, unsigned ns);

		void						twoTips(CondLikelihood & cond_like) const;
		void						oneTip(CondLikelihood & cond_like, const CondLikelihood & other_cond_like, const CountVectorType & counts) const;
		void						noTips(CondLikelihood & cond_like, const CondLikelihood & left_cond_like, const CondLikelihood & right_cond_like, const CountVectorType & counts) const;

		void						correctSiteLike(double & site_like, unsigned pat, ConstCondLikelihoodShPtr condlike_shptr) const;
		void						correctLnLike(double & ln_like, ConstCondLikelihoodShPtr condlike_shptr) const;

	protected:

		unsigned					num_rates;				/**< Number of among-site rate categories */
		unsigned					num_patterns;			/**< Number of data patterns */
		unsigned					num_states;				/**< Number of states */
		mutable std::vector<double>	underflow_work;			/**< Workspace used when correcting for underflow (will have length equal to num_patterns) */
		unsigned					underflow_num_edges;	/**< Number of edges to traverse before underflow risk is evaluated */
		double						underflow_max_value;	/**< Maximum of the `num_states' conditional likelihoods for a given rate and pattern after underflow correction */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Underflow policy for use with TreeLikelihood class that keeps track of underflow correction factors for each data
|	pattern. This policy would only be needed if individual site likelihoods are needed. It adds a vector the length of
|	which is the number of site patterns to each CondLikelihood object. In general, the default SimpleUnderflowPolicy 
|	should suffice.
*/
class PatternSpecificUnderflowPolicy : public SimpleUnderflowPolicy
	{
	public:
									PatternSpecificUnderflowPolicy() : SimpleUnderflowPolicy() {}

		void						twoTips(CondLikelihood & cond_like) const;
		void						oneTip(CondLikelihood & cond_like, const CondLikelihood & other_cond_like, const CountVectorType & counts) const;
		void						noTips(CondLikelihood & cond_like, const CondLikelihood & left_cond_like, const CondLikelihood & right_cond_like, const CountVectorType & counts) const;

		void						correctSiteLike(double & site_like, unsigned pat, ConstCondLikelihoodShPtr condlike_shptr) const;
		void						correctLnLike(double & ln_like, ConstCondLikelihoodShPtr condlike_shptr) const;
	};

} // namespace phycas

#include "phycas/src/underflow_policy.inl"

#endif
