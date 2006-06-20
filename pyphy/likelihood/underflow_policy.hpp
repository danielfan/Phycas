#if !defined(UNDERFLOW_POLICY_HPP)
#define UNDERFLOW_POLICY_HPP

#include "pyphy/common/states_patterns.hpp"
#include "pyphy/likelihood/cond_likelihood.hpp"
#include "pyphy/likelihood/cond_likelihood_storage.hpp"

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

	private:

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
	};

} // namespace phycas

#include "pyphy/likelihood/underflow_policy.inl"

#endif
