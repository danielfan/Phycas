#ifndef PHO_GIBBSPARAMETER_H
#define PHO_GIBBSPARAMETER_H

#include <boost/shared_ptr.hpp>
#include "phycas/modules/mcmc/slice_sampler.hpp"

class HKYAdHocEvaluator;
typedef boost::shared_ptr<HKYAdHocEvaluator>		HKYAdHocEvaluatorShPtr;
typedef boost::shared_ptr<ProbabilityDistribution>	ProbabilityDistributionShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	A parameter that can update itself by sampling from its full conditional distribution.
*/
class GibbsParameter : public SliceSampler
	{
	public:
						GibbsParameter(double v);
		virtual			~GibbsParameter(){}
		
		virtual double	CalcLnProb(double v);	// overrides SliceSampler pure virtual function

		void			SetModel(HKYAdHocEvaluatorShPtr m);
		void			SetPriorDistr(ProbabilityDistributionShPtr p);
		double			GetLastSampledLnL();

		void			IgnoreData(bool b = true);

		void			SetValue(double v);
		double			GetValue();

	protected:

		double							val;			/**< The numerical value of this parameter. */
		double							lnL;			/**< The log-likelihood corresponding to `val'. */
		ProbabilityDistributionShPtr	priorDistr;		/**< Shared pointer to the prior probability distribution. */
		HKYAdHocEvaluatorShPtr			model;			/**< Shared pointer to the model and analysis manager. */
		bool							ignoring_data;	/**< normally false, but if true, posterior equals prior */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Calls base class (SliceSampler) constructor, initializes `lnL' to 0.0 and `val' to the supplied value.
*/
inline GibbsParameter::GibbsParameter(
  double v)	/**< is the numerical value with which to initialize this parameter */
  : SliceSampler(), val(v), lnL(0.0)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Used to cause the slice sampler to explore the prior (if `b' is true) or the posterior (if `b' is false). Simply
|	sets the underlying data member `ignoring_data' to `b'.
*/
inline void GibbsParameter::IgnoreData(
  bool b)	/**< is the value to which the data member `ignoring_data' should be set */
	{
	ignoring_data = b;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets data member `val' to the supplied value `v'. Calls SliceSampler::SetXValue as well.
*/
inline void GibbsParameter::SetValue(
  double v)	/**< is the new parameter value */
	{
	val = v; 
	SetXValue(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of data member `val'.
*/
inline double GibbsParameter::GetValue() 
	{
	return val;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `lnL' data member, which is set each time CalcLnProb is called.
*/
inline double GibbsParameter::GetLastSampledLnL()
	{
	return lnL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `priorDistr' data member to supplied shared pointer.
*/
inline void GibbsParameter::SetPriorDistr(ProbabilityDistributionShPtr p)
	{
	assert(p);
	priorDistr = p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `model' data member to supplied shared pointer.
*/
inline void GibbsParameter::SetModel(HKYAdHocEvaluatorShPtr m)
	{
	assert(m);
	model = m;
	}

#endif
