#include "pyphy/likelihood/underflow_policy.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Handles case of a node subtending one internal node and one tip node. In this case, the number of edges traversed 
|	is set to 2 plus the number traversed edges stored by `other_cond_like'. If this new number of edges tranversed is
|	greater than or equal to `underflow_num_edges', `cond_like' is corrected for underflow.
*/
void SimpleUnderflowPolicy::oneTip(
  CondLikelihood &			cond_like,			/**< */
  const CondLikelihood &	other_cond_like, 	/**< */
  const CountVectorType	&	counts) 			/**< */
  const
	{
	unsigned nedges = 2 + other_cond_like.getUnderflowNumEdges();
	UnderflowType & uf_sum			= cond_like.getUFSumRef();
	if (nedges >= underflow_num_edges)
		{
		underflow_work.resize(num_patterns*num_states);
		underflow_work.assign(num_patterns*num_states, 0.0);
		LikeFltType * cla = cond_like.getCLA();

		// Begin by finding the largest conditional likelihood for each pattern
		// over all rates and states (store these in underflow_work vector)
		for (unsigned r = 0; r < num_rates; ++r)
			{
			for (unsigned pat = 0; pat < num_patterns; ++pat)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					double curr = *cla++;
					if (curr > underflow_work[pat])
						underflow_work[pat] = curr;
					}
				}
			}

		// Assume that x is the largest of the num_rates*num_states conditional likelihoods
		// for a given pattern. Find factor f such that f*x = underflow_max_value. Suppose
		// x = 0.05 and underflow_max_value = 10,000, f = 10,000/0.05 = 200,000 = e^{12.206}
		// In this case, we would add the integer component of the exponent (i.e. 12) to
		// the underflow correction for this pattern and adjust all conditional likelihoods 
		// by a factor f* = e^{12} = 22026.46579. So, the conditional likelihood 0.05 now
		// becomes 0.05*e^{12} = 1101.32329.
		cla = cond_like.getCLA();
		uf_sum = 0;
		unsigned condlikes_per_rate = num_patterns*num_states;
		for (unsigned pat = 0; pat < num_patterns; ++pat)
			{
			assert(underflow_work[pat] > 0.0);
			double f = floor(log(underflow_max_value/underflow_work[pat]));
			double expf = exp(f);
			LikeFltType * claPat = &cla[num_states*pat];
			for (unsigned r = 0; r < num_rates; ++r, claPat += condlikes_per_rate)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					claPat[i] *= expf;
					}
				}
			UnderflowType curr = (UnderflowType)f;
			UnderflowType cnt = (UnderflowType)counts[pat];
			uf_sum += curr*cnt;
			}
		uf_sum += other_cond_like.getUFSum();
		nedges = 0;
		}
	else
		{
		uf_sum = other_cond_like.getUFSum();
		}

	cond_like.setUnderflowNumEdges(nedges);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Handles case of a node subtending one internal node and one tip node. In this case, the number of edges traversed 
|	is set to 2 plus the number traversed edges stored by `other_cond_like'. If this new number of edges tranversed is
|	greater than or equal to `underflow_num_edges', `cond_like' is corrected for underflow.
*/
void PatternSpecificUnderflowPolicy::oneTip(
  CondLikelihood &			cond_like,			/**< */
  const CondLikelihood &	other_cond_like, 	/**< */
  const CountVectorType	&	counts) 			/**< */
  const
	{
	unsigned nedges = 2 + other_cond_like.getUnderflowNumEdges();
	if (nedges >= underflow_num_edges)
		{
		underflow_work.resize(num_patterns*num_states);
		underflow_work.assign(num_patterns*num_states, 0.0);
		LikeFltType * cla = cond_like.getCLA();

		// Begin by finding the largest conditional likelihood for each pattern
		// over all rates and states (store these in underflow_work vector)
		for (unsigned r = 0; r < num_rates; ++r)
			{
			for (unsigned pat = 0; pat < num_patterns; ++pat)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					double curr = *cla++;
					if (curr > underflow_work[pat])
						underflow_work[pat] = curr;
					}
				}
			}

		// Assume that x is the largest of the num_rates*num_states conditional likelihoods
		// for a given pattern. Find factor f such that f*x = underflow_max_value. Suppose
		// x = 0.05 and underflow_max_value = 10,000, f = 10,000/0.05 = 200,000 = e^{12.206}
		// In this case, we would add the integer component of the exponent (i.e. 12) to
		// the underflow correction for this pattern and adjust all conditional likelihoods 
		// by a factor f* = e^{12} = 22026.46579. So, the conditional likelihood 0.05 now
		// becomes 0.05*e^{12} = 1101.32329.
		cla = cond_like.getCLA();
		UnderflowType *	uf = cond_like.getUF();
		UnderflowType const * uf_right = other_cond_like.getUF();
		unsigned condlikes_per_rate = num_patterns*num_states;
		for (unsigned pat = 0; pat < num_patterns; ++pat)
			{
			assert(underflow_work[pat] > 0.0);
			double f = floor(log(underflow_max_value/underflow_work[pat]));
			double expf = exp(f);
			LikeFltType * claPat = &cla[num_states*pat];
			for (unsigned r = 0; r < num_rates; ++r, claPat += condlikes_per_rate)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					claPat[i] *= expf;
					}
				}
			*uf = *uf_right + (unsigned)f;
			++uf;
			++uf_right;
			}
		nedges = 0;
		}
	else
		{
		UnderflowType *	uf = cond_like.getUF();
		UnderflowType const * uf_right = other_cond_like.getUF();
		for (unsigned pat = 0; pat < num_patterns; ++pat, ++uf, ++uf_right)
			*uf = *uf_right;
		}
	cond_like.setUnderflowNumEdges(nedges);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Handles case of a node subtending two internal nodes. In this case, the number of edges traversed is set to 2 plus  
|	the number of traversed edges stored by `left_cond_like' plus the number of traversed edges stored by 
|	`right_cond_like'. If this new number of edges tranversed is greater than or equal to `underflow_num_edges', 
|	`cond_like' is corrected for underflow.
*/
void SimpleUnderflowPolicy::noTips(
  CondLikelihood & cond_like,				/**< */
  const CondLikelihood & left_cond_like,	/**< */ 
  const CondLikelihood & right_cond_like, 	/**< */
  const CountVectorType & counts)			/**< */
  const
	{
	unsigned nedges = 2 + left_cond_like.getUnderflowNumEdges() + right_cond_like.getUnderflowNumEdges();
	UnderflowType & uf_sum = cond_like.getUFSumRef();
	if (nedges >= underflow_num_edges)
		{
		underflow_work.resize(num_patterns*num_states);
		underflow_work.assign(num_patterns*num_states, 0.0);
		LikeFltType * cla = cond_like.getCLA();

		// Begin by finding the largest conditional likelihood for each pattern
		// over all rates and states (store these in underflow_work vector)
		for (unsigned r = 0; r < num_rates; ++r)
			{
			for (unsigned pat = 0; pat < num_patterns; ++pat)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					double curr = *cla++;
					if (curr > underflow_work[pat])
						underflow_work[pat] = curr;
					}
				}
			}

		// Assume that x is the largest of the num_rates*num_states conditional likelihoods
		// for a given pattern. Find factor f such that f*x = underflow_max_value. Suppose
		// x = 0.05 and underflow_max_value = 10,000, f = 10,000/0.05 = 200,000 = e^{12.206}
		// In this case, we would add the integer component of the exponent (i.e. 12) to
		// the underflow correction for this pattern and adjust all conditional likelihoods 
		// by a factor f* = e^{12} = 22026.46579. So, the conditional likelihood 0.05 now
		// becomes 0.05*e^{12} = 1101.32329.
		cla = cond_like.getCLA();
		uf_sum = 0;
		unsigned condlikes_per_rate = num_patterns*num_states;
		for (unsigned pat = 0; pat < num_patterns; ++pat)
			{
			assert(underflow_work[pat] > 0.0);
			double f = floor(log(underflow_max_value/underflow_work[pat]));
			double expf = exp(f);
			LikeFltType * claPat = &cla[num_states*pat];
			for (unsigned r = 0; r < num_rates; ++r, claPat += condlikes_per_rate)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					claPat[i] *= expf;
					}
				}
			UnderflowType curr = (UnderflowType)f;
			UnderflowType cnt = (UnderflowType)counts[pat];
			uf_sum += curr*cnt;
			}
		uf_sum += left_cond_like.getUFSum();
		uf_sum += right_cond_like.getUFSum();
		nedges = 0;
		}
	else
		{
		uf_sum = left_cond_like.getUFSum();
		uf_sum += right_cond_like.getUFSum();
		}
	cond_like.setUnderflowNumEdges(nedges);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Handles case of a node subtending two internal nodes. In this case, the number of edges traversed is set to 2 plus  
|	the number of traversed edges stored by `left_cond_like' plus the number of traversed edges stored by 
|	`right_cond_like'. If this new number of edges tranversed is greater than or equal to `underflow_num_edges', 
|	`cond_like' is corrected for underflow.
*/
void PatternSpecificUnderflowPolicy::noTips(
  CondLikelihood & cond_like,				/**< */
  const CondLikelihood & left_cond_like,	/**< */ 
  const CondLikelihood & right_cond_like, 	/**< */
  const CountVectorType & counts)			/**< */
  const
	{
	unsigned nedges = 2 + left_cond_like.getUnderflowNumEdges() + right_cond_like.getUnderflowNumEdges();
	if (nedges >= underflow_num_edges)
		{
		underflow_work.resize(num_patterns*num_states);
		underflow_work.assign(num_patterns*num_states, 0.0);
		LikeFltType * cla = cond_like.getCLA();

		// Begin by finding the largest conditional likelihood for each pattern
		// over all rates and states (store these in underflow_work vector)
		for (unsigned r = 0; r < num_rates; ++r)
			{
			for (unsigned pat = 0; pat < num_patterns; ++pat)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					double curr = *cla++;
					if (curr > underflow_work[pat])
						underflow_work[pat] = curr;
					}
				}
			}

		// Assume that x is the largest of the num_rates*num_states conditional likelihoods
		// for a given pattern. Find factor f such that f*x = underflow_max_value. Suppose
		// x = 0.05 and underflow_max_value = 10,000, f = 10,000/0.05 = 200,000 = e^{12.206}
		// In this case, we would add the integer component of the exponent (i.e. 12) to
		// the underflow correction for this pattern and adjust all conditional likelihoods 
		// by a factor f* = e^{12} = 22026.46579. So, the conditional likelihood 0.05 now
		// becomes 0.05*e^{12} = 1101.32329.
		cla = cond_like.getCLA();
		UnderflowType *	uf = cond_like.getUF();
		UnderflowType const *	uf_left = left_cond_like.getUF();
		UnderflowType const *	uf_right = right_cond_like.getUF();
		unsigned condlikes_per_rate = num_patterns*num_states;
		for (unsigned pat = 0; pat < num_patterns; ++pat)
			{
			assert(underflow_work[pat] > 0.0);
			double f = floor(log(underflow_max_value/underflow_work[pat]));
			double expf = exp(f);
			LikeFltType * claPat = &cla[num_states*pat];
			for (unsigned r = 0; r < num_rates; ++r, claPat += condlikes_per_rate)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					claPat[i] *= expf;
					}
				}
			*uf = *uf_left + *uf_right + (unsigned)f;
			++uf;
			++uf_left;
			++uf_right;
			}
		nedges = 0;
		}
	else
		{
		UnderflowType *	uf = cond_like.getUF();
		UnderflowType const *	uf_left = left_cond_like.getUF();
		UnderflowType const *	uf_right = right_cond_like.getUF();
		for (unsigned pat = 0; pat < num_patterns; ++pat, ++uf, ++uf_left, ++uf_right)
			*uf = *uf_left + *uf_right;
		}
	cond_like.setUnderflowNumEdges(nedges);
	}

} // namespace phycas
