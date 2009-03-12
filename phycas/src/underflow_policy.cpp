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

#include "phycas/src/underflow_policy.hpp"
#if POLPY_NEWWAY
#	include <fstream>
#	include <vector>
#endif

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Handles case of an internal node (`cond_like') with two internal node children (`left_cond_like' not equal to 
|	`right_cond_like') or the case of an internal node (`cond_like') with one tip child and one internal node child
|	(in which case `left_cond_like' will equal `right_cond_like'). The `numEdgesSinceUnderflowProtection' data member of
|	`cond_like' is set to 2 plus the number of traversed edges stored by `left_cond_like' plus (if there are two
|	internal node children) the number of traversed edges stored by `right_cond_like'. If this new number of edges
|	tranversed is greater than or equal to `underflow_num_edges', `cond_like' is corrected for underflow.
*/
void SimpleUnderflowPolicy::check(
  CondLikelihood & cond_like,				/**< the conditional likelihood array object of the focal internal node */
  const CondLikelihood & left_cond_like,	/**< the conditional likelihood array object of an internal node that is one immediate descendant of the focal node */ 
  const CondLikelihood & right_cond_like, 	/**< the conditional likelihood array object of an internal node that is the other immediate descendant of the focal node (if one descendant is a tip, left_cond_like and right_cond_like should refer to the same object) */
  const CountVectorType & counts,			/**< the pattern counts vector */
  bool polytomy)							/**< true if in the process of dealing with additional children (beyond first two) in a polytomy */
  const
	{
    if (num_patterns == 0)
        return;
        
    // Determine whether we are dealing with 
    // case 1: one tip child and one internal node child
    // case 2: two internal node children (= no_tips)
    // case 3: an extra tip (can only happen in case of polytomy)
    // case 4: an extra internal node (can only happen in case of polytomy)
    unsigned which = 0;
    const unsigned char one_tip = 1;
    const unsigned char no_tips = 2;
    const unsigned char xtra_tip = 3;
    const unsigned char xtra_internal = 4;
    if (&left_cond_like == &right_cond_like)
    	which = (polytomy ? xtra_tip : one_tip);
    else
    	which = (polytomy ? xtra_internal : no_tips);
    		
#if POLPY_NEWWAY
	std::ofstream doof("condlikes.txt", std::ios::app);
	double uf_sum_before = 0.0;
	double uf_sum_after = 0.0;
	double uf_sum_tmp = 0.0;
	std::vector<double> before;
	std::vector<double> after;
#endif

	unsigned nedges = 0;
	if (which == one_tip)
		{
		nedges += 2 + left_cond_like.getUnderflowNumEdges();
		}
	else if (which == no_tips)
		{
		nedges += 2 + left_cond_like.getUnderflowNumEdges() + right_cond_like.getUnderflowNumEdges();
		}
	else if (which == xtra_tip)
		{
		nedges += 1 + cond_like.getUnderflowNumEdges();
		}
	else	// xtra_internal
		{
		nedges += 1 + cond_like.getUnderflowNumEdges() + right_cond_like.getUnderflowNumEdges();
		}
			
	UnderflowType & uf_sum = cond_like.getUFSumRef();
	bool do_correction = (nedges >= underflow_num_edges);
	if (do_correction)
		{
		// Ok, we've traversed enough edges that it is time to take another factor out for underflow control
		
		// If we're evaluating the 3rd, 4th, etc., child in a polytomy (i.e. polytomy is true)
		// then cond_like's uf_sum data member already has information passed along from its first
		// two children that we do not want to erase
		if (!polytomy)
			{
			uf_sum = 0;
			}
			
//#error Need to take account of TreeLikelihood::pinvar_like here

		// Begin by finding the largest conditional likelihood for each pattern
		// over all rates and states (store these in underflow_work vector)
		underflow_work.resize(num_patterns*num_states);
		underflow_work.assign(num_patterns*num_states, 0.0);
		LikeFltType * cla = cond_like.getCLA();
		for (unsigned r = 0; r < num_rates; ++r)
			{
			for (unsigned pat = 0; pat < num_patterns; ++pat)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					double curr = *cla++;
#if POLPY_NEWWAY
					before.push_back(curr);
#endif
					if (curr > underflow_work[pat])
						underflow_work[pat] = curr;
					}
				}
			}

		// Assume that x is the largest of the num_rates*num_states conditional likelihoods
		// for a given pattern. Find factor g such that g*x = underflow_max_value. Suppose
		// x = 0.05 and underflow_max_value = 10000, g = 10000/0.05 = 200,000 = e^{12.206}
		// In this case, we would add the integer component of the exponent (i.e. 12) to
		// the underflow correction for this pattern and adjust all conditional likelihoods 
		// by a factor f = e^{12} = 22026.46579. So, the conditional likelihood 0.05 now
		// becomes 0.05*e^{12} = 1101.32329. The value f*count is added to uf_sum so that
		// it can be removed from the final log-likelihood.
		cla = cond_like.getCLA();
		unsigned condlikes_per_rate = num_patterns*num_states;
		for (unsigned pat = 0; pat < num_patterns; ++pat)
			{
			PHYCAS_ASSERT(underflow_work[pat] > 0.0);
			double ratio = underflow_max_value/underflow_work[pat];
			double log_ratio = std::log(ratio);
			double f = std::floor(log_ratio);
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
			
#if POLPY_NEWWAY
		cla = cond_like.getCLA();
		for (unsigned r = 0; r < num_rates; ++r)
			{
			for (unsigned pat = 0; pat < num_patterns; ++pat)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					double curr = *cla++;
					after.push_back(curr);
					}
				}
			}
		uf_sum_tmp = uf_sum;
#endif

		// uf_sum now takes account of the factor that we just removed, but doesn't yet
		// account for the factors that we've taken out previously (in nodes above this node)
		if (which == one_tip)
			{
			uf_sum += left_cond_like.getUFSum();
			}
		else if (which == no_tips)
			{
			uf_sum += left_cond_like.getUFSum();
			uf_sum += right_cond_like.getUFSum();
			}
		else if (which == xtra_tip)
			{
			// nothing above this point to take account of
			}
		else	// xtra_internal
			{
			uf_sum += right_cond_like.getUFSum();
			}

#if POLPY_NEWWAY
		uf_sum_after = uf_sum;
		uf_sum_before = uf_sum - uf_sum_tmp;
#endif
		
		nedges = 0;
		}
	else
		{
		// Have not yet traversed enough edges to trigger underflow correction, but still need
		// to carry the results of previous corrections down
		if (which == one_tip)
			{
			uf_sum = left_cond_like.getUFSum();
			}
		else if (which == no_tips)
			{
			uf_sum = left_cond_like.getUFSum();
			uf_sum += right_cond_like.getUFSum();
			}
		else if (which == xtra_tip)
			{
			// nothing above this point to take account of
			}
		else	// xtra_internal
			{
			uf_sum += right_cond_like.getUFSum();
			}
			
#if POLPY_NEWWAY
		uf_sum_after = uf_sum;
		uf_sum_before = uf_sum;
		LikeFltType * cla = cond_like.getCLA();
		for (unsigned r = 0; r < num_rates; ++r)
			{
			for (unsigned pat = 0; pat < num_patterns; ++pat)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					double curr = *cla++;
					before.push_back(curr);
					after.push_back(curr);
					}
				}
			}
		uf_sum_tmp = uf_sum;
#endif
		}
	cond_like.setUnderflowNumEdges(nedges);
	
#if POLPY_NEWWAY
	if (do_correction)
		{
		//	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
		//	|                     rate 0                    |                     rate 1                    | 
		//	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
		//	|   pattern 0   |   pattern 1   |   pattern 2   |   pattern 0   |   pattern 1   |   pattern 2   | 
		//	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
		//	| A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T |
		//	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
		// i  0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23 
		// p  0   0   0   0   1   1   1   1   2   2   2   2   0   0   0   0   1   1   1   1   2   2   2   2  
		doof << "\n----------\n";
		doof << "num_rates     = "  << num_rates << '\n';
		doof << "num_states    = "  << num_states << '\n';
		doof << "num_patterns  = "  << num_patterns << '\n';
		doof << "uf_sum_before = "  << uf_sum_before << '\n';
		doof << "uf_sum_after  = "  << uf_sum_after << '\n';
		unsigned nbefore = (unsigned)before.size();
		for (unsigned i = 0; i < nbefore; ++i)	
			{
			unsigned r = i/(num_patterns*num_states);
			unsigned p = (i/num_states) % num_patterns;
			double ratio = underflow_max_value/underflow_work[p];
			double log_ratio = std::log(ratio);
			double f = std::floor(log_ratio);
			doof << boost::str(boost::format("%6d %6d %6d %20.10f %20.10f %20.10f %20.10f\n") % i % r % p % before[i] % after[i] % underflow_work[p] % f);		
			}
		doof << '\n';
		}
	doof.close();
#endif
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns 0.0 because SimpleUnderflowPolicy does not store correction factors for individual patterns.
*/
double SimpleUnderflowPolicy::getCorrectionFactor(
  unsigned pat,								/**< is the index of the pattern for which the correction factor is desired */
  ConstCondLikelihoodShPtr condlike_shptr)	/**< is the conditional likelihood array of the likelihood root node */
  const
	{
	return 0.0;
	}	

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains a pointer to the underflow array for `cond_like', then returns the value stored in that underflow array 
|	for pattern `pat'.
*/
double PatternSpecificUnderflowPolicy::getCorrectionFactor(
  unsigned pat,								/**< is the index of the pattern for which the correction factor is desired */
  ConstCondLikelihoodShPtr condlike_shptr)	/**< is the conditional likelihood array of the likelihood root node */
  const
	{
	double site_uf_factor = 0.0;
    if (num_patterns > 0)
        {
	    PHYCAS_ASSERT(condlike_shptr);
	    UnderflowType const * uf = condlike_shptr->getUF();
	    PHYCAS_ASSERT(uf != NULL);
	    site_uf_factor = (double)uf[pat];
        }
	return site_uf_factor;
	}	

/*----------------------------------------------------------------------------------------------------------------------
|	Handles case of an internal node (`cond_like') with two internal node children (`left_cond_like' not equal to 
|	`right_cond_like') or the case of an internal node (`cond_like') with one tip child and one internal node child
|	(in which case `left_cond_like' will equal `right_cond_like'). The `numEdgesSinceUnderflowProtection' data member of
|	`cond_like' is set to 2 plus the number of traversed edges stored by `left_cond_like' plus (if there are two
|	internal node children) the number of traversed edges stored by `right_cond_like'. If this new number of edges
|	tranversed is greater than or equal to `underflow_num_edges', `cond_like' is corrected for underflow.
*/
void PatternSpecificUnderflowPolicy::check(
  CondLikelihood & cond_like,				/**< the conditional likelihood array object of the focal internal node */
  const CondLikelihood & left_cond_like,	/**< the conditional likelihood array object of an internal node that is one immediate descendant of the focal node */ 
  const CondLikelihood & right_cond_like, 	/**< the conditional likelihood array object of an internal node that is the other immediate descendant of the focal node (if one descendant is a tip, left_cond_like and right_cond_like should refer to the same object) */
  const CountVectorType & counts,			/**< the pattern counts vector */
  bool polytomy)							/**< true if in the process of dealing with additional children (beyond first two) in a polytomy */
  const
	{
    if (num_patterns == 0)
        return;

    // Determine whether we are dealing with 
    // case 1: one tip child and one internal node child
    // case 2: two internal node children (= no_tips)
    // case 3: an extra tip (can only happen in case of polytomy)
    // case 4: an extra internal node (can only happen in case of polytomy)
    unsigned which = 0;
    const unsigned char one_tip = 1;
    const unsigned char no_tips = 2;
    const unsigned char xtra_tip = 3;
    const unsigned char xtra_internal = 4;
    if (&left_cond_like == &right_cond_like)
    	which = (polytomy ? xtra_tip : one_tip);
    else
    	which = (polytomy ? xtra_internal : no_tips);
    		
	unsigned nedges = 0;
	if (which == one_tip)
		{
		nedges += 2 + left_cond_like.getUnderflowNumEdges();
		}
	else if (which == no_tips)
		{
		nedges += 2 + left_cond_like.getUnderflowNumEdges() + right_cond_like.getUnderflowNumEdges();
		}
	else if (which == xtra_tip)
		{
		nedges += 1 + cond_like.getUnderflowNumEdges();
		}
	else	// xtra_internal
		{
		nedges += 1 + cond_like.getUnderflowNumEdges() + right_cond_like.getUnderflowNumEdges();
		}
			
	bool do_correction = (nedges >= underflow_num_edges);
	if (do_correction)
		{
		// Ok, we've traversed enough edges that it is time to take another factor out for underflow control
		
//#error Need to take account of TreeLikelihood::pinvar_like here
		
		// Begin by finding the largest conditional likelihood for each pattern
		// over all rates and states (store these in underflow_work vector)
		underflow_work.resize(num_patterns*num_states);
		underflow_work.assign(num_patterns*num_states, 0.0);
		LikeFltType * cla = cond_like.getCLA();
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
		// for a given pattern. Find factor g such that g*x = underflow_max_value. Suppose
		// x = 0.05 and underflow_max_value = 10000, g = 10000/0.05 = 200,000 = e^{12.206}
		// In this case, we would add the integer component of the exponent (i.e. 12) to
		// the underflow correction for this pattern and adjust all conditional likelihoods 
		// by a factor f = e^{12} = 22026.46579. So, the conditional likelihood 0.05 now
		// becomes 0.05*e^{12} = 1101.32329. 
		cla = cond_like.getCLA();
		UnderflowType *	uf = cond_like.getUF();
		UnderflowType const * uf_left = left_cond_like.getUF();
		UnderflowType const * uf_right = right_cond_like.getUF();
		unsigned condlikes_per_rate = num_patterns*num_states;
		for (unsigned pat = 0; pat < num_patterns; ++pat)
			{
			PHYCAS_ASSERT(underflow_work[pat] > 0.0);
			double ratio = underflow_max_value/underflow_work[pat];
			double log_ratio = std::log(ratio);
			double f = std::floor(log_ratio);
			//double f = floor(log(underflow_max_value/underflow_work[pat]));
			double expf = exp(f);
			LikeFltType * claPat = &cla[num_states*pat];
			for (unsigned r = 0; r < num_rates; ++r, claPat += condlikes_per_rate)
				{
				for (unsigned i = 0; i < num_states; ++i)
					{
					claPat[i] *= expf;
					}
				}
			if (which == one_tip)
				{
				// uf_left is same object as uf_right in this case
				*uf = *uf_left + (UnderflowType)f;
				++uf_left;
				}
			else if (which == no_tips)
				{
				*uf = *uf_left + *uf_right + (UnderflowType)f;
				++uf_left;
				++uf_right;
				}
			else if (which == xtra_tip)
				{
				// nothing above this point to take account of
				}
			else	// xtra_internal
				{
				PHYCAS_ASSERT(0);	// THIS CASE NOT YET IMPLEMENTED
				}
			++uf;
			}
		nedges = 0;
		}
	else
		{
		UnderflowType *	uf = cond_like.getUF();
		UnderflowType const * uf_left = left_cond_like.getUF();
		UnderflowType const * uf_right = right_cond_like.getUF();
		for (unsigned pat = 0; pat < num_patterns; ++pat, ++uf)
			{
			if (which == one_tip)
				{
				// uf_left is same object as uf_right in this case
				*uf = *uf_left;
				++uf_left;
				}
			else if (which == no_tips)
				{
				*uf = *uf_left + *uf_right;
				++uf_left;
				++uf_right;
				}
			else if (which == xtra_tip)
				{
				// nothing above this point to take account of
				}
			else	// xtra_internal
				{
				PHYCAS_ASSERT(0);	// THIS CASE NOT YET IMPLEMENTED
				}
			}
		}
	cond_like.setUnderflowNumEdges(nedges);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Handles case of a node subtending two internal nodes. In this case, the number of edges traversed is set to 2 plus  
|	the number of traversed edges stored by `left_cond_like' plus the number of traversed edges stored by 
|	`right_cond_like'. If this new number of edges tranversed is greater than or equal to `underflow_num_edges', 
|	`cond_like' is corrected for underflow.
*/
#if 0	// this class will need some work if it is ever used again
void PatternSpecificUnderflowPolicy::noTips(
  CondLikelihood & cond_like,				/**< */
  const CondLikelihood & left_cond_like,	/**< */ 
  const CondLikelihood & right_cond_like, 	/**< */
  const CountVectorType & counts)			/**< */
  const
	{
    if (num_patterns == 0)
        return;

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
			PHYCAS_ASSERT(underflow_work[pat] > 0.0);
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
#endif

} // namespace phycas
