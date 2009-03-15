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

#include "phycas/src/underflow_manager.hpp"
#include <fstream>
#include <vector>

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor.
*/
UnderflowManager::UnderflowManager()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `underflow_max_value', which is the target value to which the largest 
|	conditional likelihood for each pattern is set.
*/
double UnderflowManager::getUnderflowMaxValue() const
	{
	return underflow_max_value;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of the data member `underflow_num_edges' to `nedges'. This is the number of edges to accumulate before 
|	correcting for underflow. Assumes `nedges' is greater than zero. Often several hundred taxa are required before
|	underflow becomes a problem, so a reasonable value for `underflow_num_edges' is 25.
*/
void UnderflowManager::setTriggerSensitivity(
  unsigned nedges)		/**< is the number of edges to traverse before correcting for underflow */
	{
	PHYCAS_ASSERT(nedges > 0);
	underflow_num_edges = nedges;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of data member `underflow_max_value' to `maxval'. The largest conditional likelihood for every pattern 
|	will be adjusted to a value close to this. Suppose the largest conditional likelihood for pattern i was x. A factor 
|	exp(f) is found such that x*exp(f) = `underflow_max_value'. The fractional component of f is then removed, yielding
|	an unsigned int k, which is what is actually stored. The largest conditional likelihood for pattern i after the
|	correction is thus x*exp(k), which will be less than or equal to `underflow_max_value'. Assumes `maxval' is greater
|	than zero. A reasonable value for `underflow_max_value' is 10000.
*/
void UnderflowManager::setCorrectToValue(
  double maxval)	/**< is the target value to which the largest conditional likelihood for any given pattern will be scaled */
	{
	PHYCAS_ASSERT(maxval > 0.0);
	underflow_max_value = maxval;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `num_rates', `num_patterns', and `num_states' data members to `nr', `np' and `ns', respectively. Assumes
|	that all three values are greater than zero.
*/
void UnderflowManager::setDimensions(
  unsigned np,	/**< is the number of patterns in the data */
  unsigned nr, 	/**< is the number of relative rates used in modeling among-site rate variation */
  unsigned ns)	/**< is the number of states */
	{
	// Note: num_patterns can legitimately be 0 if running with no data. In this case, no underflow
    // correction is ever needed, and all member functions simply return without doing anything
	PHYCAS_ASSERT(nr > 0);
	PHYCAS_ASSERT(ns > 0);
	num_patterns = np;
	num_rates = nr;
	num_states = ns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Handles case of a node subtending two tips. In this case, the number of edges traversed is just two, and thus no
|	corrective action is taken other than to set the number of underflow edges traversed to 2. The zeroUF member
|	function of the supplied `cond_like' object is called, however, to ensure that the cumulative correction factor is 
|	zero (this `cond_like' might have just been brought in from storage with some correction factor already in place).
*/
void UnderflowManager::twoTips(
  CondLikelihood & cond_like) /**< is the conditional likelihood array to correct */
  const
	{
    if (num_patterns > 0)
        {
	    cond_like.setUnderflowNumEdges(2);
	    cond_like.zeroUF();
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains a pointer to the underflow array for `cond_like', then subtracts the value stored in that underflow array 
|	for pattern `pat' from the site likelihood `site_like'.
*/
double UnderflowManager::correctSiteLike(
  double & site_like,						/**< is the log-site-likelihood to correct */
  unsigned pat,								/**< is the index of the pattern representing the site to correct */
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
	    site_like -= site_uf_factor;
        }
	return site_uf_factor;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains a pointer to the underflow array for `cond_like', then returns the value stored in that underflow array 
|	for pattern `pat'.
*/
double UnderflowManager::getCorrectionFactor(
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
void UnderflowManager::check(
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
		// cond_like, left_cond_like == right_cond_like
		nedges += 2 + left_cond_like.getUnderflowNumEdges();
		}
	else if (which == no_tips)
		{
		// cond_like, left_cond_like, right_cond_like
		nedges += 2 + left_cond_like.getUnderflowNumEdges() + right_cond_like.getUnderflowNumEdges();
		}
	else if (which == xtra_tip)
		{
		// cond_like == left_cond_like == right_cond_like
		nedges += 1 + cond_like.getUnderflowNumEdges();
		}
	else	// xtra_internal
		{
		// cond_like == left_cond_like, right_cond_like
		nedges += 1 + cond_like.getUnderflowNumEdges() + right_cond_like.getUnderflowNumEdges();
		}
			
	bool do_correction = (nedges >= underflow_num_edges);
	if (do_correction)
		{
		// Ok, we've traversed enough edges that it is time to take another factor out for underflow control
		
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
				// cond_like, left_cond_like == right_cond_like
				// uf_left is same object as uf_right in this case
				*uf = *uf_left + (UnderflowType)f;
				++uf_left;
				}
			else if (which == no_tips)
				{
				// cond_like, left_cond_like, right_cond_like
				*uf = *uf_left + *uf_right + (UnderflowType)f;
				++uf_left;
				++uf_right;
				}
			else if (which == xtra_tip)
				{
				// cond_like == left_cond_like == right_cond_like
				// nothing above this point to take account of
				*uf += (UnderflowType)f;
				}
			else	// xtra_internal
				{
				// cond_like == left_cond_like, right_cond_like
				*uf += *uf_right + (UnderflowType)f;
				++uf_right;
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
		// should move loop inside the conditionals for speed
		for (unsigned pat = 0; pat < num_patterns; ++pat, ++uf)
			{
			if (which == one_tip)
				{
				// cond_like, left_cond_like == right_cond_like
				// uf_left is same object as uf_right in this case
				*uf = *uf_left;
				++uf_left;
				}
			else if (which == no_tips)
				{
				// cond_like, left_cond_like, right_cond_like
				*uf = *uf_left + *uf_right;
				++uf_left;
				++uf_right;
				}
			else if (which == xtra_tip)
				{
				// cond_like == left_cond_like == right_cond_like
				// nothing above this point to take account of
				}
			else	// xtra_internal
				{
				// cond_like == left_cond_like, right_cond_like
				*uf += *uf_right;
				++uf_right;
				}
			}
		}
	cond_like.setUnderflowNumEdges(nedges);
	}

} // namespace phycas
