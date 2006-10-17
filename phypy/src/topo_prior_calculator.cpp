/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

//#include "phycas/force_include.h"
#include "phypy/src/topo_prior_calculator.hpp"
#include "phypy/src/basic_tree.hpp"
//#include "boost/shared_ptr.hpp"
#include <cmath>

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes `counts' vector for the supplied number of internal nodes (`n') using the method outlined by Joe 
|	Felsenstein in his 2004 book and also in Felsenstein (1978) and Felsenstein (1981). 
|	
|	Felsenstein, J. 1978. The number of evolutionary trees. Syst. Zool. 27: 27-33.
|	Felsenstein, J. 1981. Syst. Zool. 30: 122. 
|	
|	As an example, if `n' equals 4, the `counts' vector would look like this when the function returned:
|>
|	      0          1         2          3          4      
|	+----------+----------+----------+----------+----------+
|	|   236.0  |    1.0   |   25.0   |   105.0  |   105.0  |
|	+----------+----------+----------+----------+----------+
|>
|	For fully-resolved, unrooted trees, there are 4 internal nodes for 6 taxa, so `counts' needs only be computed up to 
|	`counts'[4]. `counts'[0] equals the sum of all the other elements (i.e. 236 = 1 + 25 + 105 + 105). Thus, there are 
|	105 unrooted tree topologies having 6 taxa and m = 3 interior nodes. The number of fully-resolved, unrooted trees 
|	for 6 taxa is the last entry, `counts'[4] = 105.
|	
|	If `n' equals 5, the `counts' vector would look like this when the function returned:
|>
|	      0          1         2          3          4          5 
|	+----------+----------+----------+----------+----------+----------+
|	|  2752.0  |    1.0   |   56.0   |   490.0  |  1260.0  |   945.0  |
|	+----------+----------+----------+----------+----------+----------+
|>
|	For fully-resolved, rooted trees, there are 5 internal nodes for 6 taxa, so `counts' needs only be computed up to 
|	`counts'[5]. `counts'[0] equals the sum of all the other elements (i.e. 2752 = 1 + 56 + 490 + 1260 + 945). Thus, 
|	there are 490 rooted tree topologies having 6 taxa and m = 3 interior nodes. The number of fully-resolved, rooted
|	trees for 6 taxa is the last entry, `counts'[5] = 945.
*/
void TopoPriorCalculator::RecalcCountsAndPriorsImpl(
  unsigned n) /**< is the number of internal nodes (equals ntax - 1 for rooted trees and ntax - 2 for unrooted trees) */
	{
	counts.clear();
	counts.push_back(0.0); // This element will hold sum of all other elements in the end
	counts.push_back(1.0); // m = 1

	for (unsigned m = 2; m <= n; ++m)
		{
		AddNextCount(m);
		}

	// Recalculate the `topology_prior' vector too
	RecalcPriorsImpl();

	topo_priors_dirty = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds element to `counts' vector representing `m' internal nodes given that it is already correctly computed for `m'
|	minus 1 internal nodes. On entry, assumes that the size of the `counts' vector is `m'. One new element is added to 
|	the `counts' vector by this function, so upon exit, `counts' has length `m' + 1 (the first element of `counts' is
|	reserved for holding the sum of all other elements. Assumes that `m' is greater than 1.
*/
void TopoPriorCalculator::AddNextCount(
  unsigned m)	/**< is the number of internal nodes (should equal current length of `counts' vector) */
	{
	PHYCAS_ASSERT(m > 1);
	PHYCAS_ASSERT(counts.size() == m);
	counts.push_back(0.0);

	// This will probably not make any sense at all unless you read the section on computing numbers of trees
	// in Joe Felsenstein's 2004 book
	counts[0] = 1.0;
	counts[1] = 1.0;
	double a = 1.0;
	for (unsigned k = 2; k <= m; ++k)
		{
		double b = counts[k];
		double c = a*(m + k - 1);
		if (k < m)
			c += b*k;
		counts[k] = c;
		counts[0] += c;
		a = b;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes `topology_prior' vector for the supplied number of internal nodes (`n'). The element `topology_prior'[m]
|	is the natural logarithm of the unnormalized prior probability of a (rooted/unrooted) tree with `ntax' taxa and m 
|	internal nodes. The rooted/unrooted status is determined by the state of the data member `is_rooted'. The element 
|	`topology_prior'[0] is the natural log of the normalizing constant. Normally, only unnormalized values are needed 
|	because MCMC deals in prior ratios, but if for some reason the normalized prior probability is needed, it can be 
|	computed as exp{`topology_prior'[m] - `topology_prior'[0]}. This function produces a vector having the same length 
|	as `counts', so it is important to make sure `counts' is correct before calling this function. To ensure this, a 
|	call to RecalcPriorsImpl is the last thing done by RecalcCountsImpl before it returns.
|
|	For example, consider an unrooted tree with ntax = 6. Such a tree has 4 internal nodes, and calling RecalcPriorsImpl
|	would yield the following if `is_resolution_class_prior' is false:
|>
|	     unnormalized                                     
|	m    polytomy prior     C = 2	                      
|	----------------------------------
|	1        C^3              8		     topology_prior[1] = ln(8)  = 2.079   
|	2        C^2              4		     topology_prior[2] = ln(4)  = 1.386   
|	3        C^1              2		     topology_prior[3] = ln(2)  = 0.693   
|	4        C^0              1		     topology_prior[4] = ln(1)  = 0.000
|	----------------------------------
|	       C^4 - 1         16 - 1	                      
|	       -------         ------ = 15   topology_prior[0] = ln(15) = 2.708                  
|	        C - 1           2 - 1	                      
|>
|	If instead `is_resolution_class_prior' is true, we have:
|>
|	              unnormalized
|	              resolution                                         
|	m    counts   class prior        C = 2		                     
|	---------------------------------------------
|	1         1     (C^3)/1         8/1   = 8.000   topology_prior[1] = ln(8.000) =  2.079
|	2        25     (C^2)/25        4/25  = 0.160   topology_prior[2] = ln(0.160) = -1.833
|	3       105     (C^1)/105       2/105 = 0.019   topology_prior[3] = ln(0.019) = -3.963
|	4       105     (C^0)/105       1/105 = 0.010   topology_prior[4] = ln(0.010) = -4.605 
|	---------------------------------------------
|	            (no easy formula)           8.189   topology_prior[0] = ln(8.189) =  2.103
|>												    
*/												                     
void TopoPriorCalculator::RecalcPriorsImpl()	 
	{
	topology_prior.clear();
	topology_prior.push_back(0.0);	// This will hold the normalizing constant in the end

	// Figure out the maximum possible value for m, the number of internal nodes
	unsigned maxm = ntax - (is_rooted ? 1 : 2);

	// counts vector should have length equal to maxm - 1 if everything is ok
	PHYCAS_ASSERT(maxm == (unsigned)counts.size() - 1);

	double total = 0.0;
	double logC = std::log(C);
	for (unsigned m = 1; m <= maxm; ++m)
		{
		double logCterm = (double)(maxm - m)*logC;
		double log_count_m = (is_resolution_class_prior ? std::log(counts[m]) : 0.0);
		double log_v = logCterm - log_count_m;
		total += std::exp(log_v);
		topology_prior.push_back(log_v);
		}
	topology_prior[0] = std::log(total);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of trees having `n' taxa and `m' internal nodes. Calls RecalcCountsAndPriors function if `n' is 
|	not equal to `ntax'. Assumes `m' is greater than 0. If `is_rooted' is true, assumes `m' is less than `ntax'. If 
|	`is_rooted' is false, assumes `m' less than `ntax' - 1. 
*/
double TopoPriorCalculator::GetCount(
  unsigned n,	/**< is the number of taxa */
  unsigned m)	/**< is the number of internal nodes */
	{
	PHYCAS_ASSERT((is_rooted && (m < n)) || (!is_rooted && (m < n - 1)));
	if (n != ntax)
		SetNTax(n);
	if (topo_priors_dirty)
		Reset();
	return counts[m];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of saturated (i.e. fully-resolved and thus having as many internal nodes as possible) trees 
|	of `n' taxa. Calls RecalcCountsAndPriors function if `n' is not equal to `ntax'.
*/
double TopoPriorCalculator::GetSaturatedCount(
  unsigned n)	/**< is the number of taxa */
	{
	if (n != ntax)
		SetNTax(n);
	if (topo_priors_dirty)
		Reset();
	return counts[counts.size() - 1];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the total number of trees (including all resolution classes from star to saturated trees) for `n' taxa. 
|	Calls RecalcCountsAndPriors function if `n' is not equal to `ntax'.
*/
double TopoPriorCalculator::GetTotalCount(
  unsigned n)	/**< is the number of taxa */
	{
	if (n != ntax)
		SetNTax(n);
	if (topo_priors_dirty)
		Reset();
	return counts[0];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructs a vector of realized resolution class priors from the values in the `topology_prior' vector. If 
|	`topo_priors_dirty' is true, it recomputes the `counts' and `topology_prior' vectors first. The mth element of the
|	returned vector is set to T_{n,m}*`topology_prior'[m] for m > 0. The 0th element of the returned vector holds the
|	normalization constant (sum of all other elements). This function is not efficient because it is intended only to 
|	be used for providing information to the user on request. Table 2, p. 248, in our "Polytomies and Bayesian 
|	Phylogenetic Inference" paper (Lewis, P. O., M. T. Holder and K. E. Holsinger. 2005. Systematic Biology 54(2):
|	241-253) presented (normalized) values from this vector.
*/
std::vector<double> TopoPriorCalculator::GetRealizedResClassPriorsVect()
	{
	if (topo_priors_dirty)
		Reset();

	std::vector<double> v;
	v.reserve(topology_prior.size());
	v.push_back(0.0);

	//@POL should use a version of the transform algorithm here
	double total = 0.0;
	unsigned sz = (unsigned)topology_prior.size();
	for (unsigned i = 1; i < sz; ++i)
		{
		double log_Tnm = std::log(counts[i]);
		double log_prior = log_Tnm + topology_prior[i];
		v.push_back(log_prior);
		total += std::exp(log_prior);
		}
	v[0] = std::log(total);

	return v;
	}

}	// namespace phycas
