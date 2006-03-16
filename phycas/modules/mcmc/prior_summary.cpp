#include "phycas/force_include.h"
#include <stdexcept>
#include "phycas/modules/mcmc/prior_summary.hpp"
#include "phycas/trees/tree.hpp"
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::log;
	using std::pow;
	using std::exp;
#endif
/*----------------------------------------------------------------------------------------------------------------------
|	Recalculates the vectors `vTopoPriorByResClass', `vInducedResClassPrior' and `vExpectedResClassPrior' for the case in
|	which every topology (regardless of resolution) receives equal prior probability. That is, under the flat prior,
|	the star tree is just as probable a priori as any fully-resolved tree.
|
|	Example: n = 5, lnc = 1, M ranges from 0 (star tree) to n - 3 = 2 (fully-resolved trees)
|>
|	M = m - 1   No. trees    vInducedResClassPrior[M]    vExpectedResClassPrior[M]   vTopoPriorByResClass[M]
|	--------------------------------------------------------------------------------------------------
|	  0 (star)      1              1/26 = 0.03846               1/26 = 0.03846        1/26 = 0.03846
|	  1            10             10/26 = 0.38462              10/26 = 0.38462        1/26 = 0.03846
|	  2            15             15/26 = 0.57692              15/26 = 0.57692        1/26 = 0.03846
|	--------------------------------------------------------------------------------------------------
|	               26                     1.00000                 1.00000
|<
|	vTopoPriorByResClass[M] holds the log(prior) needed for a particular tree topology having m = M + 1 internal nodes
|	vInducedResClassPrior[M] holds the log of the induced prior for the resolution class characterized by m = M + 1 
|	  internal nodes (this is the resolution class prior expected when a flat prior is applied across all tree
|	  topologies, regardless of how resolved they are
|	vExpectedResClassPrior[M] is identical to vInducedResClassPrior[M] in this case
*/
void PriorSummary::RecalcFlatTopologyPriors()
	{
	assert(nTax > 0);
	vInducedResClassPrior.clear();
	vExpectedResClassPrior.clear();
	vTopoPriorByResClass.clear();
	if (nTax < 3)
		return;
	const unsigned mmax = nTax - 2;	
	double normalizing_constant = 0.0;
	for (unsigned m = 1; m <= mmax; ++m)
		{
		// Calculate value proportional to the desired prior if all resolution classes were equal in size
		//
		const double nmTopos =  treeCounter.GetUnrootedCount(nTax, m);
		const double ln_num_trees = log(nmTopos);
		vInducedResClassPrior.push_back(ln_num_trees);
		normalizing_constant += nmTopos;
		}
	const double ln_normalizing_constant = log(normalizing_constant);
	for (unsigned M = 0; M < mmax; ++M)
		vInducedResClassPrior[M]  -= ln_normalizing_constant;
	vExpectedResClassPrior = vInducedResClassPrior;
	vTopoPriorByResClass = VecDbl(mmax, -ln_normalizing_constant);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the vector `vTopoPriorByResClass' of which the Mth element (0-offset) is the natural log of the prior for a tree
|	topology having M + 1 internal nodes. A truncated geometric prior distribution is assumed, with probability of 
|	success equal to `p'. The random variable that is geometrically distributed in this case is M, and the distribution
|	is truncated because the largest possible value for M = n - 3. The polytomy prior is proportional to q^(M), making
|	the prior on the star tree highest and fully-resolved trees least. For p near 0, the prior will be close to flat 
|	across resolution classes, and for p near 1, the prior drops steeply. Note that p = 0 is invalid, and p = 1 allows
|	only the star tree (which is nonsensical). The function thus assumes 0 < p < 1.
|	
|	Example: n = 6, p = 1 - q = 0.6, M ranges from 0 (star tree) to n - 3 = 3 (fully-resolved trees)
|>
|	M = m - 1                                                   p = 0.6  p = 0.1
|	-------------------------------------------------------------------  -------
|	  0         p q^0 / (1 - q^{n-2}) = (0.6)(1)     / 0.9744 = 0.61576  0.29078
|	  1         p q^1 / (1 - q^{n-2}) = (0.6)(0.4)   / 0.9744 = 0.24631  0.26170
|	  2         p q^2 / (1 - q^{n-2}) = (0.6)(0.16)  / 0.9744 = 0.09852  0.23553
|	  3         p q^3 / (1 - q^{n-2}) = (0.6)(0.064) / 0.9744 = 0.03941  0.21198
|	-------------------------------------------------------------------  -------
|	                                                            1.00000  1.00000
|<
|	The normalizing factor 1 - q^{n-2} is just the cumulative distribution up to M = n - 3, which can either be looked
|	up or derived as follows:
|>
|	   sum prob = p (q^0 + q^1 + q^2 + q^3)
|	            = p sum_{k=1}^{n-2} q^{k-1}
|	              p (q^{n-2} - 1)   p (1 - q^{n-2})
|	            = --------------- = --------------- = 1 - q^{n-2}
|	                   q - 1               p
|<
|	Despite the explanation above about normalization, we needn't worry about the normalizing constant. One thing we 
|	do need to worry about is the induced prior, which we must overcome	in order to have the desired actual prior 
|	prevail. To do this, we must make the prior actually used in the analysis (i.e. the one returned by this function)
|	be proportional to q^{m - 1} / pi[m], where pi[m] = T_{n,m}/T_n is the induced prior for the resolution class 
|	having m internal nodes.
*/
void PriorSummary::RecalcGeometricTopologyPriors(
  double p)	/**< is the probability of success parameter of the truncated geometric distribution; specifying the invalid value 0 for `polytomy_p' results in the use of a flat prior across all topologies (resolution classes containing relatively more trees thus have relatively greater prior weight) */
	{
	assert(nTax > 0);
	assert(p > 0.0);
	assert(p < 1.0);
	vTopoPriorByResClass.clear();
	vInducedResClassPrior.clear();
	vExpectedResClassPrior.clear();
	if (nTax < 3)
		return;
	const double q = 1.0 - p;
	const double ln_polytomy_p = log(p);
	const double ln_polytomy_q = log(q);
	const double ln_normalizing_factor = log(1.0 - pow(q, (double)(nTax - 2)));
	for (unsigned M = 0; M < nTax - 2; ++M)
		{
		// Compute actual prior
		//
		double ln_desired_prior = (double)M*ln_polytomy_q;
		double ln_induced_prior = GetLnInducedPrior(M + 1);
		vInducedResClassPrior.push_back(ln_induced_prior);
		double ln_actual_prior = ln_desired_prior - ln_induced_prior;
		vTopoPriorByResClass.push_back(ln_actual_prior);

		// Compute expected prior
		//
		double ln_numer = ln_polytomy_p + (double)M*ln_polytomy_q;
		double ln_expected_prior = ln_numer - ln_normalizing_factor;
		vExpectedResClassPrior.push_back(ln_expected_prior);
		}
	}
	
double PriorSummary::GetLnInducedPrior(UInt m) const
	{
	if (m == 0 || nTax < 2 || m > nTax - 1)
		{
		assert(false);
		std::string s;
		s << "m = " << m << ", nTax = " << nTax << " in GetLnInducedPrior()";
		throw std::domain_error(s);
		}
		//@ can easily imagine total unrooted trees being greater than UINT_MAX
	return log(treeCounter.GetUnrootedCount(nTax, m)) - log(treeCounter.GetTotalUnrooted(nTax));
	}
	
double PriorSummary::GetLnTopologyPrior(const Tree &tree)
	{
	if (tree.GetNLeaves() != nTax)
		{
		nTax = tree.GetNLeaves();
		RecalcTopologyPriors();
		}
	const UInt m = tree.GetNInternals();
	if (m == 0 || nTax < 3)
		return 0.0;  //@ throw std::domain_error("GetLnTopologyPrior() for non-tree");
	assert(vTopoPriorByResClass.size() >= m);
	return vTopoPriorByResClass[m - 1];
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Recalculates the vectors `vTopoPriorByResClass', `vInducedResClassPrior' and `vExpectedResClassPrior' for the case in
|	which a constant factor approach is used to determining the prior. Given two tree topologies, one having m internal
|	nodes and the other having m + 1 internal nodes, this prior makes the one having m internal nodes exp(lnc) times 
|	more probable a priori than the one having m + 1 internal nodes. That is, less resolved trees are given higher prior
|	probability, and the prior is proportional (on the log scale) to the difference in the number of internal nodes.
|	The value lnc = 1 works well. This value punishes trees one log-likelihood unit for each extra internal node 
|	parameter they possess.
|
|	Example: n = 5, lnc = 1, M ranges from 0 (star tree) to n - 3 = 2 (fully-resolved trees)
|>
|	M = m - 1   No. trees    vInducedResClassPrior[M]     unnormalized     vExpectedResClassPrior[M]  vTopoPriorByResClass[M]
|	-------------------------------------------------------------------------------------------------------------------
|	  0 (star)      1              1/26 = 0.03846     (1)(c^2) =  7.38906     (1)(c^2)/T = 0.14906    (c^2)/T = 0.14906
|	  1            10             10/26 = 0.38462      (10)(c) = 27.18282      (10)(c)/T = 0.54835      (c)/T = 0.05484
|	  2            15             15/26 = 0.57692      (15)(1) = 15.00000      (15)(1)/T = 0.30259      (1)/T = 0.02017
|	-------------------------------------------------------------------------------------------------------------------
|	               26                     1.00000            T = 49.57187                  1.00000
|<
|	vTopoPriorByResClass[M] holds the log(prior) needed for a particular tree topology having m = M + 1 internal nodes
|	vInducedResClassPrior[M] holds the log of the induced prior for the resolution class characterized by m = M + 1 
|	  internal nodes (this is the resolution class prior expected when a flat prior is applied across all tree
|	  topologies, regardless of how resolved they are
|	vExpectedResClassPrior[M] holds the log of the prior for the resolution class characterized by m = M + 1 
|	  internal nodes (this is the resolution class prior expected when the constant factor prior is applied)
|
|	Note: this function will fail (or be exceedingly slow) for large numbers of taxa as currently written
*/
void PriorSummary::RecalcConstTopologyPriors(
  double lnc)	/**< is the natural log of the constant factor separating trees different by one internal node */
	{
	//@ this function will fail (or be exceedingly slow) for large numbers of taxa as currently written
	assert(nTax > 0);
	vTopoPriorByResClass.clear();
	vInducedResClassPrior.clear();
	vExpectedResClassPrior.clear();
	if (nTax < 3)
		return;
	unsigned mmax = nTax - 2;	
	double normalizing_constant = 0.0;
	for (unsigned m = 1; m <= mmax; ++m)
		{
		// Calculate value proportional to the desired prior if all resolution classes were equal in size
		//
		const double ln_unit_desired_prior = (double)(mmax - m)*lnc;
		const double ln_class_desired_prior = ln_unit_desired_prior + log(treeCounter.GetUnrootedCount(nTax, m));
		vTopoPriorByResClass.push_back(ln_unit_desired_prior);
		vExpectedResClassPrior.push_back(ln_class_desired_prior);
		normalizing_constant += exp(ln_class_desired_prior);
		// Calculate and store the log of the induced prior
		//
		vInducedResClassPrior.push_back(GetLnInducedPrior(m));
		}
	const double ln_normalizing_constant = log(normalizing_constant);
	for (unsigned M = 0; M < mmax; ++M)
		{
		vTopoPriorByResClass[M] -= ln_normalizing_constant;
		vExpectedResClassPrior[M] -= ln_normalizing_constant;
		}
	}

