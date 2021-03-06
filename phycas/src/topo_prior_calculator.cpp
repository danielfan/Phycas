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

#include <cmath>
#include <limits>
#include "phycas/src/topo_prior_calculator.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/basic_lot.hpp"

namespace phycas
{


FocalTreeTopoProbCalculator::FocalTreeTopoProbCalculator(TreeShPtr t)
    :focalTree(t)
    {
    assert(bool(focalTree));
    ntips = focalTree->GetNObservables();
    focalTree->RecalcAllSplits(ntips);
    TreeNode * nd = focalTree->GetFirstPreorder();
    assert(nd->IsTipRoot());
    nd = nd->GetNextPreorder();
    while (nd)
        {
        if (!nd->IsExternalEdge()) 
            {
            double  split_prob = nd->GetEdgeLen();
            assert(split_prob > 0.0 && split_prob < 1.0);
            Split & s = nd->GetSplit();
            if (s.IsBitSet(0))
                s.InvertSplit();
            splitToProbMap[s] = split_prob;
            }
        nd = nd->GetNextPreorder();
        }
    buildScratchTree();
    }

void FocalTreeTopoProbCalculator::buildScratchTree() 
    {
    scratchTree.MirrorTopology(*focalTree);
    }

void Tree::MirrorTopology(Tree &source) 
    {
    Clear();
    TreeNode * fnd = source.GetFirstPreorder();
    assert(fnd);
    assert(fnd->IsTipRoot());
    TreeNode * tmpSNd = GetNewNode();
    SetFirstPreorder(tmpSNd);


    tmpSNd->SetCorrespondingNode(fnd);
    fnd->SetCorrespondingNode(tmpSNd);
    tmpSNd->CopyNonPointerFields(*fnd);
    tmpSNd->par = 0L;
    tmpSNd->lChild = 0L;
    tmpSNd->rSib = 0L;
    tmpSNd->nextPreorder = 0L;
    tmpSNd->prevPreorder = 0L;
    
    TreeNode * snd = GetRoot();
    assert(snd);
    assert(snd->IsTipRoot());

    TreeNode * prevS = snd;
    fnd = fnd->GetNextPreorder();
    while (fnd)
        {
        tmpSNd = GetNewNode();
        tmpSNd->SetCorrespondingNode(fnd);
        fnd->SetCorrespondingNode(tmpSNd);
        tmpSNd->CopyNonPointerFields(*fnd);
        tmpSNd->lChild = 0L;
        tmpSNd->rSib = 0L;
        tmpSNd->nextPreorder = 0L;
        prevS->nextPreorder = tmpSNd; 
        tmpSNd->prevPreorder = prevS;
        TreeNode * sPar = fnd->par->GetCorrespondingNode();
        assert(sPar);
        sPar->AddChild(tmpSNd);
        prevS = tmpSNd;
        fnd = fnd->GetNextPreorder();
        }
    
    }

// Assumes that all nodes in focalTree have a "CorrespondingNode" that is allocated in 
//  "this."
// does not alter anything except the "key" navigational pointers (par, lChild, rSib), and the split field
void Tree::RebuildTopologyFromMirror(const Tree & source)
    {
    const TreeNode * fnd = source.GetFirstPreorderConst();
    assert(fnd);
    assert(fnd->IsTipRoot());
    TreeNode * tmpSNd = fnd->GetCorrespondingNode();
    SetFirstPreorder(tmpSNd);
    tmpSNd->par = 0L;
    tmpSNd->lChild = (fnd->lChild ? fnd->lChild->GetCorrespondingNode() : 0L);
    tmpSNd->rSib = 0L;
    tmpSNd->nextPreorder = tmpSNd->lChild;
    tmpSNd->prevPreorder = 0L;
    TreeNode * prevS = tmpSNd;
    fnd = fnd->GetNextPreorderConst();
    while (fnd)
        {
        tmpSNd = fnd->GetCorrespondingNode();
        tmpSNd->lChild = (fnd->lChild ? fnd->lChild->GetCorrespondingNode() : 0L);
        tmpSNd->rSib = (fnd->rSib ? fnd->rSib->GetCorrespondingNode() : 0L);
        tmpSNd->par = (fnd->par ? fnd->par->GetCorrespondingNode() : 0L);
        tmpSNd->nextPreorder = (fnd->nextPreorder ? fnd->nextPreorder->GetCorrespondingNode() : 0L);
        tmpSNd->prevPreorder = prevS;
        tmpSNd->split = fnd->split;
        TreeNode * sPar = fnd->par->GetCorrespondingNode();
        assert(sPar);
        prevS = tmpSNd;
        fnd = fnd->GetNextPreorderConst();
        }
    
    }

// Replaces "this" with its children in the tree
// only fixes the "key" navigational pointers (par, lChild, rSib)
// does not alter "this" node's pointers
void TreeNode::CollapseEdge()
    {
    TreeNode * p = this->par;
    if (!p)
        return;
    TreeNode * ls = this->FindLeftSib();
    TreeNode * c = this->lChild;
    if (ls == 0L)
        p->lChild = c;
    else
        ls->rSib = c;
    if (!c)
        return;
    for (;;c = c->rSib) 
        {
        c->par = p;
        if (c->rSib == 0L)
            {
            c->rSib = this->rSib;
            break;
            }
        }
    }

double FocalTreeTopoProbCalculator::CalcLnNumTreesMaxDistFromTreeInSelectedRegion(const TreeNode *selectedFirstFork, unsigned numLeaves) const
{
    PHYCAS_ASSERT(selectedFirstFork);
    double lnNumTrees = 0.0;
    
    if (numLeaves < 7) 
        {
        if (numLeaves == 4)
            return log(2);
        if (numLeaves == 5)
            return log(10);
        const TreeNode * lc = selectedFirstFork->lChild;
        const TreeNode * rc = lc->rSib;
        if (lc->IsSelected())
            {
            if (rc->IsSelected())
                return log(89); // pectinate 6-leaf
            lc = lc->lChild;
            rc = lc->rSib;
            if (lc->IsSelected() && rc->IsSelected())
                return log(68); // symmetric 6-leaf
            else
                return log(89); // pectinate 6-leaf
            }
        else
            {
            PHYCAS_ASSERT(rc->IsSelected());
            lc = rc->lChild;
            rc = lc->rSib;
            if (lc->IsSelected() && rc->IsSelected())
                return log(68); // symmetric 6-leaf
            else
                return log(89); // pectinate 6-leaf
            }
        
        }
    PHYCAS_ASSERT(false); // not implemented...
    return lnNumTrees;
}

double FocalTreeTopoProbCalculator::CalcTopologyLnProb(Tree & testTree) const
    {
    testTree.RecalcAllSplits(ntips);
    const TreeID & testTreeID = testTree.getTreeID();
    const TreeID::const_iterator ttIDIt = testTreeID.end();
    scratchTree.RebuildTopologyFromMirror(*focalTree);
    TreeNode * fnd = scratchTree.GetFirstPreorder();
    fnd->SetIsSelected(false);
    fnd = fnd->GetNextPreorder();
    double lnProb = 0.0;
    
    omittedNodes.clear();
    std::map<TreeNode *, std::vector<TreeNode *> > polytomyToCollapsed;
    while (fnd)
        {
        if (fnd->IsExternalEdge()) 
            {
            if (testTreeID.find(fnd->GetSplitConst()) == ttIDIt)
                {
                lnProb += log(1 - fnd->GetEdgeLen()); // could store log(1-p) in support and log(p) in edge_len to cut down on logs
                omittedNodes.insert(fnd);
                fnd->SetIsSelected(true);
                TreeNode * p = fnd->GetParent();
                if (p->IsSelected()) 
                    fnd->prevPreorder = p->prevPreorder; // @ DANGEROUS overloading of prevPreorder to store the "deepest node" that will represent the polytomy.
                else
                    {
                    fnd->prevPreorder = p->GetParent();
                    }
                fnd->CollapseEdge();
                std::vector<TreeNode *> & vc = polytomyToCollapsed[fnd->prevPreorder->GetCorrespondingNode()];
                vc.push_back(fnd->GetCorrespondingNode());
                }
            else
                {
                lnProb += log(fnd->GetEdgeLen()); // could store log(1-p) in support and log(p) in edge_len to cut down on logs
                fnd->SetIsSelected(false);
                }
            }
        fnd = fnd->GetNextPreorder();
        }
    
    double lnDenominator = 0.0;
    for (std::map<TreeNode *, std::vector<TreeNode *> >::const_iterator pNdIt = polytomyToCollapsed.begin(); pNdIt != polytomyToCollapsed.end(); ++pNdIt)
        {
        TreeNode * polytomyNd = pNdIt->first;
        const std::vector<TreeNode *> & vc = pNdIt->second;

        // flag edges that are collapsed in the consensus
        for (std::vector<TreeNode *>::const_iterator ndIt = vc.begin(); ndIt != vc.end(); ++ndIt)
            (*ndIt)->SetIsSelected(true);

        lnDenominator += CalcLnNumTreesMaxDistFromTreeInSelectedRegion(polytomyNd, 3 + vc.size());

        // "unflag" edges that are collapsed in the consensus
        for (std::vector<TreeNode *>::const_iterator ndIt = vc.begin(); ndIt != vc.end(); ++ndIt)
            (*ndIt)->SetIsSelected(false);
        }
    
    return lnProb - lnDenominator;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes `counts' vector for the supplied number of internal nodes (`n') using the method outlined by Joe 
|	Felsenstein in his 2004 book and also in Felsenstein (1978) and Felsenstein (1981). 
|	
|	Felsenstein, J. 1978. The number of evolutionary trees. Syst. Zool. 27: 27-33.
|	Felsenstein, J. 1981. Syst. Zool. 30: 122. 
|   Felsenstein, J. 2004. Inferring Phylogeny. Sinauer, Sunderland, Massachusetts.
|   
|   Below I've reproduced the table from Felsenstein 2004 illustrating how to calculate the number of possible rooted
|   tree topologies for any number of internal nodes:
|>
|                      number of taxa in rooted tree            counts    nfactors
|             2       3       4       5        6        7         8        
|	n     +-------+-------+-------+--------+--------+--------+---------+   +---+  This table shows results for 8     
|	u  1  |   1   |   1   |   1   |   1    |    1   |    1   |    1    |   | 0 |  taxa and rooted trees. Columns 
|	m     +-------+-------+-------+--------+--------+--------+---------+   +---+  corresponding to 2-7 taxa would have
|	b  2          |   3   |   10  |   25   |   56   |   119  |   246   |   | 0 |  been replaced by the column 
|	e             +-------+-------+--------+--------+--------+---------+   +---+  corresponding to 8 taxa. This column
|	r  3                  |   15  |  105   |   490  |   1918 |   6825  |   | 0 |  is now stored in the `counts' vector.
|	                      +-------+--------+--------+--------+---------+   +---+  The scaling factor used in this case
|	o  4                          | 0.0105 | 0.1260 |  0.945 |  5.6980 |   | 1 |  was 10,000. The number of times this
|	f                             +--------+--------+--------+---------+   +---+  factor was applied is stored in the 
|	   5                                   | 0.0945 | 1.7325 | 19.0575 |   | 1 |  nfactors vector. Thus, the "count"
|	n                                      +--------+--------+---------+   +---+  for 5 internal nodes is not really
|	o  6                                            | 1.0395 | 27.0270 |   | 1 |  19.0575, but instead 190575:
|	d                                               +--------+---------+   +---+  count[4]*(scaling_factor)^nfactors[4]
|	e  7                                                     | 13.5135 |   | 1 |
|	s                                                        +---------+   +---+
|>
|   The RecalcCountsAndPriorsImpl function works from left to right, calculating each column in turn. Within a column,
|   it works down. Each cell except the first in a given column requires knowledge of two cells from the column to the
|   left (the cell to its immediate left as well as the cell above the cell to its immediate left). This is because 
|   in order to know how many tree topologies there are for N taxa, one needs to know how many places a taxon could be
|   added to a tree with one fewer taxa. The N-1 taxon tree could have the same number of internal nodes as the new
|   tree (the new taxon was added to an existing node, either creating a new polytomy or enlarging an existing one), or
|   the N-1 taxon tree could have one fewer internal nodes (in which case the new taxon inserts a new node).
*/
void PolytomyTopoPriorCalculator::RecalcCountsAndPriorsImpl(
  unsigned n) /**< is the number of internal nodes (equals ntax - 1 for rooted trees and ntax - 2 for unrooted trees) */
	{
    if (is_resolution_class_prior)
        counts_dirty = true;
    if (counts_dirty)
        {
        double scaling_factor = exp(log_scaling_factor);

	    counts.clear();
	    counts.push_back(1.0); // counts are always 1 for m = 1

        int last_factor = 0;
	    nfactors.clear();
	    nfactors.push_back(last_factor); // never need to scale this one

        // temporary variables
        double a, b, c;
        double max_log_double = log(DBL_MAX);
        double epsilon = scaling_factor/10.0; // value arbitrary, but must be larger than zero and less than scaling_factor

        // Compute the vector of counts for the number of internal nodes specified
        // This is the main loop over columns. z is the number of taxa minus 1 for rooted trees, and is the number
        // of taxa minus 2 for unrooted trees
	    for (unsigned z = 2; z <= n; ++z)
		    {
		    counts.push_back(0.0);  // this column is one element longer than the column to its left
            nfactors.push_back(last_factor);
            b = counts[0];  // counts[0] is always 1.0 because there is only one star tree topology

            // This is the loop over rows within the current column. m + 1 is the number of internal nodes.
            for (unsigned m = 1; m < z; ++m)
                {
                unsigned num_internal_nodes = m + 1;
                double diff = (double)(nfactors[m - 1] - nfactors[m]);
                double log_factor = diff*log_scaling_factor;
                if (log_factor >= max_log_double)
                    {
                    std::cerr << "Oops! log_factor >= max_log_double" << std::endl;
                    }
                PHYCAS_ASSERT(log_factor < max_log_double);
                a = b*exp(log_factor);
                b = counts[m];
                c = a*((double)(z + num_internal_nodes - 1));
                if (num_internal_nodes < z)
                    {
                    c += b*(double)num_internal_nodes;
                    }
                if (c > scaling_factor)
                    {
                    //unsigned prev_nfactors = nfactors[m];
                    double incr = floor(log(c - epsilon)/log_scaling_factor);
                    nfactors[m] += (unsigned)incr;
                    last_factor = nfactors[m];
                    counts[m] = exp(log(c) - incr*log_scaling_factor);
                    b = exp(log(b) - incr*log_scaling_factor);
                    }
                else
                    counts[m] = c;
                }
		    }

        // Now compute the log of the total number of tree topologies over all possible resolution classes
        // (i.e. number of internal nodes)
        // Begin by creating a vector of log counts and finding the largest value (this will be factored out
        // to avoid overflow)
        std::vector<double> v;
        unsigned sz = (unsigned)nfactors.size();
        PHYCAS_ASSERT(sz == counts.size());
        double max_log_count = 0.0;
        for (unsigned i = 0; i < sz; ++i)
            {
            double num_factors = (double)nfactors[i];
            double log_count = num_factors*log_scaling_factor + log(counts[i]);
            if (log_count > max_log_count)
                max_log_count = log_count;
            v.push_back(log_count);
            }

        // Compute log sum of counts by factoring out the largest count. Underflow will occur, but only for 
        // counts that are so much smaller than the dominant counts that the underflow can be ignored for
        // all practical purposes
        double sum = 0.0;
        for (std::vector<double>::const_iterator it = v.begin(); it != v.end(); ++it)
            {
            double diff = (*it) - max_log_count;
            sum += exp(diff);
            }
        PHYCAS_ASSERT(sum > 0.0);
        log_total_count = log(sum) + max_log_count;
        counts_dirty = false;
        }
    else
        {
        nfactors.clear();
        counts.clear();
        counts_dirty = true;    // ensures that counts will be calcualated if someone asks for one, say by calling GetCount
        }

	// Recalculate the `topology_prior' vector too
	RecalcPriorsImpl();

	topo_priors_dirty = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes `topology_prior' vector for the supplied number of internal nodes (`n'). The element `topology_prior'[m]
|	is the natural logarithm of the unnormalized prior probability of a (rooted/unrooted) tree with `ntax' taxa and m 
|	internal nodes. The rooted/unrooted status is determined by the state of the data member `is_rooted'. The element 
|	`topology_prior'[0] is the natural log of the normalizing constant. Normally, only unnormalized values are needed 
|	because MCMC deals in prior ratios, but if for some reason the normalized prior probability is needed, it can be 
|	computed as exp{`topology_prior'[m] - `topology_prior'[0]}. This function requires `counts' to be correct if using
|   the resolution class prior. Thus, never call RecalcPriorsImpl directly, only invoke it indirectly by calling the
|   function RecalcCountsAndPriorsImpl.
|
|	Consider an unrooted tree with ntax = 6. Such a tree has 4 internal nodes, and calling RecalcPriorsImpl
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
void PolytomyTopoPriorCalculator::RecalcPriorsImpl()	 
	{
	topology_prior.clear();
	topology_prior.push_back(0.0);	// This will hold the normalizing constant in the end

	// Figure out the maximum possible value for m, the number of internal nodes
	unsigned maxm = ntax - (is_rooted ? 1 : 2);

    if (is_resolution_class_prior)
        {
	    // counts vector should have length equal to maxm if everything is ok
	    PHYCAS_ASSERT(maxm == (unsigned)counts.size());

	    double total = 0.0;
	    double logC = std::log(C);
	    for (unsigned m = 1; m <= maxm; ++m)
		    {
		    double logCterm = (double)(maxm - m)*logC;
		    double log_count_m = std::log(counts[m - 1]) + log_scaling_factor*(double)nfactors[m - 1];
		    double log_v = logCterm - log_count_m;
		    total += std::exp(log_v);
		    topology_prior.push_back(log_v);
		    }
	    topology_prior[0] = std::log(total);
        }
    else
        {
	    double total = 0.0;
	    double logC = std::log(C);
	    for (unsigned m = 1; m <= maxm; ++m)
		    {
		    double logCterm = (double)(maxm - m)*logC;
		    total += std::exp(logCterm);
		    topology_prior.push_back(logCterm);
		    }
	    topology_prior[0] = std::log(total);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the number of trees having `n' taxa and `m' internal nodes. Calls RecalcCountsAndPriors
|   function if `n' is not equal to `ntax'. Assumes `m' is greater than 0. If `is_rooted' is true, assumes `m' is less
|   than `ntax'. If `is_rooted' is false, assumes `m' less than `ntax' - 1. 
*/
double PolytomyTopoPriorCalculator::GetLnCount(
  unsigned n,	/**< is the number of taxa */
  unsigned m)	/**< is the number of internal nodes */
	{
	PHYCAS_ASSERT((is_rooted && (m < n)) || (!is_rooted && (m < n - 1)));
	if (n != ntax)
		SetNTax(n);
	if (counts_dirty)
		Reset();
    double nf = (double)(nfactors[m - 1]);
    double log_count = nf*log_scaling_factor + log(counts[m - 1]);
	return log_count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of saturated (i.e. fully-resolved and thus having as many internal nodes as possible) trees 
|	of `n' taxa. Calls RecalcCountsAndPriors function if `n' is not equal to `ntax'.
*/
double PolytomyTopoPriorCalculator::GetLnSaturatedCount(
  unsigned n)	/**< is the number of taxa */
	{
	if (n != ntax)
		SetNTax(n);
	if (counts_dirty)
		Reset();
    unsigned last = (unsigned)(counts.size() - 1);
    double nf = (double)(nfactors[last]);
    double log_count = nf*log_scaling_factor + log(counts[last]);
	return log_count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the total number of trees for `n' taxa, including all resolution classes from the star 
|   tree to fully resolved (saturated) trees. Calls RecalcCountsAndPriors function if `n' is not equal to `ntax' or if
|   not using the resolution class prior (in which case counts have not been calculated).
*/
double PolytomyTopoPriorCalculator::GetLnTotalCount(
  unsigned n)	/**< is the number of taxa */
	{
	if (n != ntax)
		SetNTax(n);
	if (counts_dirty)
		Reset();
	return log_total_count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructs a vector of realized resolution class priors from the values in the `topology_prior' vector. If 
|	`topo_priors_dirty' is true, it recomputes the `topology_prior' vectors first. The mth element of the
|	returned vector is set to T_{n,m}*`topology_prior'[m] for m > 0. The 0th element of the returned vector holds the
|	normalization constant (sum of all other elements). This function is not efficient because it is intended only to 
|	be used for providing information to the user on request. Table 2, p. 248, in the "Polytomies and Bayesian 
|	Phylogenetic Inference" paper (Lewis, P. O., M. T. Holder and K. E. Holsinger. 2005. Systematic Biology 54(2):
|	241-253) presented (normalized) values from this vector.
*/
std::vector<double> PolytomyTopoPriorCalculator::GetRealizedResClassPriorsVect()
	{
	if (!is_resolution_class_prior)
        counts_dirty = true;
	if (topo_priors_dirty || counts_dirty)
		Reset();

	std::vector<double> v;
	v.reserve(topology_prior.size());
	v.push_back(0.0);

	unsigned sz = (unsigned)topology_prior.size();

    // First loop will be to determine largest value, which will be factored out
    // the second time through so that the total does not overflow
    double log_factored_out = 0.0;
	for (unsigned i = 1; i < sz; ++i)
		{
        double c = counts[i - 1];
        double nf = (double)nfactors[i - 1];
        double log_Tnm = log(c) + log_scaling_factor*nf;
		double log_prior = log_Tnm + topology_prior[i];
		v.push_back(log_prior);
        if (log_prior > log_factored_out)
            log_factored_out = log_prior;
		}

    // Now we can compute the total
    double total = 0.0;
    std::vector<double>::const_iterator it = v.begin();
    for (++it; it != v.end(); ++it)
		{
		total += exp((*it) - log_factored_out);
		}
	v[0] = log(total) + log_factored_out;

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructs a vector in which the element having index m (i = 0, 1, ..., max. num. internal nodes) represents the
|   natural logarithm of the number of tree topologies having m internal nodes. If `counts_dirty' is true, it recomputes
|   the `counts' vectors first. The 0th element of the returned vector holds the natural log of the total number of tree
|   topologies (log of the sum of all other elements).
*/
std::vector<double> PolytomyTopoPriorCalculator::GetLnCounts()
	{
    if (is_resolution_class_prior)
        counts_dirty = true;
	if (counts_dirty)
		Reset();

	std::vector<double> v;
	unsigned sz = ntax - (is_rooted ? 0 : 1);
	v.reserve(sz);
	v.push_back(log_total_count);

	//@POL could use a version of the transform algorithm here
	for (unsigned i = 1; i < sz; ++i)
		{
        double log_Tnm = log(counts[i - 1]) + log_scaling_factor*(double)(nfactors[i - 1]);
		v.push_back(log_Tnm);
		}

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Samples a resolution class (i.e. number of internal nodes) from the realized resolution class distribution. This
|   function is not very efficient because it calls PolytomyTopoPriorCalculator::GetRealizedResClassPriorsVect, resulting in an
|   unnecessary vector copy operation.
*/
unsigned PolytomyTopoPriorCalculator::sample(
  LotShPtr rng) /**< is the random number generator to use for sampling */
	{
    std::vector<double> v = GetRealizedResClassPriorsVect();
#if 0
    // Geez, what was I thinking when I wrote this?!
    double u = rng->Uniform(FILE_AND_LINE);
    double logu = (u > 0.0 ? log(u) : -DBL_MAX);
    double log_x = logu + v[0];
    double cum = 0.0;
    for (unsigned i = 1; i < v.size(); ++i)
        {
        cum += v[1];
        if (log_x <= cum)
            return i;
        }
#else
    double u = rng->Uniform(FILE_AND_LINE);
    //std::cerr << "PolytomyTopoPriorCalculator::sample: seed = " << rng->GetSeed() << ", u = " << u << std::endl;
    //std::cerr << "v[0]  = " << v[0] << std::endl;
    double z = v[0];
    double cum = 0.0;
    for (unsigned i = 1; i < v.size(); ++i)
        {
        //std::cerr << "v[" << i << "]  = " << v[i] << std::endl;
        cum += exp(v[i] - z);
        if (u <= cum)
            return i;
        }
#endif
    PHYCAS_ASSERT(0);
    return (unsigned)(v.size() - 1);
    }

}	// namespace phycas
