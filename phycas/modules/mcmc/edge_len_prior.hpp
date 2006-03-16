#if !defined(PHO_EDGELENPRIOR_H)
#define PHO_EDGELENPRIOR_H

class ProbabilityDistribution;
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include "phycas/trees/tree_node.hpp" //inlines only
#include "phycas/trees/tree_inl.hpp"
typedef boost::shared_ptr<ProbabilityDistribution> ProbabilityDistributionShPtr;

//#define HYPERDOOF

#if defined(HYPERDOOF)
#	include <iostream>
#	include <fstream>
#endif

class EdgeLenPriorCalculator
	{
	public:
				EdgeLenPriorCalculator();
		void	SetProbabilityDistribution(ProbabilityDistributionShPtr d);
		double	GetLnEdgeLenPrior(Tree &t);
		void	SetAllEdgeLengths(Tree &t, double common_edgelen);
		bool	CalcLnPrior(TreeNode *nd);
		bool	SetEdgeLen(TreeNode *nd);

	private:
		ProbabilityDistributionShPtr	probDist;
		double							ln_prior;
		double							edge_len;
#if defined(HYPERDOOF)
		double num_edgelens;
		double sum_edgelens;
#endif
	};

typedef boost::shared_ptr<EdgeLenPriorCalculator> EdgeLenPriorCalculatorShPtr;

inline EdgeLenPriorCalculator::EdgeLenPriorCalculator()
	{
	ln_prior = 0.0;
	edge_len = 1.0;
	}

inline void EdgeLenPriorCalculator::SetProbabilityDistribution(ProbabilityDistributionShPtr d)
	{
	probDist = d;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets every edge length in the tree provided to `common_edgelen'. Useful for debugging, but perhaps this
|	functionality really belongs in the class TreeManip.
*/
inline void EdgeLenPriorCalculator::SetAllEdgeLengths(Tree &t, double common_edgelen)
	{
	edge_len = common_edgelen;
	t.PostorderTraverse(boost::function1<bool, TreeNode*>(std::bind1st(std::mem_fun(&EdgeLenPriorCalculator::SetEdgeLen), this)));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the floating point edge length for node `nd' to `len'.
*/
inline bool EdgeLenPriorCalculator::SetEdgeLen(TreeNode *nd)
	{
	if (!nd->IsRoot())
		{
		nd->SetFltEdgeLen(edge_len);
		}
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the joint edge length prior using the prior distribution stored in `probDist'. Traverses the tree using
|	Tree::PostorderTraverse, and sums the logs of the prior calculated for each edge length encountered during the 
|	traversal.
*/
inline double EdgeLenPriorCalculator::GetLnEdgeLenPrior(Tree &t)
	{
	ln_prior = 0.0;

	//POL 29 Dec 2003 
	// PostorderTraverse used take a pointer to TreeNodeListener parameter, then call the TreeNodeListener object's
	// NodeWork virtual function. This restricted each TreeNodeListener to just one NodeWork function, and the name
	// of the function gave no clue about the job it performed. Thus, the TreeNodeListener system was replaced 
	// by a functor system. Here, a functor is constructed and passed to the PostorderTraverse function. The functor
	// object is a wrapper around a member function of this EdgeLenPriorCalculator object, namely the 
	// EdgeLenPriorCalculator::CalcLnPrior member function, which takes a TreeNode pointer parameter and returns
	// a bool. Searching for PostorderTraverse will uncover the other instances where TreeNodeListener was formerly
	// used.
	//
	// The substitution of functors for TreeNodeListener seems to slow things down a tiny bit (run times are 
	// about a third of a percent longer). In one test, the new way took 659.348 sec whereas the old way required
	// 657.165 sec.
	//

#if defined(HYPERDOOF)
	num_edgelens = 0.0;
	sum_edgelens = 0.0;
#endif

	t.PostorderTraverse(boost::function1<bool, TreeNode*>(std::bind1st(std::mem_fun(&EdgeLenPriorCalculator::CalcLnPrior), this)));

#if defined(HYPERDOOF)
	char tmps[256];
	std::ofstream hyperdoof("hyperdoof.txt", std::ios::out | std::ios::app);
	sprintf(tmps, "  edge length prior mean = %.6f", probDist->GetMean());
	hyperdoof << tmps << std::endl;
	sprintf(tmps, "  mean edge length = %.6f", (sum_edgelens/num_edgelens));
	hyperdoof << tmps << std::endl;
	hyperdoof.close();
#endif

	return ln_prior;
	}

inline bool EdgeLenPriorCalculator::CalcLnPrior(TreeNode *nd)
	{
	if (!nd->IsRoot())
		{
		double p = probDist->GetLnPDF(nd->GetEdgeLen().f);

#if defined(HYPERDOOF)
		sum_edgelens += nd->GetEdgeLen().f;
		num_edgelens += 1.0;
#endif

		ln_prior += p;
		}
	return true;
	}

#endif
