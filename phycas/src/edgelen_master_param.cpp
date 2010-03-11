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

#include <numeric>										// for std::accumulate
#include "mcmc_param.hpp"
#include "phycas/src/basic_tree.hpp"				// for Tree::begin() and Tree::end()
#include "phycas/src/mcmc_chain_manager.hpp"

// these were at the top of basic_tree.inl
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets both `has_slice_sampler' and `is_move' to false. It also sets the value of `edgeLenType' data 
|   member, which can be `internal', `external' or `both'. If `edgeLenType' is `internal', this EdgeLenMasterParam will
|   only calculate the prior for internal edges. If `edgeLenType' is `external', this EdgeLenMasterParam will only 
|   calculate the prior for external edges. If `edgeLenType' is `both', this EdgeLenMasterParam will calculate the 
|   prior for all edges in the tree.
*/
EdgeLenMasterParam::EdgeLenMasterParam(
  EdgeLenMasterParam::EdgeLenType t)    /**> is the edge length type (internal, external or both) */
  : MCMCUpdater(), edgeLenType(t)
	{
	has_slice_sampler = false;
	is_move = false;
	is_master_param = true;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor.
*/
EdgeLenMasterParam::~EdgeLenMasterParam()
	{
	//std::cerr << "\n>>>>> EdgeLenMasterParam dying..." << std::endl;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Member function that exists only to facilitate using boost::lambda::bind to be used in the getLnPrior() function. It
|	returns the log of the probability density evaluated at the current edge length associated with `nd'.
*/
double EdgeLenMasterParam::lnPriorOneEdge(TreeNode & nd) const
	{
    bool skip = (nd.IsTipRoot())
                || ((edgeLenType == EdgeLenMasterParam::internal) && (!nd.IsInternal()))
                || ((edgeLenType == EdgeLenMasterParam::external) && (nd.IsInternal()));
	if (skip)
        {
		return 0.0;
        }
    else
        {
        double v = nd.GetEdgeLen();

	    double retval = 0.0;
	    try 
		    {
		    retval = prior->GetLnPDF(v);
		    }
	    catch(XProbDist &)
		    {
		    PHYCAS_ASSERT(0);
		    }

	    return retval;
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reports the distribution being used for both the generic `working_prior' as well as the edge-specific working
|	priors for edges that were sampled before the working priors were finalized.
*/
std::string EdgeLenMasterParam::getWorkingPriorDescr() const
	{
	std::string s;
	if (prior && working_prior)
		{
		s += working_prior->GetDistributionDescription();
		s += " will be used for lengths of edges not yet seen";
		}
	else if (mv_prior && mv_working_prior)
		{
		s += mv_working_prior->GetDistributionDescription();
		s += " will be used for lengths of edges not yet seen";
		}
	else 
		{
		s += "no generic edge working prior exists at this point";
		}
	
	s += "\n    Here are the working priors that will be used for edges already seen:";
	for (WorkingPriorMapConstIter it = edge_working_prior.begin(); it != edge_working_prior.end(); ++it)
		{
		s += "\n      ";
		s += it->second.wp->GetDistributionDescription();
		s += " will be used for ";
		s += it->first.CreateNewickRepresentation();
		}

	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Member function that exists only to facilitate using boost::lambda::bind to be used in the getLnPrior() function. It
|	returns the log of the working prior probability density evaluated at the current edge length associated with `nd'.
*/
double EdgeLenMasterParam::lnWorkingPriorOneEdge(const TreeNode & nd, double v) const
	{
    bool skip = (nd.IsTipRoot())
                || ((edgeLenType == EdgeLenMasterParam::internal) && (!nd.IsInternal()))
                || ((edgeLenType == EdgeLenMasterParam::external) && (nd.IsInternal()));
	if (skip)
        {
		return 0.0;
        }
    else
        {
		double retval = 0.0;
		
#if USING_EDGE_SPECIFIC_WORKING_PRIORS
		const Split & s = nd.GetSplitConst();
		if (edge_working_prior.find(s) == edge_working_prior.end())
			{
			// There is no edge-specific working prior for this edge, so use generic one
			std::cerr << boost::str(boost::format("@@@@@@@@@@ Using generic working prior for edge having this split: %s") % s.CreateNewickRepresentation()) << std::endl;//temp
			PHYCAS_ASSERT(working_prior);
			try 
				{
				retval = working_prior->GetLnPDF(v);
				}
			catch(XProbDist &)
				{
				PHYCAS_ASSERT(0);
				}
			//std::cerr << boost::str(boost::format("%.8f <-- %.8f <-- %s <-- %s") % retval % v % working_prior->GetDistributionDescription() % getName()) << std::endl;//temp
			}
		else 
			{
			// Found edge-specific working prior for this edge 
			WorkingPriorMapConstIter ewp = edge_working_prior.find(s);
			try 
				{
				// edge_working_prior[s] retrieves EdgeWorkingPrior struct for this split, of 
				// which 'second' is the working prior distribution shared pointer
				PHYCAS_ASSERT((*ewp).second.wp);
				retval = (*ewp).second.wp->GetLnPDF(v);
				}
			catch(XProbDist &)
				{
				PHYCAS_ASSERT(0);
				}
			//std::cerr << boost::str(boost::format("%.8f <-- %.8f <-- %s <-- %s") % retval % v % (*ewp).second.wp->GetDistributionDescription() % getName()) << std::endl;//temp
			}

#else
		PHYCAS_ASSERT(working_prior);
	    try 
		    {
		    retval = working_prior->GetLnPDF(v);
		    }
	    catch(XProbDist &)
		    {
		    PHYCAS_ASSERT(0);
		    }
		//std::cerr << boost::str(boost::format("%.8f <-- %.8f <-- %s <-- %s") % retval % v % working_prior->GetDistributionDescription() % getName()) << std::endl;//temp
#endif

		return retval;
		}
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Loops through all nodes in the tree and computes the log of the working prior for each edge that it is responsible
|	for (according to its `edgeLenType'). Returns sum of these log working prior values.
*/
double EdgeLenMasterParam::recalcWorkingPrior() const
	{
	// don't need to check isPriorSteward() (as in MCMCUpdater::recalcWorkingPrior) because we know
	// that EdgeLenMasterParam is a prior steward

#if USING_EDGE_SPECIFIC_WORKING_PRIORS
	double lnwp = 0.0;
	if (!isFixed())
		{
		for (preorder_iterator nd = tree->begin(); nd != tree->end(); ++nd)
			{
			lnwp += lnWorkingPriorOneEdge(*nd, nd->GetEdgeLen());
			}
		}
	return lnwp;
#else
	if (!isFixed())
		{
	    double lnwp = std::accumulate(tree->begin(), tree->end(), 0.0,
		    boost::lambda::_1 += boost::lambda::bind(&EdgeLenMasterParam::lnWorkingPriorOneEdge, this, boost::lambda::_2));
		}
	return lnwp;
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the joint log prior over a set of edges in the associated tree and sets `curr_ln_prior'. The set of edges
|   included depends on the value of `edgeLenType', which can be `internal', `external' or `both'.
*/
double EdgeLenMasterParam::recalcPrior()
	{
	curr_ln_prior = 0.0;
	if (!chain_mgr.lock()->getSAMCLikelihoodOnly())
		{
	    curr_ln_prior = std::accumulate(tree->begin(), tree->end(), 0.0,
		    boost::lambda::_1 += boost::lambda::bind(&EdgeLenMasterParam::lnPriorOneEdge, this, boost::lambda::_2));
		}
	return curr_ln_prior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If using a hyperprior model in which the mean of the edge length prior is itself a parameter in the model, when the 
|	value of the hyperparameter changes, this change must be propagated to all EdgeLenParam objects so that they can
|	calculate their prior correctly. This function provides the means for changing the mean and variance of the prior
|	distribution. It simply calls prior->SetMeanAndVariance(), and then recalculates `curr_ln_prior' using the new prior 
|	distribution.
*/
void EdgeLenMasterParam::setPriorMeanAndVariance(double m, double v)
	{
	MCMCUpdater::setPriorMeanAndVariance(m, v);
	recalcPrior();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to add the current sum of edge lengths to the data already stored in `fitting_sample'.
|	Depending on the value of the data member `edgeLenType' (internal, external, or both), the value stored will be the
|	sum of internal, external or all edge lengths, respectively.
*/
void EdgeLenMasterParam::educateWorkingPrior()
	{
	if (!isFixed())
		{
		PHYCAS_ASSERT(isPriorSteward());			// only prior stewards should be building working priors
		PHYCAS_ASSERT(tree);	
		PHYCAS_ASSERT(tree->GetNInternals() > 0);	
		PHYCAS_ASSERT(tree->GetNTips() > 0);	
		PHYCAS_ASSERT(tree->GetNInternals() + tree->GetNTips() == tree->GetNNodes());
		
		// First, use current edge lengths to educate the generic working prior. This prior is stored in 
		// the data member `working prior' which is provided by the base class MCMCUpdater. It is used for
		// any branch that does not have a specific working prior stored in the map data member edge_working_prior
		double edgelen_sum = 0.0;
		double num_edgelens = 0.0;
		if (edgeLenType == internal)
			{
			edgelen_sum = tree->internalEdgeLenSum();
			num_edgelens = (double)(tree->GetNInternals() - 1);
			}
		else if (edgeLenType == external)
			{
			edgelen_sum = tree->externalEdgeLenSum();
			num_edgelens = (double)tree->GetNTips();
			}
		else
			{
			edgelen_sum = tree->EdgeLenSum();
			num_edgelens = (double)(tree->GetNNodes() - 1);
			}
		double edgelen_mean = edgelen_sum/num_edgelens;
		fitting_sample.push_back(edgelen_mean);
		
#if USING_EDGE_SPECIFIC_WORKING_PRIORS
		// Second, use the current edge lengths to educate the working prior specific to the edges that
		// happen to be in the current tree. If the tree topology is fixed, then all working priors will
		// be specific working priors stored in the `edge_working_prior' map, and in this case the 
		// workign prior stored in the data member `working_prior' will never be needed. If the topology
		// is not fixed, however, then after the working priors in the `edge_working_prior' map have been
		// finalized, we may (probably will) encounter splits that have never been seen before. For these
		// edges, the MCMCUpdater::working_prior will be used.
		//std::cerr << "--------- processing tree ----------" << std::endl;
		unsigned ntips = tree->GetNTips();
		tree->RecalcAllSplits(ntips);
		for (preorder_iterator nd = tree->begin(); nd != tree->end(); ++nd)
			{
			//Split & sref = nd->GetSplit();
			//std::cerr << "@@@@@@@@@@ nd = " << nd->GetNodeNumber() << ", s = " << sref.CreateNewickRepresentation() << std::endl;
			bool skip = (nd->IsAnyRoot())
						|| ((edgeLenType == EdgeLenMasterParam::internal) && (!nd->IsInternal()))
						|| ((edgeLenType == EdgeLenMasterParam::external) && (nd->IsInternal()));
			if (!skip)
				{
				double edgelen = nd->GetEdgeLen();
				Split s = nd->GetSplit();
				WorkingPriorMapIter it = edge_working_prior.find(s);
				if (it == edge_working_prior.end())
					{
					// Might need to check inverse of split - following line is designed to find out if this will be needed
					//PHYCAS_ASSERT(edge_working_prior.find(s.InvertSplit()) == edge_working_prior.end());
					
					// Split s was NOT found in the map, so add a new element
					//std::cerr << "@@@@@@@@@@ adding a new element to edge_working_prior map" << std::endl;
					EdgeWorkingPrior e;
					e.fs.push_back(edgelen);
					edge_working_prior.insert(WorkingPriorMapPair(s,e));
					}
				else 
					{
					// Split s WAS found in the map, so just add edgelen to the existing fitting sample vector
					EdgeWorkingPrior & e = (*it).second;
					e.fs.push_back(edgelen);
					//unsigned n = (unsigned)e.fs.size();
					//std::cerr << boost::str(boost::format("@@@@@@@@@@ extending existing edge_working_prior element that now has %d elements") % n)<< std::endl;
					}
				}
			}
#endif
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use samples in `fitting_sample' to parameterize `working_prior'. This function is called during s-cubed style
|	steppingstone sampling after the initial phase of sampling from the posterior so that the working prior can be
|	used for the remaining phases. Assumes `fitting_sample' has more than 1 element. Assigns a GammaDistribution object
|	to `working_prior'.
*/
void EdgeLenMasterParam::finalizeWorkingPrior()
	{
	if (!isFixed())
		{
		PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
		fitGammaWorkingPrior();
		
		// correct variance of working prior to reflect number of edge lengths 
		PHYCAS_ASSERT(working_prior);
		double m = working_prior->GetMean();
		double v = working_prior->GetVar();
		double num_edgelens = 0.0;
		if (edgeLenType == internal)
			num_edgelens = (double)(tree->GetNInternals() - 1);
		else if (edgeLenType == external)
			num_edgelens = (double)tree->GetNTips();
		else
			num_edgelens = (double)(tree->GetNNodes() - 1);
		v *= num_edgelens;
		working_prior->SetMeanAndVariance(m, v);
		
#if USING_EDGE_SPECIFIC_WORKING_PRIORS
		for (WorkingPriorMapIter it = edge_working_prior.begin(); it != edge_working_prior.end(); ++it)
			{
			PHYCAS_ASSERT((*it).second.fs.size() > 1);
			double n = (double)(*it).second.fs.size();
			double sum = 0.0;
			double sum_of_squares = 0.0;
			for (double_vect_t::iterator i = (*it).second.fs.begin(); i != (*it).second.fs.end(); ++i)
				{
				double v = (*i);
				sum += v;
				sum_of_squares += v*v;
				}
			double mean = sum/n;	// shape*scale
			double variance = (sum_of_squares - n*mean*mean)/(n - 1.0);	// shape*scale^2
			double scale = variance/mean;
			double shape = mean/scale;
			(*it).second.wp = ProbDistShPtr(new GammaDistribution(shape, scale));
			
			//std::cerr << boost::str(boost::format("@@@@@@@@@ split %s has %d elements: mean = %g, variance = %g, shape = %g, scale = %g") 
			//	% (*it).first.CreateNewickRepresentation() 
			//	% (*it).second.fs.size() 
			//	% mean
			//	% variance
			//	% shape
			//	% scale)
			//	<< std::endl;
			}
#endif
		}
	}
}
