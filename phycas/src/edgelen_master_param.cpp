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
#include <boost/format.hpp>

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
  : MCMCUpdater(), edgeLenType(t), use_edge_specific_working_priors(false), min_working_prior_sample_size(10)
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
|	Sets the minimum acceptable sample size for creating a split-specific edge length working prior. If the actual 
|	sample size is less than the supplied value `n', a generic working prior will be used for that particular split.
|	Assumes that `n' is greater than 1.
*/
void EdgeLenMasterParam::setMinWorkingPriorSampleSize(
  unsigned n)	/**< is the minimum sample size to set */
	{
	PHYCAS_ASSERT(n > 1);
	min_working_prior_sample_size = n;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	If `use_it' is true, then split-specific edge length working priors will be constructed for all edges seen during
|	an MCMC analysis. If `use_it' is false, then a single generic edge length working prior will be used for all edges.
*/
void EdgeLenMasterParam::useEdgeSpecificWorkingPriors(
  bool use_it)	/**< should be true if you wish to use split-specific edge length working priors */
	{
	use_edge_specific_working_priors = use_it;
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
			//std::cerr << "+++++ " << nd.GetNodeNumber() << " --> " << v << " --> " << retval << std::endl;	//@@@
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
	std::string typestr = "all";
	if (edgeLenType == internal)
		typestr = "internal";
	else if (edgeLenType == external)
		typestr = "external";
		
	if (prior && working_prior)
		{
		s += boost::str(boost::format("%s will be used for lengths of %s edges") % working_prior->GetDistributionDescription() % typestr);
		}
	else if (mv_prior && mv_working_prior)
		{
		s += boost::str(boost::format("%s will be used for lengths of %s edges") % mv_working_prior->GetDistributionDescription() % typestr);
		}
	else 
		{
		s += "no generic edge working prior exists at this point";
		}

	if (use_edge_specific_working_priors)
		{
		s += "\n    Here are the working priors that will be used for edges already seen:";
		for (WorkingPriorMapConstIter it = edge_working_prior.begin(); it != edge_working_prior.end(); ++it)
			{
			if (it->second.wp)
				{
				s += "\n      ";
				s += it->second.wp->GetDistributionDescription();
				s += " will be used for ";
				s += it->first.CreateNewickRepresentation();
				}
			}
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
		
		bool use_generic = true;
		if (use_edge_specific_working_priors)
			{
			const Split & s = nd.GetSplitConst();
			if (edge_working_prior.find(s) != edge_working_prior.end())
				{
				// Found edge-specific working prior for this edge 
				WorkingPriorMapConstIter ewp = edge_working_prior.find(s);
				if ((*ewp).second.wp)
					{
					try 
						{
						// edge_working_prior[s] retrieves EdgeWorkingPrior struct for this split, of 
						// which 'second' is the working prior distribution shared pointer
						retval = (*ewp).second.wp->GetLnPDF(v);
						}
					catch(XProbDist &)
						{
						PHYCAS_ASSERT(0);
						}
					use_generic = false;
					}
				}
			}
		if (use_generic)
			{
			// There is no edge-specific working prior for this edge, so use generic one
			PHYCAS_ASSERT(working_prior);
			try 
				{
				retval = working_prior->GetLnPDF(v);
				}
			catch(XProbDist &)
				{
				PHYCAS_ASSERT(0);
				}
			}

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

	double lnwp = 0.0;
	if (!isFixed())
		{
		for (preorder_iterator nd = tree->begin(); nd != tree->end(); ++nd)
			{
			lnwp += lnWorkingPriorOneEdge(*nd, nd->GetEdgeLen());
			}
		}
	return lnwp;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the joint log prior over a set of edges in the associated tree and sets `curr_ln_prior'. The set of edges
|   included depends on the value of `edgeLenType', which can be `internal', `external' or `both'.
*/
double EdgeLenMasterParam::recalcPrior()
	{
	curr_ln_prior = std::accumulate(tree->begin(), tree->end(), 0.0,
		boost::lambda::_1 += boost::lambda::bind(&EdgeLenMasterParam::lnPriorOneEdge, this, boost::lambda::_2));
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
		
		if (use_edge_specific_working_priors)
			{
			// Second, use the current edge lengths to educate the working prior specific to the edges that
			// happen to be in the current tree. If the tree topology is fixed, then all working priors will
			// be specific working priors stored in the `edge_working_prior' map, and in this case the 
			// working prior stored in the data member `working_prior' will never be needed. If the topology
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
			}
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
		
		if (use_edge_specific_working_priors)
			{
			for (WorkingPriorMapIter it = edge_working_prior.begin(); it != edge_working_prior.end(); ++it)
				{
				if ((*it).second.fs.size() < min_working_prior_sample_size)
					{
					// Sample size is not great enough to create a split-specific edge length working prior
					// Reset wp to indicate that no working prior actually exists for this split. If wp does not point to anything, 
					// then a generic working prior will be used in EdgeLenMasterParam::lnWorkingPriorOneEdge
					(*it).second.wp.reset();
					}
				else 
					{
					// We have a sufficient sample size to create a split-specific edge length working prior
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
				}
			}
		}
	}
}
