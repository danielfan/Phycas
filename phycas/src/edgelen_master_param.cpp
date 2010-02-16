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
|	Member function that exists only to facilitate using boost::lambda::bind to be used in the getLnPrior() function. It
|	returns the log of the working prior probability density evaluated at the current edge length associated with `nd'.
*/
double EdgeLenMasterParam::lnWorkingPriorOneEdge(bool temp_debugging, TreeNode & nd) const
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
		    retval = working_prior->GetLnPDF(v);
		    }
	    catch(XProbDist &)
		    {
		    PHYCAS_ASSERT(0);
		    }
		if (temp_debugging)
			std::cerr << boost::str(boost::format("%.8f <-- %.8f <-- %s <-- %s") % retval % v % working_prior->GetDistributionDescription() % getName()) << std::endl;//temp
	    return retval;
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Loops through all nodes in the tree and computes the log of the working prior for each edge that it is responsible
|	for (according to its `edgeLenType'). Returns sum of these log working prior values.
*/
double EdgeLenMasterParam::recalcWorkingPrior(bool temp_debug) const
	{
    double lnwp = std::accumulate(tree->begin(), tree->end(), 0.0,
	    boost::lambda::_1 += boost::lambda::bind(&EdgeLenMasterParam::lnWorkingPriorOneEdge, this, temp_debug, boost::lambda::_2));
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
	PHYCAS_ASSERT(isPriorSteward());			// only prior stewards should be building working priors
	PHYCAS_ASSERT(tree);	
	PHYCAS_ASSERT(tree->GetNInternals() > 0);	
	PHYCAS_ASSERT(tree->GetNTips() > 0);	
	PHYCAS_ASSERT(tree->GetNInternals() + tree->GetNTips() == tree->GetNNodes());	
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
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Use samples in `fitting_sample' to parameterize `working_prior'. This function is called during s-cubed style
|	steppingstone sampling after the initial phase of sampling from the posterior so that the working prior can be
|	used for the remaining phases. Assumes `fitting_sample' has more than 1 element. Assigns a GammaDistribution object
|	to `working_prior'.
*/
void EdgeLenMasterParam::finalizeWorkingPrior()
	{
	PHYCAS_ASSERT(isPriorSteward());	// only prior stewards should be building working priors
	fitGammaWorkingPrior();
	
	// correct variance of working prior to reflect number of edge lengths 
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
	
	// std::cerr << boost::str(boost::format("@@@@@@@@@ --> working prior of %s adjusted to have mean = %g, variance = %g, based on %g %s edge lengths") % getName() % m % v % num_edgelens % (edgeLenType == internal ? "internal" : (edgeLenType == external ? "external" : "total") ) ) << std::endl;
	}
}
