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
	std::cerr << "\n>>>>> EdgeLenMasterParam dying..." << std::endl;
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
		    retval = prior->GetRelativeLnPDF(v);
		    }

	    return retval;
        }
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

}
