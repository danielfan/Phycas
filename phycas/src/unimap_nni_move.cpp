/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
|																			  |
|  This program is free software; you can redistribute it and/or modify		  |
|  it under the terms of the GNU General Public License as published by		  |
|  the Free Software Foundation; either version 2 of the License, or		  |
|  (at your option) any later version.										  |
|																			  |
|  This program is distributed in the hope that it will be useful,			  |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of			  |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the			  |
|  GNU General Public License for more details.								  |
|																			  |
|  You should have received a copy of the GNU General Public License along	  |
|  with this program; if not, write to the Free Software Foundation, Inc.,	  |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.				  |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/likelihood_models.hpp"
#if POLPY_NEWWAY
#	include "phycas/src/basic_tree_node.hpp"
#endif
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/mcmc_chain_manager.hpp"
#include "phycas/src/unimap_nni_move.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/internal_data.hpp"
#include "phycas/src/tip_data.hpp"

#include "boost/format.hpp"

namespace phycas
{

UniventManager * GetUniventManager(TreeNode * nd)
	{
	if (nd->IsTip())
		{
		TipData * td = nd->GetTipData();
		return (td ? td->getUniventManager() : NULL);
		}
	else
		{
		InternalData * id = nd->GetInternalData();
		return (id ? id->getUniventManager() : NULL);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool UnimapNNIMove::update()
	{
	// The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed UnimapNNIMove cannot be accepted without changing edge lengths, so it is best to bail out now.
	if (is_fixed)
		return false;

	tree->renumberInternalNodes(tree->GetNTips()); //@POL this should be somewhere else

	proposeNewState();
	ChainManagerShPtr p = chain_mgr.lock();

	double prev_posterior = 0.0;
	double curr_posterior = 0.0;
	if (is_standard_heating)
		{
		prev_posterior = heating_power*(prev_ln_like + prev_ln_prior);
		curr_posterior = heating_power*(curr_ln_like + curr_ln_prior);
		}
	else
		{
		prev_posterior = heating_power*prev_ln_like + prev_ln_prior;
		curr_posterior = heating_power*curr_ln_like + curr_ln_prior;
		}

	double ln_accept_ratio = curr_posterior - prev_posterior + getLnHastingsRatio() + getLnJacobian();
	//double lnu = std::log(rng->Uniform(FILE_AND_LINE));
	//bool accepted = (ln_accept_ratio >= 0.0 || lnu <= ln_accept_ratio);
	//bool accepted = (ln_accept_ratio >= 0.0 || std::log(rng->Uniform(FILE_AND_LINE)) <= ln_accept_ratio);
	double lnu = DBL_MAX;
	bool accepted = (ln_accept_ratio >= 0.0);
	if (!accepted)
		{
		double u = rng->Uniform(FILE_AND_LINE);
		lnu = std::log(u);
		accepted = (lnu <= ln_accept_ratio);
		}

	if (accepted)
		{
		p->setLastLnPrior(curr_ln_prior);
		p->setLastLnLike(curr_ln_like);
		accept();
		return true;
		}
	else
		{
		curr_ln_like	= p->getLastLnLike();
		curr_ln_prior	= p->getLastLnPrior();
		revert();
		return false;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Selects an internal node at random from a discrete uniform distribution with the constraint that the returned node
|   is not equal to the subroot (the sole child of the tip node serving as the root).
*/
TreeNode * UnimapNNIMove::randomInternalAboveSubroot()
    {
	// Avoiding the "subroot" node (only child of the tip serving as the root), so the number of 
	// acceptable nodes is one fewer than the number of internal nodes
	unsigned numAcceptableNodes = tree->GetNInternals() - 1;

	unsigned ypos = rng->SampleUInt(numAcceptableNodes);
	unsigned i = 0;
    TreeNode * nd = tree->GetFirstPreorder();
	for (; nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsInternal() && !nd->GetParentConst()->IsTipRoot())
			{
			if (i == ypos)
				break;
			++i;
			}
		}
	PHYCAS_ASSERT(nd->GetLeftChild() != NULL);
	PHYCAS_ASSERT(nd->GetParentConst() != NULL);
	PHYCAS_ASSERT(!nd->GetParent()->IsTipRoot());
    return nd;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapNNIMove::proposeNewState()
	{
	TreeNode * nd  = randomInternalAboveSubroot();
	PHYCAS_ASSERT(nd);
	TreeNode * ndP = nd->GetParent();
	PHYCAS_ASSERT(ndP);
	z = nd->FindNextSib();
	PHYCAS_ASSERT(z);
	PHYCAS_ASSERT(nd->CountChildren() == 2); // we haven't figured this out for polytomies
	PHYCAS_ASSERT(ndP->CountChildren() == 2); // we haven't figured this out for polytomies
	x = nd->GetLeftChild();
	PHYCAS_ASSERT(x);
	y = x->GetRightSib();
	PHYCAS_ASSERT(y);
	x_is_left = rng->Boolean();
	if (!x_is_left)
		std::swap(x,y);
	
	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);

    //temporary!
    likelihood->storeAllCLAs(tree);

	prev_ln_like = FourTaxonLnL(nd);
	prev_ln_prior = (x->IsInternal() ? p->calcInternalEdgeLenPriorUnnorm(x->GetEdgeLen()) : p->calcExternalEdgeLenPriorUnnorm(x->GetEdgeLen()));
	prev_ln_prior += (y->IsInternal() ? p->calcInternalEdgeLenPriorUnnorm(y->GetEdgeLen()) : p->calcExternalEdgeLenPriorUnnorm(y->GetEdgeLen()));
	prev_ln_prior += (z->IsInternal() ? p->calcInternalEdgeLenPriorUnnorm(z->GetEdgeLen()) : p->calcExternalEdgeLenPriorUnnorm(z->GetEdgeLen()));
	prev_ln_prior += p->calcInternalEdgeLenPriorUnnorm(nd->GetEdgeLen());
	prev_ln_prior += p->calcInternalEdgeLenPriorUnnorm(ndP->GetEdgeLen());

    x->SelectNode();
    y->SelectNode();
    z->SelectNode();
    nd->SelectNode();

    //likelihood->startTreeViewer(tree, "before");
	
	tree_manipulator.NNISwap(z, x);
	
    //likelihood->startTreeViewer(tree, "after");
	
    x->UnselectNode();
    y->UnselectNode();
    z->UnselectNode();
    nd->UnselectNode();

    curr_ln_like = FourTaxonLnL(nd);
	curr_ln_prior = (x->IsInternal() ? p->calcInternalEdgeLenPriorUnnorm(x->GetEdgeLen()) : p->calcExternalEdgeLenPriorUnnorm(x->GetEdgeLen()));
	curr_ln_prior += (y->IsInternal() ? p->calcInternalEdgeLenPriorUnnorm(y->GetEdgeLen()) : p->calcExternalEdgeLenPriorUnnorm(y->GetEdgeLen()));
	curr_ln_prior += (z->IsInternal() ? p->calcInternalEdgeLenPriorUnnorm(z->GetEdgeLen()) : p->calcExternalEdgeLenPriorUnnorm(z->GetEdgeLen()));
	curr_ln_prior += p->calcInternalEdgeLenPriorUnnorm(nd->GetEdgeLen());
	curr_ln_prior += p->calcInternalEdgeLenPriorUnnorm(ndP->GetEdgeLen());
	
	}

double UnimapNNIMove::FourTaxonLnL(TreeNode * nd)
	{
	PHYCAS_ASSERT(nd);
	TreeNode * ndP = nd->GetParent();
	PHYCAS_ASSERT(ndP);
	TreeNode * parent = ndP->GetParent();
	PHYCAS_ASSERT(parent);
	TreeNode * lower = ndP->FindNextSib();
	PHYCAS_ASSERT(lower);
	PHYCAS_ASSERT(nd->CountChildren() == 2); // we haven't figured this out for polytomies
	PHYCAS_ASSERT(ndP->CountChildren() == 2); // we haven't figured this out for polytomies
	TreeNode * upLeftNd = nd->GetLeftChild();
	PHYCAS_ASSERT(upLeftNd);
	TreeNode * upRightNd = upLeftNd->GetRightSib();
	PHYCAS_ASSERT(upRightNd);

	xTipData = allocTipDataFromUnivents(upLeftNd, true);  /* LEAK !!!!*/
	yTipData = allocTipDataFromUnivents(upRightNd, true);
	if (!x_is_left)
		std::swap(xTipData, yTipData);
	zTipData = allocTipDataFromUnivents(lower, true);
	wTipData = allocTipDataFromUnivents(parent, false);

	double lnlike = FourTaxonLnLFromCorrectTipDataMembers(nd);

    //DebugSaveNexusFile(xTipData, yTipData, zTipData, wTipData, lnlike);

    return lnlike;
	}

void UnimapNNIMove::DebugSaveNexusFile(TipData * xtd, TipData * ytd, TipData * ztd, TipData * wtd, double lnlike)
    {
    typedef boost::shared_array<const int8_t> StateArr;
    StateArr xdata = xtd->getTipStatesArray();
    StateArr ydata = ytd->getTipStatesArray();
    StateArr zdata = ztd->getTipStatesArray();
    StateArr wdata = wtd->getTipStatesArray();
    unsigned nchar = likelihood->getNPatterns();
    unsigned i;

    std::ofstream nxsf("tmp.nex");
    nxsf << "#nexus\n" << std::endl;
    nxsf << "begin data;" << std::endl;
    nxsf << "  dimensions ntax=4 nchar=" << nchar << ";" << std::endl;
    nxsf << "  format datatype=dna gap=- missing=?;" << std::endl;
    nxsf << "  matrix" << std::endl;
    nxsf << "    x ";
    for (i = 0; i < nchar; ++i)
        {
        if (xdata[i] == 0)
            nxsf << "a";
        else if (xdata[i] == 1)
            nxsf << "c";
        else if (xdata[i] == 2)
            nxsf << "g";
        else
            nxsf << "t";
        }
    nxsf << std::endl;
    nxsf << "    y ";
    for (i = 0; i < nchar; ++i)
        {
        if (ydata[i] == 0)
            nxsf << "a";
        else if (ydata[i] == 1)
            nxsf << "c";
        else if (ydata[i] == 2)
            nxsf << "g";
        else
            nxsf << "t";
        }
    nxsf << std::endl;
    nxsf << "    z ";
    for (i = 0; i < nchar; ++i)
        {
        if (zdata[i] == 0)
            nxsf << "a";
        else if (zdata[i] == 1)
            nxsf << "c";
        else if (zdata[i] == 2)
            nxsf << "g";
        else
            nxsf << "t";
        }
    nxsf << std::endl;
    nxsf << "    w ";
    for (i = 0; i < nchar; ++i)
        {
        if (wdata[i] == 0)
            nxsf << "a";
        else if (wdata[i] == 1)
            nxsf << "c";
        else if (wdata[i] == 2)
            nxsf << "g";
        else
            nxsf << "t";
        }
    nxsf << std::endl;
    nxsf << "  ;\nend;" << std::endl;

    nxsf << "\nbegin paup;" << std::endl;
    nxsf << "  log file=tmp.log start replace;" << std::endl;
    nxsf << "  set criterion=likelihood storebrlens;" << std::endl;
    nxsf << "end;" << std::endl;

    nxsf << "\nbegin trees;" << std::endl;
    nxsf << "  translate" << std::endl;
    nxsf << "    1 x," << std::endl;
    nxsf << "    2 y," << std::endl;
    nxsf << "    3 z," << std::endl;
    nxsf << "    4 w" << std::endl;
    nxsf << "    ;" << std::endl;
    TreeNode * xpar = x->GetParent();
    TreeNode * zpar = z->GetParent();
    if (xpar->GetParent() == zpar)
        nxsf << boost::str(boost::format("  utree curr = (x:%.8f, y:%.8f, (z:%.8f, w:%.8f):%.8f);") % x->GetEdgeLen() % y->GetEdgeLen() % z->GetEdgeLen() % zpar->GetEdgeLen() % xpar->GetEdgeLen()) << std::endl;
    else
        nxsf << boost::str(boost::format("  utree curr = (z:%.8f, y:%.8f, (x:%.8f, w:%.8f):%.8f);") % z->GetEdgeLen() % y->GetEdgeLen() % x->GetEdgeLen() % xpar->GetEdgeLen() % zpar->GetEdgeLen()) << std::endl;
    nxsf << "end;" << std::endl;

    nxsf << "\nbegin paup;" << std::endl;
    nxsf << boost::str(boost::format("  [!***** phycas lnL = %.8f *****]") % lnlike) << std::endl;
    if (model->getModelName().compare("JC69") == 0)
        {
        JC * p = dynamic_cast<JC *>(model.get());
        nxsf << "  lset nst=1 basefreq=equal;" << std::endl;
        }
    else if (model->getModelName().compare("HKY85") == 0)
        {
        HKY * p = dynamic_cast<HKY *>(model.get());
		const std::vector<double> freq = model->getStateFreqs();
        nxsf << boost::str(boost::format("  lset nst=2 variant=hky basefreq=(%.8f %.8f %.8f) tratio=%.8f;") % freq[0] % freq[1] % freq[2] % p->calcTRatio()) << std::endl;
        }
    else if (model->getModelName().compare("GTR") == 0)
        {
        GTR * p = dynamic_cast<GTR *>(model.get());
		const std::vector<double> freq = model->getStateFreqs();
        std::vector<double> rmat = p->getRelRates();
        nxsf << boost::str(boost::format("  lset nst=6 basefreq=(%.8f %.8f %.8f) rmatrix=(%.8f %.8f %.8f %.8f %.8f);") % freq[0] % freq[1] % freq[2] % rmat[0] % rmat[1] % rmat[2] % rmat[3] % rmat[4]) << std::endl;
        }
    nxsf << "  lscores 1 / userbrlens;" << std::endl;
    nxsf << "  log stop;" << std::endl;
    nxsf << "end;" << std::endl;

    nxsf.close();
    std::exit(0);
    }
	
double UnimapNNIMove::FourTaxonLnLFromCorrectTipDataMembers(TreeNode * nd)
	{
	InternalData * nd_internal_data = nd->GetInternalData();
	PHYCAS_ASSERT(nd_internal_data);
	CondLikelihoodShPtr nd_childCLPtr =  nd_internal_data->getChildCondLikePtr();
	CondLikelihoodShPtr nd_parentCLPtr =  nd_internal_data->getParentalCondLikePtr();

    likelihood->calcPMatTranspose(xTipData->getTransposedPMatrices(), xTipData->getConstStateListPos(), x->GetEdgeLen());
    likelihood->calcPMatTranspose(yTipData->getTransposedPMatrices(), yTipData->getConstStateListPos(), y->GetEdgeLen());
    likelihood->calcPMatTranspose(zTipData->getTransposedPMatrices(), zTipData->getConstStateListPos(), z->GetEdgeLen());
    likelihood->calcPMatTranspose(wTipData->getTransposedPMatrices(), wTipData->getConstStateListPos(), nd->GetParent()->GetEdgeLen());
	likelihood->calcCLATwoTips(*nd_childCLPtr, *xTipData, *yTipData);
	likelihood->calcCLATwoTips(*nd_parentCLPtr, *zTipData, *wTipData);

	likelihood->calcPMat(nd_internal_data->getMutablePMatrices(), nd->GetEdgeLen());
	const double * const *  childPMatrix = nd_internal_data->getConstPMatrices()[0];
	return HarvestLnLikeFromCondLikePar(nd_childCLPtr, nd_parentCLPtr, childPMatrix);
	}
	
double UnimapNNIMove::HarvestLnLikeFromCondLikePar(
  ConstCondLikelihoodShPtr focalCondLike, 
  ConstCondLikelihoodShPtr neighborCondLike, 
  const double * const * childPMatrix)
	{
	const LikeFltType * focalNdCLAPtr = focalCondLike->getCLA(); //PELIGROSO
	PHYCAS_ASSERT(focalNdCLAPtr);
	const unsigned num_patterns = likelihood->getNPatterns();
	const unsigned num_states = likelihood->getNStates();
		// Get state frequencies from model and alias rate category probability array for speed
	const double * stateFreq = &model->getStateFreqs()[0]; //PELIGROSO
	double lnLikelihood = 0.0;

	const double * focalNeighborCLAPtr = neighborCondLike->getCLA(); //PELIGROSO
	for (unsigned pat = 0; pat < num_patterns; ++pat)
		{
		double siteLike = 0.0;
		for (unsigned i = 0; i < num_states; ++i)
			{
			double neigborLike = 0.0;
			const double * childP_i = childPMatrix[i];
			for (unsigned j = 0; j < num_states; ++j)
				neigborLike += childP_i[j]*focalNeighborCLAPtr[j];
			siteLike += stateFreq[i]*focalNdCLAPtr[i]*neigborLike;
			}
		focalNdCLAPtr += num_states;
		focalNeighborCLAPtr += num_states;
		double site_lnL = std::log(siteLike);

#if defined(DO_UNDERFLOW_POLICY)
		underflow_policy.correctSiteLike(site_lnL, pat, focalCondLike);
		underflow_policy.correctSiteLike(site_lnL, pat, neighborCondLike);
#endif
		lnLikelihood += site_lnL;
		}
#if defined(DO_UNDERFLOW_POLICY)
	underflow_policy.correctLnLike(lnLikelihood, focalCondLike);
	underflow_policy.correctLnLike(lnLikelihood, neighborCondLike);
#endif

	return lnLikelihood;
	}

TipData * UnimapNNIMove::allocTipDataFromUnivents(TreeNode * nd , bool use_last)
	{
	PHYCAS_ASSERT(nd);
	const unsigned num_patterns = likelihood->getNPatterns();
	const UniventManager * um = GetUniventManager(nd);
	PHYCAS_ASSERT(um);
	PHYCAS_ASSERT(num_patterns == um->size());
	int8_t * tipSpecificStateCode = new int8_t[num_patterns];
	const int8_t num_states = (int8_t) likelihood->getNStates();

	for (unsigned site = 0; site < num_patterns; ++site)
		{
		const StateTimeList & state_time_pair_vec =  (*um)[site];
        PHYCAS_ASSERT(!state_time_pair_vec.empty());
		const StateTimePair & state_time_pair  = (use_last ? *state_time_pair_vec.rbegin() : *state_time_pair_vec.begin());
		int8_t globalStateCode = state_time_pair.first;
		PHYCAS_ASSERT (globalStateCode >= 0 && globalStateCode < (int8_t)likelihood->getNStates());
		tipSpecificStateCode[site] = globalStateCode;
		}

	std::vector<unsigned int> emptyStateListVec;
    TipData * retval = new TipData(	true,
					num_patterns,
					emptyStateListVec,												// stateListPosVec
					boost::shared_array<const int8_t>(tipSpecificStateCode),	// stateCodesShPtr
					1,													// number of relative rate categories
					likelihood->getNStates(),													// number of states in the model
					NULL,
					true,														// managePMatrices
					likelihood->getCondLikelihoodStorage());
    return retval;
	}
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapNNIMove::revert()
	{
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void UnimapNNIMove::accept()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The default constructor sets `lambda' to the default value (0.2), sets `topol_changed' to false, and `m' and `mstar'
|	to 0.0. All other data members are automatically initialized (shared pointers) or are initialized via a call to 
|	reset().
*/
UnimapNNIMove::UnimapNNIMove() : MCMCUpdater(),
  x(0),
  y(0),
  z(0),
  xTipData(0),
  yTipData(0),
  zTipData(0),
  wTipData(0)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move.
*/
double UnimapNNIMove::getLnHastingsRatio() const
	{
	PHYCAS_ASSERT(false);
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Jacobian for this move.
*/
double UnimapNNIMove::getLnJacobian() const
	{
	PHYCAS_ASSERT(false);
	return 0.0;
	}

}	// namespace phycas
