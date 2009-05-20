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
#include <fstream>
#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/basic_tree_node.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/unimap_nni_move.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/internal_data.hpp"
#include "phycas/src/tip_data.hpp"

#include "boost/format.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapNNIMove::setLot(LotShPtr p)
	{
	MCMCUpdater::setLot(p);
	gammaDist.SetLot(rng.get());
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

	std::cerr << "****** UnimapNNIMove::update" << std::endl;

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

	double ln_accept_ratio = curr_posterior - prev_posterior + getLnHastingsRatio();
	//std::cerr << "ln_accept_ratio = curr_posterior - prev_posterior + ln_density_reverse_move - ln_density_forward_move;\n";
	//std::cerr << ln_accept_ratio << "   " << curr_posterior << "   " <<prev_posterior << "   " << ln_density_reverse_move << "   " << ln_density_forward_move << "\n\n";
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
	//std::cerr << "lnu    ln_accept_ratio = curr_posterior - prev_posterior + ln_density_reverse_move - ln_density_forward_move;\n";
	//std::cerr << lnu << "   " << ln_accept_ratio << "   " << curr_posterior << "   " <<prev_posterior << "   " << ln_density_reverse_move << "   " << ln_density_forward_move << "\n\n";

    if (save_debug_info)
        debug_info = str(boost::format("swapping %d <-> %d (%s, lnR = %.5f)") % x->GetNodeNumber() % z->GetNodeNumber() % (accepted ? "accepted" : "rejected") % ln_accept_ratio);
    
    if (accepted)
		accept();
	else
		revert();

	return accepted;
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

void UnimapNNIMove::calculatePairwiseDistances()
	{
	dXY = dWX =  dXZ =  dWY =  dYZ =  dWZ = 0.0;
	/* This is called before the swap so "x" is "ySis" and "z" is "wSis" */
	const int8_t * xStates = ySisTipData->getTipStatesArray().get();
	const int8_t * yStates = yTipData->getTipStatesArray().get();
	const int8_t * wStates = wTipData->getTipStatesArray().get();
	const int8_t * zStates = wSisTipData->getTipStatesArray().get();
	const unsigned num_patterns = likelihood->getNPatterns();
	for (unsigned i = 0; i < num_patterns; ++i)
		{
		const int8_t x = *xStates++;
		const int8_t y = *yStates++;
		const int8_t w = *wStates++;
		const int8_t z = *zStates++;
		if (x != y)
			{
			dXY += 1.0;
			if (x != w)
				{
				dWX += 1.0;
				if (x != z)
					{
					dXZ += 1.0;
					if (y != w)
						dWY += 1.0;
					if (y != z)
						dYZ += 1.0;
					if (w != z)
						dWZ += 1.0;
					}
				else
					{
					dYZ += 1.0;
					dWZ += 1.0;
					if (y != w)
						dWY += 1.0;
					}
				}
			else
				{
				dWY += 1.0;
				if (x != z)
					{
					dXZ += 1.0;
					dWZ += 1.0;
					}
				if (y != z)
					dYZ += 1.0;
				}
			
			}
		else
			{
			if (x != w)
				{
				dWX += 1.0;
				dWY += 1.0;
				if (x != z)
					{
					dXZ += 1.0;
					dYZ += 1.0;
					}
				if (w != z)
					dWZ += 1.0;
				}
			else
				{			
				if (x != z)
					{
					dXZ += 1.0;
					dYZ += 1.0;
					dWZ += 1.0;
					}
				}
}
		}
	double dnp = (double) num_patterns;
	dXY /= dnp;
	dWX /= dnp;
	dXZ /= dnp;
	dWY /= dnp;
	dYZ /= dnp;
	dWZ /= dnp;
	}

void UnimapNNIMove::calculateProposalDist(bool before_swap)
	{
	if (!before_swap)
		{
		propMeanX = std::max(min_edge_len_mean, (2*dXY + dXZ + dWX - dYZ - dWY)/4.0);
		propMeanY = std::max(min_edge_len_mean, (2*dXY + dYZ + dWY - dXZ - dWX)/4.0);
		propMeanZ = std::max(min_edge_len_mean, (2*dWZ + dXZ + dYZ - dWX - dWY)/4.0); 
		propMeanW = std::max(min_edge_len_mean, (2*dWZ + dWX + dWY - dXZ - dXY)/4.0);
		propMeanInternal = std::max(min_edge_len_mean, (dXZ + dWX + dYZ + dWY - 2*dXY - 2*dWZ)/4.0);
		}
	else
		{
		propMeanX = std::max(min_edge_len_mean, (2*dWX + dXZ + dXY - dWZ - dWY)/4.0); 
		propMeanY = std::max(min_edge_len_mean, (2*dYZ + dXY + dWY - dXZ - dWZ)/4.0);
		propMeanW = std::max(min_edge_len_mean, (2*dWX + dWZ + dWY - dXZ - dXY)/4.0);
		propMeanZ = std::max(min_edge_len_mean, (2*dYZ + dXZ + dWZ - dXY - dWY)/4.0);
		propMeanInternal = std::max(min_edge_len_mean, (dXZ + dWZ + dXY + dWY - 2*dYZ - 2*dWX)/4.0);
		}
	}

double UnimapNNIMove::calcProposalLnDensity(double mean, double x)
	{
	const double variance = (edge_len_prop_cv*edge_len_prop_cv*mean*mean);
	gammaDist.SetMeanAndVariance(mean, variance);
	const double d = gammaDist.GetLnPDF(x);
	//std::cerr << x << " from Gamma(" << mean << ", " << variance << ") = " << d << "\n";
	return d;
	}

double UnimapNNIMove::proposeEdgeLen(double mean)
	{
	const double variance = (edge_len_prop_cv*edge_len_prop_cv*mean*mean);
	gammaDist.SetMeanAndVariance(mean, variance);
	return gammaDist.Sample();
	}

void UnimapNNIMove::AlterBranchLengths(ChainManagerShPtr & p)
	{
	prev_x_len = x->GetEdgeLen();
	prev_y_len = y->GetEdgeLen();
	prev_z_len = z->GetEdgeLen();
	prev_nd_len = nd->GetEdgeLen();
	prev_ndP_len = ndP->GetEdgeLen();
	
	prev_ln_prior = calcEdgeLenLnPrior(*x, prev_x_len, p);
	prev_ln_prior += calcEdgeLenLnPrior(*y, prev_y_len, p);
	prev_ln_prior += calcEdgeLenLnPrior(*z, prev_z_len, p);
	prev_ln_prior += calcEdgeLenLnPrior(*nd, prev_nd_len, p);
	prev_ln_prior += calcEdgeLenLnPrior(*ndP, prev_ndP_len, p);
	
	calculatePairwiseDistances();
	calculateProposalDist(true);
	ln_density_reverse_move = calcProposalLnDensity(propMeanX, prev_x_len);
	ln_density_reverse_move += calcProposalLnDensity(propMeanY, prev_y_len);
	ln_density_reverse_move += calcProposalLnDensity(propMeanW, prev_ndP_len);
	ln_density_reverse_move += calcProposalLnDensity(propMeanZ, prev_z_len);
	ln_density_reverse_move += calcProposalLnDensity(propMeanInternal, prev_nd_len);
	
	// code that followed the swap:
		
	calculateProposalDist(false);
	double xLen = proposeEdgeLen(propMeanX);
	double yLen = proposeEdgeLen(propMeanY);
	double zLen = proposeEdgeLen(propMeanZ);
	double ndPLen = proposeEdgeLen(propMeanW);
	double ndLen = proposeEdgeLen(propMeanInternal);
	
	//std::cerr << boost::str(boost::format("tree before [%.5f] = (x:%.5f,y:%.5f,(z:%.5f,w:%.5f):%.5f);\n") % prev_ln_like % prev_x_len % prev_y_len % prev_z_len % prev_ndP_len % prev_nd_len);
	
	x->SetEdgeLen(xLen);
	y->SetEdgeLen(yLen);
	z->SetEdgeLen(zLen);
	nd->SetEdgeLen(ndLen);
	ndP->SetEdgeLen(ndPLen);

	ln_density_forward_move = calcProposalLnDensity(propMeanX, xLen);
	ln_density_forward_move += calcProposalLnDensity(propMeanY, yLen);
	ln_density_forward_move += calcProposalLnDensity(propMeanW, ndPLen);
	ln_density_forward_move += calcProposalLnDensity(propMeanZ, zLen);
	ln_density_forward_move += calcProposalLnDensity(propMeanInternal, ndLen);

	//std::cerr << boost::str(boost::format("tree after = (z:%.5f,y:%.5f,(x:%.5f,w:%.5f):%.5f);\n") % zLen % yLen % xLen % ndPLen % ndLen);
        
	curr_ln_prior 	= calcEdgeLenLnPrior(*x, xLen, p)
					+ calcEdgeLenLnPrior(*y, yLen, p)
					+ calcEdgeLenLnPrior(*z, zLen, p)
					+ calcEdgeLenLnPrior(*nd, ndLen, p)
					+ calcEdgeLenLnPrior(*ndP, ndPLen, p);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapNNIMove::proposeNewState()
	{
    // Choose random internal node nd and randomly choose one of nd's children to call x (the other is y)
    // Here is the actual layout of the tree in memory
    //
    //      x  y        x = child that is swapped
    //       \/         y = x's sibling
    //   z   nd         z = nd's sibling
    //    \ /           w = nd's parent
    //    ndP
    //    /
    //   w
    //
    // Here is the layout that conforms to the terminology used
    //
    // x     z  <-- swapping x and z
    //  \___/ 
    //  /   \   Before move, x = ySis and z = wSis
    // y     w  After move,  x = wSis and z = ySis
    //
	scoringBeforeMove = true;
	nd  = randomInternalAboveSubroot();
	PHYCAS_ASSERT(nd);
	ndP = nd->GetParent();
	PHYCAS_ASSERT(ndP);
	wSis = z = nd->FindNextSib();
	PHYCAS_ASSERT(z);
	PHYCAS_ASSERT(nd->CountChildren() == 2); // we haven't figured this out for polytomies
	PHYCAS_ASSERT(ndP->CountChildren() == 2); // we haven't figured this out for polytomies

    x = nd->GetLeftChild();
	PHYCAS_ASSERT(x);
	y = x->GetRightSib();
	PHYCAS_ASSERT(y);

    // Decide which of the two children of nd to involve in the NNI move
	x_is_left = rng->Boolean();
	if (!x_is_left)
		std::swap(x,y);
    ySis = x;

	InternalData * nd_internal_data = nd->GetInternalData();
	PHYCAS_ASSERT(nd_internal_data);
	pre_root_posterior = nd_internal_data->getParentalCondLikePtr();
	pre_cla = nd_internal_data->getChildCondLikePtr();
	pre_p_mat = nd_internal_data->getMutablePMatrices();
	
	InternalData * ndp_internal_data = ndP->GetInternalData();
	PHYCAS_ASSERT(ndp_internal_data);
	post_root_posterior = ndp_internal_data->getParentalCondLikePtr();
	post_cla = ndp_internal_data->getChildCondLikePtr();
	post_p_mat = ndp_internal_data->getMutablePMatrices();


	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);

    //temporary!
    likelihood->storeAllCLAs(tree);

	prev_ln_like = FourTaxonLnLBeforeMove(nd);
	
	AlterBranchLengths(p);
	
	
	/* This swap is the equivalent of an NNI swap of the the nodes that are closest to y and w */
	scoringBeforeMove = false;
	std::swap(ySisTipData, wSisTipData);
	std::swap(ySis, wSis);


    curr_ln_like = FourTaxonLnLFromCorrectTipDataMembers(nd);
	//std::cerr << "curr_ln_like = " << curr_ln_like << '\n';
	}

double calcEdgeLenLnPrior(const TreeNode &x, double edge_len, ChainManagerShPtr & chain_mgr) 
	{
	if (x.IsTip())
		return chain_mgr->calcExternalEdgeLenPriorUnnorm(edge_len);
	return chain_mgr->calcInternalEdgeLenPriorUnnorm(edge_len);
	}


UnimapNNIMove::~UnimapNNIMove()
	{
	delete ySisTipData;
	delete yTipData;
	delete wSisTipData;
	delete wTipData;
	DeleteTwoDArray<double> (pre_w_pmat_transposed);
	DeleteTwoDArray<double> (pre_x_pmat_transposed);
	DeleteTwoDArray<double> (pre_y_pmat_transposed);
	DeleteTwoDArray<double> (pre_z_pmat_transposed);
	}

double UnimapNNIMove::FourTaxonLnLBeforeMove(TreeNode * nodePointer)
	{
	PHYCAS_ASSERT(nd);
	PHYCAS_ASSERT(nd == nodePointer);
	PHYCAS_ASSERT(ndP);
	TreeNode * lower = nd->FindNextSib();
	PHYCAS_ASSERT(lower);
	PHYCAS_ASSERT(nd->CountChildren() == 2); // we haven't figured this out for polytomies
	PHYCAS_ASSERT(ndP->CountChildren() == 2); // we haven't figured this out for polytomies
	TreeNode * upLeftNd = nd->GetLeftChild();
	PHYCAS_ASSERT(upLeftNd);
	TreeNode * upRightNd = upLeftNd->GetRightSib();
	PHYCAS_ASSERT(upRightNd);

	ySisTipData = createTipDataFromUnivents(getUniventsConstRef(*upLeftNd), ySisTipData);
	yTipData = createTipDataFromUnivents(getUniventsConstRef(*upRightNd), yTipData);
	if (!x_is_left)
		std::swap(ySisTipData, yTipData);
	wSisTipData = createTipDataFromUnivents(getUniventsConstRef(*lower), wSisTipData);		
	TreeNode * ndGp = ndP->GetParent();
	PHYCAS_ASSERT(ndGp);
	wTipData = createTipDataFromUnivents(getUniventsConstRef(*ndGp), wTipData);

	double lnlike = FourTaxonLnLFromCorrectTipDataMembers(nd);
	
	storePMatTransposed(pre_x_pmat_transposed, (const double ***) ySisTipData->getTransposedPMatrices());
	storePMatTransposed(pre_y_pmat_transposed, (const double ***) yTipData->getTransposedPMatrices());
	storePMatTransposed(pre_z_pmat_transposed, (const double ***) wSisTipData->getTransposedPMatrices());
	storePMatTransposed(pre_w_pmat_transposed, (const double ***) wTipData->getTransposedPMatrices());


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

    std::ofstream nxsf("tmp.nex", std::ios::app);
    //nxsf << "#nexus\n" << std::endl;
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
    	{
        nxsf << boost::str(boost::format("  utree curr = (x:%.8f, y:%.8f, (z:%.8f, w:%.8f):%.8f);") % ySis->GetEdgeLen() % y->GetEdgeLen() % wSis->GetEdgeLen() % nd->GetParent()->GetEdgeLen() % nd->GetEdgeLen()) << std::endl;
        nxsf << "  [y->blen    = " << y->GetEdgeLen() << '\n';
        nxsf << "   ySis->blen = " << ySis->GetEdgeLen() << '\n';
        nxsf << "   nd->GetParent()->blen    = " << nd->GetParent()->GetEdgeLen() << '\n';
        nxsf << "   wSis->blen = " << wSis->GetEdgeLen() << '\n';
        nxsf << "   ]\n";
        }
    else
    	{
    	PHYCAS_ASSERT(false);
        nxsf << boost::str(boost::format("  utree curr = (z:%.8f, y:%.8f, (x:%.8f, w:%.8f):%.8f);") % z->GetEdgeLen() % y->GetEdgeLen() % x->GetEdgeLen() % xpar->GetEdgeLen() % zpar->GetEdgeLen()) << std::endl;
        }
    nxsf << "end;" << std::endl;

    nxsf << "\nbegin paup;" << std::endl;
    nxsf << boost::str(boost::format("  [!***** phycas lnL = %.8f *****]") % lnlike) << std::endl;
    if (model->getModelName().compare("JC69") == 0)
        nxsf << "  lset nst=1 basefreq=equal;" << std::endl;
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
    }

void UnimapNNIMove::storePMatTransposed(double **& cached, const double *** p_mat_array)
	{
	const unsigned nStates = likelihood->getNStates();
	if (!cached)
		cached = NewTwoDArray<double>(nStates + 1, nStates + 1);
	for (unsigned i = 0; i < nStates; ++i)
		{
		for (unsigned j = 0; j < nStates; ++j)
			cached[i][j] = p_mat_array[0][i][j];
		}
	}
	
double UnimapNNIMove::FourTaxonLnLFromCorrectTipDataMembers(TreeNode * nodePointer)
	{
	PHYCAS_ASSERT(nd == nodePointer);
	CondLikelihoodShPtr nd_childCLPtr, nd_parentCLPtr;
	double *** p_mat;
	if (scoringBeforeMove)
		{
		nd_childCLPtr = pre_root_posterior;
		nd_parentCLPtr = pre_cla;
		p_mat = pre_p_mat;
		}
	else
		{
		nd_childCLPtr = post_root_posterior;
		nd_parentCLPtr = post_cla;
		p_mat = post_p_mat;
		}

    likelihood->calcPMatTranspose(ySisTipData->getTransposedPMatrices(), ySisTipData->getConstStateListPos(), ySis->GetEdgeLen());
    likelihood->calcPMatTranspose(yTipData->getTransposedPMatrices(), yTipData->getConstStateListPos(), y->GetEdgeLen());
    likelihood->calcPMatTranspose(wSisTipData->getTransposedPMatrices(), wSisTipData->getConstStateListPos(), wSis->GetEdgeLen());
    likelihood->calcPMatTranspose(wTipData->getTransposedPMatrices(), wTipData->getConstStateListPos(), nd->GetParent()->GetEdgeLen());

	
	likelihood->calcPMat(p_mat, nd->GetEdgeLen());
	const double * const *  childPMatrix = p_mat[0];
	likelihood->calcCLATwoTips(*nd_childCLPtr, *ySisTipData, *yTipData);
	likelihood->calcCLATwoTips(*nd_parentCLPtr, *wSisTipData, *wTipData);

	double lnl =  HarvestLnLikeFromCondLikePar(nd_childCLPtr, nd_parentCLPtr, childPMatrix);
	DebugSaveNexusFile(ySisTipData, yTipData, wSisTipData, wTipData, lnl);
	return lnl;
	}
	
double UnimapNNIMove::HarvestLnLikeFromCondLikePar(
  CondLikelihoodShPtr focalCondLike, 
  ConstCondLikelihoodShPtr neighborCondLike, 
  const double * const * childPMatrix)
	{
	LikeFltType * focalNdCLAPtr = focalCondLike->getCLA(); //PELIGROSO
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
			focalNdCLAPtr[i] *= stateFreq[i]*neigborLike;
			siteLike += focalNdCLAPtr[i];
			}
		focalNdCLAPtr += num_states;
		focalNeighborCLAPtr += num_states;
		double site_lnL = std::log(siteLike);

		lnLikelihood += site_lnL;
		}
	return lnLikelihood;
	}

	
TipData * UnimapNNIMove::createTipDataFromUnivents(const Univents & u, TipData *td)
	{
	if (td)
		{
		/* this is the one place in which we overwrite the state codes */
		int8_t * stateCodes = const_cast<int8_t *>(td->getTipStatesArray().get());
		u.fillStateCodeArray(stateCodes);
		}
	else
		{
		const unsigned num_patterns = likelihood->getNPatterns();
		int8_t * tipSpecificStateCode = new int8_t[num_patterns];
		u.fillStateCodeArray(tipSpecificStateCode);
		std::vector<unsigned int> emptyStateListVec;
		td = new TipData(	true,
						num_patterns,
						emptyStateListVec,												// stateListPosVec
						boost::shared_array<const int8_t>(tipSpecificStateCode),	// stateCodesShPtr
						1,													// number of relative rate categories
						likelihood->getNStates(),													// number of states in the model
						NULL,
						true,														// managePMatrices
						likelihood->getCondLikelihoodStorage());
		}
	return td;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapNNIMove::revert()
	{
	//std::cerr << "REVERTED" << std::endl;
	
	x->SetEdgeLen(prev_x_len);
	y->SetEdgeLen(prev_y_len);
	z->SetEdgeLen(prev_z_len);
	TreeNode * nd =  y->GetParent();
	assert(nd);
	assert(nd->GetParent());
	nd->SetEdgeLen(prev_nd_len);
	TreeNode * ndP = nd->GetParent();
	ndP->SetEdgeLen(prev_ndP_len);
	

	const UniventProbMgr & upm = likelihood->GetUniventProbMgrConstRef();
	Lot & rngRef = *rng;
	/* using the statecode arrays at  ySisTipData  and yTipData
		as storage for the new samples of nd and ndP sequences.
	*/
	
	LikeFltType * root_state_posterior = pre_root_posterior->getCLA(); //PELIGROSO
	Univents & ndU = getUniventsRef(*nd);
	upm.sampleRootStates(ndU, root_state_posterior, rngRef, false, NULL);

	Univents & ndPU = getUniventsRef(*ndP);
	const int8_t * nd_states = &(ndU.getEndStatesVecConstRef()[0]);
	const LikeFltType * des_cla = pre_cla->getCLA(); //PELIGROSO
	upm.sampleDescendantStates(ndPU, (const double **) pre_p_mat[0], des_cla, nd_states, rngRef);

	if (doSampleUnivents)
		{
		const int8_t * ndP_states = &(ndPU.getEndStatesVecConstRef()[0]);
		upm.sampleUniventsKeepEndStates(getUniventsRef(*x), prev_x_len, nd_states, (const double **) pre_x_pmat_transposed, rngRef);
		upm.sampleUniventsKeepEndStates(getUniventsRef(*y), prev_y_len, nd_states, (const double **) pre_y_pmat_transposed, rngRef);
		upm.sampleUniventsKeepEndStates(getUniventsRef(*z), prev_z_len, ndP_states, (const double **) pre_z_pmat_transposed, rngRef);
		upm.sampleUnivents(ndU, prev_nd_len, ndP_states, const_cast<const double *const *>(pre_p_mat[0]), rngRef, NULL);
		
		const TreeNode * ndGP = ndP->GetParent();
		PHYCAS_ASSERT(ndGP);
		const Univents & ndGPU = getUniventsConstRef(*ndGP);
		const int8_t * ndGP_states = &(ndGPU.getEndStatesVecConstRef()[0]);

		/*	We want the w to ndP branch's pmat, but we don't want it transposed
			since we are done with pre_z_pmat_transposed for this move, we can 
			clobber the stored values 
		*/
		fillTranspose(pre_z_pmat_transposed, const_cast<const double *const *>(pre_w_pmat_transposed), likelihood->getNStates());
		upm.sampleUnivents(ndPU, prev_ndP_len, ndGP_states, const_cast<const double *const *>(pre_z_pmat_transposed), rngRef, NULL);

		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void UnimapNNIMove::accept()
	{
	//std::cerr << "ACCEPTED" << std::endl;
	
	/*x->SelectNode();
    y->SelectNode();
    z->SelectNode();
    nd->SelectNode();

	likelihood->startTreeViewer(tree, "before");
	*/
	tree_manipulator.NNISwap(z, x);
    /*
    likelihood->startTreeViewer(tree, "after");
    
    x->UnselectNode();
    y->UnselectNode();
    z->UnselectNode();
    nd->UnselectNode();
    */
    
   // post_child_cla post_parent_cla
	/* using the statecode arrays at  ySisTipData  and yTipData
		as storage for the new samples of nd and ndP sequences.
	*/
	TreeNode * nd =  y->GetParent();
	assert(nd);
	assert(nd->GetParent());
	TreeNode * ndP = nd->GetParent();

	const UniventProbMgr & upm = likelihood->GetUniventProbMgrConstRef();
	Lot & rngRef = *rng;

	LikeFltType * root_state_posterior = post_root_posterior->getCLA(); //PELIGROSO
	Univents & ndU = getUniventsRef(*nd);
	upm.sampleRootStates(ndU, root_state_posterior, rngRef, false, NULL);

	Univents & ndPU = getUniventsRef(*ndP);
	const int8_t * nd_states = &(ndU.getEndStatesVecConstRef()[0]);
	const LikeFltType * des_cla = post_cla->getCLA(); //PELIGROSO
	upm.sampleDescendantStates(ndPU, (const double **) post_p_mat[0], des_cla, nd_states, rngRef);
	if (doSampleUnivents)
		{
		const int8_t * ndP_states = &(ndPU.getEndStatesVecConstRef()[0]);
		upm.sampleUniventsKeepEndStates(getUniventsRef(*x), x->GetEdgeLen(), ndP_states, (const double **) wSisTipData->getTransposedPMatrices()[0], rngRef);
		upm.sampleUniventsKeepEndStates(getUniventsRef(*y), x->GetEdgeLen(), nd_states, (const double **) yTipData->getTransposedPMatrices()[0], rngRef);
		upm.sampleUniventsKeepEndStates(getUniventsRef(*z), x->GetEdgeLen(), nd_states, (const double **) ySisTipData->getTransposedPMatrices()[0], rngRef);
		upm.sampleUnivents(ndU, nd->GetEdgeLen(), ndP_states, (const double **) post_p_mat[0], rngRef, NULL);
		
		const TreeNode * ndGP = ndP->GetParent();
		PHYCAS_ASSERT(ndGP);
		const Univents & ndGPU = getUniventsConstRef(*ndGP);
		const int8_t * ndGP_states = &(ndGPU.getEndStatesVecConstRef()[0]);
		/*	We want the w to ndP branch's pmat, but we don't want it transposed
			since we are done with pre_z_pmat_transposed for this move, we can 
			clobber the stored values 
		*/
		fillTranspose(pre_z_pmat_transposed, const_cast<const double *const *>(wTipData->getTransposedPMatrices()[0]), likelihood->getNStates());
		upm.sampleUnivents(ndPU, ndP->GetEdgeLen(), ndGP_states, pre_z_pmat_transposed, rngRef, NULL);
		}
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
  ySisTipData(0),
  yTipData(0),
  wSisTipData(0),
  wTipData(0),
  doSampleUnivents(false)
	{
	min_edge_len_mean = 0.02;
	edge_len_prop_cv = 1;
	pre_x_pmat_transposed = 0L;
	pre_y_pmat_transposed = 0L;
	pre_w_pmat_transposed = 0L;
	pre_z_pmat_transposed = 0L;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move.
*/
double UnimapNNIMove::getLnHastingsRatio() const
	{
	return ln_density_reverse_move - ln_density_forward_move;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Jacobian for this move.
*/
double UnimapNNIMove::getLnJacobian() const
	{
	return 0.0;
	}

}	// namespace phycas
