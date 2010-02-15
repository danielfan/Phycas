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
#include "phycas/src/gtr_model.hpp"

#include "boost/format.hpp"

namespace phycas
{
bool modify_terminal_edges = false;
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapTopoMove::setLot(LotShPtr p)
	{
	MCMCUpdater::setLot(p);
	}

void UnimapNNIMove::setLot(LotShPtr p)
	{
	UnimapTopoMove::setLot(p);
	gammaDist.SetLot(rng.get());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool UnimapTopoMove::update()
	{
#if POLPY_NEWWAY
	PHYCAS_ASSERT(0);	// need to take account of working_prior
#endif
	// The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed UnimapNNIMove cannot be accepted without changing edge lengths, so it is best to bail out now.
	if (is_fixed)
		return false;

	//std::cerr << "****** UnimapNNIMove::update" << std::endl;

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
TreeNode * UnimapTopoMove::randomInternalAboveSubroot()
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
#if DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
	dXY = dWX =  dXZ =  dWY =  dYZ =  dWZ = 0.0;
	/* This is called before the swap so "x" is "b" and "z" is "d" */
	const int8_t * xStates = bTipData->getTipStatesArray().get();
	const int8_t * yStates = aTipData->getTipStatesArray().get();
	const int8_t * wStates = cTipData->getTipStatesArray().get();
	const int8_t * zStates = dTipData->getTipStatesArray().get();
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
#endif
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
	return gammaDist.GetLnPDF(x);
	}

double UnimapNNIMove::proposeEdgeLen(double mean)
	{
	const double variance = (edge_len_prop_cv*edge_len_prop_cv*mean*mean);
	gammaDist.SetMeanAndVariance(mean, variance);
	return gammaDist.Sample();
	}



/*******************************************************************************
*	Alter branch lengths is called in the context:

##########################################################
	prev_ln_like = FourTaxonLnLBeforeMove();
	ProposeStateWithTemporaries(p);
    curr_ln_like = FourTaxonLnLFromCorrectTipDataMembers();
##########################################################

	so it should:
		1. calculate the temporaries needed to calc the Hastings ratio,
		2. calculate the curr_ln_prior
		3. set the branches to the proposed branch lengths,
		4. swap ySis with wSis and ySisTipData with wSisTipData if there is an NNI move
*/
void UnimapNNIMove::ProposeStateWithTemporaries(ChainManagerShPtr & p)
	{
	
	calculatePairwiseDistances();
	calculateProposalDist(true);
	ln_density_reverse_move = calcProposalLnDensity(propMeanX, prev_x_len);
	ln_density_reverse_move += calcProposalLnDensity(propMeanY, prev_y_len);
	ln_density_reverse_move += calcProposalLnDensity(propMeanW, prev_ndP_len);
	ln_density_reverse_move += calcProposalLnDensity(propMeanZ, prev_z_len);
	ln_density_reverse_move += calcProposalLnDensity(propMeanInternal, prev_nd_len);
	
		
	calculateProposalDist(false);
	double xLen, yLen, zLen, ndLen, ndPLen;
	if (modify_terminal_edges)
		{
		xLen= proposeEdgeLen(propMeanX);
		yLen = proposeEdgeLen(propMeanY);
		zLen = proposeEdgeLen(propMeanZ);
		ndPLen = proposeEdgeLen(propMeanW);
		ndLen = proposeEdgeLen(propMeanInternal);
		}
	else
		{
		xLen = prev_x_len;
		yLen = prev_y_len;
		zLen = prev_z_len;
		ndPLen = prev_ndP_len;
		ndLen = proposeEdgeLen(0.1);
		}
	
	//std::cerr << boost::str(boost::format("tree before [%.5f] = (x:%.5f,y:%.5f,(z:%.5f,w:%.5f):%.5f);\n") % prev_ln_like % prev_x_len % prev_y_len % prev_z_len % prev_ndP_len % prev_nd_len);
	
	x->SetEdgeLen(xLen);
	y->SetEdgeLen(yLen);
	z->SetEdgeLen(zLen);
	origNode->SetEdgeLen(ndLen);
	origNodePar->SetEdgeLen(ndPLen);

	ln_density_forward_move = calcProposalLnDensity(propMeanX, xLen);
	ln_density_forward_move += calcProposalLnDensity(propMeanY, yLen);
	ln_density_forward_move += calcProposalLnDensity(propMeanW, ndPLen);
	ln_density_forward_move += calcProposalLnDensity(propMeanZ, zLen);
	ln_density_forward_move += calcProposalLnDensity(propMeanInternal, ndLen);

	//std::cerr << boost::str(boost::format("tree after = (z:%.5f,y:%.5f,(x:%.5f,w:%.5f):%.5f);\n") % zLen % yLen % xLen % ndPLen % ndLen);
        
	curr_ln_prior 	= calcEdgeLenLnPrior(*x, xLen, p)
					+ calcEdgeLenLnPrior(*y, yLen, p)
					+ calcEdgeLenLnPrior(*z, zLen, p)
					+ calcEdgeLenLnPrior(*origNode, ndLen, p)
					+ calcEdgeLenLnPrior(*origNodePar, ndPLen, p);

	
	/* This swap is the equivalent of an NNI swap of the the nodes that are closest to y and w */
	std::swap(bTipData, dTipData);
	std::swap(b, d);
	std::swap(bLenNd, dLenNd);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapTopoMove::proposeNewState()
	{
#if DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
    // Choose random internal node origNode and randomly choose one of origNode's children to call x (the other is y)
    // Here is the actual layout of the tree in memory
    //
    //      x  y        x = child that is swapped
    //       \/         y = x's sibling
    //   z   origNode	z = origNode's sibling
    //    \ /           w = origNodePar's parent
    //    origNodePar
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
    
    // When calculating likelihoods we refer to the tip data as 
    // aTipData, bTipData, cTipData, and dTipData where the tree is ((a,b), (c,d))
    // with the internal length from origNd

	origNode = randomInternalAboveSubroot();
	PHYCAS_ASSERT(origNode);
	PHYCAS_ASSERT(origNode->CountChildren() == 2); // we haven't figured this out for polytomies
	origNodePar = origNode->GetParent();
	PHYCAS_ASSERT(origNodePar);
	PHYCAS_ASSERT(origNodePar->CountChildren() == 2); // we haven't figured this out for polytomies
	w = origNodePar->GetParent();
	PHYCAS_ASSERT(w);
	z = origNode->FindNextSib();
	PHYCAS_ASSERT(z);
    x = origNode->GetLeftChild();
	PHYCAS_ASSERT(x);
	y = x->GetRightSib();
	PHYCAS_ASSERT(y);

    // Decide which of the two children of origNode to involve in the NNI move
	x_is_left = rng->Boolean();
	if (!x_is_left)
		std::swap(x,y);

	a = y;
	b = x;
	c = w;
	d = z;
	aLenNd = a;
	bLenNd = b;
	cLenNd = origNodePar;
	dLenNd = d;
	
	InternalData * nd_internal_data = origNode->GetInternalData();
	PHYCAS_ASSERT(nd_internal_data);
	pre_root_posterior = nd_internal_data->getParentalCondLikePtr();
	pre_cla = nd_internal_data->getChildCondLikePtr();
	pre_p_mat = nd_internal_data->getMutablePMatrices();
	
	InternalData * ndp_internal_data = origNodePar->GetInternalData();
	PHYCAS_ASSERT(ndp_internal_data);
	post_root_posterior = ndp_internal_data->getParentalCondLikePtr();
	post_cla = ndp_internal_data->getChildCondLikePtr();
	post_p_mat = ndp_internal_data->getMutablePMatrices();


	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);

    //temporary!
    scoringBeforeMove = true;
	likelihood->storeAllCLAs(tree);

	prev_ln_like = FourTaxonLnLBeforeMove();


	prev_x_len = x->GetEdgeLen();
	prev_y_len = y->GetEdgeLen();
	prev_z_len = z->GetEdgeLen();
	prev_nd_len = origNode->GetEdgeLen();
	prev_ndP_len = origNodePar->GetEdgeLen();
	
	prev_ln_prior = calcEdgeLenLnPrior(*x, prev_x_len, p);
	prev_ln_prior += calcEdgeLenLnPrior(*y, prev_y_len, p);
	prev_ln_prior += calcEdgeLenLnPrior(*z, prev_z_len, p);
	prev_ln_prior += calcEdgeLenLnPrior(*origNode, prev_nd_len, p);
	prev_ln_prior += calcEdgeLenLnPrior(*origNodePar, prev_ndP_len, p);
	
	ProposeStateWithTemporaries(p);

	scoringBeforeMove = false;
    curr_ln_like = FourTaxonLnLFromCorrectTipDataMembers();
#endif
	}

double calcEdgeLenLnPrior(const TreeNode &x, double edge_len, ChainManagerShPtr & chain_mgr) 
	{
	if (x.IsTip())
		return chain_mgr->calcExternalEdgeLenPriorUnnorm(edge_len);
	return chain_mgr->calcInternalEdgeLenPriorUnnorm(edge_len);
	}


UnimapTopoMove::~UnimapTopoMove()
	{
	delete aTipData;
	delete bTipData;
	delete cTipData;
	delete dTipData;
	DeleteTwoDArray<double> (pre_w_pmat_transposed);
	DeleteTwoDArray<double> (pre_x_pmat_transposed);
	DeleteTwoDArray<double> (pre_y_pmat_transposed);
	DeleteTwoDArray<double> (pre_z_pmat_transposed);
	}

double UnimapTopoMove::FourTaxonLnLBeforeMove()
	{
#if DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING	
	aTipData = createTipDataFromUnivents(getUniventsConstRef(*a), aTipData);
	bTipData = createTipDataFromUnivents(getUniventsConstRef(*b), bTipData);
	cTipData = createTipDataFromUnivents(getUniventsConstRef(*c), cTipData);
	dTipData = createTipDataFromUnivents(getUniventsConstRef(*d), dTipData);		

	double lnlike = FourTaxonLnLFromCorrectTipDataMembers();
	
	storePMatTransposed(pre_y_pmat_transposed, (const double ***) aTipData->getTransposedPMatrices());
	storePMatTransposed(pre_x_pmat_transposed, (const double ***) bTipData->getTransposedPMatrices());
	storePMatTransposed(pre_w_pmat_transposed, (const double ***) cTipData->getTransposedPMatrices());
	storePMatTransposed(pre_z_pmat_transposed, (const double ***) dTipData->getTransposedPMatrices());


    return lnlike;
#else
	return 0.0;
#endif
	}

bool UnimapTopoMove::CheckWithPaup(double lnlike)
{
	std::ofstream nxsf;
	nxsf.open("debug-with-paup.nex");
	nxsf << "#NEXUS\n";
	DebugSaveNexusFile(nxsf, lnlike);
	nxsf.close();
	nxsf << "begin PAUP;\n quit ;\nend;\n";
	return (system("./checkWithPAUP.py") == 0);
}

void UnimapTopoMove::DebugSaveNexusFile(std::ostream & nxsf, double lnlike)
    {
#if DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
    typedef boost::shared_array<const int8_t> StateArr;
    StateArr adata = aTipData->getTipStatesArray();
    StateArr bdata = bTipData->getTipStatesArray();
    StateArr cdata = cTipData->getTipStatesArray();
    StateArr ddata = dTipData->getTipStatesArray();
    unsigned nchar = likelihood->getNPatterns();
    unsigned i;
	const char * alphabet = "ACGT";
    //nxsf << "#nexus\n" << std::endl;
    nxsf << "begin data;" << std::endl;
    nxsf << "  dimensions ntax=4 nchar=" << nchar << ";" << std::endl;
    nxsf << "  format datatype=dna gap=- missing=?;" << std::endl;
    nxsf << "  matrix" << std::endl;
    nxsf << "    a ";
    for (i = 0; i < nchar; ++i)
        nxsf << alphabet[adata[i]];
    nxsf << "\n    b ";
    for (i = 0; i < nchar; ++i)
        nxsf << alphabet[bdata[i]];
    nxsf << "\n    c ";
    for (i = 0; i < nchar; ++i)
        nxsf << alphabet[cdata[i]];
    nxsf << "\n    d ";
    for (i = 0; i < nchar; ++i)
        nxsf << alphabet[ddata[i]];
    nxsf << "\n  ;\nend;" << std::endl;

    nxsf << "\nbegin paup;" << std::endl;
    nxsf << "  log file=tmp.log start replace;" << std::endl;
    nxsf << "  set criterion=likelihood storebrlens;" << std::endl;
    nxsf << "end;" << std::endl;

    nxsf << "\nbegin trees;" << std::endl;
	nxsf << boost::str(boost::format("  utree curr = (a:%.8f, b:%.8f, (c:%.8f, d:%.8f):%.8f);") % aLenNd->GetEdgeLen() % bLenNd->GetEdgeLen() % cLenNd->GetEdgeLen() % dLenNd->GetEdgeLen() % origNode->GetEdgeLen()) << std::endl;
	nxsf << " \nend;\n\nbegin paup;\n";
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
#endif
    }

void UnimapTopoMove::storePMatTransposed(double **& cached, const double *** p_mat_array)
	{
#if DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
	const unsigned nStates = likelihood->getNStates();
	if (!cached)
		cached = NewTwoDArray<double>(nStates + 1, nStates + 1);
	for (unsigned i = 0; i < nStates; ++i)
		{
		for (unsigned j = 0; j < nStates; ++j)
			cached[i][j] = p_mat_array[0][i][j];
		}
#endif
	}
	
double UnimapTopoMove::FourTaxonLnLFromCorrectTipDataMembers()
	{
#if DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
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

    likelihood->calcPMatTranspose(aTipData->getTransposedPMatrices(), aTipData->getConstStateListPos(), aLenNd->GetEdgeLen());
    likelihood->calcPMatTranspose(bTipData->getTransposedPMatrices(), bTipData->getConstStateListPos(), bLenNd->GetEdgeLen());
    likelihood->calcPMatTranspose(cTipData->getTransposedPMatrices(), cTipData->getConstStateListPos(), cLenNd->GetEdgeLen());
    likelihood->calcPMatTranspose(dTipData->getTransposedPMatrices(), dTipData->getConstStateListPos(), dLenNd->GetEdgeLen());

	
	likelihood->calcPMat(p_mat, origNode->GetEdgeLen());
	const double * const *  childPMatrix = p_mat[0];
	likelihood->calcCLATwoTips(*nd_childCLPtr, *bTipData, *aTipData);
	likelihood->calcCLATwoTips(*nd_parentCLPtr, *dTipData, *cTipData);

	double lnl =  HarvestLnLikeFromCondLikePar(nd_childCLPtr, nd_parentCLPtr, childPMatrix);
#	if defined CHECK_EACH_CALC_AGAINST_PAUP
		std::cerr << "Scoring " << (scoringBeforeMove ? "before" : "after") << " the proposal\n";
		PHYCAS_ASSERT(CheckWithPaup(lnl));
#	endif
	return lnl;
#else // old way
	return 0.0;
#endif
	}
	
double UnimapTopoMove::HarvestLnLikeFromCondLikePar(
  CondLikelihoodShPtr focalCondLike, 
  ConstCondLikelihoodShPtr neighborCondLike, 
  const double * const * childPMatrix)
	{
#if DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
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
#else
	return 0.0;
#endif
	}

	
TipData * UnimapTopoMove::createTipDataFromUnivents(const Univents & u, TipData *td)
	{
#if DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
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
#else
	return NULL;
#endif
	}


/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapTopoMove::revert()
	{
	MCMCUpdater::revert();
//	std::cerr << "REVERTED" << std::endl;
	
	x->SetEdgeLen(prev_x_len);
	y->SetEdgeLen(prev_y_len);
	z->SetEdgeLen(prev_z_len);
	PHYCAS_ASSERT(origNode ==  y->GetParent());
	PHYCAS_ASSERT(origNode);
	PHYCAS_ASSERT(origNodePar == origNode->GetParent());
	origNode->SetEdgeLen(prev_nd_len);
	origNodePar->SetEdgeLen(prev_ndP_len);
	
	if (doSampleInternalStates)
		resampleInternalNodeStates(pre_root_posterior->getCLA(), pre_cla->getCLA());
	}

void UnimapTopoMove::resampleInternalNodeStates(const LikeFltType * root_state_posterior, const LikeFltType * des_cla)
{
	const UniventProbMgr & upm = likelihood->GetUniventProbMgrConstRef();
	Lot & rngRef = *rng;
	
	TreeNode * aPar = a->GetParent();
	PHYCAS_ASSERT(aPar && (aPar == origNode || aPar == origNodePar));
	Univents & ndU = getUniventsRef(*aPar);
	upm.sampleRootStates(ndU, root_state_posterior, rngRef, false, NULL);
	TreeNode * otherInternal = (aPar == origNode ? origNodePar : origNode);
	PHYCAS_ASSERT(otherInternal);
	Univents & ndPU = getUniventsRef(*otherInternal);
	const int8_t * nd_states = &(ndU.getEndStatesVecConstRef()[0]);
	upm.sampleDescendantStates(ndPU, (const double **) pre_p_mat[0], des_cla, nd_states, rngRef);
	
	likelihood->flagNodeWithInvalidUnivents(x);
	getUniventsRef(*x).setValid(false);
	likelihood->flagNodeWithInvalidUnivents(y);
	getUniventsRef(*y).setValid(false);
	likelihood->flagNodeWithInvalidUnivents(z);
	getUniventsRef(*z).setValid(false);
	likelihood->flagNodeWithInvalidUnivents(origNode);
	getUniventsRef(*origNode).setValid(false);
	likelihood->flagNodeWithInvalidUnivents(origNodePar);
	getUniventsRef(*origNodePar).setValid(false);
}
/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void UnimapNNIMove::accept()
	{
	MCMCUpdater::accept();
	//std::cerr << "ACCEPTED" << std::endl;
	tree_manipulator.NNISwap(z, x);
    
	PHYCAS_ASSERT(origNode);
	PHYCAS_ASSERT(origNode->GetParent() == origNodePar);

	if (doSampleInternalStates)
		resampleInternalNodeStates(post_root_posterior->getCLA(), post_cla->getCLA());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
UnimapTopoMove::UnimapTopoMove() : MCMCUpdater(),
  x(0),
  y(0),
  z(0),
  w(0),
  a(0),
  b(0),
  c(0),
  d(0),
  aTipData(0),
  bTipData(0),
  cTipData(0),
  dTipData(0),
  pre_x_pmat_transposed(0L),
  pre_y_pmat_transposed(0L),
  pre_w_pmat_transposed(0L),
  pre_z_pmat_transposed(0L),
  doSampleInternalStates(true)
  	{
	}
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
UnimapNNIMove::UnimapNNIMove() : UnimapTopoMove()
	{
	min_edge_len_mean = 0.02;
	edge_len_prop_cv = 1;
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
