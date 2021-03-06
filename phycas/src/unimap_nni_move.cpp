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
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>


#include "phycas/src/unimap_nni_move.hpp"


#include "phycas/src/probability_distribution.hpp"
#include "phycas/src/likelihood_models.hpp"
#include "phycas/src/basic_tree_node.hpp"
#include "phycas/src/tree_likelihood.hpp"
#include "phycas/src/xlikelihood.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/internal_data.hpp"
#include "phycas/src/tip_data.hpp"
#include "phycas/src/gtr_model.hpp"
#include "phycas/src/cond_likelihood_storage.hpp"
#include "boost/format.hpp"

namespace phycas
{
bool modify_terminal_edges = true;
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapTopoMove::setLot(LotShPtr p)
	{
	MCMCUpdater::setLot(p);
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool UnimapTopoMove::update()
	{
	// The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed UnimapNNIMove cannot be accepted without changing edge lengths, so it is best to bail out now.
	if (is_fixed)
		return false;

	//std::cerr << "****** UnimapNNIMove::update" << std::endl;

	if (!ndQ.empty() && ndQ.front() == 0L) // a 0L in the ndQ is a signal to do nothing.
		{
		ndQ.pop();
		PHYCAS_ASSERT(ndQ.empty());
		return false;
		}
	bool accepted = false;
	
#if 1 ||  DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
	//tree->renumberInternalNodes(tree->GetNTips()); //@POL this should be somewhere else

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
	const double ln_hastings =  getLnHastingsRatio();
	double ln_accept_ratio = curr_posterior - prev_posterior + ln_hastings;
	//double lnu = std::log(rng->Uniform(FILE_AND_LINE));
	//bool accepted = (ln_accept_ratio >= 0.0 || lnu <= ln_accept_ratio);
	//bool accepted = (ln_accept_ratio >= 0.0 || std::log(rng->Uniform(FILE_AND_LINE)) <= ln_accept_ratio);
	double lnu = DBL_MAX;
	accepted = (ln_accept_ratio >= 0.0);
	if (!accepted)
		{
		double u = rng->Uniform(FILE_AND_LINE);
		lnu = std::log(u);
		accepted = (lnu <= ln_accept_ratio);
		}
	//std::cerr << "lnu    ln_accept_ratio = curr_posterior - prev_posterior + ln_density_reverse_move - ln_density_forward_move;\n";
	//std::cerr << lnu << "   " << ln_accept_ratio << "   " << curr_posterior << "   " <<prev_posterior << "   " << ln_density_reverse_move << "   " << ln_density_forward_move << "\n\n";

    if (save_debug_info)
        {
		//std::cerr << "ln_accept_ratio = curr_posterior - prev_posterior + ln_hastings \n";
		//std::cerr << ln_accept_ratio << "   " << curr_posterior << "   " <<prev_posterior << ' '<< ln_hastings << "\n\n";
		
		debug_info = str(boost::format("swapping %d <-> %d (%s, lnR = %.5f)") % x->GetNodeNumber() % z->GetNodeNumber() % (accepted ? "accepted" : "rejected") % ln_accept_ratio);
    	}
    if (accepted)
		accept();
	else
		revert();

#	if FOUR_TAXON_TEST
		if (isFirstTime)
			isFirstTime = false;
		else
			PHYCAS_ASSERT(fabs(cached_posterior - prev_posterior) < 0.0001);

		if (accepted)
			cached_posterior = curr_posterior;
		else 
			cached_posterior = prev_posterior;
#	endif
#endif

	return accepted;
	}

TreeNode * randomInternalAboveSubroot(Tree & tree, Lot &rng)
    {
	// Avoiding the "subroot" node (only child of the tip serving as the root), so the number of 
	// acceptable nodes is one fewer than the number of internal nodes
	unsigned numAcceptableNodes = tree.GetNInternals() - 1;

	unsigned ypos = rng.SampleUInt(numAcceptableNodes);
	unsigned i = 0;
    TreeNode * nd = tree.GetFirstPreorder();
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


void  UnimapTopoMove::releaseCondLikeArrays()
	{
	if (cla_pool)
		{
		if (pre_root_posterior)
			cla_pool->putCondLikelihood(pre_root_posterior);
		pre_root_posterior = CondLikelihoodShPtr();
		if (pre_cla)
			cla_pool->putCondLikelihood(pre_cla);
		pre_cla = CondLikelihoodShPtr();
		if (post_root_posterior)
			cla_pool->putCondLikelihood(post_root_posterior);
		post_root_posterior = CondLikelihoodShPtr();
		if (post_cla)
			cla_pool->putCondLikelihood(post_cla);
		pre_root_posterior = CondLikelihoodShPtr();
		}
	}
void  UnimapTopoMove::acquireCondLikeArrays(InternalData & id)
	{
	releaseCondLikeArrays();
	cla_pool = id.getCondLikeStoragePool();
	pre_root_posterior = cla_pool->getCondLikelihood();
	pre_cla = cla_pool->getCondLikelihood();
	post_root_posterior = cla_pool->getCondLikelihood();
	post_cla = cla_pool->getCondLikelihood();
	}
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapTopoMove::proposeNewState()
	{
#if 1 || DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
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
	//std::cerr << "UnimapTopoMove::proposeNewState " << this->getName() << '\n';
	if (ndQ.empty())
		{
		origNode = randomInternalAboveSubroot(*tree, *rng);
		PHYCAS_ASSERT(false); //temporary
		}
	else
		{
		origNode = ndQ.front();
		ndQ.pop();
		PHYCAS_ASSERT(ndQ.empty());
		}
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
	InternalData * ndp_internal_data = origNodePar->GetInternalData();
	PHYCAS_ASSERT(ndp_internal_data);
	const unsigned numSubsets = likelihood->getNumSubsets();

#	if defined (BORROW_COND_LIKE_PTR)
		pre_root_posterior = nd_internal_data->getParentalCondLikePtr();
		PHYCAS_ASSERT(pre_root_posterior);
		pre_cla = nd_internal_data->getChildCondLikePtr();
		PHYCAS_ASSERT(pre_cla);
		post_root_posterior = ndp_internal_data->getParentalCondLikePtr();
		PHYCAS_ASSERT(post_root_posterior);
		post_cla = ndp_internal_data->getChildCondLikePtr();
		PHYCAS_ASSERT(post_cla);
#	else
		if (!pre_root_posterior)
			{
			acquireCondLikeArrays(*nd_internal_data);
			}
#	endif

	pre_p_mat.resize(numSubsets);
	post_p_mat.resize(numSubsets);
	
	for (unsigned i = 0; i < numSubsets; ++i)
		{
		pre_p_mat[i] = nd_internal_data->getMutablePMatrices(i);
		post_p_mat[i] = ndp_internal_data->getMutablePMatrices(i);
		}

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);

    //temporary!
    scoringBeforeMove = true;
	//likelihood->storeAllCLAs(tree);

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
	
	//origNode->SelectNode();
	//likelihood->startTreeViewer(tree, "Just after ProposeStateWithTemporaries called in UnimapTopoMove::proposeNewState");
	//origNode->UnselectNode();

	scoringBeforeMove = false;
    curr_ln_like = 0.0;
	for (unsigned i = 0; i < numSubsets; ++i)
	    curr_ln_like += FourTaxonLnLFromCorrectTipDataMembers(i);
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
	
	const unsigned numSubsets = likelihood->getNumSubsets();
	for (unsigned i = 0; i < numSubsets; ++i)
		{
		delete aTipData;
		delete bTipData;
		delete cTipData;
		delete dTipData;
		DeleteTwoDArray<double> (pre_w_pmat_transposed[i]);
		DeleteTwoDArray<double> (pre_x_pmat_transposed[i]);
		DeleteTwoDArray<double> (pre_y_pmat_transposed[i]);
		DeleteTwoDArray<double> (pre_z_pmat_transposed[i]);
		}
	releaseCondLikeArrays();

	}

double UnimapTopoMove::FourTaxonLnLBeforeMove()
	{
#if 1 || DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING	
	const unsigned numSubsets = likelihood->getNumSubsets();
	aTipData = createTipDataFromUnivents(getUniventsVectorConstRef(*a), aTipData);
	bTipData = createTipDataFromUnivents(getUniventsVectorConstRef(*b), bTipData);
	cTipData = createTipDataFromUnivents(getUniventsVectorConstRef(*c), cTipData);
	dTipData = createTipDataFromUnivents(getUniventsVectorConstRef(*d), dTipData);		

	double lnlike = 0.0;
	for (unsigned i = 0; i < numSubsets; ++i)
		{
		lnlike += FourTaxonLnLFromCorrectTipDataMembers(i);
		
		storePMatTransposed(pre_y_pmat_transposed[i], (const double ***) aTipData->getTransposedPMatrices(i), i); // we cache this, but never use it!!!!
		storePMatTransposed(pre_x_pmat_transposed[i], (const double ***) bTipData->getTransposedPMatrices(i), i); 
		storePMatTransposed(pre_w_pmat_transposed[i], (const double ***) cTipData->getTransposedPMatrices(i), i);
		storePMatTransposed(pre_z_pmat_transposed[i], (const double ***) dTipData->getTransposedPMatrices(i), i);

		}
    return lnlike;
#else
	return 0.0;
#endif
	}

bool UnimapTopoMove::CheckWithPaup(double lnlike, unsigned subsetIndex)
{
	std::ofstream nxsf;
	nxsf.open("debug-with-paup.nex");
	nxsf << "#NEXUS\n";
	DebugSaveNexusFile(nxsf, lnlike, subsetIndex);
	nxsf.close();
	nxsf << "begin PAUP;\n quit ;\nend;\n";
	return (system("./checkWithPAUP.py") == 0);
}

void UnimapTopoMove::DebugSaveNexusFile(std::ostream & nxsf, double lnlike, unsigned subsetIndex)
    {
#if 1 || DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
    typedef boost::shared_array<const int8_t> StateArr;
    state_list_t & adata = aTipData->getTipStatesArray(subsetIndex);
    state_list_t & bdata = bTipData->getTipStatesArray(subsetIndex);
    state_list_t & cdata = cTipData->getTipStatesArray(subsetIndex);
    state_list_t & ddata = dTipData->getTipStatesArray(subsetIndex);
    unsigned nchar = likelihood->getNumPatterns();
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
    nxsf << boost::str(boost::format("  [!***** phycas lnL = %.8f (model = %s) *****]") % lnlike % model->getModelName()) << std::endl;
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

void UnimapTopoMove::storePMatTransposed(double **& cached, const double *** p_mat_array, unsigned subsetIndex)
	{
#if 1 || DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
	const unsigned nStates = likelihood->getNumStates(subsetIndex);
	if (!cached)
		cached = NewTwoDArray<double>(nStates + 1, nStates + 1);
	for (unsigned i = 0; i < nStates; ++i)
		{
		for (unsigned j = 0; j < nStates; ++j)
			cached[i][j] = p_mat_array[0][i][j];
		}
#endif
	}


/*
| uses the aliases a, b, c, and d to calculate a 4-taxon likelihood.  The likelihood is always calculated as if the tree were:
    //      a  b        
    //       \/         
    //   c   nd_childCLPtr	
    //    \ /             <=== length is always origNode->GetEdgeLen()
    //    nd_parentCLPtr
    //    /
    //   d
*/
double UnimapTopoMove::FourTaxonLnLFromCorrectTipDataMembers(unsigned subsetIndex)
	{
#if 1 || DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
	CondLikelihoodShPtr nd_childCLPtr, nd_parentCLPtr;
	double *** p_mat;
	if (scoringBeforeMove)
		{
		nd_childCLPtr = pre_root_posterior;
		nd_parentCLPtr = pre_cla;
		p_mat = pre_p_mat[subsetIndex];
		}
	else
		{
		nd_childCLPtr = post_root_posterior;
		nd_parentCLPtr = post_cla;
		p_mat = post_p_mat[subsetIndex];
		}
	PHYCAS_ASSERT(dLenNd == d);
	PHYCAS_ASSERT(bLenNd == b);

    likelihood->calcPMatTranspose(subsetIndex, aTipData->getTransposedPMatrices(subsetIndex), aTipData->getConstStateListPos(subsetIndex), aLenNd->GetEdgeLen());
    likelihood->calcPMatTranspose(subsetIndex, bTipData->getTransposedPMatrices(subsetIndex), bTipData->getConstStateListPos(subsetIndex), bLenNd->GetEdgeLen());
    likelihood->calcPMatTranspose(subsetIndex, cTipData->getTransposedPMatrices(subsetIndex), cTipData->getConstStateListPos(subsetIndex), cLenNd->GetEdgeLen());
    likelihood->calcPMatTranspose(subsetIndex, dTipData->getTransposedPMatrices(subsetIndex), dTipData->getConstStateListPos(subsetIndex), dLenNd->GetEdgeLen());

	
	likelihood->calcPMat(subsetIndex, p_mat, origNode->GetEdgeLen());
	const double * const *  childPMatrix = p_mat[0];
	likelihood->calcCLATwoTips(*nd_childCLPtr, *bTipData, *aTipData);
	likelihood->calcCLATwoTips(*nd_parentCLPtr, *dTipData, *cTipData);

	double lnl =  HarvestLnLikeFromCondLikePar(nd_childCLPtr, nd_parentCLPtr, childPMatrix, subsetIndex);
#	if defined CHECK_EACH_CALC_AGAINST_PAUP
		std::cerr << "Scoring " << (scoringBeforeMove ? "before" : "after") << " the proposal\n";
		if (!CheckWithPaup(lnl, subsetIndex))
			{
			PHYCAS_ASSERT(false);
			}
#	endif
	return lnl;
#else // old way
	return 0.0;
#endif
	}
	
double UnimapTopoMove::HarvestLnLikeFromCondLikePar(
  CondLikelihoodShPtr focalCondLike, 
  ConstCondLikelihoodShPtr neighborCondLike, 
  const double * const * childPMatrix,
  unsigned subsetIndex)
	{
#if 1 || DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
	double lnLikelihood = 0.0;
	for (unsigned rep = 0; rep < 1; ++rep){
	LikeFltType * focalNdCLAPtr = focalCondLike->getCLA(); //PELIGROSO
	PHYCAS_ASSERT(focalNdCLAPtr);
	const unsigned num_patterns = likelihood->getPartitionModel()->getNumPatterns(subsetIndex);
	const unsigned num_states = likelihood->getNumStates(subsetIndex);
		// Get state frequencies from model and alias rate category probability array for speed
	const double * stateFreq = &model->getStateFreqs()[0]; //PELIGROSO
	
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
	}
	return lnLikelihood;
#else
	return 0.0;
#endif
	}

	
TipData * UnimapTopoMove::createTipDataFromUnivents(const std::vector<Univents> & uv, TipData *td)
	{
#if 1 || DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
	if (td)
		{
		unsigned subsetIndex = 0;
		for (std::vector<Univents>::const_iterator uvIt = uv.begin(); uvIt != uv.end(); ++uvIt, ++subsetIndex)
			{
			const Univents & u = *uvIt;
			/* this is the one place in which we overwrite the state codes */
			int8_t * stateCodes = const_cast<int8_t *>(&(td->getTipStatesArray(subsetIndex)[0]));
			u.fillStateCodeArray(stateCodes);
			}
		}
	else
		{
		state_list_vect_t tipSpecificStateCode;
		const unsigned numSubsets = likelihood->getNumSubsets();
		tipSpecificStateCode.resize(numSubsets);
		PartitionModelShPtr partModPtr = likelihood->getPartitionModel();
		unsigned subsetIndex = 0;
		state_list_pos_vect_t emptyStateListPosVec(numSubsets);
		for (std::vector<Univents>::const_iterator uvIt = uv.begin(); uvIt != uv.end(); ++uvIt, ++subsetIndex)
			{
			const Univents & u = *uvIt;
			state_list_t & tssc = tipSpecificStateCode[subsetIndex];
			tssc.resize(partModPtr->getNumPatterns(subsetIndex));
			u.fillStateCodeArray(&tssc[0]);
			}
		td =new TipData(	true,
						partModPtr,
						emptyStateListPosVec,												// stateListPosVec
						tipSpecificStateCode,	// stateCodesShPtr
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
		{
		unsigned numModelSubsets = getUniventsVectorConstRef(*x).size();

		for (unsigned i = 0 ; i < numModelSubsets; ++i)
			resampleInternalNodeStates(pre_root_posterior->getCLA(), pre_cla->getCLA(), i);
		}
	}

void UnimapTopoMove::resampleInternalNodeStates(const LikeFltType * root_state_posterior, const LikeFltType * des_cla, unsigned subsetIndex)
{
	const UniventProbMgr & upm = likelihood->GetUniventProbMgrConstRef(subsetIndex);
	Lot & rngRef = *rng;
	
	TreeNode * aPar = a->GetParent();
	PHYCAS_ASSERT(aPar && (aPar == origNode || aPar == origNodePar));
	Univents & ndU = getUniventsRef(*aPar, subsetIndex);
	upm.sampleRootStates(ndU, root_state_posterior, rngRef, false, NULL);
	TreeNode * otherInternal = (aPar == origNode ? origNodePar : origNode);
	PHYCAS_ASSERT(otherInternal);
	Univents & ndPU = getUniventsRef(*otherInternal, subsetIndex);
	const int8_t * nd_states = &(ndU.getEndStatesVecConstRef()[0]);
	upm.sampleDescendantStates(ndPU, (const double **) pre_p_mat[subsetIndex][0], des_cla, nd_states, rngRef);
	
	likelihood->flagNodeWithInvalidUnivents(x);
	getUniventsRef(*x, subsetIndex).setValid(false);
	likelihood->flagNodeWithInvalidUnivents(y);
	getUniventsRef(*y, subsetIndex).setValid(false);
	likelihood->flagNodeWithInvalidUnivents(z);
	getUniventsRef(*z, subsetIndex).setValid(false);
	likelihood->flagNodeWithInvalidUnivents(origNode);
	getUniventsRef(*origNode, subsetIndex).setValid(false);
	likelihood->flagNodeWithInvalidUnivents(origNodePar);
	getUniventsRef(*origNodePar, subsetIndex).setValid(false);
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
		{
		unsigned numModelSubsets = getUniventsVectorConstRef(*a).size();

		for (unsigned i = 0 ; i < numModelSubsets; ++i)
			resampleInternalNodeStates(post_root_posterior->getCLA(), post_cla->getCLA(), i);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
UnimapTopoMove::UnimapTopoMove(TreeLikeShPtr treeLikePtr) : MCMCUpdater(),
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
  pre_x_pmat_transposed(),
  pre_y_pmat_transposed(),
  pre_w_pmat_transposed(),
  pre_z_pmat_transposed(),
  doSampleInternalStates(true)
  	{
	is_move = true;
	this->setTreeLikelihood(treeLikePtr);
	isFirstTime = true;
	}

void UnimapTopoMove::setTreeLikelihood(TreeLikeShPtr treeLike)
	{
	PartitionModelShPtr partModel = treeLike->getPartitionModel();
	const unsigned numSubsets = partModel->getNumSubsets();
	pre_x_pmat_transposed.clear();
	pre_y_pmat_transposed.clear();
	pre_w_pmat_transposed.clear();
	pre_z_pmat_transposed.clear();
	pre_p_mat.clear();
	post_p_mat.clear();

	pre_x_pmat_transposed.resize(numSubsets);
	pre_y_pmat_transposed.resize(numSubsets);
	pre_w_pmat_transposed.resize(numSubsets);
	pre_z_pmat_transposed.resize(numSubsets);
	pre_p_mat.resize(numSubsets);
	post_p_mat.resize(numSubsets);
	MCMCUpdater::setTreeLikelihood(treeLike);
	}


void UnimapNNIMove::setLot(LotShPtr p)
	{
	UnimapTopoMove::setLot(p);
	gammaDist.SetLot(rng.get());
	}

void UnimapNNIMove::calculatePairwiseDistances()
	{
	double tmpXY = 0.0;
	double tmpWX = 0.0;
	double tmpXZ = 0.0;
	double tmpWY = 0.0;
	double tmpYZ = 0.0;
	double tmpWZ = 0.0;
	double totalNumPatterns = 0.0;
	PartitionModelShPtr partModel = likelihood->getPartitionModel();
	for (unsigned subsetIndex = 0; subsetIndex < likelihood->getNumSubsets(); ++subsetIndex)
		{
		const double num_patterns = (double) partModel->getNumPatterns(subsetIndex);
		calculatePairwiseDistancesForSubset(subsetIndex);
		tmpXY += dXY*num_patterns;
		tmpWX += dWX*num_patterns;
		tmpXZ += dXZ*num_patterns;
		tmpWY += dWY*num_patterns;
		tmpWZ += dWZ*num_patterns;
		tmpYZ += dYZ*num_patterns;
		totalNumPatterns += num_patterns;
		}
	dXY = tmpXY /totalNumPatterns;
	dWX = tmpWX /totalNumPatterns;
	dXZ = tmpXZ /totalNumPatterns;
	dWY = tmpWY /totalNumPatterns;
	dWZ = tmpWZ /totalNumPatterns;
	dYZ = tmpYZ /totalNumPatterns;
	
//	std::cerr << "pairwiseDiffs: " << dXY << ' ' << dWX << ' ' << dXZ << ' ' << dWY << ' ' << dWZ << ' ' << dYZ << '\n';
	}

void UnimapNNIMove::calculatePairwiseDistancesForSubset(unsigned subsetIndex)
	{
#if 1 || DISABLED_UNTIL_UNIMAP_WORKING_WITH_PARTITIONING
	dXY = dWX =  dXZ =  dWY =  dYZ =  dWZ = 0.0;
	/* This is called before the swap so "x" is "b" and "z" is "d" */
	const int8_t * xStates = &(bTipData->getTipStatesArray(subsetIndex)[0]);
	const int8_t * yStates = &(aTipData->getTipStatesArray(subsetIndex)[0]);
	const int8_t * wStates = &(cTipData->getTipStatesArray(subsetIndex)[0]);
	const int8_t * zStates = &(dTipData->getTipStatesArray(subsetIndex)[0]);
	PartitionModelShPtr partModel = likelihood->getPartitionModel();
	const unsigned num_patterns = partModel->getNumPatterns(subsetIndex);
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
	if (before_swap)
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
		xLen = proposeEdgeLen(propMeanX);
		yLen = proposeEdgeLen(propMeanY);
		zLen = proposeEdgeLen(propMeanZ);
		ndPLen = proposeEdgeLen(propMeanW);
		ndLen = proposeEdgeLen(propMeanInternal);
		}
	else
		{
		xLen = proposeEdgeLen(propMeanX);
		yLen = prev_y_len; // proposeEdgeLen(propMeanY);
		zLen = proposeEdgeLen(propMeanZ);
		ndPLen = proposeEdgeLen(propMeanW);
		ndLen = proposeEdgeLen(propMeanInternal);
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
UnimapNNIMove::UnimapNNIMove(TreeLikeShPtr treeLike) : UnimapTopoMove(treeLike)
	{
	min_edge_len_mean = 0.02;
	edge_len_prop_cv = 1;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move.
*/
double UnimapNNIMove::getLnHastingsRatio() const
	{
	if (save_debug_info)
		{
		std::cerr << "ln_density_reverse_move - ln_density_forward_move;\n";
		std::cerr << ln_density_reverse_move << "   " << ln_density_forward_move << "\n\n";
		}

	return ln_density_reverse_move - ln_density_forward_move;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Jacobian for this move.
*/
double UnimapNNIMove::getLnJacobian() const
	{
	return 0.0;
	}




void UnimapLargetSimonMove::setLot(LotShPtr p)
	{
	UnimapTopoMove::setLot(p);
	}


/*******************************************************************************
*	the move selects a three-edge path: X ---- origNode --- origNodePar ----- (W or Z)
*	 the edge connecting the last selected node (W or Z) will be referred to as the lowerEdge
*/
void UnimapLargetSimonMove::ProposeStateWithTemporaries(ChainManagerShPtr & p)
	{
	wOnThreeEdgePath = rng->Boolean();
	TreeNode * lowerEdgeNd = (wOnThreeEdgePath ? origNodePar : z );
	const double lowerEdgeLen =  lowerEdgeNd->GetEdgeLen();
	const double m = prev_x_len + prev_nd_len + lowerEdgeLen;
	expandContractFactor = exp(lambda*(rng->Uniform(FILE_AND_LINE) - 0.5));
	const double mstar = m*expandContractFactor;
	detachUpperNode = rng->Boolean();
	const double newPlacement = mstar*rng->Uniform(FILE_AND_LINE);
	topoChanging = false;
	if (detachUpperNode)
		{
		const double lowerEdgeLenStar = lowerEdgeLen*expandContractFactor;
		if (newPlacement > lowerEdgeLenStar)
			{// no topo change
			x->SetEdgeLen(mstar - newPlacement);
			origNode->SetEdgeLen(newPlacement - lowerEdgeLenStar);
			lowerEdgeNd->SetEdgeLen(lowerEdgeLenStar);
			}
		else
			{// topo change
			x->SetEdgeLen(mstar - lowerEdgeLenStar);
			origNode->SetEdgeLen(lowerEdgeLenStar - newPlacement);
			lowerEdgeNd->SetEdgeLen(newPlacement);
			topoChanging = true;
			}
		}
	else
		{
		const double upperEdgeLenStar = (prev_nd_len + lowerEdgeLen)*expandContractFactor;
		if (newPlacement < upperEdgeLenStar)
			{// no topo change
			x->SetEdgeLen(mstar - upperEdgeLenStar);
			origNode->SetEdgeLen(upperEdgeLenStar - newPlacement);
			lowerEdgeNd->SetEdgeLen(newPlacement);
			}
		else
			{// topo change
			x->SetEdgeLen(mstar - newPlacement);
			origNode->SetEdgeLen(newPlacement - upperEdgeLenStar);
			lowerEdgeNd->SetEdgeLen(upperEdgeLenStar);
			topoChanging = true;
			}
		}
	if (topoChanging)
		{
		if (wOnThreeEdgePath)
			{ // w and y should be "sister".  So swap the a (y's alias) and d (z's alias)
			std::swap(aTipData, dTipData);
			std::swap(a, d);
			std::swap(aLenNd, dLenNd);
			swapYwithZ = true;
			}
		else
			{ // w and x should be "sister".  So swap the b (x's alias) and d (z's alias)
			std::swap(bTipData, dTipData);
			std::swap(b, d);
			std::swap(bLenNd, dLenNd);
			swapYwithZ = false;
			}
		}
	curr_ln_prior 	= calcEdgeLenLnPrior(*x, x->GetEdgeLen(), p)
					+ calcEdgeLenLnPrior(*y, y->GetEdgeLen(), p)
					+ calcEdgeLenLnPrior(*z, z->GetEdgeLen(), p)
					+ calcEdgeLenLnPrior(*origNode, origNode->GetEdgeLen(), p)
					+ calcEdgeLenLnPrior(*origNodePar, origNodePar->GetEdgeLen(), p);

	
	}
	

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void UnimapLargetSimonMove::accept()
	{
	MCMCUpdater::accept();
	if (topoChanging)
		{
		if (swapYwithZ)
			tree_manipulator.NNISwap(y, z);
		else
			tree_manipulator.NNISwap(z, x);
		}
    
	PHYCAS_ASSERT(origNode);
	PHYCAS_ASSERT(origNode->GetParent() == origNodePar);

	if (doSampleInternalStates)
		{
		unsigned numModelSubsets = getUniventsVectorConstRef(*a).size();

		for (unsigned i = 0 ; i < numModelSubsets; ++i)
			resampleInternalNodeStates(post_root_posterior->getCLA(), post_cla->getCLA(), i);
		}

	}

	
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
UnimapLargetSimonMove::UnimapLargetSimonMove(TreeLikeShPtr treeLike) : UnimapTopoMove(treeLike)
	{
	lambda = 0.2;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move.
*/
double UnimapLargetSimonMove::getLnHastingsRatio() const
	{
	return 3.0*(std::log(expandContractFactor));	
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Jacobian for this move.
*/
double UnimapLargetSimonMove::getLnJacobian() const
	{
	return 0.0;
	}
	

void FillSetWithNdsInQuartet(TreeNode *nd, std::set<TreeNode*> & ndForbidden)
{
	ndForbidden.insert(nd);
	
	TreeNode * par = nd->GetParent();
	PHYCAS_ASSERT(par);
	ndForbidden.insert(par);
	
	TreeNode * scratch = par->GetLeftChild(); 
	PHYCAS_ASSERT(scratch);
	ndForbidden.insert(scratch);

	scratch = scratch->GetRightSib(); 
	PHYCAS_ASSERT(scratch);
	ndForbidden.insert(scratch);
	
	PHYCAS_ASSERT(scratch->GetRightSib() == 0L); // does not work with polytomies
	
	scratch = par->GetParent(); 
	PHYCAS_ASSERT(scratch);
	ndForbidden.insert(scratch);

	scratch = nd->GetLeftChild(); 
	PHYCAS_ASSERT(scratch);
	ndForbidden.insert(scratch);

	scratch = scratch->GetRightSib(); 
	PHYCAS_ASSERT(scratch);
	ndForbidden.insert(scratch);

	PHYCAS_ASSERT(scratch->GetRightSib() == 0L); // does not work with polytomies
}

/// assumes nd is from randomInternalAboveSubroot
bool UnimapTopoMoveSpreader::conflictsWithPrevious(TreeNode *nd) const 
	{
	std::set<TreeNode*> ndForbidden;
	FillSetWithNdsInQuartet(nd, ndForbidden);
	
	for(std::set<TreeNode*>::const_iterator nIt = selectedNodes.begin(); nIt != selectedNodes.end(); ++nIt)
		{
		TreeNode * prev = *nIt;
		if(prev == 0L)
			continue;
		std::set<TreeNode*> prevForbidden;
		FillSetWithNdsInQuartet(prev, prevForbidden);
		
		if (ndForbidden.find(prev) != ndForbidden.end())
			return true;
		if (ndForbidden.find(prev->GetParent()) != ndForbidden.end())
			return true;
		if (prevForbidden.find(nd) != prevForbidden.end())
			return true;
		if (prevForbidden.find(nd->GetParent()) != prevForbidden.end())
			return true;
		}
	return false;
	}

/*struct LargetSimonCallable
	{
		LargetSimonCallable(UnimapTopoMove *m) :move (m) {}
		
		void operator()();
		
	private:
		UnimapTopoMove * move; 
	}
	
void LargetSimonCallable::operator()()
	{
	move->update();
	}
*/	
#if defined(USING_THREAD_BARRIERS) && USING_THREAD_BARRIERS
	class DualBarrierUpdater 
		{
		public:
			DualBarrierUpdater(MCMCUpdater * m, boost::barrier * startBarrier, boost::barrier * afterBarrier, boost::mutex & mx)
				:number(0),
				theStartUpdateBarrier(startBarrier),
				theAfterUpdateBarrier(afterBarrier),
				io_mutex(mx),
				move(m)
				//,zzz(0.0)
				{}
	
			void operator()() 
				{
				while (true)
					{
					//std::cerr << "DBU # " << number << " start =" << (long)theStartUpdateBarrier << "\n";
					if (theStartUpdateBarrier)
						{
						theStartUpdateBarrier->wait();
						//if (theStartUpdateBarrier->wait())
						//	{
						//	boost::mutex::scoped_lock scoped_lock(io_mutex);
						//	std::cerr << boost::str(boost::format("before >> Thread %d received TRUE from theStartUpdateBarrier->wait()") % number) << "\n";
						//	}
						//else 
						//	{
						//	boost::mutex::scoped_lock scoped_lock(io_mutex);
						//	std::cerr << boost::str(boost::format("before >> Thread %d received FALSE from theStartUpdateBarrier->wait()") % number) << "\n";
						//	}
						}

					if (move)
						move->update();
						
					//for (unsigned z = 0; z < 2000000; z++)
					//	{
					//	for (unsigned zz = 0; zz < 2000000; zz++)
					//		{
					//		zzz += log(z);
					//		zzz += log(zz);
					//		}	
					//	}	
						
					//std::cerr << "DBU # " << number << " after =" << (long)theAfterUpdateBarrier << "\n";
					// while this wait blocks, the "start" barrier will be installed
					if (theAfterUpdateBarrier) 
						{
						theAfterUpdateBarrier->wait();
						//if (theAfterUpdateBarrier->wait())
						//	{
						//	boost::mutex::scoped_lock scoped_lock(io_mutex);
						//	std::cerr << boost::str(boost::format("after >> Thread %d received TRUE from theAfterUpdateBarrier->wait()") % number) << "\n";
						//	}
						//else 
						//	{
						//	boost::mutex::scoped_lock scoped_lock(io_mutex);
						//	std::cerr << boost::str(boost::format("after >> Thread %d received FALSE from theAfterUpdateBarrier->wait()") % number) << "\n";
						//	}
						}
					else
						{
						PHYCAS_ASSERT(0);
						//std::cerr << "DBU # " << number << " exiting for lack of afterBarrier\n";
						break;
						}
					}
				}

			unsigned number;
			
		private:
			boost::mutex	& io_mutex; 
			boost::barrier	* theStartUpdateBarrier;
			boost::barrier	* theAfterUpdateBarrier;
			MCMCUpdater		* move;
			//double			  zzz;
	};

	class SimplePtrCaller {
		public:
			SimplePtrCaller(DualBarrierUpdater * m): dbu(m){}
	
			void operator()() 
			{
				if (dbu)
					(*dbu)();
			}
		private:
			DualBarrierUpdater * dbu;
	};
#else 
	class CU {
		public:
			CU(UnimapTopoMove * m): move(m){}
	
			void operator()() 
			{
				if (move)
					move->update();
			}
		UnimapTopoMove * move;
	};
#endif

// Not currently working for some reason...
//
// void UnimapTopoMoveSpreader::debugShowSelectedSubtrees()
// 	{
// 	std::set<TreeNode *>::iterator ndIt;
// 	for (ndIt = selectedNodes.begin(); ndIt != selectedNodes.end(); ++ndIt)
// 		{
// 		(*ndIt)->SelectNode();
// 		}
// 	likelihood->startTreeViewer(tree, "Nodes selected by UnimapTopoMoveSpreader");	
// 	for (ndIt = selectedNodes.begin(); ndIt != selectedNodes.end(); ++ndIt)
// 		{
// 		(*ndIt)->UnselectNode();
// 		}
// 	}

bool UnimapTopoMoveSpreader::update()
	{
	//std::cerr << "UnimapTopoMoveSpreader::update " << this->getName() << '\n';
	PHYCAS_ASSERT(!topoMoves.empty());
	const unsigned numUpdaters = topoMoves.size();
	
	selectedNodes.clear();
	queued.clear();

	for(unsigned i = 0; i < numUpdaters; ++i)
		{
		PHYCAS_ASSERT(tree);
		PHYCAS_ASSERT(rng);
		TreeNode * nextNd = randomInternalAboveSubroot(*tree, *rng);
		if (conflictsWithPrevious(nextNd))
			{ // give it one more shot...
			nextNd = randomInternalAboveSubroot(*tree, *rng);
			if (conflictsWithPrevious(nextNd))
				nextNd = 0L;
			}
		queued.push_back(nextNd);
		selectedNodes.insert(nextNd);
		}
	PHYCAS_ASSERT(queued.size() == numUpdaters);
	
#	if defined(USING_THREAD_BARRIERS) && USING_THREAD_BARRIERS
		if (threadVec.empty())
			{
			PHYCAS_ASSERT(topoMovesVec.empty());
			PHYCAS_ASSERT(updaterWithBarriers.empty());
			PHYCAS_ASSERT(afterUpdateBarrier == 0L);
			startUpdateBarrier = new boost::barrier(1 + numUpdaters);
			afterUpdateBarrier = new boost::barrier(1 + numUpdaters);
			unsigned i = 0;
			unsigned prevSeed = rng->GetSeed();
			for (std::set<UnimapTopoMove *>::iterator upIt = topoMoves.begin(); upIt != topoMoves.end(); ++upIt, ++i)
				{
				PHYCAS_ASSERT(*upIt);
				UnimapTopoMove *move = *upIt;
				
				// give the move a new Lot object, so that the threads will not compete for random numbers from the same Lot
				//	this should make runs repeatable...
				LotShPtr nl(new Lot());
				unsigned seedForMove;
				unsigned seedOffset = 10*(1+i);
				if (prevSeed > seedOffset + 1)
					seedForMove = prevSeed - seedOffset;
				else
					seedForMove = UINT_MAX - seedOffset;
				nl->SetSeed(seedForMove);
				move->setLot(nl);
				
				TreeNode * nodeForMove = queued[i];
				PHYCAS_ASSERT(nodeForMove);
				InternalData * nd_internal_data = nodeForMove->GetInternalData();
				PHYCAS_ASSERT(nd_internal_data);
				move->acquireCondLikeArrays(*nd_internal_data);
				
				topoMovesVec.push_back(move);
				DualBarrierUpdater * dbu = new DualBarrierUpdater(move, startUpdateBarrier, afterUpdateBarrier, io_mutex);
				dbu->number = i;
				updaterWithBarriers.push_back(dbu);
				SimplePtrCaller x(dbu);
				threadVec.push_back(new boost::thread(x));
				}
			}
		
		PHYCAS_ASSERT(threadVec.size() == numUpdaters);
		PHYCAS_ASSERT(updaterWithBarriers.size() == numUpdaters);
		PHYCAS_ASSERT(topoMovesVec.size() == numUpdaters);
		PHYCAS_ASSERT(startUpdateBarrier);
		PHYCAS_ASSERT(afterUpdateBarrier);
		
		for (unsigned threadInd = 0; threadInd < numUpdaters; ++threadInd)
			{
			UnimapTopoMove * move = topoMovesVec[threadInd];
			PHYCAS_ASSERT(move);
			TreeNode * nodeForMove = queued[threadInd];
			move->queueNode(nodeForMove);
			DualBarrierUpdater * dbu = updaterWithBarriers[threadInd];
			PHYCAS_ASSERT(dbu);
			}
		
		//{
		//	boost::mutex::scoped_lock scoped_lock(io_mutex);
		//	std::cerr << "Threads about to wait on start barrier" << (long)startUpdateBarrier << "\n";
		//}
		// trigger the realease of the other threads
		startUpdateBarrier->wait(); // this should return instantly
		//if (startUpdateBarrier->wait()) // this should return instantly
		//	{
		//	boost::mutex::scoped_lock scoped_lock(io_mutex);
		//	std::cerr << "before >> Master received TRUE from startUpdateBarrier->wait()" << "\n";
		//	}
		//else 
		//	{
		//	boost::mutex::scoped_lock scoped_lock(io_mutex);
		//	std::cerr << "before >> Master received FALSE from startUpdateBarrier->wait()" << "\n";
		//	}

		//std::cerr << "Master about to wait on afterUpdateBarrier barrier" << (long)afterUpdateBarrier << "\n";
		afterUpdateBarrier->wait();
		//if (afterUpdateBarrier->wait())
		//	{
		//	boost::mutex::scoped_lock scoped_lock(io_mutex);
		//	std::cerr << "after >> Master received TRUE from afterUpdateBarrier->wait()" << "\n";
		//	}
		//else 
		//	{
		//	boost::mutex::scoped_lock scoped_lock(io_mutex);
		//	std::cerr << "after >> Master received FALSE from afterUpdateBarrier->wait()" << "\n";
		//	}

		return true;
			
#	else
		unsigned i = 0;
		
		std::vector<CU> cuVec;
		cuVec.reserve(topoMoves.size());
		for (std::set<UnimapTopoMove *>::iterator upIt = topoMoves.begin(); upIt != topoMoves.end(); ++upIt, ++i)
			{
			(*upIt)->queueNode(queued[i]);
			cuVec.push_back(CU(*upIt));
			}
			
		std::vector<boost::thread *> threadVec;
		threadVec.reserve(topoMoves.size());
		for (std::vector<CU>::iterator cuIt = cuVec.begin(); cuIt != cuVec.end(); ++cuIt)
			threadVec.push_back(new boost::thread(*cuIt));
	
		//debugShowSelectedSubtrees();
	
		try
			{
			for (unsigned threadInd = 0; threadInd < threadVec.size(); ++threadInd)
				{
				boost::thread * t = threadVec[threadInd];
				t->join();
				threadVec[threadInd] = 0L;
				delete t;
				}
			}
		catch(...)
			{
			for (std::vector<boost::thread *>::iterator thIt = threadVec.begin(); thIt != threadVec.end(); ++thIt)
				{
				boost::thread * t = *thIt;
				if (t)
					delete t;
				}
			}
#	endif
	return true;
	}

}	// namespace phycas
