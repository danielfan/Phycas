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
#include "phycas/src/unimap_fast_nni_move.hpp"
#include "phycas/src/basic_tree.hpp"
#include "phycas/src/internal_data.hpp"
#include "phycas/src/tip_data.hpp"

#include "boost/format.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
UnimapFastNNIMove::UnimapFastNNIMove()
  : x(0), y(0), a(0), b(0), c(0), d(0), num_states(0), num_sites(0), log_umat(0), tuning_factor(0.5)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
UnimapFastNNIMove::~UnimapFastNNIMove()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or 
|	revert(), whichever is appropriate.
*/
bool UnimapFastNNIMove::update()
	{
	bool accepted = false;
#if DISABLED_UNTIL_WORKING_PRIOR_ACCOMMODATED
	// If this is the first time update has been called, need to initialize num_states, num_sites and log_umat
	if (num_states == 0)
		{
		num_states = likelihood->getNStates();
		PHYCAS_ASSERT(num_states > 0);
		}
		
	if (log_umat == 0)
		{
		log_umat_caretaker.CreateMatrix(num_states, 0.0);
		log_umat = log_umat_caretaker.GetMatrix();
		}
	
	// Make sure number of states hasn't changed since we instantiated this move
	PHYCAS_ASSERT(num_states == likelihood->getNStates());
	
	if (num_sites == 0)
		{
		num_sites = likelihood->getNPatterns();
		PHYCAS_ASSERT(num_sites > 0);
		}
	
	// The only case in which is_fixed is true occurs when the user decides to fix the edge lengths.
	// A proposed UnimapFastNNIMove cannot be accepted without changing edge lengths, so it is best to bail out now.
	if (is_fixed)
		return false;

	++num_moves_attempted;

	proposeNewState();
	
	// 	ChainManagerShPtr p = chain_mgr.lock();
	// 
	// 	double prev_posterior = 0.0;
	// 	double curr_posterior = 0.0;
	// 	if (is_standard_heating)
	// 		{
	// 		prev_posterior = heating_power*(prev_ln_like + prev_ln_prior);
	// 		curr_posterior = heating_power*(curr_ln_like + curr_ln_prior);
	// 		}
	// 	else
	// 		{
	// 		prev_posterior = heating_power*prev_ln_like + prev_ln_prior;
	// 		curr_posterior = heating_power*curr_ln_like + curr_ln_prior;
	// 		}
	// 
	// 	double ln_accept_ratio = curr_posterior - prev_posterior + getLnHastingsRatio();
	double lnu = DBL_MAX;
	//double ln_accept_ratio = 1.1;	//temporary!
	accepted = (ln_accept_ratio >= 0.0);
	if (!accepted)
		{
		double u = rng->Uniform(FILE_AND_LINE);
		lnu = std::log(u);
		accepted = (lnu <= ln_accept_ratio);
		}

	//     if (save_debug_info)
	//         debug_info = str(boost::format("swapping %d <-> %d (%s, lnR = %.5f)") % x->GetNodeNumber() % z->GetNodeNumber() % (accepted ? "accepted" : "rejected") % ln_accept_ratio);
    
    if (accepted)
    	{
		++num_moves_accepted;
		accept();
		}
	else
		revert();
	//std::cerr << num_moves_accepted << " accept decisions out of " << num_moves_attempted << " attempts.\n";
#endif
	return accepted;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Selects an internal node at random from a discrete uniform distribution with the constraint that the returned node
|   is not equal to the subroot (the sole child of the tip node serving as the root).
*/
TreeNode * UnimapFastNNIMove::randomInternalAboveSubroot()
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
|	Tallies up univents for each of the 5 involved in this move and ensures that the sMat structure stored in the node's
|	tip data or internal data structure is correct. Also checks the mdot count stored in each node's univent structure.
*/
void UnimapFastNNIMove::debugCheckUnivents()
	{
	std::cerr << "Fast NNI move (before):" << std::endl;
	std::cerr << "  Edge lengths\n";
	std::cerr << "    a = " << alen_before << "\n";
	std::cerr << "    b = " << blen_before << "\n";
	std::cerr << "    c = " << clen_before << "\n";
	std::cerr << "    x = " << xlen_before << "\n";
	std::cerr << "    y = " << ylen_before << "\n";
	std::cerr << std::endl;
	std::cerr << "  Number of univents per edge:\n";
	std::cerr << "    a = " << mdota_before << "\n";
	std::cerr << "    b = " << mdotb_before << "\n";
	std::cerr << "    c = " << mdotc_before << "\n";
	std::cerr << "    x = " << mdotx_before << "\n";
	std::cerr << "    y = " << mdoty_before << "\n";
	std::cerr << std::endl;

	std::cerr << "  Sum of nodeSMat matrices stored in internal or tip data structures:\n";
	std::cerr << boost::str(boost::format("    %6.0f %6.0f %6.0f %6.0f\n") % smat_before[0][0] % smat_before[0][1] % smat_before[0][2] % smat_before[0][3] );
	std::cerr << boost::str(boost::format("    %6.0f %6.0f %6.0f %6.0f\n") % smat_before[1][0] % smat_before[1][1] % smat_before[1][2] % smat_before[1][3] );
	std::cerr << boost::str(boost::format("    %6.0f %6.0f %6.0f %6.0f\n") % smat_before[2][0] % smat_before[2][1] % smat_before[2][2] % smat_before[2][3] );
	std::cerr << boost::str(boost::format("    %6.0f %6.0f %6.0f %6.0f\n") % smat_before[3][0] % smat_before[3][1] % smat_before[3][2] % smat_before[3][3] );
	std::cerr << std::endl;

	// Tally up univents along the 5 edges as a check of smat_before
	SquareMatrix tmp2(num_states, 0.0);
	addUniventsOneEdge(tmp2, a);
	addUniventsOneEdge(tmp2, b);
	addUniventsOneEdge(tmp2, c);
	addUniventsOneEdge(tmp2, x);
	addUniventsOneEdge(tmp2, y);

	std::cerr << "  Tally of univents along the 5 edges:\n";
	std::cerr << boost::str(boost::format("    %6.0f %6.0f %6.0f %6.0f\n") % tmp2[0][0] % tmp2[0][1] % tmp2[0][2] % tmp2[0][3] );
	std::cerr << boost::str(boost::format("    %6.0f %6.0f %6.0f %6.0f\n") % tmp2[1][0] % tmp2[1][1] % tmp2[1][2] % tmp2[1][3] );
	std::cerr << boost::str(boost::format("    %6.0f %6.0f %6.0f %6.0f\n") % tmp2[2][0] % tmp2[2][1] % tmp2[2][2] % tmp2[2][3] );
	std::cerr << boost::str(boost::format("    %6.0f %6.0f %6.0f %6.0f\n") % tmp2[3][0] % tmp2[3][1] % tmp2[3][2] % tmp2[3][3] );
	std::cerr << std::endl;

	std::cerr << "Aborting because Unimap Fast NNI move not yet completely implemented" << std::endl;
	std::exit(0);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Stores the information needed to compute the log-likelihood for the four-taxon subtree being considered. Assumes
|	that x, y, a, b, c and d have all been assigned.
*/
void UnimapFastNNIMove::storeOrigEdgeInfo()
	{
	double uni_lambda = model->calcLMat(log_umat);
	
	// Save edge lengths as they were before the move is proposed
	alen_before = a->GetEdgeLen();
	blen_before = b->GetEdgeLen();
	clen_before = c->GetEdgeLen();
	xlen_before = x->GetEdgeLen();
	ylen_before = y->GetEdgeLen();

	// The sum of edge lengths is used in the likelihood calculation
	tlen_before = 0.0;
	tlen_before += alen_before;
	tlen_before += blen_before;
	tlen_before += clen_before;
	tlen_before += xlen_before;
	tlen_before += ylen_before;
		
	// Get reference to the univents structure attached to each node
	const Univents & univents_a = getUniventsConstRef(*a);
	const Univents & univents_b = getUniventsConstRef(*b);
	const Univents & univents_c = getUniventsConstRef(*c);
	const Univents & univents_x = getUniventsConstRef(*x);
	const Univents & univents_y = getUniventsConstRef(*y);
	
	// Save mdot for each edge before the move is proposed
	mdota_before = univents_a.getMDot();
	mdotb_before = univents_b.getMDot();
	mdotc_before = univents_c.getMDot();
	mdotx_before = univents_x.getMDot();
	mdoty_before = univents_y.getMDot();
	
	// Initialize smat_before and smat_after to all zeros
	smat_before.CreateMatrix(num_states, 0.0);	//this involves reallocating the matrix - should just zero it out
	smat_after.CreateMatrix(num_states, 0.0);	//this involves reallocating the matrix - should just zero it out

	// get pointers to the sMat data member for each node
	// *** NEED TO WRITE getNodeSMat FUNCTION ***
	unsigned * * smata = NULL; //likelihood->getNodeSMat(a);
	unsigned * * smatb = NULL; //likelihood->getNodeSMat(b);
	unsigned * * smatc = NULL; //likelihood->getNodeSMat(c);
	unsigned * * smatx = NULL; //likelihood->getNodeSMat(x);
	unsigned * * smaty = NULL; //likelihood->getNodeSMat(y);
	
	// add all five sMat matrices to initially empty smat_before
	for (unsigned i = 0; i < num_states; ++i)
		for (unsigned j = 0; j < num_states; ++j)
			{
			double tmp = (double)smata[i][j]
				+ (double)smatb[i][j]
				+ (double)smatc[i][j]
				+ (double)smatx[i][j]
				+ (double)smaty[i][j];
			smat_before[i][j] += tmp;
			}

	double lnL_before	= (double)num_sites*uni_lambda*tlen_before
						+ (double)mdota_before*log(alen_before)
						+ (double)mdotb_before*log(blen_before)
						+ (double)mdotc_before*log(clen_before)
						+ (double)mdotx_before*log(xlen_before)
						+ (double)mdoty_before*log(ylen_before);
							
	// the following function exits after displaying debug information
	debugCheckUnivents();
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Adds the univents recorded in the univents structure of node `nd' to the corresponding elements of the num_states 
|	by num_states matrix `smat'.
*/
void UnimapFastNNIMove::addUniventsOneEdge(SquareMatrix & smat, TreeNode * nd)
	{
	// code taken from TreeLikelihood::debugCheckSMatrix

	// get reference to the univents structure for nd
	const Univents & u = getUniventsConstRef(*nd);
	PHYCAS_ASSERT(u.isValid());
	
	// get reference to the 2-d vector of univents: 
	// univents[i][j] holds a state for univent j at site i 
	const std::vector<StateMapping> &  univents	= u.getVecEventsVecConstRef();
	
	// get reference for the vector of starting states
	const int8_t * starting_state = NULL;
	if (nd == tree->GetSubroot())
		{
		// if the current node is the subroot, the starting states are the 
		// observed states at the tip root
		TreeNode * root_tip = tree->GetRoot();
		const TipData &	root_tip_data = *(root_tip->GetTipData());
		starting_state = root_tip_data.getConstStateCodes();
		}
	else
		{
		// if the current node is not the subroot, then the starting states
		// are the ending state of the node's parent
		TreeNode * nd_par = nd->GetParent();
		PHYCAS_ASSERT(nd_par);
		const Univents & upar = getUniventsConstRef(*nd_par);
		PHYCAS_ASSERT(upar.isValid());
		const std::vector<int8_t> & states_vec	= upar.getEndStatesVecConstRef();
		starting_state = &states_vec[0];
		}
		
	// loop over sites, adding the univents for that site on the focal branch to smat
	for (std::vector<StateMapping>::const_iterator sit = univents.begin(); sit != univents.end(); ++sit, ++starting_state)
		{
		// *sit is a vector of states at each univent for the current site
		const StateMapping & stlist = (*sit);
		if (!stlist.empty())
			{
			//unsigned k = 0;
			int8_t prev_state = *starting_state;
			const StateMapping::const_iterator endIt = stlist.end();
			for (StateMapping::const_iterator it = stlist.begin(); it != endIt; ++it)
				{
				const int8_t new_state = *it;
				smat[prev_state][new_state] += 1.0;
				//std::cerr << "--->  | s[" << (k++) << "] = " << (int)prev_state << " -> " << (int)new_state << '\n';
				prev_state = new_state;
				}
			}
		}
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
double UnimapFastNNIMove::calcFourTaxonLogLikelihood()
	{
//#error not yet finished
	return 0.0;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Proposes a Larget-Simon LOCAL move. 
*/
void UnimapFastNNIMove::proposeNewState()
	{
	// Find a random 4-taxon subtree within the larger tree
	//
	//      a  b
	//       \/
	//   c   x <- randomly-chosen non-subroot internal node
	//    \ /
	//     y
	//    /
	//   d
	x = randomInternalAboveSubroot();
	PHYCAS_ASSERT(x);
	PHYCAS_ASSERT(x->CountChildren() == 2); // we haven't figured this out for polytomies

	y = x->GetParent();
	PHYCAS_ASSERT(y);
	PHYCAS_ASSERT(y->CountChildren() == 2); // we haven't figured this out for polytomies
	
	d = y->GetParent();
	PHYCAS_ASSERT(d);
	
	c = x->FindNextSib();
	PHYCAS_ASSERT(c);
	
    a = x->GetLeftChild();
	PHYCAS_ASSERT(a);
	
	b = a->GetRightSib();
	PHYCAS_ASSERT(b);
	
	storeOrigEdgeInfo();
	lnL_before = calcFourTaxonLogLikelihood();
	
	// Decide whether upper terminus should be a or b
	double u1 = rng->Uniform(FILE_AND_LINE);
	a_is_top = (u1 <= 0.5 ? true : false);
	
	// Decide whether lower terminus should be c or d
	double u2 = rng->Uniform(FILE_AND_LINE);
	c_is_bottom = (u2 <= 0.5 ? true : false);
	
	// Decide whether x or y should do the sliding
	double u3 = rng->Uniform(FILE_AND_LINE);
	sliding_x = (u3 <= 0.5 ? true : false);
	
	// Draw r, the multiplicative factor by which the chosen three-edge segment is expanded or contracted
	double u4 = rng->Uniform(FILE_AND_LINE);
	three_edge_length_ratio = exp(tuning_factor*(u4 - 0.5));
	
	// Deal with the 8 possible cases (could make this more compact, but the long way makes it easier to diagnose bugs)
	// Note that the tree and branch lengths are not actually modified by this proposed move unless the move is accepted.
	// Everything needed for computing the likelihood is stored in temporaries (e.g. alen_before, alen_after, etc.)
	which_case = 0;
	alen_after = alen_before;
	blen_after = blen_before;
	clen_after = clen_before;
	xlen_after = xlen_before;
	ylen_after = ylen_before;
	if (sliding_x)
		{
		if (a_is_top)
			{
			if (c_is_bottom)
				{
				// case 1: sliding x, a is top, c is bottom
				//       a  b
				//       \\/
				//   c    x* 
				//    \\ //
				//     y
				//    /
				//   d
				which_case = 1;
				
				// modify edge lengths
				alen_after = alen_before*three_edge_length_ratio;
				xlen_after = xlen_before*three_edge_length_ratio;
				clen_after = clen_before*three_edge_length_ratio;
#				error begin again here
				}
			else 
				{
				// case 2: sliding x, a is top, d is bottom
				//       a  b
				//       \\/
				//   c    x* 
				//    \  //
				//     y
				//    //
				//   d
				which_case = 2;
				
				// modify edge lengths
				alen_after = alen_before*three_edge_length_ratio;
				xlen_after = xlen_before*three_edge_length_ratio;
				ylen_after = ylen_before*three_edge_length_ratio;
				}
			}
		else
			{
			if (c_is_bottom)
				{
				// case 3: sliding x, b is top, c is bottom
				//       a  b
				//        \//
				//   c    x*
				//    \\ //
				//     y
				//    /
				//   d
				which_case = 3;
				
				// modify edge lengths
				blen_after = blen_before*three_edge_length_ratio;
				xlen_after = xlen_before*three_edge_length_ratio;
				clen_after = clen_before*three_edge_length_ratio;
				}
			else 
				{
				// case 4: sliding x, b is top, d is bottom
				//       a  b
				//        \//
				//   c    x*
				//    \ //
				//     y
				//    //
				//   d
				which_case = 4;
				
				// modify edge lengths
				blen_after = blen_before*three_edge_length_ratio;
				xlen_after = xlen_before*three_edge_length_ratio;
				ylen_after = ylen_before*three_edge_length_ratio;
				}
			}
		}
	else
		{
		if (a_is_top)
			{
			if (c_is_bottom)
				{
				// case 5: sliding y, a is top, c is bottom
				//       a  b
				//       \\/
				//   c    x 
				//    \\ //
				//     y*
				//    /
				//   d
				which_case = 5;
				
				// modify edge lengths
				alen_after = alen_before*three_edge_length_ratio;
				xlen_after = xlen_before*three_edge_length_ratio;
				clen_after = clen_before*three_edge_length_ratio;
				}
			else 
				{
				// case 6: sliding y, a is top, d is bottom
				//       a  b
				//       \\/
				//   c    x 
				//    \  //
				//     y*
				//    //
				//   d
				which_case = 6;
				
				// modify edge lengths
				alen_after = alen_before*three_edge_length_ratio;
				xlen_after = xlen_before*three_edge_length_ratio;
				ylen_after = ylen_before*three_edge_length_ratio;
				}
			}
		else
			{
			if (c_is_bottom)
				{
				// case 7: sliding y, b is top, c is bottom
				//       a  b
				//        \//
				//   c    x
				//    \\ //
				//     y*
				//    /
				//   d
				which_case = 7;
				
				// modify edge lengths
				blen_after = blen_before*three_edge_length_ratio;
				xlen_after = xlen_before*three_edge_length_ratio;
				clen_after = clen_before*three_edge_length_ratio;
				}
			else 
				{
				// case 8: sliding y, b is top, d is bottom
				//       a  b
				//        \//
				//   c    x
				//    \  //
				//     y*
				//    //
				//   d
				which_case = 8;
				
				// modify edge lengths
				blen_after = blen_before*three_edge_length_ratio;
				xlen_after = xlen_before*three_edge_length_ratio;
				ylen_after = ylen_before*three_edge_length_ratio;
				}
			}
		}
	}

double UnimapFastNNIMove::calcEdgeLenLnPrior(const TreeNode &x, double edge_len, ChainManagerShPtr & chain_mgr) 
	{
	if (x.IsTip())
		return chain_mgr->calcExternalEdgeLenPriorUnnorm(edge_len);
	return chain_mgr->calcInternalEdgeLenPriorUnnorm(edge_len);
	}

TipData * UnimapFastNNIMove::createTipDataFromUnivents(const Univents & u, TipData *td)
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
void UnimapFastNNIMove::revert()
	{
	MCMCUpdater::revert();
//	std::cerr << "REVERTED" << std::endl;
	
// 	x->SetEdgeLen(prev_x_len);
// 	y->SetEdgeLen(prev_y_len);
// 	z->SetEdgeLen(prev_z_len);
// 	PHYCAS_ASSERT(origNode ==  y->GetParent());
// 	PHYCAS_ASSERT(origNode);
// 	PHYCAS_ASSERT(origNodePar == origNode->GetParent());
// 	origNode->SetEdgeLen(prev_nd_len);
// 	origNodePar->SetEdgeLen(prev_ndP_len);
// 	
// 	if (doSampleInternalStates)
// 		resampleInternalNodeStates(pre_root_posterior->getCLA(), pre_cla->getCLA());
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
void UnimapFastNNIMove::accept()
	{
	MCMCUpdater::accept();
// 	tree_manipulator.NNISwap(z, x);
//     
// 	PHYCAS_ASSERT(origNode);
// 	PHYCAS_ASSERT(origNode->GetParent() == origNodePar);
// 
// 	if (doSampleInternalStates)
// 		resampleInternalNodeStates(post_root_posterior->getCLA(), post_cla->getCLA());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move.
*/
double UnimapFastNNIMove::getLnHastingsRatio() const
	{
	return ln_density_reverse_move - ln_density_forward_move;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Jacobian for this move.
*/
double UnimapFastNNIMove::getLnJacobian() const
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
void UnimapFastNNIMove::setLot(LotShPtr p)
	{
	MCMCUpdater::setLot(p);
	}

void UnimapFastNNIMove::DebugSaveNexusFile(TipData * xtd, TipData * ytd, TipData * ztd, TipData * wtd, double lnlike)
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
    //TODO
	//nxsf << boost::str(boost::format("  utree curr = (a:%.8f, b:%.8f, (c:%.8f, d:%.8f):%.8f);") % aLenNd->GetEdgeLen() % bLenNd->GetEdgeLen() % cLenNd->GetEdgeLen() % dLenNd->GetParent()->GetEdgeLen() % origNode->GetEdgeLen()) << std::endl;
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

    nxsf.close();
    }

}	// namespace phycas
