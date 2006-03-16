#include "phycas/force_include.h"
//#define RECORD_MOVELOG

#include <vector>
#include "phycas/modules/mcmc/bush_master.hpp"
#include "phycas/rand/lot.hpp"
#include "phycas/trees/tree_node.hpp"
using std::vector;
using std::pair;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::log;
	using std::pow;
	using std::exp;
#endif

/*--------------------------------------------------------------------------------------------------------------------------
|	Passes along supplied pointer to a fully-resolved Tree to the TreeMove constructor. Ascertains `numTaxa' from the
|	supplied tree and sets the data member `numNodesInFullyResolvedTree' accordingly. Also supplies `numTaxa' to 
|	ComputePolytomyDistribution. Sets `edgeLenMean' and instantiates edge length proposal distribution `expDist'. Assumes 
|	`t' is non-NULL. 
*/
BushMaster::BushMaster(
  LotShPtr r,	/**< is the pseudorandom number generator, which is passed along to the base class member function MCMCMove::SetRandomNumberGenerator */
  Tree *t)	/**< is a pointer to the Tree object, passed along to base class TreeMove constructor */
  : TreeMove(t)
	{
	assert(t != NULL);
	assert(r);
	MCMCMove::SetRandomNumberGenerator(r);

	numTaxa = t->GetNLeaves();
	numNodesInFullyResolvedTree = 2*(numTaxa - 1);

	ComputePolytomyDistribution(numTaxa);

	edgeLenMean = 1.0; //@POL should let user choose this
	expDist.SetMeanAndVariance(edgeLenMean, edgeLenMean);
	Lot &lotref = *r;
	expDist.SetLot(&lotref);	//POL 22Jun2005

	Clear();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Nothing to be done.
*/
BushMaster::~BushMaster()
	{
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is accepted. Simply calls Clear function.
*/
void BushMaster::Accept()
	{
#if defined(RECORD_MOVELOG)
	ofstream tmpf("movelog.txt", ios::out | ios::app);

	if (birthMoveProposed)
		tmpf << "   *** BIRTH MOVE ACCEPTED ***" << endl;
	else
		tmpf << "   *** DEATH MOVE ACCEPTED ***" << endl;

	tmpf.close();
#endif

	Clear();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Causes tree to be returned to its state just prior to proposing the move.
*/
void BushMaster::Revert()
	{
#if defined(RECORD_MOVELOG)
	ofstream tmpf("movelog.txt", ios::out | ios::app);
#endif

	if (birthMoveProposed)
		{
		// A birth move was previously proposed, which means a new internal node (origLChild) and its edge were created. 
		// origLChild's parent is known as origPar. The children of origLChild represent a subset of the original children 
		// of origPar. To reverse the birth move, we need only transfer all children of origLChild to origPar, then store 
		// origLChild for use in some later birth move proposal
		//
#if defined(RECORD_MOVELOG)
		tmpf << "   *** BIRTH MOVE REJECTED ***" << endl;
#endif

		// Transfer all children of origLChild to origPar
		//
		while (origLChild->GetLeftChild() != NULL)
			{
			TreeNode *s = origLChild->GetLeftChild();
			DetachSubtree(s);
			InsertSubtree(s, origPar, TreeManip::kOnRight);
			}

		// Now delete origLChild
		//
		DetachSubtree(origLChild);
		tree->StoreTreeNode(origLChild);
		//origPar->InvalidateAttrDown(true, origPar->GetParent());	// no need to invalidate transition probability matrix
		origPar->InvalidateCondLikeArrays();	// invalidates conditional likelihood arrays only

		tree->InvalidateNodeCounts();
		}
	else
		{
		// A death move was proposed, which means an internal node and its edge have been deleted. This deleted
		// node's parent is known as origPar. The children of the deleted node have been added as children of origPar,
		// with the first being identified now as origLChild and the last in the series now pointed to by origRChild.
		// To revert the death move, create a new child of origPar having edgelen equal to origEdgelen, then transfer
		// the nodes from origLChild to origRChild (following rSib pointers, and including both origLChild and 
		// origRChild) to the new node.
		//
#if defined(RECORD_MOVELOG)
		tmpf << "   *** DEATH MOVE REJECTED ***" << endl;
#endif
		TreeNode *u = tree->CreateTreeNode(UINT_MAX, origEdgelen.f, false);
		InsertSubtree(u, origPar, TreeManip::kOnRight);
		u->InvalidateAttrDown(true, origPar);	// invalidates transition probability matrix only
		u->InvalidateCondLikeArrays(); 	// invalidates conditional likelihood arrays only

		TreeNode *nd = origLChild;
		for (;;)
			{
			TreeNode *s = nd;
			assert(s != NULL);
			nd = nd->GetRightSib();
			SibToChild(u, s, TreeManip::kOnRight);
			if (s == origRChild)
				break;
			}
		tree->InvalidateNodeCounts();
		}

#if defined(RECORD_MOVELOG)
	tmpf.close();
#endif

	Clear();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ln_hastings', which is computed in both BushMaster::ProposeBirthMove and 
|	BushMaster::ProposeDeathMove.
*/
double BushMaster::GetLnHastingsRatio() const
	{
	return ln_hastings;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ln_jacobian', which is computed in both BushMaster::ProposeBirthMove and 
|	BushMaster::ProposeDeathMove.
*/
double BushMaster::GetLnJacobian() const
	{
	return ln_jacobian;
	}

void BushMaster::RefreshPolytomies()
	{
	//@POL should keep polytomies list up to date rather than building anew
	polytomies.clear();
	for (TreeNode *nd = tree->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		if (nd->IsRoot())
			continue;
		
		unsigned s = nd->CountChildren();
		if (s > 2)
			{
			polytomies.push_back(nd);
			}
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Proposes either a birth move (addition of an edge to a polytomy) or a death move (deletion of an existing edge) with
|	equal probability. If tree is currently fully-resolved, proposes only a death move. If tree is currently the star tree,
|	proposes only a birth move.
*/
void BushMaster::ProposeNewState()
	{
	// Make sure MCMCMove::SetRandomNumberGenerator has been called to set rng
	//
	assert(rng != NULL); 

	unsigned nnodes = tree->GetNNodes();
	const bool fullyResolvedBefore = (numNodesInFullyResolvedTree == nnodes);	
	const bool starTreeBefore = (nnodes == numTaxa + 1);
	birthMoveProposed = (!fullyResolvedBefore) && (starTreeBefore || rng->Uniform() < 0.5);
	
	RefreshPolytomies();
	numPolytomies = (unsigned) polytomies.size();

	const unsigned numInternalEdgesBefore = tree->GetNInternals() - 1;
	TreeNode * nd;
	if (birthMoveProposed)
		{
#		if defined(RECORD_MOVELOG)
			ofstream tmpf("movelog.txt", ios::out | ios::app);
			tmpf << "\n\nBIRTH MOVE PROPOSED" << endl;
			tmpf << "  Before move, tree has " << numInternalEdgesBefore << " internal edges and " << numPolytomies << " polytomies" << endl;
#		endif

		// Choose a polytomy at random to split
		//
		NXS_ASSERT(numPolytomies > 0);
		unsigned i = rng->SampleUInt(numPolytomies);
		nd = polytomies[i];
		ProposeBirthMove(nd);

		// Compute Jacobian by working backwards to obtain the uniform used to sample new_edgelen
		// Remember that U = 1 - e^{-theta*v} where v is the new edge length, and theta = 1/edgeLenMean
		// is the hazard parameter of the exponential distribution used for sampling new edge lengths.
		// log(1 - U) = -theta*v = -v/edgeLenMean. The Jacobian for this move is edgeLenMean/(1 - U), so
		// ln_jacobian = log(edgeLenMean) - log(1 - U) = log(edgeLenMean) + v/edgeLenMean
		//
		ln_jacobian = log(edgeLenMean) + new_edgelen/edgeLenMean;

#		if defined(RECORD_MOVELOG)
			tmpf << "  Jacobian = " << exp(ln_jacobian) << endl;
#		endif

		// Hastings ratio is probability of the reverse move divided by probability of the forward move. The birth move
		// involves the following 4 steps:
		//   1) choosing a polytomy to break (probability is inverse of the number of polytomies)
		//   2) choosing the number of spokes to move over to the new node
		//   3) selecting exactly which spokes to move
		//   4) choosing U ~ Uniform(0,1) to determine the new edge length (probability density is 1).
		// The probability of steps 2 and 3 are easier to compute jointly because of some convenient cancellation:
		// For step 2, assuming x spokes out of n total are moved, the probability is:
		//
		//              {n choose x} + {n choose n-x}
		// Pr(step 2) = ------------------------------
		//                       2^n - 2(n+1)
		//
		// If x is exactly half of n, the numerator will have only the one term (n choose x).
		// For step 3, the probability is simply 1 over the number of distinct trees in which an n-spoke
		// polytomy has been divided into a node with x+1 spokes and another with n - x + 1 spokes:
		//
		//                              2
		// Pr(step 3) = -------------------------------
		//               {n choose x} + {n choose n-x}
		//
		// Again, if x is exactly half of n, the denominator will have only the one term {n choose x}.
		// Fortunately, the hard part of both calculations cancels, and the probability of steps 2 and 3 is:
		// 
		//                                2                1
		// Pr(step 2 and step 3) = -------------- = ---------------
		//                           2^n - 2(n+1)   2^(n-1) - n - 1
		//
		// The probability of proposing a birth move is thus:
		//
		//                       1                  1            1                    1
		// Pr(birth move) = ------------- X --------------- X ------- = ---------------------------------
		//                  numPolytomies   2^(n-1) - n - 1   (1 - 0)   numPolytomies * [2^(n-1) - n - 1]
		//
		// The death move is simpler, involving only the choice of the edge that must be deleted from the tree to
		// reverse the birth move (probability is inverse of the number of edges in the post-birth-move tree). 
		//
		// Thus, the Hastings ratio for a birth move is:
		//
		//                                                      1
		//                                              -----------------
		//  Pr(death of proposed birth move)             currNumEdges + 1               currNumPolytomies * [2^(n-1) - n - 1]
		//  -------------------------------- = ------------------------------------- =  -------------------------------------
		//      Pr(proposed birth move)                         1                           numInternalEdgesBefore + 1
		//                                     -------------------------------------
		//                                     currNumPolytomies * [2^(n-1) - n - 1]
		//
		//
		// If proposed state is the fully-resolved tree, or if the reverse move would generate a fully-resolved tree, 
		// then the Hastings ratio must account for the fact that one of the moves is only attempted half the time 
		// whereas the other move is the only move possible. If the proposal generates a fully-resolved tree,
		// multiply the Hastings ratio as calculated above by 2.0. If the proposal takes us away from a fully-resolved
		// tree, then divide the Hastings ratio as calculated above by 2.0.
		//
		double nspokes = (double)polytomySize;
		ln_hastings = log((double)numPolytomies);
		ln_hastings += log(pow(2.0, nspokes - 1.0) - nspokes - 1.0);
		ln_hastings -= log((double)numInternalEdgesBefore + 1.0);
		nnodes = tree->GetNNodes();
		const bool fullyResolvedAfter = (nnodes == numNodesInFullyResolvedTree);
		if (starTreeBefore && !fullyResolvedAfter)
			ln_hastings -= log(2.0);
		else if (fullyResolvedAfter && !starTreeBefore)
			ln_hastings += log(2.0);
		//@POL the pow will be problematic for large polytomies

#		if defined(RECORD_MOVELOG)
			tmpf << "  Hastings = " << exp(ln_hastings) << endl;
			tmpf << "    numPolytomies:              " << numPolytomies << endl;
			tmpf << "    polytomySize:               " << polytomySize << endl;
			tmpf << "    numInternalEdgesBefore:       " << numInternalEdgesBefore << endl;
			tmpf << "    nnodes:                     " << nnodes << endl;
			tmpf << "    numNodesInFullyResolvedTree:" << numNodesInFullyResolvedTree << endl;
			if (starTreeBefore)
				tmpf << "    leaving star tree, Hastings ratio corrected" << endl;
			else if (fullyResolvedAfter)
				tmpf << "    going to fully-resolved tree, Hastings ratio corrected" << endl;
			else
				{
				tmpf << "    birth move not leaving star tree and not going to fully-resolved tree" << endl;
				tmpf << "      no correction to Hastings ratio made" << endl;
				}
			tmpf.close();
#		endif
		}
	else
		{
#		if defined(RECORD_MOVELOG)
			ofstream tmpf("movelog.txt", ios::out | ios::app);
			tmpf << "\n\nDEATH MOVE PROPOSED" << endl;
			tmpf << "  Before move, tree has " << numInternalEdgesBefore << " internal edges and " << numPolytomies << " polytomies" << endl;
#		endif

		// Choose an internal node at random (but not the only child of the root node)
		// and delete its edge to create a polytomy (or a bigger polytomy if there is 
		// already a polytomy)
		//
		unsigned numInternals = nnodes - numTaxa - 1;
		NXS_ASSERT(numInternals > 0);
		NXS_ASSERT(numInternals < numTaxa);
		unsigned i = rng->SampleUInt(numInternals);

		for (nd = tree->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
			{
			if (nd->IsLeafOrRoot() || nd->GetParent()->IsRoot())
				continue;

			if (i == 0)
				break;
			else
				--i;
			}
		ProposeDeathMove(nd);

		// The Jacobian for this death move is theta*exp(-theta*v), where v is the length of the edge
		// being deleted and theta = 1/edgeLenMean is the hazard parameter of the exponential distribution
		// used for sampling new edge lengths. Thus,
		//
		// ln_jacobian = log[theta*exp(-theta*v)]
		//             = log[exp(-v/edgeLenMean)/edgeLenMean]
		//             = -(v/edgeLenMean) - log(edgeLenMean)
		//
		// Note that v = origEdgelen.f; origEdgelen was set in ProposeDeathMove.
		//
		ln_jacobian = -origEdgelen.f/edgeLenMean - log(edgeLenMean);

#		if defined(RECORD_MOVELOG)
			tmpf << "  Jacobian = " << exp(ln_jacobian) << endl;
#		endif

		// Hastings ratio is inverse of that for the birth move (see extensive notes above), but be
		// careful to use number of polytomies in tree *after* death move and number of internal edges
		// before death move in the formula. Here, n is the number of spokes in the polytomy created
		// by the death move. Both polytomySize and numPolytomies are correctly computed by ProposeDeathMove.
		//
		//  Pr(birth move that reverts proposed death move)           numInternalEdgesBefore  
		//  ----------------------------------------------- = ---------------------------------------
		//           Pr(proposed death move)                   numPolytomiesAfter * [2^(n-1) - n - 1]   
		//
		// If the proposed state is the star tree, then the Hastings ratio must account for the fact that
		// the forward move is only attempted half the time whereas the reverse move is the only move possible
		// for the star tree. Thus, if the number of nodes in the tree after the proposal is 1 more than the 
		// number of taxa, the Hastings ratio as computed above must be multiplied by 2.0. If the proposal takes
		// us away from the star tree, then the Hastings ratio must be divided by 2.0.
		//
		double nspokes = (double)polytomySize;
		ln_hastings = log((double)numInternalEdgesBefore);
		ln_hastings -= log((double)numPolytomies);
		ln_hastings -= log(pow(2.0, nspokes - 1.0) - nspokes - 1.0);
		nnodes = tree->GetNNodes();
		const bool starTreeAfter = (nnodes == numTaxa + 1);
		if (fullyResolvedBefore && !starTreeAfter)
			{
			ln_hastings -= log(2.0);
			}
		else if (starTreeAfter && !fullyResolvedBefore)
			{
			ln_hastings += log(2.0);
			}
		//@POL the pow will be problematic for large polytomies

#		if defined(RECORD_MOVELOG)
			tmpf << "  Hastings = " << exp(ln_hastings) << endl;
			tmpf << "    numPolytomies:        " << numPolytomies << endl;
			tmpf << "    polytomySize:         " << polytomySize << endl;
			tmpf << "    numInternalEdgesBefore: " << numInternalEdgesBefore << endl;
			tmpf << "    nnodes:                     " << nnodes << endl;
			tmpf << "    numNodesInFullyResolvedTree:" << numNodesInFullyResolvedTree << endl;
			if (fullyResolvedBefore)
				tmpf << "    leaving fully-resolved tree, Hastings ratio corrected" << endl;
			else if (starTreeAfter)
				tmpf << "    going to star tree, Hastings ratio corrected" << endl;
			else
				{
				tmpf << "    death move not leaving fully-resolved tree and not going to star tree" << endl;
				tmpf << "      no correction to Hastings ratio made" << endl;
				}
			tmpf.close();
#		endif
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Determines distribution of x given nspokes, where x is the number spokes assigned to the newly created node in a birth
|	move. The number of ways of choosing x spokes to move out of nspokes total is {nspokes \choose x}. We are not 
|	interested in the values x = 0, x = 1, x = nspokes - 1, and x = n because these lead either to a non-move (x = 0 and 
|	x = nspokes) or to a tree that has an invalid structure (x = 1 ad x = nspokes - 1). Thus, the total number of possible
|	trees considered is {nspokes \choose 2} + {nspokes \choose 3} + ... + {n \choose nspokes - 2} = 
|	2^nspokes - 2*(nspokes + 1). The 2^nspokes comes from the fact that 2^nspokes is the sum of all binomial coefficients
|	for a sample size of nspokes. The subtracted term 2(nspokes + 1) comes from the fact that the first and last binomial
|	coefficients - {nspokes \choose 0} and {nspokes \choose nspokes} - are always 1 and the second and penultimate binomial
|	coefficients - {nspokes \choose 1} and {nspokes \choose nspokes - 1} - always equal nspokes. Thus, if one wishes to 
|	choose randomly from all possible ways of splitting the polytomy into two groups of spokes, select x with probability:
|>
|	                    (nspokes \choose x}
|	  Pr(X = x) = -----------------------------, x = 2, 3, ..., nspokes - 2
|	                2^nspokes - 2*(nspokes + 1)
|>
*/
const VecDbl &BushMaster::ComputePolytomyDistribution(unsigned nspokes)
	{
	assert(nspokes > 2);
	pair<PolytomyDistrMap::const_iterator, bool> retval;
	PolytomyDistrMap::const_iterator i = polyProb.find(nspokes);
	if (i == polyProb.end())
		{
		// There is no existing probability distribution vector corresponding to nspokes
		// Need to calcuate it and insert into polyProb map.
		//
		//cerr << "\nRecomputing polytomy distribution for " << nspokes << " spokes\n" << endl;

		double ln_nfact = rng->LotLnGamma((double)(nspokes + 1));
		double denom = exp((double)nspokes * log(2.0)) - 2.0*nspokes - 2.0; //@POL may need to factor if n large
		double ln_denom = log(denom);
		VecPolytomyDistr v;
		for (unsigned x = 2; x <= nspokes - 2; ++x)
			{
			double ln_numer = ln_nfact - rng->LotLnGamma((double)(x + 1)) - rng->LotLnGamma((double)(nspokes - x + 1));
			double prob_x = exp(ln_numer - ln_denom);
			v.push_back(prob_x);
			}
		retval = polyProb.insert(PolytomyDistrMap::value_type(nspokes, v));
		assert(retval.second == true);
		i = retval.first;
		}
	return (i->second);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Split up the polytomy at `u' by creating a new internal node v and a new edge connecting u with v. Node u is saved as
|	`origPar' and node v is saved as `origLChild' in case we need to revert the proposed move.
*/
void BushMaster::ProposeBirthMove(TreeNode *u)
	{
	assert(u != NULL);
	polytomySize = 1 + u->CountChildren();

	const vector<double> &prob_n = ComputePolytomyDistribution(polytomySize);

	// Select number of spokes to move over to new node
	//
	unsigned x;
	double p = rng->Uniform();
	double cum = 0.0;
	for (unsigned k = 0; k <= polytomySize - 4; ++k)
		{
		x = k + 2;
		double prob_k_given_n = prob_n[k];
		cum += prob_k_given_n;
		if (p < cum)
			break;
		}
	NXS_ASSERT(x < polytomySize - 1);

	// Create the new node that will receive the x randomly-chosen spokes
	//
	//pol, to get this to compile I had to change this from
	//  new_edgelen = expDist.Sample(*rng);
	//	this is currently dangerous because we are getting a bare Lot pointer from 
	//  rng, so the reference to rng stored in expDist may be a dead reference
	//	if all of the shared references to rng go away.
	//	Personally, I'm happy passing a reference to a Lot obj into Sample.
	expDist.SetLot(rng.get());  
	new_edgelen = expDist.Sample();//POL 22Jun2005
	TreeNode *v = tree->CreateTreeNode(UINT_MAX, new_edgelen, false);
	InsertSubtree(v, u, TreeManip::kOnLeft);
	v->InvalidateAttrDown(true, u);	// invalidates transition probability matrix only
	v->InvalidateCondLikeArrays();	// invalidates conditional likelihood arrays only

	// Save u and v. If revert is necessary, all of origLChild's nodes will be returned
	// to origPar, and origLChild will be deleted.
	//
	origPar = u;
	origLChild = v;

	// After the move, either v or u should have x spokes and the other node polytomySize - x spokes (u and v will 
	// each have 1 additional connector spoke).Choose x spokes randomly out of the polytomySize available. 
	// If u->par is included, let u retain the x spokes and move polytomySize - x spokes to v. Otherwise, move the
	// x spokes to v leaving polytomySize - x spokes behind.
	//
	vector<TreeNode *> uspokes;
	uspokes.push_back(u->GetParent());
	for (TreeNode *uchild = u->GetLeftChild(); uchild != NULL; uchild = uchild->GetRightSib())
		{
		if (uchild != v)
			uspokes.push_back(uchild);
		}
	assert(uspokes.size() == polytomySize);

	bool reverse_polarity = false;
	vector<TreeNode *> vspokes;
#	if defined (__MWERKS__)
		typedef vector<TreeNode *>::iterator::difference_type vec_it_diff;
#	else
		typedef unsigned vec_it_diff;
#	endif
	for (unsigned k = 0; k < x; ++k)
		{
		unsigned numUSpokes = (unsigned)uspokes.size();
		NXS_ASSERT(numUSpokes > 0);
		unsigned j = rng->SampleUInt(numUSpokes);
		TreeNode *s = uspokes[j];
		if (s == u->GetParent())
			reverse_polarity = true;
		vspokes.push_back(s);
		uspokes.erase(uspokes.begin() + (vec_it_diff) j);
		} 
	assert(uspokes.size() + vspokes.size() == polytomySize);

	if (reverse_polarity)
		{
		// transfer nodes in uspokes to v 
		//
		vector<TreeNode *>::iterator s;
		for (s = uspokes.begin(); s != uspokes.end(); ++s)
			{
			DetachSubtree(*s);
			InsertSubtree(*s, v, TreeManip::kOnRight);
			}
		}
	else
		{
		// transfer nodes in vspokes to v 
		//
		vector<TreeNode *>::iterator s;
		for (s = vspokes.begin(); s != vspokes.end(); ++s)
			{
			DetachSubtree(*s);
			InsertSubtree(*s, v, TreeManip::kOnRight);
			}
		}

	tree->InvalidateNodeCounts();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Delete the edge associated with `u' to create a polytomy (or a bigger polytomy if `u->par' was already a polytomy).
|	The supplied node u should not be the only child of the root node.
|>
|	      b       c
|	       \     /
|	        \   / 
|	         \ /
|	  a       u		   a   b   c
|	   \     /		    \  |  /
|	    \   /		     \ | /
|	     \ /		      \|/
|         v                v
|	     /			      /
|	
|	    Before           After
|>
|	Returns the number of polytomies in the tree after the proposed death move. The return value will be incorrect if 
|	the polytomies vector is not up-to-date.
*/
void BushMaster::ProposeDeathMove(TreeNode *u)
	{
	// Save nd's edge length in case we need to revert
	//
	origEdgelen = u->GetEdgeLen();

	// This operation should not leave the root node (which is a tip) with more than
	// one child, so check to make sure that the supplied node is not the root nor a
	// child of root
	//
	origPar = u->GetParent();
	assert(origPar != NULL);
	assert(!origPar->IsRoot());

	numPolytomies = (unsigned)polytomies.size();

	// Compute size of polytomy after the death move, a quantity that is needed for computing the Hastings ratio.
	// Note that one of v's children (i.e. u) is deleted but this is made up for by considering v->par, which is 
	// also a spoke that counts.
	//
	unsigned u_children = u->CountChildren();
	unsigned v_children = origPar->CountChildren();
	polytomySize = v_children + u_children; 

	bool u_polytomy_before = (u_children > 2);
	bool v_polytomy_before = (v_children > 2);
	if (u_polytomy_before && v_polytomy_before)
		{
		// No. polytomies will decrease by one as a result of this death move
		//
		--numPolytomies;
		}
	else if (!u_polytomy_before && !v_polytomy_before)
		{
		// No. polytomies will increase by one as a result of this death move
		//
		++numPolytomies;
		}

	// Make all of u's children left siblings (i.e. children of u->par)
	//
	origLChild = u->GetLeftChild();
	while (u->GetLeftChild() != NULL)
		{
		origRChild = u->GetLeftChild();
		DetachSubtree(origRChild);
		InsertSubtree(origRChild, origPar, TreeManip::kOnRight);
		}
	//origPar->InvalidateAttrDown(true, origPar->GetParent());	// no need to invalidate transition probability matrix
	origPar->InvalidateCondLikeArrays();	// invalidates conditional likelihood arrays only
	DeleteLeaf(u);

	tree->InvalidateNodeCounts();
	}

