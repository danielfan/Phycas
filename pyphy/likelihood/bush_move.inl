#if ! defined(BUSH_MOVE_INL)
#define BUSH_MOVE_INL

namespace phycas
{

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns shared pointer copy of the `topo_prior_calculator' object. This allows an external Python program to modify
|	the object and thus change the topology prior used by this BushMove object.
*/
inline TopoPriorCalculatorShPtr BushMove::getTopoPriorCalculator()
	{
	return topo_prior_calculator;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `yes' is true, subsequent calls to BushMove::update will pop up a graphical tree viewer to show the edges 
|	affected by the move.
*/
inline void BushMove::viewProposedMove(bool yes)
	{
	view_proposed_move = yes;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Setup `edgelen_dist' data member.
*/
inline void BushMove::setEdgeLenDistMean(
  double mean)
	{
	assert(mean > 0.0);
	if (rng)
		{
		edgelen_mean = mean;
		edgelen_dist = ExponentialDistributionShPtr(new ExponentialDistribution(1.0/edgelen_mean));
		edgelen_dist->SetMeanAndVariance(edgelen_mean, edgelen_mean);
		edgelen_dist->SetLot(rng.get());  //@POL should just be able to pass in the shared_ptr rather than a raw pointer
		}
	else
		{
		throw XLikelihood("Must set the pseudorandom number generator before calling setEdgeLenDistMean");
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the move is accepted.
*/
inline void BushMove::accept()
	{
	if (add_edge_move_proposed)
		{
		// Keeping added edge, where orig_par was original polytomous node and orig_lchild was the added node
		likelihood->useAsLikelihoodRoot(orig_lchild);
		likelihood->discardCacheAwayFromNode(*orig_lchild);
		likelihood->discardCacheBothEnds(orig_lchild);

		if (view_proposed_move)
			likelihood->startTreeViewer(tree, "Add edge move ACCEPTED");
		}
	else
		{
		// Keeping edge deletion, orig_par is the new polytomous node
		likelihood->useAsLikelihoodRoot(orig_par);
		likelihood->discardCacheAwayFromNode(*orig_par);

		if (view_proposed_move)
			likelihood->startTreeViewer(tree, "Delete edge move ACCEPTED");
		}

	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of `add_edge_move_proposed' data member, which is true if the last move proposed was an
|	add-edge move and false if last move proposed was a delete-edge move.
*/
inline bool BushMove::addEdgeMoveProposed() const
	{
	return add_edge_move_proposed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Re-initializes all data members that are recomputed each time ProposeNewState is called.
*/
inline void BushMove::reset()
	{
	// Workspace used for computing edge length prior
	// Should always be length one.
	if (one_edgelen.empty())
        one_edgelen.push_back(0.0);

	polytomies.clear();

	add_edge_move_proposed	= false;
	orig_edgelen			= 0.0;
	orig_lchild				= NULL;
	orig_rchild				= NULL;
	orig_par				= NULL;
	ln_jacobian				= 0.0;
	ln_hastings				= 0.0;
	new_edgelen				= 0.0;
	polytomy_size			= 0;
	num_polytomies			= 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ln_hastings', which is the natural log of the Hastings ratio for this move and which is 
|	computed in both BushMove::ProposeAddEdgeMove and BushMove::ProposeDeleteEdgeMove.
*/
inline double BushMove::getLnHastingsRatio() const
	{
	return ln_hastings;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `ln_jacobian', which is the natural log of the Jacobian for this move and which is computed in
|	both BushMove::ProposeAddEdgeMove and BushMove::ProposeDeleteEdgeMove.
*/
inline double BushMove::getLnJacobian() const
	{
	return ln_jacobian;
	}

} // namespace phycas
#endif
