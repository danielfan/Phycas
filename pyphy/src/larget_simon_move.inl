#if ! defined(LARGET_SIMON_MOVE_INL)
#define LARGET_SIMON_MOVE_INL

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	The default constructor sets `lambda' to the default value (0.2), sets `topol_changed' to false, and `m' and `mstar'
|	to 0.0. All other data members are automatically initialized (shared pointers) or are initialized via a call to 
|	reset().
*/
inline LargetSimonMove::LargetSimonMove() : MCMCUpdater()
	{
	topol_changed	= false;
	lambda			= 0.2;
	m				= 0.0;
	mstar			= 0.0;
	three_edgelens.reserve(3);
	view_proposed_move = false;
	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `yes' is true, subsequent calls to LargetSimonMove::update will pop up a graphical tree viewer to show the edges
|	affected by the move.
*/
inline void LargetSimonMove::viewProposedMove(bool yes)
	{
	view_proposed_move = yes;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Forgets information saved to enable reverting a proposed move.
*/
inline void LargetSimonMove::reset()
	{
	swap1		= NULL;
	swap2		= NULL;
	ndX			= NULL;
	ndY			= NULL;
	ndZ			= NULL;
	ndBase		= NULL;
	origX		= 0.0;
	origY		= 0.0;
	origZ		= 0.0;

	// three_edgelens should have 3 elements, used for computing the edge length prior for default proposal in update
	if (three_edgelens.size() != 3)
		three_edgelens.resize(3, 0.0);

	// these are related to the star tree exception
	orig_edge_len	= 0.0;
	orig_node		= NULL;

	// one_edgelen should have 1 element, used for computing the edge length prior for star tree proposal in update
	if (one_edgelen.size() != 1)
		one_edgelen.resize(1, 0.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides a way of setting the value for the private data member `lambda'.
*/
inline void LargetSimonMove::setLambda(double x)
	{
	lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the private data member `lambda'.
*/
inline double LargetSimonMove::getLambda() const
	{
	return lambda;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move. The Hastings ratio is (`mstar'/`m')^3, where `mstar' is
|	the new length of the modified three-edge segment and `m' is the length before the move is proposed. If the tree is
|	the star tree (only one internal node), then the exponent is 1 rather than 3 because only one edge is involved.
*/
inline double LargetSimonMove::getLnHastingsRatio() const
	{
	if (star_tree_proposal)
		return std::log(orig_node->GetEdgeLen()) - std::log(orig_edge_len);
	else
		return 3.0*(std::log(mstar) - std::log(m));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Jacobian for this move.
*/
inline double LargetSimonMove::getLnJacobian() const
	{
	return 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the current value of private variable `which_case', which is the index of whichever of the eight possible 
|	cases was used last.
*/
inline unsigned LargetSimonMove::getWhichCase() const
	{
	return which_case;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the private data member `topol_changed'.
*/
inline bool	LargetSimonMove::topologyChanged() const
	{
	return topol_changed;
	}

} // namespace phycas
#endif
