#if ! defined(EDGE_MOVE_INL)
#define EDGE_MOVE_INL

namespace phycas
{

/*--------------------------------------------------------------------------------------------------------------------------
|	Sets `origNode' to NULL and `origEdgelen' to 0.0. These variables are used to save quantities required for reverting
|	a proposed move that was not accepted, so this function is called at the end of both accept() and revert() to reset the
|	object for the next proposal. It is also called by the constructor to initialize these variables.
*/
inline void EdgeMove::reset()
	{
	origEdgelen	= 0.0;
	origNode	= NULL;

	// one_edgelen should have one element, used for computing the edge length prior in update
	if (one_edgelen.empty())
		one_edgelen.push_back(0.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member lambda, which is the tuning parameter for this move.
*/
inline void EdgeMove::setLambda(double x)
	{
	lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member lambda, which is the tuning parameter for this move.
*/
inline double EdgeMove::getLambda() const
	{
	return lambda;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the natural log of the Hastings ratio for this move. The Hastings ratio is:
|>
|	Pr(new state -> old state)   [1/(lambda m)]
|	-------------------------- = --------------- = m* / m
|	Pr(old state -> new state)   [1/(lambda m*)]
|>
|	Here, m* is the new length of the modified edge and m is the length before the move is proposed.
*/
inline double EdgeMove::getLnHastingsRatio() const
	{
	return std::log(origNode->GetEdgeLen()) - std::log(origEdgelen);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	This move does not change the model dimension, so the Jacobian is irrelevant.
*/
inline double EdgeMove::getLnJacobian() const
	{
	return 0.0;
	}

} // namespace phycas

#endif
