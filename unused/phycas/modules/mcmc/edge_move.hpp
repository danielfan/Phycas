#ifndef PHO_EDGEMOVE_H
#define PHO_EDGEMOVE_H

#include "phycas/trees/tree_move.hpp"

class TreeNode;
class Tree;

/*----------------------------------------------------------------------------------------------------------------------
|	An EdgeMove changes the length of just one randomly-chosen edge in the tree. An edge chosen at random is set to the
|	value Y = m*exp(`lambda'*(u - 0.5)), where m is the original length and u is a Uniform(0,1) random deviate. Under
|	this proposal scheme, the random variable Y has the following properties:
|>
|	density         f(Y) = 1/(lambda*Y)
|	cdf             F(Y) = 0.5 + (1/lambda) log(Y/m)
|	minimum         m*exp(-lambda/2)
|	maximum         m*exp(lambda/2)
|	mean            (m/lambda)[exp(lambda/2) - exp(-lambda/2)]
|	variance        (m/lambda)^2 [(lambda/2)(exp(lambda) - exp(-lambda)) - (exp(lambda/2) - exp(-lambda/2))^2]
|>
|	With a starting edge length of 1.0, the proposed edge lengths have increasing mean and variance with increasing
|	lambda, but values of lambda in the range 0.1 to 2.0 appear to be reasonable.
|>
|	lambda       mean        s.d.
|	-----------------------------
|	 0.1      1.00042     0.02888
|	 0.5      1.01045     0.14554
|	 1.0      1.04219     0.29840
|	 2.0      1.17520     0.65752
|	10.0     14.84064    29.68297
|>
*/
class EdgeMove : public TreeMove 
	{
	public:
							EdgeMove(Tree *t);
		virtual				~EdgeMove();

		void				SetLambda(double x);
		double				GetLambda() const;

		// These are pure virtual functions in the MCMCMove base class
		//
		double				GetLnHastingsRatio() const;
		double				GetLnJacobian() const;
		void				ProposeNewState();
		void				Revert();
		void				Accept();

	private:
		double				lambda;			/* larger values result in changes of greater magnitude: CV = sqrt[(lambda/2) - 1] */
		Length				origEdgelen;	/* length of modified edge saved (in case revert is necessary) */
		TreeNode			*origNode;		/* node owning the modified edge (in case revert is necessary) */

	private:
		EdgeMove			&operator=(const EdgeMove &);	// never use - don't define
	};
typedef boost::shared_ptr<EdgeMove> EdgeMoveShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the value for the data member lambda, which is the tuning parameter for this move.
*/
inline void EdgeMove::SetLambda(double x)
	{
	lambda = x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Provides read access to the data member lambda, which is the tuning parameter for this move.
*/
inline double EdgeMove::GetLambda() const
	{
	return lambda;
	}

#endif
