#if ! defined(LARGET_SIMON_MOVE_HPP)
#define LARGET_SIMON_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
//#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "pyphy/src/mcmc_updater.hpp"		// for base class MCMCUpdater

//struct CIPRES_Matrix;

//class ProbabilityDistribution;

namespace phycas
{

//class TreeNode;

//class Tree;
//typedef boost::shared_ptr<Tree>					TreeShPtr;

//class Model;
//typedef boost::shared_ptr<Model>					ModelShPtr;

//class TreeLikelihood;
//typedef boost::shared_ptr<TreeLikelihood>			TreeLikeShPtr;

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

//class HKY;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the Larget-Simon local move. The LS move is the default move, but an exception must be made in the case
|	of a star tree, which does not have three contiguous edges required for the LS move. For star trees, the proposal 
|	changes the length of just one randomly-chosen edge in the tree. An edge chosen at random is set to the value 
|	Y = m*exp(`lambda'*(u - 0.5)), where m is the original length and u is a Uniform(0,1) random deviate. Under this 
|	proposal scheme, the random variable Y has the following properties:
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
class LargetSimonMove : public MCMCUpdater
	{
	public:
						LargetSimonMove();
						virtual ~LargetSimonMove() 
							{
							//std::cerr << "LargetSimonMove dying..." << std::endl;
							}

		//virtual void	setModel(ModelShPtr p);

		unsigned		getWhichCase() const;
		void			setLambda(double x);
		double			getLambda() const;
		bool			topologyChanged() const;
		void			defaultProposeNewState();
		void			starTreeProposeNewState();

		void			viewProposedMove(bool yes);

		// These are virtual functions in the MCMCUpdater base class
		//
		virtual void	update();
		virtual double	getLnHastingsRatio() const;
		virtual double	getLnJacobian() const;
		virtual void	proposeNewState();
		virtual void	revert();
		virtual void	accept();

	private:

		double			lambda;						/**< Factor used in modifying backbone length */

		TreeNode *		ndBase;						/**< Most ancestral node involved in the move, used as the center of the likelihood calcuation (and in revert) */
		TreeNode *		ndX;						/**< Node at one end of segment involved in move; used by Revert to undo a move */
		TreeNode *		ndY;						/**< One of two node in middle of segment involved in move; used by Revert to undo a move */
		TreeNode *		ndZ;						/**< Node at other end (from ndX) of segment involved in move; used by Revert to undo a move */
		double			origX;						/**< Original length of ndX's branch; used by Revert to undo a move */
		double			origY;						/**< Original length of ndX's branch; used by Revert to undo a move */
		double			origZ;						/**< Original length of ndX's branch; used by Revert to undo a move */
		TreeNode *		swap1;						/**< First of the two nodes involved in an NNI swap; NULL if no swap was performed; used by Revert to undo a move */
		TreeNode *		swap2;						/**< Second of the two nodes involved in an NNI swap; NULL if no swap was performed; used by Revert to undo a move */

		double			m;							/**< Original 3-segment length; needed for computing Hastings ratio */
		double			mstar;						/**< Modified 3-segment length; needed for computing Hastings ratio */

		bool			topol_changed;				/**< If true, last proposal changed topology */
		unsigned		which_case;					/**< Which of the eight possible cases was tried last */

		std::vector<double> three_edgelens;			/**< workspace declared here to avoid unnecessary allocs/deallocs */

		void			reset();					/**< Returns variables involved with reversing a proposed move to the state needed for the start of another proposal */

		bool			star_tree_proposal;			/**< True if last proposed move was on a star tree (only one randomly-chosen edge changed); False if last proposed move was not on a star tree */

		bool			view_proposed_move;			/**< If set to true, graphical tree viewer will pop up showing edges affected by the next proposed Larget-Simon move */

		// These are needed for the star tree exception
		double						orig_edge_len;	/**< Length of modified edge saved (in case revert is necessary) */
		TreeNode *					orig_node;		/**< Node owning the modified edge (in case revert is necessary) */
		std::vector<double>			one_edgelen;	/**< workspace declared here to avoid unnecessary allocs/deallocs */
	};

} // namespace phycas

#include "pyphy/src/larget_simon_move.inl"

#endif
