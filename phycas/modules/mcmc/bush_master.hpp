#ifndef PHO_BUSHMASTER_H
#define PHO_BUSHMASTER_H
#include <cmath>
#include "phycas/trees/tree_move.hpp"
#include "phycas/trees/tree_counter.hpp"
#include "phycas/rand/probability_distribution.hpp" //ExponentialDistribution data member

typedef std::map<unsigned, std::vector<double> > PolytomyDistrMap;
typedef std::vector<double> VecPolytomyDistr;

class TreeNode;
class Tree;

/*----------------------------------------------------------------------------------------------------------------------
|	A BushMaster proposes tree moves that change the dimension of the tree by either adding or subtracting edges.
*/
class BushMaster : public TreeMove 
	{
	public:
								BushMaster(LotShPtr r, Tree *t);
		virtual					~BushMaster();

		// Accessors
		//
		bool					IsBirthMoveProposed();

		// These are pure virtual functions in the MCMCMove base class
		//
		double					GetLnHastingsRatio() const;
		double					GetLnJacobian() const;
		void					ProposeNewState();
		void					Revert();
		void					Accept();

	private:
		unsigned				numTaxa;						/* number of taxa */
		unsigned				numNodesInFullyResolvedTree;	/* total number of nodes in a fully resolved tree having the same number of taxa, N (= 2N - 2) */
		bool					birthMoveProposed;				/* true if last move proposed was a birth move, and false otherwise */

		double					edgeLenMean;					/* the mean of the exponential distribution from which new edge lengths are chosen in birth moves */
		ExponentialDistribution	expDist;						/* the distribution from which new edge lengths are chosen */

		PolytomyDistrMap		polyProb;						/* map of probability distributions for x where x is the number of spokes split off in a birth move that divides the spokes of a polytomy of size n into two groups. A polytomy with n spokes is split into two nodes, with x spokes going to one node and n - x remaining with the other node. polyProb[n][x] is prob(x|n). */

		std::vector<TreeNode *>		polytomies;						/* vector containing pointers to TreeNode objects representing polytomies */

		double					new_edgelen;					/* stores new edge length created by a birth move */
		unsigned				polytomySize;					/* stores number of spokes in polytomy broken by last proposed birth move */
		unsigned				numPolytomies;					/* stores number of polytomies in tree */

		Length					origEdgelen;					/* length of deleted edge saved (in case revert of death move is necessary) */
		TreeNode				*origPar;						/* parent of deleted node (in case revert of death move is necessary) */
		TreeNode				*origLChild;					/* leftmost child of deleted node (in case revert of death move is necessary) */
		TreeNode				*origRChild;					/* rightmost child of deleted node (in case revert of death move is necessary) */
		double					ln_jacobian;					/* the natural log of the jacobian for the move last proposed */
		double					ln_hastings;					/* the natural log of the hastings ratio for the move last proposed */

	private:
		void						Clear();
		void						ProposeBirthMove(TreeNode *nd);
		void						ProposeDeathMove(TreeNode *nd);
		const VecPolytomyDistr	  & ComputePolytomyDistribution(unsigned nspokes);
		void						RefreshPolytomies();
		BushMaster				  & operator=(const BushMaster &);	// never use - don't define
	};
typedef boost::shared_ptr<BushMaster> BushMasterShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of `birthMoveProposed' data member, which is true if the last move proposed was a birth move
|	and false if last move proposed was a death move.
*/
inline bool BushMaster::IsBirthMoveProposed()
	{
	return birthMoveProposed;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Re-initializes all data members that are recomputed each time ProposeNewState is called.
*/
inline void BushMaster::Clear()
	{
	polytomies.clear();

	birthMoveProposed	= false;
	origEdgelen.f		= 0.0;
	origLChild			= NULL;
	origRChild			= NULL;
	origPar				= NULL;
	ln_jacobian			= 0.0;
	ln_hastings			= 0.0;
	new_edgelen			= 0.0;
	polytomySize		= 0;
	numPolytomies		= 0;
	}

#endif
