#ifndef PHO_MCMCMOVE_H
#define PHO_MCMCMOVE_H

class TreeNode;
class MoveInitializer;

#include <boost/shared_ptr.hpp>
typedef boost::shared_ptr<Lot>	LotShPtr;

class MCMCMove
	{
	protected:
		LotShPtr			rng;

	public:
							MCMCMove();
							virtual ~MCMCMove();

		// Pure virtual functions that must be overridden in derived classes
		//
		virtual	bool 		IsMetropolisHastings() const	= 0;
		virtual double		GetLnHastingsRatio() const		= 0;
		virtual double		GetLnJacobian() const			= 0;
		virtual void 		ProposeNewState()				= 0;
		virtual void		Accept()						= 0;
		virtual void		Revert()						= 0;

		void SetRandomNumberGenerator(LotShPtr r);
	};

inline MCMCMove::MCMCMove()
	{
	}

inline MCMCMove::~MCMCMove()
	{
	}

inline void MCMCMove::SetRandomNumberGenerator(LotShPtr r)
	{
	assert(r);
	rng = r;
	}

#endif
