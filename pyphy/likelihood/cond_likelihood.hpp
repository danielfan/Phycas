#if ! defined(COND_LIKELIHOOD_HPP)
#define COND_LIKELIHOOD_HPP

#include <vector>
namespace phycas
{
typedef double LikeFltType;

class CondLikelihood
	{
	public:
		CondLikelihood(unsigned len)
			:underflowExpon(0),
			numEdgesSinceUnderflowProtection(UINT_MAX),
			claVec(len)
			{
			cla = &claVec[0];
			}

		/* enum CLADirection {CLA_DIR_UP, CLA_DIR_DOWN}; 
		
		CLADirection getDirection() const
			{
			return direction;
			}
		*/
		LikeFltType *	getCLA()
			{
			return cla;
			}
		LikeFltType *	getCLA() const
			{
			return cla;
			}
		unsigned	getCLASize() const
			{
			return claVec.size();
			}
	private:
			//CLADirection 	direction;
			LikeFltType *	cla;
			unsigned 		underflowExpon;
			unsigned 		numEdgesSinceUnderflowProtection;
			std::vector<LikeFltType> 	claVec;
	};

} //namespace phycas
#endif
