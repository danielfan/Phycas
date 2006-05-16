#if ! defined(BASIC_LOT_HPP)
#define BASIC_LOT_HPP

#include <cassert>
#include <boost/shared_ptr.hpp>
//#include "phycas/rand/cdf.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	This class was called Lot because the noun lot is defined as "an object used in deciding something by chance" 
|	according to The New Merriam-Webster Dictionary.
*/
class Lot
	{
	public:
								Lot();
								Lot(unsigned);
								~Lot();

		// Accessors
		unsigned				GetSeed() const;
		unsigned 				GetInitSeed() const;

		// Modifiers
		void 					UseClockToSeed();
		void 					SetSeed(unsigned s);

		// Utilities
		unsigned 				SampleUInt(unsigned);
		unsigned				GetRandBits(unsigned nbits);

#if POLPY_NEWWAY
#	if defined(NDEBUG)
		double 					Uniform();
#	else
		double 					Uniform(const char * file = "", int line = 0);
#	endif
#else
		double 					Uniform();
#endif

	private:    	

		unsigned 				last_seed_setting;
		unsigned				curr_seed;
		//CDF *					cdf_converter;

#if POLPY_NEWWAY
		unsigned				num_seeds_generated;
#endif
	};

typedef boost::shared_ptr<Lot> LotShPtr;

} // namespace phycas

#include "pyphy/prob_dist/basic_lot.inl"

#endif

