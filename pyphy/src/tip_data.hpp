#if ! defined(TIP_DATA_HPP)
#define TIP_DATA_HPP

#include <vector>
#include <cassert>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include "pyphy/src/cipres/AllocateMatrix.hpp"
#include "pyphy/src/cipres/ConfigDependentHeaders.h"

struct CIPRES_Matrix;

namespace CipresNative
{
class DiscreteMatrix;
}
	
namespace phycas
{
class CondLikelihood;
class CondLikelihoodStorage;
typedef boost::shared_ptr<CondLikelihood> CondLikelihoodShPtr;
typedef boost::shared_ptr<const CondLikelihood> ConstCondLikelihoodShPtr;
class Tree;
typedef std::vector<unsigned int> StateListPos;

/*----------------------------------------------------------------------------------------------------------------------
|	The TipData class stores the data associated with each tip in the tree. The description below assumes the data was
|	read from the following NEXUS data file (ambig.nex):
|>
|	#nexus
|	
|	begin data;
|	  dimensions ntax=4 nchar=5;
|	  format datatype=dna missing=? gap=-;
|	  matrix
|	    taxon1 A ? T  G    {CGT}
|	    taxon2 A C - {ACG}  T 
|	    taxon3 A C -  N    (AG)
|	    taxon4 A C T  R     Y
|	  ;
|	end;
|>
|	The above NEXUS file would be stored internally (i.e., CipresNative::DiscreteMatrix) in the following form, where 
|	states or state combinations have been replaced with integers ranging from -1 to 9:
|>
|	taxon1  0 4  3 2 6
|	taxon2  0 1 -1 7 3
|	taxon3  0 1 -1 5 8
|	taxon4  0 1  3 8 9
|>
|	The global state list position vector (i.e. `CipresNative::DiscreteMatrix::getStateListPos()') is 
|	  0 2 4 6 8 14 19 23 27 30
|	Each of the above values represents an index into the global state list, which looks like this:
|	  1 0 1 1 1 2 1 3 5 -1 0 1 2 3 4 0 1 2 3 3 1 2 3 3 0 2 3 2 0 2 4 1 3
|	The table below explains this enigmatic global state list. The asterisks in the first column mark the elements of 
|	the global state list position vector given above. The second column is the global state list, which is also given 
|	above. The third column contains the codes used for states internally. Finally, the last column relates	these global 
|	state codes to the representation in the original nexus file. Horizonatal lines mark the boundaries between states. 
|	The first state list element in each section holds the number of subsequent state list elements that need to be read
|	in order to fully characterize the state or state combination. The unambiguous states come first, followed by the 
|	state code representing complete ambiguity (corresponding to the missing data symbol - usually ? - in the nexus data
|	file). After that, ambiguous states are added as they are encountered in the data file (e.g., the next global state
|	code is 5, corresponding to 'N' in the nexus file, which means "A, C, G or T"). The gap (or inapplicable) state is 
|	always represented by -1 and has no separate entry in the global state list.
|>
|             global   global     representation
|	          state    state      in nexus file
|	index     list     code  -->  ambig.nex
|	---------------------------------------------
|	  0*        1       0    -->  A
|	  1         0
|	---------------------------------------------
|	  2*        1       1    -->  C
|	  3         1
|	---------------------------------------------
|	  4*        1       2    -->  G
|	  5         2
|	---------------------------------------------
|	  6*        1       3    -->  T
|	  7         3
|	---------------------------------------------
|	  8*        5       4    -->  ?
|	  9        -1                 
|	 10         0       Note that ? allows gaps
|	 11         1       in addition to A, C, G or
|	 12         2       T, so it is even more
|	 13         3       ambiguous than N
|	---------------------------------------------
|	 14*        4       5    -->  N
|	 15         0
|	 16         1
|	 17         2
|	 18         3
|	---------------------------------------------
|	 19*        3       6    -->  {CGT}
|	 20         1
|	 21         2
|	 22         3
|	---------------------------------------------
|	 23*        3       7    -->  {ACG}
|	 24         0
|	 25         2
|	 26         3
|	---------------------------------------------
|	 27*        2       8    -->  R, {AG}
|	 28         0
|	 29         2
|	---------------------------------------------
|	 30*        4       9    -->  Y
|	 31         1
|	 32         3
|	---------------------------------------------
|>	
|	For large datasets with many different ambiguity combinations, the global state list could grow quite large. Each 
|	state code ends up being a row in the augmented transposed transition probability matrix used in computing 
|	conditional likelihood arrays, and it could be quite inefficient if these augmented matrices had many unnecessary 
|	rows. Thus, when a row in the global matrix is copied to the tip of a tree, where it is used to compute likelihoods,
|	a translation is done into local state codes. The local codes are identical to the global codes for the primary states
|	and the state repreenting complete ambiguity, but may differ from the global codes for codes representing partial
|	ambiguities. To illustrate, consider taxon4 in the example. It has these states in the NEXUS data file:
|	  A  C  T  R  Y
|	which translate to these global state codes:
|	  0  1  3  8  9
|	When copied to a tip node in a tree, however, these codes become translated to
|	  0  1  3  5  6
|	Global codes 8 and 9 have become 5 and 6, respectively.
|	
|	Note that the constructor is private and thus TipData objects can only be created by the friend function
|	allocateTipData().
*/
class TipData
  : boost::noncopyable
	{
	friend class TreeLikelihood;
	friend class SimData;

	public:
											//~TipData()
											//	{
											//	std::cerr << "TipData dying..." << std::endl;
											//	}

		const double * const * const *		getConstTransposedPMatrices() const;
		double * * *						getTransposedPMatrices();
		double * * *						getMutableTransposedPMatrices() const;
		const int8_t *						getConstStateCodes() const;

		CondLikelihoodShPtr					getParentalCondLikePtr();
		ConstCondLikelihoodShPtr			getValidParentalCondLikePtr() const;

		bool								parentalCLAValid() const;
		bool								parentalCLACached() const;

	private:

											TipData(unsigned nRates, unsigned nStates, CondLikelihoodStorage & cla_storage);
											TipData(const std::vector<unsigned int> & stateListPosVec, boost::shared_array<const int8_t> stateCodesShPtr, unsigned nRates, unsigned nStates, double * * * pMatTranspose, bool managePMatrices, CondLikelihoodStorage & cla_storage);
		const StateListPos &				getConstStateListPos() const;

		friend void							calcPMatTranspose(const TreeLikelihood & treeLikeInfo, const TipData & tipData, double edgeLength);
	
	private:
											// conditional likelihood of the rest of the tree
		//bool								parCLAValid;
		CondLikelihoodShPtr					parWorkingCLA; 	/**< conditional likelihood array for parent and beyond (valid if it points to something, invalid otherwise) */
		CondLikelihoodShPtr					parCachedCLA; 	/**< parental conditional likelihood array is stored here to make reverting MCMC moves cheap */
		

		int8_t								state;				/**< Used in simulation to temporarily store the state for one character */
		StateListPos						state_list_pos;		/**< Vector of indices into the tip-specific `state_codes' array */
		boost::shared_array<const int8_t>	state_codes; 		/**< Array of tip-specific state codes */
		mutable double * * *				pMatrixTranspose;	/**< The (rate category) x (stateCode) x (ancestor state) augmented transposed transition probability matrices (points to `ownedPMatrices' if the constructor parameter `managePMatrices' parameter is true) */
		ScopedThreeDMatrix<double>			ownedPMatrices;		/**< Vector of transposed transition matrices */
		CondLikelihoodStorage & 			cla_pool;			/**< Source of CondLikelihood objects if needed */
	};
	
typedef boost::shared_ptr<TipData> TipDataShPtr;
typedef std::vector<TipDataShPtr> VecTipDataShPtr;

// *********************************************************************************
// ***** TipData inlines ***********************************************************
// *********************************************************************************

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `parWorkingCLA' data member as a shared pointer to const CondLikelihood. If `parWorkingCLA' does not 
|	currently point to anything, a CondLikelihood object is first retrieved from `cla_pool', so this function always 
|	returns a shared pointer that actually points to something.
*/
inline ConstCondLikelihoodShPtr	TipData::getValidParentalCondLikePtr() const
	{
	//assert(parCLAValid);
	//TipData * t = const_cast<TipData *>(this);
	//return t->getParentalCondLikePtr();
	assert(parWorkingCLA);
	return ConstCondLikelihoodShPtr(parWorkingCLA);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `pMatrixTranspose' (const).
*/
inline const double * const * const * TipData::getConstTransposedPMatrices() const
	{
	return pMatrixTranspose;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `pMatrixTranspose' (non-const).
*/
inline double * * * TipData::getTransposedPMatrices()
	{
	return pMatrixTranspose;
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `pMatrixTranspose' (non-const).
*/
inline double * * * TipData::getMutableTransposedPMatrices() const
	{
	return pMatrixTranspose;
	}	

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that calls the get() function of the `state_codes' data member. The returned array contains the 
|	state codes for this particular tip.
*/
inline const int8_t * TipData::getConstStateCodes() const
	{
	return state_codes.get();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that calls returns the `state_list_pos' data member. The returned vector of unsigned values holds
|	the position of each particular coded state for this tip in the global state list. (correct?)
*/
inline const StateListPos & TipData::getConstStateListPos() const
	{
	return state_list_pos;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `parWorkingCLA' data member actually points to something, which means that the parental CLA is 
|	valid. Returns false if `parWorkingCLA' does not point to a CondLikelihood object.
*/
inline bool TipData::parentalCLAValid() const
	{
	return parWorkingCLA;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `parCachedCLA' data member actually points to something, which means the working CLA has been cached
|	in case a revert of the current move is needed. Returns false if `parCachedCLA' does not point to a CondLikelihood 
|	object.
*/
inline bool TipData::parentalCLACached() const
	{
	return parCachedCLA;
	}

} // namespace phycas

#endif
