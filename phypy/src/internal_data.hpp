#if ! defined(INTERNAL_DATA_HPP)
#define INTERNAL_DATA_HPP

#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include "phypy/src/cipres/AllocateMatrix.hpp"
#include "phypy/src/cipres/ConfigDependentHeaders.h"

struct CIPRES_Matrix;

namespace CipresNative
{
class DiscreteMatrix;
}
	
namespace phycas
{
class CondLikelihood;
class CondLikelihoodStorage;

/*----------------------------------------------------------------------------------------------------------------------
|	The InternalData class stores the data structures needed for computing likelihoods on trees. This includes both
|	condition likelihood arrays and transition probability matrices. A Conditional Likelihood Array (CLA) is stored
|	internally as a std::vector<double>, but the accessor function getCLA() returns a pointer to the first element so
|	the array can be efficiently traversed during update loops. The CLAs are laid out as follows (for 3 relative rate
|	categores and DNA data):
|>
|	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---
|	|                    pattern 1                  |                    pattern 2                  |                   
|	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---
|	|     rate 1    |     rate 2    |     rate 3    |     rate 1    |     rate 2    |     rate 3    |     rate 1    |   
|	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---
|	| A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A 
|	+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---
|	                  ^
|>
|	The element above indicated by the symbol ^ would represent the probability of that part of data pattern 1 above 
|	this node given that the relative rate was rate 2 and this node had state A. The `ownedPMatrices' data member stores
|	an array of two-dimensional transition matrices similar to the following (for 3 relative rate categories, DNA data).
|	All elements of all of these transition matrices are laid out to occupy one contiguous block of memory; the index of 
|	each element in memory is indicated by the numbers inside the cells below. Note that this order is compatible with
|	the layout of the CLA to minimize cache faults as the CLAs are computed.
|>
|	            to state               to state               to state     
|	        A    C    G    T       A    C    G    T       A    C    G    T
|	f    +----+----+----+----+  +----+----+----+----+  +----+----+----+----+
|	r  A |  0 |  1 |  2 |  3 |  | 16 | 17 | 18 | 19 |  | 32 | 33 | 34 | 35 |
|	o    +----+----+----+----+  +----+----+----+----+  +----+----+----+----+
|	m  C |  4 |  5 |  6 |  7 |  | 20 | 21 | 22 | 23 |  | 36 | 37 | 38 | 39 |
|	     +----+----+----+----+  +----+----+----+----+  +----+----+----+----+
|	s  G |  8 |  9 | 10 | 11 |  | 24 | 25 | 26 | 27 |  | 40 | 41 | 42 | 43 |
|	t    +----+----+----+----+  +----+----+----+----+  +----+----+----+----+
|	a  T | 12 | 13 | 14 | 15 |  | 28 | 29 | 30 | 31 |  | 44 | 45 | 46 | 47 |
|	t    +----+----+----+----+  +----+----+----+----+  +----+----+----+----+
|	e            rate 1                 rate 2                 rate 3     
|>
*/
class InternalData
  : boost::noncopyable
	{
	friend class TreeLikelihood;

	public:

										//~InternalData()
										//	{
										//	std::cerr << "InternalData dying..." << std::endl;
										//	}

		unsigned						getCLASize() const;
		double * * *					getPMatrices();
		double * * *					getMutablePMatrices() const;
		const double * const * const *	getConstPMatrices() const;
		
		CondLikelihoodShPtr				getChildCondLikePtr();
		CondLikelihoodShPtr				getParentalCondLikePtr();
		ConstCondLikelihoodShPtr		getValidChildCondLikePtr() const;
		ConstCondLikelihoodShPtr		getValidParentalCondLikePtr() const;
		
		bool							filialCLAValid() const;
		bool							filialCLACached() const;
		bool							parentalCLAValid() const;
		bool							parentalCLACached() const;

	private:
										InternalData(unsigned nPatterns, unsigned nRates, unsigned nStates, double * * * pMatrices, bool managePMatrices, CondLikelihoodStorage & cla_storage);

		//CLA's for an edge from a node to its parent are stored in the node's InternalData (or TipData).
		//bool							parCLAValid;	/**< true if parWorkingCLA is valid */
		CondLikelihoodShPtr				parWorkingCLA; 	/**< conditional likelihood array for the parent node and beyond (valid if it points to something, invalid otherwise) */
		CondLikelihoodShPtr				parCachedCLA; 	/**< parental conditional likelihood array is stored here to make reverting MCMC moves cheap */
		//bool							childCLAValid;	/**< true if childWorkingCLA is valid. */
		CondLikelihoodShPtr				childWorkingCLA;/**< conditional likelihood array for this node and above (valid if it points to something, invalid otherwise) */ 
		CondLikelihoodShPtr				childCachedCLA; /**< filial conditional likelihood array is stored here to make reverting MCMC moves cheap */
		
		int8_t							state;			/**< Used in simulation to temporarily store the state for one character */
		mutable double * * *			pMatrices; 		/**< Either an alias to a pMatrix array or an alias to ownedPMatrices.ptr */
		ScopedThreeDMatrix<double>		ownedPMatrices; /**< Transition probability matrices for this interior node  */
		CondLikelihoodStorage &			cla_pool;		/**< CondLikelihood object storage facility */
	};
	
typedef boost::shared_ptr<InternalData>			InternalDataShPtr;
typedef std::vector<InternalDataShPtr>			VecInternalDataShPtr;

// **************************************************************************************
// ***** InternalData inlines ***********************************************************
// **************************************************************************************

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `pMatrices'.
*/
inline double * * * InternalData::getPMatrices()
	{
	return pMatrices;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `pMatrices'.
*/
inline double * * * InternalData::getMutablePMatrices() const
	{
	return pMatrices;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the data member `pMatrices'.
*/
inline const double * const * const * InternalData::getConstPMatrices() const
	{
	return pMatrices;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `childWorkingCLA' data member. If `childWorkingCLA' does not currently point to anything, a CondLikelihood 
|	object is first retrieved from `cla_pool', so this function always returns a shared pointer that actually points to
|	something.
*/
inline CondLikelihoodShPtr InternalData::getChildCondLikePtr()
	{
	if (!childWorkingCLA)
		childWorkingCLA = cla_pool.getCondLikelihood();
	return childWorkingCLA;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns `parWorkingCLA' data member. If `parWorkingCLA' does not currently point to anything, a CondLikelihood 
|	object is first retrieved from `cla_pool', so this function always returns a shared pointer that actually points to
|	something.
*/
inline CondLikelihoodShPtr InternalData::getParentalCondLikePtr()
	{
	if (!parWorkingCLA)
		parWorkingCLA = cla_pool.getCondLikelihood();
	return parWorkingCLA;
	}

inline ConstCondLikelihoodShPtr InternalData::getValidChildCondLikePtr() const
	{
	//@POL-NESCENT Mark, I don't understand this - why not just assert that childWorkingCLA actually 
	// points to a CondLikelihood object, then return childWorkingCLA directly?
	//PHYCAS_ASSERT(childCLAValid);
	//InternalData * t = const_cast<InternalData *>(this);
	//return t->getChildCondLikePtr();
	PHYCAS_ASSERT(childWorkingCLA);
	return ConstCondLikelihoodShPtr(childWorkingCLA);
	}

inline ConstCondLikelihoodShPtr InternalData::getValidParentalCondLikePtr() const
	{
	//PHYCAS_ASSERT(parCLAValid);
	//InternalData * t = const_cast<InternalData *>(this);
	//return t->getParentalCondLikePtr();
	PHYCAS_ASSERT(parWorkingCLA);
	return ConstCondLikelihoodShPtr(parWorkingCLA);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `childWorkingCLA' data member actually points to something, which means that the parental CLA is 
|	valid. Returns false if `childWorkingCLA' does not point to a CondLikelihood object.
*/
inline bool InternalData::filialCLAValid() const
	{
	return childWorkingCLA;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `childCachedCLA' data member actually points to something, which means the working CLA has been 
|	cached in case a revert of the current move is needed. Returns false if `childCachedCLA' does not point to a 
|	CondLikelihood object.
*/
inline bool InternalData::filialCLACached() const
	{
	return childCachedCLA;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `parWorkingCLA' data member actually points to something, which means that the parental CLA is 
|	valid. Returns false if `parWorkingCLA' does not point to a CondLikelihood object.
*/
inline bool InternalData::parentalCLAValid() const
	{
	return parWorkingCLA;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if `parCachedCLA' data member actually points to something, which means the working CLA has been cached
|	in case a revert of the current move is needed. Returns false if `parCachedCLA' does not point to a CondLikelihood 
|	object.
*/
inline bool InternalData::parentalCLACached() const
	{
	return parCachedCLA;
	}

} // namespace phycas

#endif