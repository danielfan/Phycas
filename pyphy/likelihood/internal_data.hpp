#if ! defined(INTERNAL_DATA_HPP)
#define INTERNAL_DATA_HPP

#include <vector>
#include <cassert>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include "CipresCommlib/AllocateMatrix.hpp"
#include "CipresCommlib/ConfigDependentHeaders.h"

struct CIPRES_Matrix;

namespace CipresNative
{
class DiscreteMatrix;
}
	
namespace phycas
{

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
		double *						getCLA();
		const double * const			getConstCLA() const;
		double * * *					getPMatrices();
		const double * const * const *	getConstPMatrices() const;

	private:

										InternalData(unsigned nPatterns, unsigned nRates, unsigned nStates, double * * * pMatrices = NULL, bool managePMatrices = false);
	
	private:

		int8_t							state;					/**< Used in simulation to temporarily store the state for one character */
		std::vector<double>				cla;					/**< Conditional likelihood array */
		double * * *					pMatrices; 				/**< Either an alias to a pMatrix array or an alias to ownedPMatrices.ptr */
		ScopedThreeDMatrix<double>		ownedPMatrices; 		/**< Transition probability matrices for this interior node  */
	};
	
typedef boost::shared_ptr<InternalData>			InternalDataShPtr;
typedef std::vector<InternalDataShPtr>			VecInternalDataShPtr;

// **************************************************************************************
// ***** InternalData inlines ***********************************************************
// **************************************************************************************

/*----------------------------------------------------------------------------------------------------------------------
|	Returns size of the vector of conditional likelihoods `cla'. Useful for debugging to check length of `cla' before
|	playing dangerously with pointers using the getCLA or getConstCLA functions.
*/
inline unsigned InternalData::getCLASize() const
	{
	return (unsigned)cla.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns address of first element of the vector conditional likelihoods `cla'. Assumes `cla' 
|	is not empty.
*/
inline double * InternalData::getCLA()
	{
	assert(!cla.empty());
	return &cla[0]; //PELIGROSO
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns address of first element of the vector conditional likelihoods `cla'. Assumes `cla' 
|	is not empty.
*/
inline const double * const InternalData::getConstCLA() const
	{
	assert(!cla.empty());
	return &cla[0]; //PELIGROSO
	}

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
inline const double * const * const * InternalData::getConstPMatrices() const
	{
	return pMatrices;
	}

} // namespace phycas

#endif
