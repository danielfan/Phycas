/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __NXS_UTILITIES_
#define __NXS_UTILITIES_
#include "phycas/src/ncl/nxs_defs.hpp"
#include <vector>

#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::strchr;
#endif

/*--------------------------------------------------------------------------------------------------------------------------
| Compares two patterns as the < operator would.   Each taxon's bitCode is evaluated by using the numerical < operator. 
| If more than one bit Code Unit is required for each taxon the lower the bit Code  Unit indices are less important:  
| So if each taxon needs two bitCodeUnits and  
| a[0] > b[0] but  a[1] < b[1]
| Then a is less than b.  
|
| As expected, the first taxon is the most important in sorting.
*/
template<typename T> 
bool BitCodeLessThan(
  const T     * f, 		/* "left hand" or first pattern	*/
  const T     * s, 		/* "right hand" or second pattern	*/
  unsigned		nWords,			/* the number of words in each array (ntaxa for data patterns)	*/
  unsigned char wordLen)		/* the number of bitCodeUnits needed to encode the states for one taxon	*/
	{
	if (wordLen == 1)
		{	//faster version for types with few states.
		for (unsigned taxInd = 0 ; taxInd< nWords ; ++taxInd)
			{
			if (*f < *s)
				return true;
			if (*f > *s)
				return false;
			++f;
			++s;
			}
		}
	else
		{
		// mth bug fixed here 3-18-01
		
			//	Need to compare the last bitCodeUnits first 
			// (ie having state 8 > having state 1, but in our encoding they are stored 0x0 0x1 and 0x1 0x0, respectively)
			//	Which forces an odd stride through the arrays  The alternative(less intuitive) way
			// 	would be sorting with the last taxon's states having precedence over the first
			//	in which case we could hop to the end and straight iterate back through 
		const unsigned twoWords = 2U * wordLen;
			// jump to the last BitCode of the first taxon
		f += wordLen - 1;
		s += wordLen - 1;
		for (unsigned taxInd = 0 ; taxInd< nWords ; ++taxInd)
			{	//jump to the last bitcodeunit of the next taxon
			for (unsigned char currBUnit = 0 ; currBUnit < wordLen ; ++currBUnit)
				{
				if (*f < *s)
					return true;
				if (*f > *s)
					return false;
				--f;
				--s;
				}
			//	now jump to the last bitCode of the next taxon
			//
			f += twoWords; 
			s += twoWords;
			}
		}
	return false;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Compares two patterns as the > operator would.   Each taxon's bitCode is evaluated by using the numerical > operator. 
| If more than one bit Code Unit is required for each taxon the lower the bit Code  Unit indices are less important:  
| So if each taxon needs two bitCodeUnits and  
| a[0] > b[0] but  a[1] < b[1]
| Then a is less than b.  
|
| As expected, the first taxon is the most important in sorting.
*/
template <typename T>
bool BitCodeGreaterThan(
  const T *f, 			/* "left hand" or first pattern	*/
  const T *s, 			/* "right hand" or second pattern	*/
  unsigned nWords,			/* the number of words in each array (ntaxa for data patterns)	*/
  unsigned char	wordLen)		/* the number of bitCodeUnits needed to encode the states for one taxon	*/
	{
	if (wordLen==1)
		{	//faster version for types with few states.
		for (unsigned i = 0 ; i< nWords ; ++i)
			{
			if (f[i] > s[i])
				return true;
			if (f[i] < s[i])
				return false;
			}
		}
	else
		{
		// mth bug fixed here 3-18-01
		
			//	Need to compare the last bitCodeUnits first 
			// (ie having state 8 > having state 1, but in our encoding they are stored 0x0 0x1 and 0x1 0x0, respectively)
			//	Which forces an odd stride through the arrays  The alternative(less intuitive) way
			// 	would be sorting with the last taxon's states having precedence over the first
			//	in which case we could hop to the end and straight iterate back through 
		const unsigned twoWords = 2 * wordLen;
			// jump to the last BitCode of the first taxon
		f += wordLen - 1;
		s += wordLen - 1;
		for (unsigned taxInd = 0 ; taxInd< nWords ; ++taxInd)
			{	//jump to the last bitcodeunit of the next taxon
			for (unsigned char currBUnit = 0 ; currBUnit < wordLen ; ++currBUnit)
				{
				if (*f > *s)
					return true;
				if (*f < *s)
					return false;
				--f;
				--s;
				}
			//	now jump to the last bitCode of the next taxon
			//
			f += twoWords; 
			s += twoWords;
			}
		}
	return false;
	}




template<class Key,class Val>std::vector<Key> GetKeys(
  const std::map<Key,Val> &m)
	{
  	std::vector<Key> retVec;
  	typedef typename std::map<Key,Val>::const_iterator MIterator;
  	MIterator mIt = m.begin();
	for (; mIt != m.end(); ++mIt)
		retVec.push_back(mIt->first);
	return retVec;
	}
	
/*-------------------------------------------------------------------------------------------------------------------------- 
|	returns 1 if the leftOp is > rightOp, -1 if leftOp < rightOp, and 0 if the arrays are equal
|	The last elements in the array [dim -1] have highest priority in the comparison
|
*/
template<class T> int ArrayCompare(const T * leftOp, const T *rightOp, unsigned dim)
	{
#	if defined(C_FUNCS_IN_STD_NAMESPACE)
		using std::memcmp;
#	endif
	return memcmp(leftOp, rightOp, dim*sizeof(T));
	}

/*-------------------------------------------------------------------------------------------------------------------------- 
|	returns 1 if the leftOp is > rightOp, -1 if leftOp < rightOp, and 0 if the arrays are equal
|	The last elements in the array [dim -1] have highest priority in the comparison
|
*/
template<class T> int ArrayCompareReversePrecedence(const T * leftOp, const T *rightOp, unsigned dim)
	{
	NXS_ASSERT(dim>0);
	unsigned i = dim-1;
	for (;;)
		{
		if (leftOp[i] < rightOp[i])
			return -1;
		if (rightOp[i] < leftOp[i])
			return 1;
		if (i == 0)
			return 0;
		--i;
		}
	}

/*-------------------------------------------------------------------------------------------------------------------------- 
|	returns 1 if the leftOp is > rightOp, -1 if leftOp < rightOp, and 0 if the arrays are equal
|	The last elements in the array [dim -1] have highest priority in the comparison
|
*/
template<class T> int ArrayLessThanReversePrecedence(const T * leftOp, const T *rightOp, unsigned dim)
	{
	NXS_ASSERT(dim>0);
	unsigned i = dim-1;
	for (;;)
		{
		if (leftOp[i] < rightOp[i])
			return -1;
		if (rightOp[i] < leftOp[i])
			return 1;
		if (i == 0)
			return 0;
		--i;
		}
	}


template<class InputIterator> inline void DeleteRange(InputIterator first, InputIterator last)
	{
	for (; first != last; ++first)
		delete (*first);
	}

template<class T> void DeletePtr(
  boost::shared_ptr<T> )	
	{
	}

template<class T> void DeletePtr(
  T *ptr)	/* reference to the vector to be cleared */
	{
	delete ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	walks through a vector of pointers and deletes each one then clears the vector.  This is necessary whenever a 
|	vector is being used to hold objecet allocated with new.
|	NOTE THIS VERSION only works with single object pointers NOT ARRAYS!!!
*/
template<class T> void DeleteAndClearVecOfPtrs(
  std::vector<T*> &vec)	/* reference to the vector to be cleared */
	{
	DeleteRange(vec.begin(),vec.end());
	vec.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	We don't need to delete shared pointers (in fact we shouldn't) just clear.
|	This form of the function lets us change between using dumb pointers to shared_ptrs without leaking memory (when 
|	using dumb pointers) or double deleting (when using smart pointers)
*/
template<class T> void DeleteAndClearVecOfPtrs(
  std::vector< boost::shared_ptr<T> > &vec)	/* reference to the vector to be cleared */
	{
	vec.clear();
	}

template <class T> unsigned GetMaxSize(
  T begIt,
  T endIt)
  	{
  	unsigned maxS = 0;
  	for (; begIt != endIt; ++begIt)
  		{
  		const unsigned tempS = (unsigned)begIt->size();
  		if (tempS > maxS)
  			maxS = tempS;
  		}
  	return maxS;
  	}

#if 0
/*----------------------------------------------------------------------------------------------------------------------
|	walks through a vector of pointers and deletes each one then clears the vector.  This is necessary whenever a 
|	vector is being used to hold objecet allocated with new.
|	NOTE THIS VERSION only works with ARRAYS THE delete [] operator is called NOT just delete!!!
*/
template<class T> inline void DeleteAndClearVecOfArrPtrs(
  std::vector<T *> &vec)	/* reference to the vector to be cleared */
	{
	class std::vector<T *>::iterator vIt;
	vIt = vec.begin();
	for (; vIt != vec.end() ; ++vIt)
		delete [] *vIt;
	vec.clear();
	}
/*--------------------------------------------------------------------------------------------------------------------------
| Allocates a three dimensional array of T as one contiguous block of memory
| the dimensions are f two dimensional arrays that are s by t.  
| the array is set up so that 
| for(i = 0 ; i < f ; ++i)
|	for (j = 0 ; j < s ; ++j)
|		for (k = 0 ; k < t; ++k)
|			{
|			array[i][j][k];
|			}
| would be the same order of access as: 
| 
|	double *temp = **array;
|	for (i = 0 ; i < f*s*t ; ++i)
|		{
|		*temp++;
|		}
*/
template <class T> T ***NewThreeDArray(
  unsigned f, 
  unsigned s, 
  unsigned t)
	{
	NXS_ASSERT(f > 0 && s > 0 && t> 0);
	T ***temp;
	temp = new T **[f];
	temp[0] = new T *[f * s];
	temp[0][0] = new T[f * s * t];
	for (unsigned sIt = 1 ; sIt < s ; ++sIt)
		temp[0][sIt] = temp[0][sIt-1] + t ;
	for (unsigned fIt = 1 ; fIt < f ; fIt ++)
		{
		temp[fIt] = temp[fIt -1] +  s ;
		temp[fIt][0] = temp[fIt -1][0] + (s*t);
		for (int sIt = 1 ; sIt < s ; ++mIt)
			temp[fIt][sIt] = temp[fIt][sIt-1] + t ;
		}
	return temp;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Delete a Thwo Dimensional Array that has been allocated using NewThreeDArray().
|	Note that this function takes a *** reference so that the pointer in the calling function can be set to NULL
*/
template <class T> void DeleteThreeDArray(
  T ***&temp)
	{
	if (temp != NULL)
		{
		if (*temp != NULL)
			{
			delete [] **temp;
			delete [] * temp;
			}
		delete [] temp;
		temp = NULL;
		}
	}


template<class Key,class Val>std::set<Key> MapToSet(
  const std::map<Key,Val> &m)
	{
  	std::set<Key> retSet;
  	class std::map<Key,Val>::const_iterator mIt = m.begin();
	for (; mIt != m.end(); ++mIt)
		retSet.insert(mIt->first);
	return retSet;
	}

/*-------------------------------------------------------------------------------------------------------------------------- 
|	returns leftOp == rightOp by comparing the first dim elements using class T's == operator
*/
template<class T> bool ArrayEq(const T * leftOp, const T *rightOp, unsigned dim)
	{
	for (unsigned i = 0; i < dim; ++i)
		{
		if (leftOp[i] != rightOp[i] )
			return false;
		}
	return true;
	}
	
/*-------------------------------------------------------------------------------------------------------------------------- 
|	returns 1 if the leftOp is > rightOp, -1 if leftOp < rightOp, and 0 if the arrays are equal
|	The last elements in the array [dim -1] have highest priority in the comparison
|
*/
template<class T> int ArrayCompReversePrecedence(const T * leftOp, const T *rightOp, unsigned dim)
	{
	NXS_ASSERT(dim>0);
	for (unsigned i = dim-1; i >= 0; --i)
		{
		if (leftOp[i] < rightOp[i])
			return -1;
		if (rightOp[i] < leftOp[i])
			return 1;
		}
	return 0;
	}
/*-------------------------------------------------------------------------------------------------------------------------- 
|	returns 1 if the leftOp is > rightOp, -1 if leftOp < rightOp, and 0 if the arrays are equal
|	The first elements in the array [0] have highest priority in the comparison
|
*/
template<class T> int ArrayComp(const T * leftOp, const T *rightOp, unsigned dim)
	{
	NXS_ASSERT(dim>0);
	for (unsigned i = 0; i < dim; ++i)
		{
		if (leftOp[i] < rightOp[i])
			return -1;
		if (rightOp[i] < leftOp[i])
			return 1;
		}
	return 0;
	}
#endif
#endif

