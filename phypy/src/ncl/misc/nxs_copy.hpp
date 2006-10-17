/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

#if !defined NXS_COPY_H
#define NXS_COPY_H
#include <algorithm>
#include <cstring>
#include "phypy/src/ncl/misc/generic_type_mapping.hpp"

#if 0
template<typename T> struct SupportsBitwiseCopy { enum {kResult = false};	};
template<typename T> struct SupportsBitwiseCopy<T*> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<short int> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<int> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<char> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<long int> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<double> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<unsigned short int> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<unsigned int> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<unsigned char> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<unsigned long int> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<bool> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<wchar_t> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<float> {	enum {kResult = true}; 	};
template<> struct SupportsBitwiseCopy<long double> {	enum {kResult = true}; 	};


//	This file uses tricks discussed in Andrei Alexandrescu's book to
//	implement a call to memcpy for primitive types or any class which
//  has the statement template<> struct SupportsBitwiseCopy<CLASS> {	enum {kResult = true}; 	};
//	because of potential portability issues with TypeList, primitive types are
//	have SupportsBitwiseCopy specialized here by brute force enumeration
	

class NullType {};

template <typename T>
class TypeTraits
	{
	private :
		template <class U> struct PointerTraits
			{
				enum {kResult = false};
				enum {kCopyWithMemCopy = false};	// only allowing memcpy on bare pointers
				enum {kSizeOfPointee = 0};	// only allowing memcpy on bare pointers
			};
		template <class U> struct PointerTraits<U*>
			{
				enum {kResult = true};
				enum {kCopyWithMemCopy = SupportsBitwiseCopy<U>::kResult};	
				enum {kSizeOfPointee = sizeof(U)};	
			};
		template <class U> struct PointerTraits<const U*>
			{
				enum {kResult = true};
				enum {kCopyWithMemCopy = SupportsBitwiseCopy<U>::kResult};	
				enum {kSizeOfPointee = sizeof(U)};	
			};
	public:
		enum {kIsPointer = PointerTraits<T>::kResult};
		enum {kCanUseMemCpyOnPointee = PointerTraits<T>::kCopyWithMemCopy};
		enum {kPointeeSize = PointerTraits<T>::kSizeOfPointee}; //only valid if kIsPointer !!
	//	typedef PointerTraits<T>::PointeeType  PointeeType;
	};
	
template<class T, class U> 
class Conversion
	{
	public:
		enum {kSameType = false};
	};

template<class T> 
class Conversion<T,T>
	{
	public:
		enum {kSameType = true};
	};

template<class T> 
class Conversion<const T*,T*>
	{
	public:
		enum {kSameType = true};
	};

template<class T> 
class Conversion<T*, const T*>
	{
	public:
		enum {kSameType = true};
	};

enum CopyAlgoSeclector 
	{
		kConservative, 
		kFast
	};

template <typename InIt, typename OutIt>
inline OutIt Copy_n_Impl(InIt first, OutIt resultP, unsigned n, Int2Type<kConservative>)
	{
	for (unsigned i = 0; i < n; ++i)
		{
		*resultP = *first;
		++first;
		++resultP;
		}
	return resultP; 
	}

template <typename InIt, typename OutIt>
inline OutIt Copy_n_Impl(InIt first, OutIt resultP, unsigned n, Int2Type<kFast>)
	{
#	if !defined (NDEBUG)
		InIt temp = resultP; //this statement is used to generate a warning if InIt and OutIt, don't point to the sametype
#	endif
	return (OutIt) std::memcpy(resultP, first,  n * sizeof(*first));
	}
	
template <typename InIt, typename OutIt>
inline OutIt CopyImpl(InIt first, InIt last, OutIt resultP, Int2Type<kConservative>)
	{
	return std::copy(first, last, resultP);
	}
	
template <typename InIt, typename OutIt>
inline OutIt CopyImpl(InIt first, InIt last, OutIt resultP, Int2Type<kFast>)
	{
#	if !defined (NDEBUG)
		InIt temp = resultP; //this statement is used to generate a warning if InIt and OutIt, don't point to the sametype
#	endif
	return (OutIt) std::memcpy(resultP, first,  ((std::size_t) (last - first)) * sizeof(*first));
	}
	
//	In Phorest and NCL we are explicitly casting safe casts and maximizing compiler cast warnings.
//	because of this the 
//	
template <typename InIt, typename OutIt>
OutIt nxs_copy(InIt first, InIt last, OutIt resultP)
	{
		enum { kUseMemCpy =(TypeTraits<InIt>::kIsPointer && 
							TypeTraits<OutIt>::kIsPointer &&
							TypeTraits<InIt>::kCanUseMemCpyOnPointee &&
							TypeTraits<OutIt>::kCanUseMemCpyOnPointee &&
							TypeTraits<InIt>::kPointeeSize == TypeTraits<OutIt>::kPointeeSize) ? kFast : kConservative};
		return CopyImpl(first, last, resultP, Int2Type<kUseMemCpy>());
	}

//	In Phorest and NCL we are explicitly casting safe casts and maximizing compiler cast warnings.
//	because of this the 
//	
template <typename InIt, typename OutIt>
OutIt nxs_copy_n(InIt first, OutIt resultP, unsigned n)
	{
		enum { kUseMemCpy =(TypeTraits<InIt>::kIsPointer && 
							TypeTraits<OutIt>::kIsPointer &&
							TypeTraits<InIt>::kCanUseMemCpyOnPointee &&
							TypeTraits<OutIt>::kCanUseMemCpyOnPointee &&
							TypeTraits<InIt>::kPointeeSize == TypeTraits<OutIt>::kPointeeSize) ? kFast : kConservative};
		return Copy_n_Impl(first, resultP, n, Int2Type<kUseMemCpy>());
	}

#else
	
template <typename InIt, typename OutIt>
inline OutIt nxs_copy(InIt first, InIt last, OutIt resultP)
	{
	return std::copy(first, last, resultP);
	}

template <typename InIt, typename OutIt>
OutIt nxs_copy_n(InIt first, OutIt resultP, unsigned n)
	{
	for (unsigned i = 0; i < n; ++i)
		*resultP++ = *first++;
	return resultP;
	}
#endif

#endif

