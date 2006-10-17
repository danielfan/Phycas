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

#if !defined NXS_STL_ALGORITHM_EXTENSIONS
#define NXS_STL_ALGORITHM_EXTENSIONS

#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <string>
// This file contains simple extensions to stl algorithms that 
//	are listed in boost's unofficial wiki site's STLAlgorithmExtensions site
//	(or are very similar to the algorithm's suggested there.

#include <list>
#include "phypy/src/ncl/misc/nxs_copy.hpp"


/*----------------------------------------------------------------------------------------------------------------------
|	copy_if
|	from Meyers, Effective STL p 156
*/
template <	typename InputIterator,
			typename OutputIterator,
			typename Predicate>
OutputIterator copy_if(InputIterator beginIt, InputIterator endIt, OutputIterator destIt, Predicate p)
	{
	while (beginIt != endIt)
		{
		if (p(*beginIt))
			{
			*destIt = *beginIt;
			++destIt;
			}
		++beginIt;
		}
	return destIt;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	copy_if_ptr is like copy_if, but dereferences ptr before calling the predicate
*/
template <	typename InputIterator,
			typename OutputIterator,
			typename Predicate>
OutputIterator copy_if_ptr(InputIterator beginIt, InputIterator endIt, OutputIterator destIt, Predicate p)
	{
	while (beginIt != endIt)
		{
		if (p(**beginIt))
			{
			*destIt = *beginIt;
			++destIt;
			}
		++beginIt;
		}
	return destIt;
	}


template <typename T> 
void AppendUnique(std::vector<T> *dest, const std::vector<T> toAppend)
	{
	typedef typename std::vector<T>::const_iterator VecIt;
	for (VecIt newIt = toAppend.begin(); newIt != toAppend.end(); ++newIt)
		{
		VecIt dIt = find(dest->begin(), dest->end(), *newIt);
		if (dIt == dest->end())
			dest->push_back(*newIt);
		}
	}
	

/*----------------------------------------------------------------------------------------------------------------------
|	Simple structure that stores the value of a variable at creation and restores it on destruction.
|	Makes it easier to write exception safe code in cases when a variable needs to be returned to it original value 
|	before leaving a function.
*/
template <typename T> 
class RestoreOnExit
	{
		T	&ref;
		T	orig;
	public:
		RestoreOnExit(T &r) : ref(r), orig(r){}
		~RestoreOnExit() 
			{
			ref = orig;
			}
	};
	
/*----------------------------------------------------------------------------------------------------------------------
|	Similar to the Perl split function.  Copies the original container of T objects into the list outList as 
|	smaller containers broken at the specified value (which is omitted).
|	e.g. if the incoming string "usr/home/mine" and the splitAtVal was the char '/', then the resulting list would
|	have three strings "usr" "home" and "mine"
*/
template <typename T, class ORIG_CONTAINER>
void split(const ORIG_CONTAINER &origContainer, T splitAtVal, std::list<ORIG_CONTAINER> *outList)
	{
	typedef typename ORIG_CONTAINER::const_iterator OrigCIt;
	OrigCIt begIt = origContainer.begin();
	const OrigCIt endIt = origContainer.end();
	if (begIt == endIt)
		return;	//empty container sent
	OrigCIt copyToIt;
	do	{
		copyToIt = find(begIt, endIt, splitAtVal);
		if (begIt == copyToIt)	
			outList->push_back(ORIG_CONTAINER());
		else
			{
			outList->push_back(ORIG_CONTAINER(begIt,copyToIt));
			begIt = copyToIt;
			}
		++begIt;
		}
	while (copyToIt != endIt);
	}

template <typename InputIterator, typename Predicate>
inline bool all(InputIterator first, InputIterator last, Predicate p)
	{
	for (; first != last; ++first)
		{
		if (!p(*first))
			return false;
		}
	return true;
	} 
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the sum of squares of all elements in the range.
*/
template <typename InputIterator, typename T>
T sum_sq(InputIterator start, const InputIterator finish, T initialV)
	{
	for (;start != finish; ++start)
		{
		T x = (*start)*(*start);  // Comment from similar operation in praxis code DLS 22jun00: Use of the temporary variable here seems to be important for numerical accuracy (not sure why, but different code is generated in both PPC Metrowerks and Dec Alpha compilers)
		initialV += x;
		}
	return initialV;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the sum of squares of all elements in the range.
*/
template <typename InputIterator, typename T>
void scale(InputIterator start, const InputIterator finish, const T multiplier)
	{
	for (;start != finish; ++start)
		*start *= multiplier;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Swaps the references to min_v and max_v if min_v < max_v
*/
template <typename T>
inline void pho_sort(T &min_v, T &max_v)
	{
	if (max_v < min_v)
		std::swap<T>(min_v, max_v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Swaps the references to min_v and max_v if min_v < max_v
*/
template <typename T>
inline void sum_sq(T &min_v, T &max_v)
	{
	if (max_v < min_v)
		std::swap<T>(min_v, max_v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Swaps the references to min_v and max_v if min_v < max_v
*/
template <typename T>
inline std::pair<T, T> &swap_pair_order(std::pair<T, T> &p)
	{
	std::swap<T>(p.first, p.second);
	return p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	return Sum of fir[i]*sec[i] for i < n
*/
template <typename T>
inline T sum_of_array_product_n(const T *fir, const T *sec, unsigned n)
	{
	T temp = T(0);
	for (unsigned i = 0; i< n; ++i)
		temp += fir[i]*sec[i];
	return temp; 
	}

template <typename InputIterator, typename Predicate>
inline bool any(InputIterator first, InputIterator last, Predicate p)
	{
	for (; first != last; ++first)
		{
		if (p(*first))
			return true;
		}
	return false;
	} 
	
template <typename OutputIterator, typename T>
inline void replace_all(OutputIterator first, OutputIterator last, const T &new_val)
	{
	for (; first != last; ++first) 
		*first = new_val;
	} 

//deals with const pointers
template <typename InputIterator>
inline bool equals_n(InputIterator l, InputIterator r, unsigned nElements)
	{
	for (unsigned i = 0; i < nElements; ++i) 
		if (!(*l == *r))
			return false;
	return true;
	} 

template<typename InputIterator, typename T>
inline bool exact_count(InputIterator first, InputIterator last, T const&v, unsigned count = 1)
        {
        first = find(first, last, v);
        if (first == last) 
        	return (count == 0); 
        unsigned nFound = 1;
        if (nFound <= count)
        	{
        	first = find(++first,last, v);
        	if (first == last)
        		return (nFound == count);
        	else
        		++nFound;
        	}
        return false;	
        }

template<typename InputIterator, typename Predicate>
inline bool exact_count_if(InputIterator first, InputIterator last, Predicate pred, int count = 1)
        {
        first = find_if(first, last, pred);
        if (first == last) 
        	return (count == 0); 
        unsigned nFound = 1;
        if (nFound <= count)
        	{
        	first = find_if(++first,last, pred);
        	if (first == last)
        		return (nFound == count);
        	else
        		++nFound;
        	}
        return false;	
        }
        
template<typename InputIterator, typename T>
inline bool exists_and_only(InputIterator first, InputIterator last, T const&v)
	{
	first = find(first, last, v);
	return (first != last && find(++first,last, v) == last);
	}

template<typename InputIterator, typename Predicate>
inline bool exists_and_only_if(InputIterator first, InputIterator last, Predicate pred)
	{
	first = find_if(first, last, pred);
	return (first != last && find_if(++first,last, pred) == last);
	}

//Performs f() on any element elem in [first, last) where pred(elem) is true
template<typename InputIterator, typename Predicate, typename UnaryFunction>
 UnaryFunction for_each_if(
   InputIterator first, 
   InputIterator last,
   Predicate pred, 
   UnaryFunction f)
	{
	for(;first != last; ++first)
		{
		if (pred(*first))
			f(*first);
		}
	return f;
 	}
namespace ncl
{

class XEmptyContainer: public std::invalid_argument 
	{
	public:
		XEmptyContainer(const std::string & msg):invalid_argument(msg){}
	};
#if 1
		/// Calculates the mean in double precision.  
	template<typename T>
	double calc_mean_dbl_unchecked(const std::vector<T> & container)
		{
		return (double) std::accumulate(container.begin(), container.end(), T(0)) /(double) container.size();
		}
		/// Calculates the mean in double precision.  
	template<typename T>
	double calc_mean_dbl(const std::vector<T> & container)
		{
		if (container.empty())
			throw XEmptyContainer("empty container in calc_mean_dbl");
		return calc_mean_dbl_unchecked<T>(container);
		}

	template<typename T>
	T sort_to_get_median_unchecked(std::vector<T> & container) /// should be specialized for sorted containers
		{
		std::sort(container.begin(), container.end());
		return container[(unsigned) container.size()/2];
		}

	template<typename T>
	T sort_to_get_median(std::vector<T> & container) /// should be specialized for sorted containers
		{
		if (container.empty())
			throw XEmptyContainer("empty container in sort_to_get_median");
		return sort_to_get_median_unchecked<T>(container);
		}
#else 
		/// Calculates the mean in double precision.  
	template<template <typename> class C, typename T>
	double calc_mean_dbl_unchecked(const C<T> & container)
		{
		return (double) std::accumulate(container.begin(), container.end(), T(0)) /(double) container.size();
		}
		/// Calculates the mean in double precision.  
	template<template <typename> class C, typename T>
	double calc_mean_dbl(const C<T> & container)
		{
		if (container.empty())
			throw XEmptyContainer("empty container in calc_mean_dbl");
		return calc_mean_dbl_unchecked<C, T>(container);
		}

	template<template <typename> class C, typename T>
	T sort_to_get_median_unchecked(C<T> & container) /// should be specialized for sorted containers
		{
		std::sort(container.begin(), container.end());
		return container[(unsigned) container.size()/2];
		}

	template<template <typename> class C, typename T>
	T sort_to_get_median(C<T> & container) /// should be specialized for sorted containers
		{
		if (container.empty())
			throw XEmptyContainer("empty container in sort_to_get_median");
		return sort_to_get_median_unchecked<C, T>(container);
		}
#endif
}// namespace ncl
#	if 0 // simple stats to be added later
		template<typename AccumulatorType>
		class order_2_accumulator
			{
			public:
				typedef AccumulatorType value_type;

				order_2_accumulator()
				  :Count(), 
				  Sum(), 
				  Sum2() 
				  {}
				  
				order_2_accumulator(unsigned int count_, const value_type& sum_, const value_type& sum2_)
				  :Count(count_), Sum(sum_), Sum2(sum2_)
					{}

			unsigned int	count() const     	{ return Count; }
			value_type 		sum() const         { return Sum; }
			value_type 		sum_squares() const { return Sum2; }
			value_type 		mean() const        { return Sum/Count; }
			value_type 		variance() const    { return Sum2/Count - (Sum*Sum)/(Count*Count); }

			// Scott Kirkwood: added these - would require sqrt() though.
			value_type std() const         			{ return std::sqrt(variance); }
			value_type std_error_of_mean() 			{ return std() / std::sqrt(Count); }
			value_type root_mean_square()  			{ return std::sqrt(Sum2 / Count); }
			value_type coefficient_of_variation() 	{ return 100 * std() / mean(); }

			template<typename T>
			order_2_accumulator<value_type> bump(const T& value_)
				{
				++Count;
				Sum  += value_;
				Sum2 += value_*value_;
				return *this;
				}
			private:
				unsigned int Count;
				value_type   Sum;
				value_type   Sum2;
			};


		// Default operator used by std::accumulate
		template<typename AccumulatorType, typename T>
		AccumulatorType operator+(const AccumulatorType& init_, const T& value_)
			{
			AccumulatorType accum(init_);
			return accum.bump(value_);
			} 

		typedef order_2_accumulator<double> DblStatArray;

		/* can be used like this:
			int seq[] ={1, 2, 3, 4};
			DblStatArray sum=std::accumulate(seq, seq+4, DblStatArray());
			out << sum.count() << ' ' << sum.sum() << ' ' << sum.sum_squares() << std::endl;
			out << sum.mean() << ' ' << sum.variance() << std::endl;

			double seq2[] = {1., 2., 3., 4.};
			sum=std::accumulate(seq2, seq2+4, sum);
		 	out << sum.count() << ' ' << sum.sum() << ' ' << sum.sum_squares() << std::endl;
			out << sum.mean() << ' ' << sum.variance() << std::endl; 
		*/


		template<class InputIterator>
		inline typename std::iterator_traits<InputIterator>::value_type mean(InputIterator begin, InputIterator end, typename std::iterator_traits<InputIterator>::value_type zero = 0)     
			{
			unsigned int count = 0;
			typename std::iterator_traits<InputIterator>::value_type sum = zero;
			for(InputIterator i=begin; i < end; ++i) 
				{
				sum += *i;
				++count;
				}
			return sum/count;
			}


		template<class InputIterator>
		inline typename std::iterator_traits<InputIterator>::value_type variance(InputIterator begin, InputIterator end, typename std::iterator_traits<InputIterator>::value_type zero = 0)
			{
			typename std::iterator_traits<InputIterator>::value_type mn = mean(begin,end);
			unsigned int count = 0;
			typename std::iterator_traits<InputIterator>::value_type sum = zero;
			for(InputIterator i=begin; i < end; ++i) 
				{
				sum += std::pow(*i - mn, 2);
				++count;
				}
			return sum/(count-1);
			} 
#	endif // 0 // simple stats to be added later
#endif
