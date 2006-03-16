#ifndef NCL_INDEXSET_H
#define NCL_INDEXSET_H

#include "ncl/nxs_defs.hpp"

class NxsIndexSet 
	{
	public:
		class const_iterator
			{
			friend class NxsIndexSet;
					std::set<unsigned>::const_iterator includedIterator;
			public	:
				bool operator==(const const_iterator r) const
					{
					return (includedIterator == r.includedIterator);
					}
				bool operator!=(const const_iterator r)	const
					{
					return !((*this) == r);
					}
				unsigned operator*() const
					{
					return *includedIterator;
					}
				const_iterator operator++() {++includedIterator; return *this;}	//prefix
				const_iterator operator--() {--includedIterator; return *this;}	//prefix
				const_iterator operator++(int) {const_iterator ret(*this); ++includedIterator; return ret;}	//postfix
				const_iterator operator--(int) {const_iterator ret(*this); --includedIterator; return ret;}	//postfix
				const_iterator(const std::set<unsigned>::const_iterator &r)
				  : includedIterator(r)
					{
					}					
			};

	private:
	
		friend class const_iterator;
		std::set<unsigned> included;
	
	public:
			// perhaps the name field should be removed ?, no, we're storing sets by a capitalized version of their name, this field should hold the original name
		std::string name;	
		
		NxsIndexSet() {}
		explicit NxsIndexSet(unsigned, unsigned);
		explicit NxsIndexSet(std::string pname);
		explicit NxsIndexSet(std::string pname, unsigned, unsigned);
		explicit NxsIndexSet(std::string pname,const NxsIndexSet&);
		explicit NxsIndexSet(const std::set<unsigned> &stdIS);
		
		void 		insert(const NxsIndexSet &);
		void 		insert(unsigned c) {included.insert(c);}
		void 		InsertRange(unsigned, unsigned);
		void		AddToAllIndices(int n);
		void		clear();
		void		Crop(unsigned lowest, unsigned highest);
		const std::string & GetName() const {return name;}
		std::string	GetNexusDescription(bool addOne = true) const;
		unsigned 	size() const;
		unsigned	GetLast() const;
		unsigned	GetFirst() const;
		bool		empty() const;
		void 		erase(const NxsIndexSet &);
		void 		erase(unsigned);
		void		EraseRange(unsigned, unsigned);
		void		Mask(const NxsIndexSet &);
		void 		SetToIntersection(const NxsIndexSet &, const NxsIndexSet&);
		void		SetName(const std::string &n);
		bool 		operator==(const NxsIndexSet &r) const;
		bool 		operator!=(const NxsIndexSet &r) const;
		bool 		IsAMember(unsigned c) const ;
		
		unsigned 	FindNthElement(unsigned n) const;
		unsigned 	CountIndicesBeforeMember(unsigned i) const;
		const_iterator begin() const ; 
		const_iterator end() const;
		
		NxsIndexSet &operator=(const std::set<unsigned> &stdIS);
	};

inline void NxsIndexSet::EraseRange(
  unsigned lowest, 
  unsigned highest)
  	{
  	const std::set<unsigned>::iterator begEraseIt = included.lower_bound(lowest);
  	if (begEraseIt == included.end())
  		return;
  	const std::set<unsigned>::iterator endEraseIt = included.upper_bound(highest);
  	if (endEraseIt != begEraseIt)
  		included.erase(begEraseIt, endEraseIt);
  	}
  	
inline void NxsIndexSet::SetName(const std::string &n)	
	{
	name = n;
	}
	
inline void NxsIndexSet::insert(const NxsIndexSet &r)
	{
	included.insert(r.begin(), r.end());
	}

inline void NxsIndexSet::InsertRange(unsigned f, unsigned t)
	{
	assert(f <= t);
	for (;f <= t; ++f)
		included.insert(f);
	}
		
	
inline void NxsIndexSet::clear()	// does not erase the name
	{
	included.clear();
	}

inline NxsIndexSet::NxsIndexSet(
  unsigned f, 
  unsigned t)
	{
	InsertRange(f,t);
	}
	
inline NxsIndexSet::NxsIndexSet(
  std::string pname, 
  unsigned f, 
  unsigned t)
  	:name(pname)
  	{
  	InsertRange(f,t);
  	}

inline NxsIndexSet::NxsIndexSet(
  std::string pname)
  	:name(pname)
  	{
  	}

inline NxsIndexSet::NxsIndexSet(
  std::string pname,
  const NxsIndexSet &c)
  	:included(c.included),
        name(pname)
  	{
  	}

inline NxsIndexSet::NxsIndexSet(
  const std::set<unsigned> &c)
  	{
  	for (std::set<unsigned>::const_iterator cIt = c.begin(); cIt != c.end(); ++cIt)
  		included.insert(*cIt);
  	}

inline bool NxsIndexSet::operator!=(
 const NxsIndexSet &c) const
  	{
	return !(*this == c);
	}
	
inline bool NxsIndexSet::operator==(const NxsIndexSet &r) const
	{
	return included == r.included;
	}
	
inline bool NxsIndexSet::IsAMember(unsigned c)	const
	{
	return (included.find(c) != included.end());
	}


inline	unsigned 	NxsIndexSet::size() const
	{
	return (unsigned)included.size();
	}
	
inline	unsigned	NxsIndexSet::GetLast() const	
		{
		if(empty())
			return UINT_MAX;
		return *included.rbegin();
		}
		
inline	unsigned	NxsIndexSet::GetFirst() const
		{
		if(empty())
			return UINT_MAX;
		return *included.begin();
		}
inline	bool NxsIndexSet::empty() const	{return (included.empty());}

inline	void NxsIndexSet::erase(
  const NxsIndexSet &c)
	{
	std::set<unsigned> temp;
	set_difference(included.begin(), included.end(), c.included.begin(), c.included.end(), std::insert_iterator< std::set<unsigned> >(temp, temp.begin()));
	included = temp;
	}
	
inline	void NxsIndexSet::erase(
  unsigned c)
	{
	std::set<unsigned>::iterator cIt = included.find(c);
	if (cIt != included.end())
		included.erase(cIt);
	}
	
inline	void NxsIndexSet::SetToIntersection(const NxsIndexSet &f, const NxsIndexSet&s)
	{
	assert( &f != this && &s != this);
	clear();
	set_intersection(f.included.begin(),f.included.end(), s.included.begin(),s.included.end(), std::insert_iterator< std::set<unsigned> >(included, included.begin()));
	}
	
inline NxsIndexSet::const_iterator NxsIndexSet::begin() const
	{ 
	return NxsIndexSet::const_iterator(included.begin());
	}
	 
inline NxsIndexSet::const_iterator NxsIndexSet::end() const
	{ 
	return NxsIndexSet::const_iterator(included.end());
	}

inline NxsIndexSet &NxsIndexSet::operator=(const std::set<unsigned> &stdIS)
	{
	for (std::set<unsigned>::const_iterator cIt = stdIS.begin(); cIt != stdIS.end(); ++cIt)
  		included.insert(*cIt);
  	return *this;
	}


#endif
