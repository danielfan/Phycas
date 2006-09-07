#if !defined(PHO_TREEID_H)
#define PHO_TREEID_H

#include "phypy/src/oldphycas/split.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	A wrapper around a set of splits
*/
class TreeID
	{
	public:
		bool 			ContainsSplit(const Split& s) const;
		std::string 		CreateIdRepresentation() const;
		bool 			Equals( const TreeID& other ) const;
		unsigned 		GetNumSplitsInID() const;
		SplitSet const &GetSplitSet() const;
		bool 			IsLessThan( const TreeID& other ) const;
		bool 			IsValid() const;
		double 			PctIdentity( const TreeID& other ) const;
		void 			ShowComparison( std::string& msg, const char* this_name, const char* other_name, const TreeID& other ) const;
		
		bool 			operator<( const TreeID& other ) const;
		bool 			operator==( const TreeID& other) const;
		bool 			operator!=( const TreeID& other) const;
        		
		void 			AddSplit(const Split& s);
		void 			Clear();
	private:
		
		void 	GetIntersection(const TreeID &otherSet, SplitSet *intersect) const;
		SplitSet 				 split_set;
		friend class SplitManager;

		friend std::istream& operator>>(std::istream& in, TreeID& id );
		friend std::string& operator<<(std::string& out, const TreeID& id );
	};
std::string& operator<<(std::string& out, const TreeID& id );

inline SplitSet const &TreeID::GetSplitSet() const
	{
	return split_set;
	}

inline void TreeID::Clear() 
	{
	split_set.clear();
	}

inline bool TreeID::ContainsSplit(const Split& s) const
	{
	return (split_set.find(s) != split_set.end());
	}

inline void TreeID::AddSplit(const Split& s)
	{
	split_set.insert(s);
	}
	
inline bool TreeID::IsLessThan( const TreeID& other ) const
	{
	return (split_set < other.split_set);
	}

inline unsigned TreeID::GetNumSplitsInID() const
	{ 
	return ((unsigned)split_set.size()); 
	}

inline bool TreeID::IsValid() const
	{ 
	return !(split_set.empty()); 
	}

inline bool TreeID::operator<( const TreeID& other ) const
	{ 
	return IsLessThan( other ); 
	}
	
inline bool TreeID::operator==( const TreeID& other) const
	{ 
	return Equals( other ); 
	}
	
inline bool TreeID::operator!=( const TreeID& other) const
	{ 
	return !Equals( other ); 
	}
			
inline bool TreeID::Equals( const TreeID& other ) const
	{
	return (split_set == other.split_set);
	}


#endif
