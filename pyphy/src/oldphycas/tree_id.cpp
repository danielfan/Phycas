#include "phycas/force_include.h"
#include "ncl/misc/string_extensions.hpp"
#include "phycas/trees/split.hpp"
#include "phycas/trees/tree_id.hpp"
#include "ncl/misc/algorithm_extensions.hpp"
using std::max;
using std::istream;
using std::string;
string TreeID::CreateIdRepresentation() const
	{
	string s;
	for(SplitSet::const_iterator i = split_set.begin(); i != split_set.end(); ++i )
		{
		if (i != split_set.begin())
			s << ':';
		s << i->CreateIdRepresentation();
		}
	return s;
	}


inline void TreeID::GetIntersection(const TreeID &otherSet, SplitSet *intersect) const
	{
	PHYCAS_ASSERT(intersect != NULL);
	set_intersection( split_set.begin(), split_set.end(), otherSet.split_set.begin(), otherSet.split_set.end(), inserter( *intersect, intersect->begin() ) );
	}
	
double TreeID::PctIdentity( const TreeID& other ) const
	{
	SplitSet intersectionSet;
	GetIntersection(other, &intersectionSet);
	double x = (double) intersectionSet.size();
	double n = (double) max(split_set.size(), other.split_set.size());
	return ( n > 0.0 ? (100.0 * x / n) : 0.0);
	}

void TreeID::ShowComparison( 
  string& msg, 
  const char* this_name, 
  const char* other_name, 
  const TreeID& other ) const
	{
	// Get intersection
	//
	SplitSet intersectionSet;
	GetIntersection(other, &intersectionSet);
	unsigned x = (unsigned)intersectionSet.size();
	unsigned n = (unsigned)split_set.size();
	msg << "There are " << x << " (out of " << n << ") splits in common\n";
	for(SplitSet::const_iterator i = intersectionSet.begin(); i != intersectionSet.end(); ++i ) 
		msg << "  " << (*i) << '\n';
	
	// Get splits found in this and not found in other
	//
	SplitSet thisUniqueSet;

	set_difference( split_set.begin(), split_set.end(), other.split_set.begin(), other.split_set.end(),
	inserter( thisUniqueSet, thisUniqueSet.begin() ) );

	x = (unsigned)thisUniqueSet.size();
	msg << "There are " << x << " splits found in " << this_name << " that are not in " << other_name << '\n';
	for(SplitSet::const_iterator i = thisUniqueSet.begin(); i != thisUniqueSet.end(); ++i )
		msg << "  " << (*i) << '\n';
	
	// Get splits found in other and not found in this
	//
	SplitSet otherUniqueSet;

	set_difference( other.split_set.begin(), other.split_set.end(), split_set.begin(), split_set.end(),
	inserter( otherUniqueSet, otherUniqueSet.begin() ) );

	x = (unsigned) otherUniqueSet.size();
	msg << "There are " << x << " splits found in " << other_name << " that are not in " << this_name << '\n';
	for(SplitSet::const_iterator i = otherUniqueSet.begin(); i != otherUniqueSet.end(); ++i ) 
		msg << "  " << (*i) << '\n';
	}


string& operator<<(string& out, const TreeID& id )
	{
	out << (unsigned)id.split_set.size() << ' ' << (unsigned)Split::nunits;
	for(SplitSet::const_iterator i = id.split_set.begin(); i != id.split_set.end(); ++i )
		out << (*i);
	return out;
	}

istream& operator>>( istream& in, TreeID& id )
	{
	id.Clear();	
	unsigned sz;
	unsigned nu;
	in >> sz >> nu;
	PHYCAS_ASSERT( nu == Split::nunits );
	for (unsigned i = 0; i < (unsigned)sz; ++i)
		{
		Split x;
		in >> x;
		id.split_set.insert(x);
		}
	
	return in;
	}

