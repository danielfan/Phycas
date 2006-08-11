//#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/nxs_basic_manager.hpp"

using std::pair;
using std::map;
using std::max;
using std::string;
using std::vector;

/// SLOW (constructs a vector of all labels up to GetSize() by callying GetLabels()
vector<string> NxsIndexInfo::GetAllLabels() const
	{
	const UInt s = GetSize();
	vector<string> l(s);
	for (UInt i = 0; i < s; ++i)
		l[i] = GetLabel(i);
	return l;
	}
		
const string NxsIndexInfo::GetLabel(unsigned index) const
	{
	const string & usl = GetUserSuppliedLabel(index);
	if (usl.empty())
		{
		string tempLabelScratchSpace;
		tempLabelScratchSpace << (index + 1);
		return tempLabelScratchSpace;
		}
	return usl;
	}

NxsBasicListManager::NxsBasicListManager()
	:activeSet(NULL),
	allSet(NULL)
	{   //create the all and active set and put them in the list of 
	allSet = new NxsIndexSet("All");
	activeSet = new NxsIndexSet("Active");
	knownSets["ALL"] = allSet;	
	knownSets["ACTIVE"] = activeSet;	
	}
	
void NxsBasicListManager::ClearSets()
	{
		// we remove the All and Active sets, so that they don't get deleted, and then add them
		// back to the knownSets.
	knownSets.erase("ALL"); 
	knownSets.erase("ACTIVE");
	DeleteAndClearSets(knownSets);
	allSet->clear();
	activeSet->clear();
	knownSets["ALL"] = allSet;	
	knownSets["ACTIVE"] = activeSet;	
	}
	
void NxsBasicListManager::IndicesIncreased(unsigned newSize, bool makeNewIndicesActive)
	{
	if (newSize == 0)
		return;
	unsigned prevMax = (allSet->empty() ? 0 : allSet->GetLast());
	NXS_ASSERT(newSize > prevMax);
	if (allSet->empty() || newSize-1 >= prevMax)
		{
		NxsIndexSet temp(prevMax, newSize-1);
		allSet->insert(temp);
		if (makeNewIndicesActive)
			activeSet->insert(temp);
		}
	}
	
void NxsBasicListManager::IndicesDecreased(unsigned newSize)
	{
	assert(newSize <= (allSet->empty() ? 0 : allSet->GetLast()));
	if (newSize == 0)
		{
		allSet->clear();
		activeSet->clear();
		ClearSets();
		}
	else
		for (IndexSetContainer::iterator sIt = knownSets.begin(); sIt != knownSets.end(); ++sIt)
			sIt->second->Crop(0, newSize);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Adds the set c to the container of sets, mappedSets.
|	If there is already a set with the same name (case-insenstive) that NxsIndexSet is made = to c.
|	otherwise a new NxsIndexSet is allocated and put into the mappedSets.
|	(To append indices to a set that is already in the map call GetSetNonConst() to get a pointer to the set and add to it
|	directly)
*/
void NxsSetManager::AddSet(
  const NxsIndexSet   &c,
  IndexSetContainer	   &mappedSets)
	{
	string s(c.GetName());
	ToUpper(s);
	IndexSetContainer::iterator dcsIt = mappedSets.find(s);
	if (dcsIt == mappedSets.end())
		mappedSets[s] = new NxsIndexSet(c);	
	else
		{
		NxsIndexSet * tempCS = dcsIt->second;
		*tempCS = c; // replace the old char set with the new one
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	called to empty the container holding taxsets or charsets (whichever is sent as an argument) by deleting the NxsIndexSet
|	that they hold and clearing the container
*/
void NxsSetManager::DeleteAndClearSets(
  IndexSetContainer &setsToDel)
	{
	IndexSetContainer::iterator dcsIt = setsToDel.begin();
	for (;  dcsIt != setsToDel.end() ; ++dcsIt)
		delete dcsIt->second;
	setsToDel.clear();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Searches through the container ofindex sets, setToSearch, and returns a pointer to the one that has a name that matches 
|	the argument (string test is NOT case sensitive);
|	returns NULL if a match isn't found.
|	
*/
NxsIndexSet *NxsSetManager::GetSetNonConst(
  const string 	& csN,			/* name of the desired NxsIndexSet */ 
  const IndexSetContainer  &setToSearch)
  	{
  	string n(csN);
	ToUpper(n);
	IndexSetContainer::const_iterator csIt = setToSearch.find(n);
	if (csIt != setToSearch.end())
		return csIt->second;
	return NULL;
	}
	

unsigned GetIndexFromNexusName(
  const VecString &v, 
  const string &label, 
  bool allowNumbers,
  bool allowNewNumericNames) //@ revisit this after Vancouver, March 2005 true if you allow "17" when you only have 10 labels in v ... should lead to extending the labels up to the index of the number sent in label
	{
	NStrCaseInsensitiveEquals compF(label);
	VecString::const_iterator iter = find_if(v.begin(), v.end(), compF);
	if (iter != v.end())
		return ((unsigned) (iter - v.begin()) );
	return allowNumbers ? GetIndexByTreatingLabelAsNumber(label, allowNewNumericNames ? UINT_MAX : (unsigned)v.size()) : UINT_MAX ;
	}

unsigned GetIndexByTreatingLabelAsNumber(const string &label, unsigned maxNumb)
	{
	UInt u;
	if (IsAnUnsigned(label, &u) &&  u > 0 && u <= maxNumb)
		return u - 1;
	return UINT_MAX;
	}

unsigned GetMaxNexusLabelLength(const VecString &v)
	{
	string s;
	s << (unsigned)v.size();
	unsigned numbLen = (unsigned)s.length();
	unsigned labLen = GetMaxSize<VecString::const_iterator>(v.begin(), v.end());
	return (numbLen > labLen ? numbLen : labLen);
	}

unsigned GetMaxNexusLabelLength(const map<unsigned, string> &v, unsigned maxNum)
	{
	string s;
	s << (unsigned) maxNum;
	unsigned numbLen = (unsigned)s.length();
	unsigned labLen = 0;
	for (map<unsigned, string>::const_iterator vIt = v.begin(); vIt != v.end(); ++vIt)
		labLen = (vIt->second.size() > labLen ? (unsigned)vIt->second.size() : labLen);
	return max(numbLen, labLen);
	}

unsigned GetIndexFromNexusName(
  const std::map<unsigned, std::string> & m,
  const std::string & label, 
  bool allowNumbers, 
  unsigned maxNum,  
  bool allowNewNumericNames) // true if you allow "17" when you only have 10 labels in v ... should lead to extending the labels up to the index of the number sent in label
	{
	for (map<unsigned, string>::const_iterator mIt = m.begin(); mIt != m.end(); ++mIt)
		{
		if (EqualsCaseInsensitive(mIt->second, label))
			return mIt->first;
		}
	return allowNumbers ? GetIndexByTreatingLabelAsNumber(label, allowNewNumericNames ? UINT_MAX : maxNum) : UINT_MAX ;
	}


	
