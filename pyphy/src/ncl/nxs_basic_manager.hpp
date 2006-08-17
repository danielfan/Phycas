#ifndef NCL_NXS_BASIC_DATA_INTERFACE_H
#define NCL_NXS_BASIC_DATA_INTERFACE_H

#include <boost/bind.hpp>

#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/misc/nxs_index_set.hpp"
#include "pyphy/src/ncl/misc/utilities.hpp"
#include "pyphy/src/ncl/misc/string_extensions.hpp"

class NxsBlock;
class NxsTaxaBlock;
class NxsCharactersBlock;
class NxsTreesBlock;
typedef std::map<std::string,NxsIndexSet *> IndexSetContainer;

	
class NxsX_IllegalLabel
	{
	public:
		NxsX_IllegalLabel(const std::string & errorMsg)
			:msg(errorMsg)
			{}
		std::string msg;
	}; // exception potentially thrown by OrderedCaseInsensitiveLabels modifiers
class NxsX_IndexNotFound {}; // exception potentially thrown by FindCharOrThrow

unsigned GetIndexByTreatingLabelAsNumber(const std::string & label, unsigned maxNumb);
unsigned GetIndexFromNexusName(const VecString & v, const std::string & n, bool allowNumbers, bool allowNewNumericNames = false);
unsigned GetIndexFromNexusName(const std::map<unsigned, std::string> & m, const std::string & label, bool allowNumbers, unsigned maxNum = UINT_MAX, bool allowNewNumericNames = false);
const std::string & GetLabelFromVec(const VecString & v, unsigned i);
const std::string & GetLabelFromMap(const LabelMap & v, unsigned i);
unsigned GetMaxNexusLabelLength(const VecString & v);
unsigned GetMaxNexusLabelLength(const std::map<unsigned, std::string> & v, unsigned m);

/*----------------------------------------------------------------------------------------------------------------------
|	A data-less class that provides functions for storing and accessing IndexSets in a std::map<std::string,NxsIndexSet *> 
|	Where the key is a Capitalized name of the set.
|	Important: The map "owns" (allocates and deletes) NxsIndexSet, so if the AddSet() function is used DeleteAndClearSets()
|		should also be used (or ownership of the set should be transferred to another object that will delete the set
|		when it is done.
*/
class NxsSetManager
	{
	public :
		STATELESS_FUNC void 				AddSet(const NxsIndexSet &, IndexSetContainer &);
		STATELESS_FUNC void					DeleteAndClearSets(IndexSetContainer & setsToDel);
		STATELESS_FUNC const NxsIndexSet  * GetSet(const std::string & n, const IndexSetContainer &);	
		STATELESS_FUNC VecString			GetSetNames(const IndexSetContainer & setToScan) ;
	protected:
		STATELESS_FUNC NxsIndexSet		  * GetSetNonConst(const std::string & n, const IndexSetContainer &);	
	};


/*----------------------------------------------------------------------------------------------------------------------
|	Stores a numbered list of labels.  
|	When finding the index for a label, the comparison is done in a case-insensitive manner.
*/
class OrderedCaseInsensitiveLabels
	{
	typedef std::map<std::string, unsigned int> CapLabelToIndex;
		
	public:
		void clear()
			{
			origLabels.clear();
			capLabelToIndex.clear();
			}
			
		bool empty() const
			{
			return this->origLabels.empty();
			}
			///	\returns the number of labels stored
		unsigned size() const
			{
			return (unsigned) this->origLabels.size();
			}
		
		void push_back(const std::string & label)
			{
			const unsigned prevSize = (unsigned)this->origLabels.size();
			unsigned int u;
			if (IsAnUnsigned(label, &u) && (u != (prevSize + 1)) )
				{
				std::string errorMsg;
				StrPrintF(errorMsg, "%s is not valid label (a label or %d)", label.c_str(), (int)(1 + prevSize));
				throw NxsX_IllegalLabel(errorMsg);
				}
			this->origLabels.push_back(label);
			this->capLabelToIndex[GetCapitalized(label)] = (unsigned)(prevSize);
			}
		
		void AppendLabels(const VecString & v)
			{
			for_each(v.begin(),  v.end(), boost::bind(&OrderedCaseInsensitiveLabels::push_back, this, _1));
			}
	
		const VecString & GetLabels() const
			{
			return this->origLabels;
			}
			/// does NOT verify uniqueness of label
		void SetLabel(unsigned index, const std::string &label)
			{
			NXS_ASSERT(index < this->size());
			this->origLabels[index] = label;
			this->capLabelToIndex[GetCapitalized(label)] = index;
			}
			///	\returns the index for the label supplied (returns UINT_MAX if the label is not recognized).
		unsigned FindIndex(const std::string &s, const bool allowNumbers = false, const bool allowNewNumericLabels = false) const
			{
			const std::string cap(GetCapitalized(s));
			return FindIndexFromCapitalized(cap, allowNumbers, allowNewNumericLabels);
			}
		
			///	\returns the index for the label supplied (returns UINT_MAX if the label is not recognized).
		unsigned FindIndexFromCapitalized(const std::string & s, const bool allowNumbers = false, const bool allowNewNumericLabels = false) const
			{
			const CapLabelToIndex::const_iterator toIndexIterator = this->capLabelToIndex.find(s);
			if (toIndexIterator != this->capLabelToIndex.end())
			return toIndexIterator->second;
			return (allowNumbers ? GetIndexByTreatingLabelAsNumber(s, allowNewNumericLabels ? UINT_MAX : this->size()) : UINT_MAX);
			}
		
			/// returns the label stored in the origLabels vector (warning: if you called FindIndex with allowNewNumericLabels, 
			/// then this may throw NxsX_IndexNotFound exception).
		const std::string & GetLabel(unsigned i) const 
			{
			return GetLabelFromVec(this->origLabels, i);
			}
			
		unsigned GetMaxLabelLength() const
			{
			return GetMaxNexusLabelLength(this->origLabels);
			}

			
	private:
		
		VecString		origLabels;
		CapLabelToIndex capLabelToIndex;
	
	};
/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector with the (Capitalized) names of all of the sets in setsToScan
*/
inline VecString NxsSetManager::GetSetNames(
  const IndexSetContainer & setToScan)
	{
	return GetKeys<std::string, NxsIndexSet *>(setToScan);
	}




template<class Obj> unsigned GetIndexFromNexusNameOnlyGeneric(const std::vector<Obj> & v, const std::string  &label)
	{
	typedef typename std::vector<Obj>::const_iterator VectIt; 
	VectIt fromIt  = v.begin();
	const VectIt toIt  = v.end();
	NamedObjEqualsStored<Obj> compF(label);
	VectIt resIt = std::find_if(fromIt, toIt, compF);
	if (resIt != toIt)
		return (unsigned) (resIt - fromIt);
	return UINT_MAX;
	}

template<class Obj> unsigned GetIndexFromNexusNameGeneric(const std::vector<Obj> & v, const std::string & label)
	{
	unsigned retVal = GetIndexFromNexusNameOnlyGeneric<Obj>(v, label);
	if (retVal == UINT_MAX)
		return (unsigned)GetIndexByTreatingLabelAsNumber(label, (unsigned)v.size());
	return retVal;
	}

template<class Obj> unsigned GetMaxNexusLabelGeneric(const std::vector<Obj> & v)
	{
	std::string s;
	AppendNumber<unsigned>(s, (unsigned)v.size());
	unsigned maxLen = (unsigned)s.length();
	typedef typename std::vector<Obj>::const_iterator VectIt;
	VectIt fromIt  = v.begin();
	const VectIt toIt  = v.end();
	for (; fromIt != toIt; ++fromIt)
		{
		const unsigned len = (unsigned)fromIt->GetName().length();
		if (len> maxLen)
			maxLen = len;
		}
	return maxLen;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The interface of a class that stores names/numbers for nxsobjects (e.g. taxa, characters, trees)
|	Written to provide consistent interface for adding labels, translating from label to index and vice versa
*/
class NxsIndexInfo
	{
	public :
		virtual ~NxsIndexInfo(){}
		virtual void				Clear() = 0;
		virtual unsigned			FindIndex(const std::string &label ) const = 0;	// returns UINT_MAX if not found
		virtual unsigned			FindIndexFromNamesOnly(const std::string &label ) const = 0;	// returns UINT_MAX if not found
		unsigned					FindIndexOrThrow(const std::string &label ) const;	// throws nosuchchar if not found
		typedef std::vector<std::string> VecString;
		virtual VecString			GetAllLabels() const;
		virtual const std::string   GetLabel( unsigned index ) const;	//ASSUMES that argument is a valid index
		virtual const std::string & GetUserSuppliedLabel( unsigned index ) const = 0;	//ASSUMES that argument is a valid index. returns an empty string if
		virtual unsigned			GetMaxLabelLength() const = 0;
		virtual unsigned			GetSize() const = 0;
	};

typedef NxsIndexInfo TaxaIndexInfo;
typedef NxsIndexInfo CharIndexInfo;
typedef NxsIndexInfo TreeIndexInfo;

inline const std::string &GetLabelFromVec(const  VecString & v, unsigned i)
	{
	if (i >= v.size())
		throw NxsX_IndexNotFound();
	return v[i];
	}

inline const std::string & GetLabelFromMap(const LabelMap & labelMap, UInt i)
	{
	STATIC_CONST const std::string emptyString;
	LabelMap::const_iterator mIt = labelMap.find(i);
	if (mIt != labelMap.end())
		return mIt->second;
	return emptyString;
	}
	
template<class Obj> inline const std::string &GetLabelFromVecGeneric(const std::vector<Obj> & v, unsigned i)
	{
	if (i >= v.size())
		throw NxsX_IndexNotFound();
	return v[i].GetName();
	}

// increases the vector's length to lastNumber by pushing on each positions number (index+1) in the form of a string
void	PadVectorWithNumbers(VecString *, unsigned lastNumber);

inline void	PadVectorWithNumbers(VecString * v, unsigned lastNumber)
	{
	if (v->size() < lastNumber)
		FillVectorWithNumbers(*v, (unsigned)v->size() + 1, lastNumber + 1);
	}
/*----------------------------------------------------------------------------------------------------------------------
|	The interface of a class that stores named sets and partitions of numbered for nxsobjects (e.g. taxa, characters, 
|	trees)
*/	
class NxsSetInfo
	{
	public :
		virtual ~NxsSetInfo(){}
		virtual void				AddSet(const NxsIndexSet &) = 0;	// copies the TaxSet
		virtual void				AddPartition(const std::string &, const PartitionDescription&) = 0;
		virtual const NxsIndexSet * GetSet(const std::string &n) const = 0; // returns NULL if there is no TaxSet of the given name
	};
	
typedef NxsSetInfo TaxaSetInfo;
typedef NxsSetInfo CharSetInfo;
typedef NxsSetInfo TreeSetInfo;
	
class NxsNamedListManager : public NxsSetInfo, public NxsIndexInfo
	{
	public:
		virtual ~NxsNamedListManager(){}
		virtual unsigned 	GetNumActive() const = 0;
		unsigned 			GetNumInactive() const {return GetSize() - GetNumActive();}
		virtual bool 		IsActive(unsigned) const = 0;
	protected:
		virtual void		SetActiveStatus(unsigned ind, bool newStatus) = 0;
		virtual void		SetActiveStatus(const NxsIndexSet &, bool newStatus) = 0;
		
	};

class NxsBlock;
/*----------------------------------------------------------------------------------------------------------------------
|	Supplies much of the basic functionality of a managing sets
|	it is up to the derived class to call ClearSets(), IndicesIncreased(), and IndicesDecreased()
*/
class NxsBasicListManager : public NxsNamedListManager , public NxsSetManager
	{
	public: 
		NxsBasicListManager();
		virtual ~NxsBasicListManager();
			const 	NxsIndexSet   & GetActiveSet() const;
					unsigned 		GetNumActive() const;
			const 	NxsIndexSet	  * GetSet(const std::string &n) const;
					VecString 		GetSetNames() const;
					bool 			IsActive(unsigned) const;
	protected:
		
					void			AddPartition(const std::string &, const PartitionDescription&);
					void			AddSet(const NxsIndexSet &);
					void			ClearSets(); //gets rid of all sets, and then creates new ALL and ACTIVE sets
					void			IndicesIncreased(unsigned newSize, bool makeNewIndicesActive); //adds indices to the "ALL" Set - to be called when indices are added (NOT REINCLUDED - call SetActiveStatus then)
					void			IndicesDecreased(unsigned newSize); //Calling with 0 results in ClearSets, Removes indices from the "ALL" and "ACTIVE" sets - to be called when indices are deleted (NOT REINCLUDED - call SetActiveStatus then)
					void			SetActiveStatus(unsigned ind, bool newStatus);
					void			SetActiveStatus(const NxsIndexSet &, bool newStatus);
			const 	NxsIndexSet   & GetFullSet() const;
			
	private:
					NxsIndexSet   * GetActiveSetPtr();

			NamedPartitions 	knownPartitions;
			IndexSetContainer	knownSets;
			NxsIndexSet		  * activeSet;	// should be an alias to the active set that is stored in knownSets !!
			NxsIndexSet		  * allSet;	// should be an alias to the all set that is stored in knownSets !!
		 
	};

class NxsBlockManager
	{
	public:
		virtual ~NxsBlockManager(){}
		virtual		CmdResult		NewBlockRead(NxsBlock *) = 0;
	};

class NxsListAndBlockManager : public NxsBlockManager, public NxsBasicListManager
	{
	};
	
inline VecString NxsBasicListManager::GetSetNames() const
	{
	return NxsSetManager::GetSetNames(knownSets);
	}
	
inline NxsIndexSet *NxsBasicListManager::GetActiveSetPtr()
	{
	return activeSet;
	}

inline const NxsIndexSet & NxsBasicListManager::GetActiveSet() const
	{
	return *activeSet;
	}

inline const NxsIndexSet & NxsBasicListManager::GetFullSet() const
	{
	return *allSet;
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	Searches through the container ofindex sets, setToSearch, and returns a pointer to the one that has a name that matches 
|	the argument (string test is NOT case sensitive);
|	returns NULL if a match isn't found.
*/
inline const NxsIndexSet *NxsSetManager::GetSet(
  const std::string 	&csN,			/* name of the desired NxsIndexSet */ 
  const IndexSetContainer  &setToSearch)
  	{
  	return GetSetNonConst(csN, setToSearch);
  	}

	
inline NxsBasicListManager::~NxsBasicListManager()
	{
	DeleteAndClearSets(knownSets);
	}
	
inline void NxsBasicListManager::AddSet(const NxsIndexSet &newSet) 
	{
	NxsSetManager::AddSet(newSet, knownSets);
	}
	
inline void NxsBasicListManager::AddPartition(const std::string &partName, const PartitionDescription& partitionD)
	{
	knownPartitions[partName] = partitionD;
	}

inline const NxsIndexSet * NxsBasicListManager::GetSet(const std::string &n) const
	{
	return NxsSetManager::GetSet(n, knownSets);
	}

inline unsigned NxsBasicListManager::GetNumActive() const
	{
	return GetActiveSet().size();
	}
	
inline bool NxsBasicListManager::IsActive(unsigned n) const
	{
	return GetActiveSet().IsAMember(n);
	}
	
inline unsigned NxsIndexInfo::FindIndexOrThrow(
  const std::string &label ) const
  	{
	unsigned ret = FindIndex(label);
	if (ret == UINT_MAX)
		throw NxsX_IndexNotFound();
	return ret;
	}
inline void  NxsBasicListManager::SetActiveStatus(unsigned s, bool makeActive)
	{
	if (makeActive)
		GetActiveSetPtr()->insert(s);
	else
		GetActiveSetPtr()->erase(s);
	}
	
inline void  NxsBasicListManager::SetActiveStatus(const NxsIndexSet &s, bool makeActive)
	{
	if (makeActive)
		GetActiveSetPtr()->insert(s);
	else
		GetActiveSetPtr()->erase(s);
	}

#endif
