#ifndef NXS_DATA_PATTERN_H
#define NXS_DATA_PATTERN_H

#include "pyphy/src/ncl/misc/nxs_data_type.hpp"
#include "pyphy/src/ncl/misc/utilities.hpp"
#include "pyphy/src/ncl/misc/algorithm_extensions.hpp"

typedef DataStorageType *DataPattern;

/*----------------------------------------------------------------------------------------------------------------------
|	DataPatternInfo stores the number of original characters that had this data pattern (`count'), a vector of original
|	character indices (`origIndex'), and the sum of character weights (`sumOfWeights') of all the original characters
|	that had this data pattern.
*/	
class DataPatternInfo
	{
	public:

								DataPatternInfo();
								DataPatternInfo(unsigned i, double w = 0.0);
		
		void					AddIndex(unsigned i, double w = 1.0);
		double					GetWeight() const;
		unsigned				GetCount() const;
		unsigned				GetOrigIndex(unsigned i) const;

	protected:
			
		unsigned 				count;			/**< the number of original characters having this pattern */
		std::vector<unsigned>	origIndex;		/**< list of indices of all original characters having this pattern */
		double					sumOfWeights;	/**< sum of weights of all original characters having this pattern */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Default constructor simply initializes `count' to 0 and `sumOfWeights' to 0.0
*/	
inline DataPatternInfo::DataPatternInfo() : count(0), sumOfWeights(0.0)
	{}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor initializes `count' to 1, `origIndex' to (1, `i') and `sumOfWeights' to `w'.
*/	
inline DataPatternInfo::DataPatternInfo(
  unsigned i,	/**< is the original index */
  double w)		/**< is the sum of weights */
  : count(1), origIndex(1, i), sumOfWeights(w)
	{}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds an original character index and sum of weights to this DataPatternInfo object, incrementing `count' by 1.
*/	
inline void DataPatternInfo::AddIndex(
  unsigned i,	/**< is the original character index */
  double w)		/**< is the sum of weights */
	{ 
	origIndex.push_back(i); 
	++count;
	sumOfWeights += w;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor for `sumOfWeights'.
*/	
inline double DataPatternInfo::GetWeight() const
	{
	return sumOfWeights;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor for `count'.
*/	
inline unsigned DataPatternInfo::GetCount() const
	{
	return count;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor for `origIndex'. Returns `origIndex'[`i']. Assumes `i' is less than `count'.
*/	
inline unsigned DataPatternInfo::GetOrigIndex(
  unsigned i)	/**< is the (0-based) index of the desired value in the `origIndex' vector */
  const
	{
	assert(i < count); 
	return origIndex[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function object used in sorting lists of data patterns.  needed because each pattern doesn't know it length
*/	
class DataPatternLess : public std::less<DataPattern>
	{
	public:
		unsigned nInclTaxa;		/* nIncludedTaxa */ 
		unsigned char nBCU;				/* numOfBitCodeUnits */
		DataPatternLess(unsigned it,unsigned char  bc) : nInclTaxa(it), nBCU(bc) {}
		bool operator() (const DataStorageType *lPat,const DataStorageType *rPat) const 
			{
			return  BitCodeLessThan<DataStorageType>(lPat, rPat, nInclTaxa, nBCU);
			}
	};

typedef std::map<DataPattern, DataPatternInfo, DataPatternLess> SortedPatterns; 
  		
/*--------------------------------------------------------------------------------------------------------------------------
|	Function object used in sorting lists of data patterns.  needed because each pattern doesn't know it length
|
|	because of Microsoft VC's implementation of stl we had to derive this from greater to do a list sort
*/	
class DPatternGreater : public std::greater<DataPattern>
	{
	public:
		unsigned nInclTaxa;	/* nIncludedTaxa */
		unsigned char nBCU;				/* numOfBitCodeUnits per character */
		DPatternGreater(unsigned it,unsigned char bc): nInclTaxa(it), nBCU(bc) {}

		bool operator()(const DataStorageType *lPat, const DataStorageType *rPat) const
			{
			return BitCodeGreaterThan<DataStorageType>(lPat, rPat, nInclTaxa, nBCU);
			}
	};

/*--------------------------------------------------------------------------------------------------------------------------
|	Function object used in comparingof data patterns.  needed because each pattern doesn't know it length
|
*/	
class DPatternEq : public std::binary_function<DataPattern, DataPattern, bool>
	{
	public:
		unsigned arrLen;
		DPatternEq(unsigned it,unsigned char bc): arrLen((unsigned)(it * bc))	{}
		bool operator() (const DataStorageType *lPat,const DataStorageType *rPat) const 
			{
			return  equals_n(lPat, rPat, arrLen);
			}
	};
	
#endif

