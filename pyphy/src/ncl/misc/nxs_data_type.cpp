//#include "phycas/force_include.h"
#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/misc/nxs_data_type.hpp"
#include "pyphy/src/ncl/nxs_exception.hpp"
#include "pyphy/src/ncl/misc/algorithm_extensions.hpp"

//@TEMP 
#include <iostream>


using std::vector;
using std::map;
using std::set;
using std::pair;
using std::make_pair;
using std::replace;
using std::string;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::islower;
	using std::isupper;
	using std::tolower;
	using std::toupper;
#endif
const unsigned char nOnesSet[] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
  4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
 
 const unsigned char firstOnBit[] ={ 0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0};

const unsigned char lastOnBit[] ={ 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,  
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,  
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,  
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,  
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,  
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,  
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,  
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};

//const unsigned char bitsAbove[] = {0xFF , 0xFE , 0xFC , 0xF8 , 0xF0 , 0xE0 , 0xC0 , 0x80};
//const unsigned char bitsBelow[] = {0x01 , 0x03 , 0x07 , 0x0F , 0x1F , 0x3F , 0x7F , 0xFF};
const unsigned char bitVersion[] = {0x01 , 0x02, 0x04 , 0x08 , 0x10 , 0x20 , 0x40 , 0x80};
const DataStorageType kAllDSTBits = (DataStorageType) ~0;



namespace ncl
{
bool CheckTriangeInequality(
  const DoubleMatrix &mat, 
  unsigned * const fromN, 
  unsigned * const toN, 
  unsigned * const intermed)
	{
	unsigned dimen = (unsigned)mat.size();
	VerifyMatrixIsSquare(mat);
	for (*fromN = 0; *fromN < dimen; *fromN += 1)
		{
		for (*toN = 0; *toN < dimen; *toN += 1)
			{
			if (*toN == *fromN)
				{
				if (mat[*fromN][*toN] != 0.0)
					return false;
				}
			else
				{
				for (*intermed = 0; *intermed < dimen; *intermed += 1)
					{
					if ((*toN != *intermed) && (*fromN != *intermed))
						{
						if (mat[*fromN][*toN] > (mat[*fromN][*intermed] + mat[*intermed][*toN]))
							return false;
						}
					}
				}
			}
		}
	return true;
	}

void VerifyMatrixIsSquare(
  const DoubleMatrix &mat) 
	{
	unsigned dimen = (unsigned)mat.size();
	for (unsigned i = 0; i < dimen; ++i)
		{
		if (mat[i].size() != dimen)
			throw XNotSquareMatrix();
		}
	}
	
}//namespace ncl
/*----------------------------------------------------------------------------------------------------------------------
|	Converts the bit storage to an array of double 0.0's and 1.0's
*/
void NxsDataType::FillInCondLikeArray(
  const DataStorageType *s, /* pointer to the code of bits that represent the state */
  double *arr) /* conditional likelihood array to fill */ const
	{
	replace_all(arr, arr+GetNumStates(), 0.0);
	vector<unsigned> statesPresent = ConvertBitArrToVecUInt(s, GetNumElementsPerCode());
	for (vector<unsigned>::const_iterator spIt = statesPresent.begin(); spIt != statesPresent.end(); ++spIt)
		arr[*spIt] = 1.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
| Public function used in determining whether all of the codes in an array have an intersection of states (for instance
| to determine if a pattern can be explained as an invariant site).  The return value is true if there is overlapping 
| state.
*/
bool NxsDataType::FillWithIntersection(
  const DataStorageType	*arrayOfChars, 		/* the array of bitCodeUnits that are being checked for an intersection */
  DataStorageType		*intersectionCode,		/* allocated space for the intersection to be stored */
  unsigned		 		nCodes) const			/* the length of the array in number of coded characters */
	{
	if (nCodeUnitsPerCharacter == 1)
		{
		*intersectionCode = *missingCode;
		for (unsigned codNum = 0; codNum < nCodes; ++codNum)
			{
			intersectionCode[0] &= *arrayOfChars++;
			if (!intersectionCode[0])
				return false;
			}	
		}
	else
		{
		nxs_copy_n(missingCode, intersectionCode, nCodeUnitsPerCharacter);
		unsigned firstCodeToCheck = 0;
		unsigned lastCodeToCheck = nCodeUnitsPerCharacter - 1;
		for (unsigned codNum = 0; codNum < nCodes; ++codNum)
			{
			bool inter =false;
			for (unsigned b = firstCodeToCheck; b <= lastCodeToCheck; ++b)
				{
				intersectionCode[b] &= arrayOfChars[codNum*nCodeUnitsPerCharacter + b];
				if (intersectionCode[b])
					inter =true;
				}
			if (inter)
				{
				while (!intersectionCode[firstCodeToCheck])
					++firstCodeToCheck;
				while (!intersectionCode[lastCodeToCheck])
					--lastCodeToCheck;
				}
			else
				return false;
			}
		}
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	read a word of length wordLen from the array s, interpretting 1's as numbers that should be added to the vector.
|	eg. 06010000 00000001 with a word length of 2 would return vector of length 3 with 4, 6, and 8
*/
vector<unsigned> ConvertBitArrToVecUInt(
  const DataStorageType *s, 
  unsigned wordLen) /* number of DataStorageType per BitCode */
  	{
	VecUInt retVec;
  	unsigned currOffset = 0;
  	for (unsigned w = 0; w < wordLen; ++w)
		{
		DataStorageType temp = *s++;

		for (unsigned i = 0; i < sizeof(DataStorageType); ++i, currOffset += 8)
			{
			const unsigned char key = (DataStorageType) (temp & kAllDSTBits);
			unsigned char fb = firstOnBit[key];

			if (nOnesSet[key] == 1)
				{
				retVec.push_back(currOffset + fb);
				}
			else if (key != 0)
				{
				const unsigned char lb = lastOnBit[key];
				retVec.push_back(currOffset + fb);
				retVec.push_back(currOffset + lb);
				
				for (fb++; fb <lb; ++fb)
					{
					if (bitVersion[fb] & key)
						retVec.push_back(currOffset + fb);
					}
				}
			temp >>= 8;
			}
		}

	return retVec;
  	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Returns a map of an equate's character code to the string that describes its expansion
*/
map<char, string> NxsDataType::GetEquates() const
	{
	map<char, string> retMap;
	string s;
	for (EquateMap::const_iterator emIt = equateCodes.begin(); emIt != equateCodes.end(); ++emIt)
		retMap[emIt->first] = GetStatesAsNexusString(emIt->second.first, emIt->second.second, true);
	return retMap;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
| 	Creates a datatype with the symbols as they occur in the sym string.
*/
NxsDataType::NxsDataType(
  const NxsDatatypesEnum dataTypeCode,
  bool respectCase)
    :nStates(0),
    nCodeUnitsPerCharacter(1),
    missingCode(NULL),
	missingSymbol('?'),
	matchSymbol('\0'),
	gapSymbol('\0'),
	scratchSpace(NULL),
#	if !CODED_STATES_AS_VEC
		allCodedStates(NULL),
		nCodedStatesStored(0),
#	endif
	respectingCase( respectCase),
	dataTypeCategory(dataTypeCode)
  	{
	for (unsigned i = 0; i < MAX_NSTATES_ALLOWED; ++i)
		symbolToCode[i]  = NULL;
	string sym;  /*list of symbols in desired order (the order in this string determines the indexing of the states) */
	switch(dataTypeCategory)
		{
		case kDNA:			//fall through
		case kNucleotide:
			sym = "ACGT";
			break;
		case kRNA:
			sym = "ACGU";
			break;
		case kProtein:
			sym = "ACDEFGHIKLMNPQRSTVWY*";
			break;
		default:
			sym = "01";
		}
	ReplaceSymbols(sym, respectCase);
	if(dataTypeCategory == kDNA || dataTypeCategory == kRNA || dataTypeCategory == kNucleotide)
		{
		CreateNewEquateEntry('R', "AG", false, false);
		CreateNewEquateEntry('M', "AC", false, false);
		CreateNewEquateEntry('S', "CG", false, false);
		CreateNewEquateEntry('V', "ACG", false, false);
		if(dataTypeCategory == kDNA || dataTypeCategory == kNucleotide)
			{
			CreateNewEquateEntry('Y', "CT", false, false);
			CreateNewEquateEntry('K', "GT", false, false);
			CreateNewEquateEntry('W', "AT", false, false);
			CreateNewEquateEntry('B', "CGT", false, false);
			CreateNewEquateEntry('D', "AGT", false, false);
			CreateNewEquateEntry('H', "ACT", false, false);
			CreateNewEquateEntry('N', "ACGT", false, false);
			CreateNewEquateEntry('X', "ACGT", false, false);
			if (dataTypeCategory == kNucleotide)
				CreateNewEquateEntry('U', "T", false, false);
			}
		else if(dataTypeCategory == kRNA)
			{
			CreateNewEquateEntry('Y', "CU", false, false);
			CreateNewEquateEntry('K', "GU", false, false);
			CreateNewEquateEntry('W', "AU", false, false);
			CreateNewEquateEntry('B', "CGU", false, false);
			CreateNewEquateEntry('D', "AGU", false, false);
			CreateNewEquateEntry('H', "ACU", false, false);
			CreateNewEquateEntry('N', "ACGU", false, false);
			CreateNewEquateEntry('X', "ACGU", false, false);
			}
		}
	else if(dataTypeCategory == kProtein)
		{
	  	CreateNewEquateEntry('B', "DN", false, false);
		CreateNewEquateEntry('Z', "EQ", false, false);
		CreateNewEquateEntry('X', "ACDEFGHIKLMNPQRSTVWY", false, false);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Normal Copy constructor (expcept that the allCodedStates is not copied)
*/
NxsDataType::NxsDataType(const NxsDataType &r)
#	if !CODED_STATES_AS_VEC
		:allCodedStates(NULL),
		nCodedStatesStored(0)
#	endif
	{
	Copy(r);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
| 	Does the work of the copy constructor
*/
void NxsDataType::Copy(const NxsDataType &r)
    {
  	Clear();
  	if (r.GetNumStates() > 0)
  		{
		dataTypeCategory = r.dataTypeCategory;
	  	const string symbList(r.GetSymbols());
	  	ReplaceSymbols(symbList, r.respectingCase);
	  	SetGapChar(r.GetGapChar());
		SetMatchChar(r.GetMatchChar());
		SetMissingChar(r.GetMissingChar());
		for (EquateMap::const_iterator ecIt = r.equateCodes.begin(); ecIt !=  r.equateCodes.end(); ++ecIt)
			{
			char c;
			const string s = r.GetStatesAsNexusStringNoPunctuation(ecIt->second.first, ecIt->second.second, &c, true);
			//std::cout << "NxsDataType::Copy, CreateNewEquateEntry(" << ecIt->first << ", " << s  << ", " << ecIt->second.second << ")"<< std::endl;
			CreateNewEquateEntry(ecIt->first, s, ecIt->second.second);
			}
		}
  	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Calls Clear
*/
NxsDataType::~NxsDataType()
	{
 	Clear();
 	}
 
/*----------------------------------------------------------------------------------------------------------------------
| 	maps the character c and it's lower/upper case equivalent (if not respecting case) to the pointer p that holds the 
|	state binary code.
*/
void NxsDataType::SetSymbolToCodePtr(char c, DataStorageTypePtr p)
	{
	symbolToCode[(int) c] = p;
	if (!respectingCase)
		{
		if (islower(c))
			symbolToCode[toupper(c)] = p;
		else if (isupper(c))
			symbolToCode[tolower(c)] = p;
		}		
	}
 	
/*----------------------------------------------------------------------------------------------------------------------
| 	Replaces the symbols with those in `s' (which are assumed to be legal).
*/
void NxsDataType::ReplaceWithCheckedSymbols(  
  const string s, 
  bool respectCase)
	{
	Clear();
	nStates = (unsigned) s.length();
	if (nStates  == 0)
		{
		NXS_ASSERT(0);	//I can't think of case when we'd want to initialize with empty symbols list, comment out the NXS_ASSERT if you can
		return;
		}
	respectingCase = respectCase;
	NXS_ASSERT(!symbolsList.empty() && symbolsList[symbolsList.size() - 1] == '\0');
	symbolsList.pop_back();
	nxs_copy(s.begin(), s.end(), back_inserter(symbolsList));
	symbolsList.push_back('\0');	// so that array will be null terminated
	const unsigned bitsPerUnit = (sizeof(DataStorageType) * 8);
	nCodeUnitsPerCharacter = 1 + (nStates -  1)/ bitsPerUnit;
	//	allocate enough room to store all of the state bit codes contiguously + 2 more 
	//	one for the missing data code and one for the scratch space
	//
#	if CODED_STATES_AS_VEC
		codedStatesVec = vector<DataStorageType>(nCodeUnitsPerCharacter * (2 + nStates), DataStorageType(0));
#	else
		nCodedStatesStored = (2 + nStates);
		allCodedStates = new DataStorageType [ nCodeUnitsPerCharacter * nCodedStatesStored];
#	endif
	//	Initialize the codes, and set the pointers in the symbolToCode lookup table
	//	to the appropriate entry
	//
	const DataStorageType lastStateInEl = kAllDSTBits - (kAllDSTBits >> 1);
	DataStorageType currBit = 1;
	unsigned currEl = 0;
	unsigned j = 0;
	for (; j < nStates; ++j)
		{
		//	fill the state code for state j with its bit pattern
		//
		unsigned firstElement = j*nCodeUnitsPerCharacter; 
		for (unsigned k = 0; k < nCodeUnitsPerCharacter; ++k)
			{
			if (k == currEl)
				{
#				if CODED_STATES_AS_VEC
					codedStatesVec[firstElement + k] = currBit;
#				else
					allCodedStates[firstElement + k] = currBit;
#				endif
				if (currBit != lastStateInEl)
					currBit <<= 1;
				else
					{
					currBit = 1;
					++currEl;
					// don't want to set the next element too, so iterate k and set the next
					//	element to zero
					//
					++k;
					if (k < nCodeUnitsPerCharacter)
#						if CODED_STATES_AS_VEC
							codedStatesVec[firstElement + k] = 0;
#						else
							allCodedStates[firstElement + k] = 0;
#						endif
					}
				}
			else
#				if CODED_STATES_AS_VEC
					codedStatesVec[firstElement + k] = 0;
#				else
					allCodedStates[firstElement + k] = 0;
#				endif
			}
		//	alias this element in the contigCodedStates array to the element in 
		//	symbolToCode that is associated with the jth symbol
		//
		NXS_ASSERT(symbolToCode[(int) s[j]] == NULL);
#		if CODED_STATES_AS_VEC
			SetSymbolToCodePtr(s[j], &codedStatesVec[0] + firstElement);
#		else
			SetSymbolToCodePtr(s[j], allCodedStates + firstElement);
#		endif
		}
	//	Now set up the missing data code after the other state codes
	//
#	if CODED_STATES_AS_VEC
		missingCode = &codedStatesVec[0] + (j*nCodeUnitsPerCharacter);
#	else
		missingCode = allCodedStates + (j*nCodeUnitsPerCharacter);
#	endif
	unsigned k = 0; 
	for (; k < nCodeUnitsPerCharacter - 1; ++k)
		missingCode[k] = kAllDSTBits;
	unsigned nStatesInLastCodeUnit = nStates % bitsPerUnit;
	if (nStatesInLastCodeUnit == 0)
		missingCode[k] = kAllDSTBits;
	else
		{
		missingCode[k] = 1;
		for (unsigned mb = 1; mb < nStatesInLastCodeUnit; ++mb)
			missingCode[k] = (DataStorageType) ((missingCode[k] << 1) + 1);
		}
	scratchSpace = missingCode + nCodeUnitsPerCharacter;
	NXS_ASSERT(symbolToCode[(int) missingSymbol] == NULL);
	SetMissingChar(missingSymbol);
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Specifies the missing data symbol
*/
void NxsDataType::SetMissingChar( char m) 
	{
	SetSymbolToCodePtr(missingSymbol, NULL);
	missingSymbol = m;
	SetSymbolToCodePtr(missingSymbol, NULL);
	//SetSymbolToCodePtr(missingSymbol, missingCode); // we commented this out because we now check for the missing symbol
		//explicitly when reading a matrix (so that we can flag sites that are missing).
	}


/*----------------------------------------------------------------------------------------------------------------------
| 	Throws a NxsException if duplicate or illegal symbols are found
*/
void NxsDataType::CheckSymbols(
  const string &rawSymbols, 
  bool respectCase) const // throws NxsException
  	{
  	set<char> uniqSymbs;
  	if (!symbolsList.empty())
  		{
  		SymbolsListType::const_iterator endIt = symbolsList.end();
  		--endIt;
  		NXS_ASSERT(*endIt == '\0');
  		nxs_copy(symbolsList.begin(), endIt, inserter(uniqSymbs, uniqSymbs.begin()));
  		}
	for (unsigned sl = 0; sl < rawSymbols.length(); ++sl)
		{
		char nextSymb = rawSymbols[sl];
		if (strchr("{}()\';", nextSymb) != NULL)
			throw NxsException("the SYMBOLS list cannot include any of the following {}()\'; ");
		CheckNewSymbol(nextSymb, uniqSymbs);
		if (!respectCase)
			{
			if (islower(nextSymb))
				CheckNewSymbol((char) toupper(nextSymb), uniqSymbs);
			else if (isupper(nextSymb))
				CheckNewSymbol((char) tolower(nextSymb), uniqSymbs);
			}
		}
  	}

bool NxsDataType::CheckNewSymbol(
  char nextSymb, 
  set<char> &uniqSymbs) const
	{
	if (equateCodes.find(nextSymb) != equateCodes.end())
		{
		string estring;
		estring << "the SYMBOL " << nextSymb << " is already a default equate code for this DATATYPE";
		throw NxsException(estring);
		}
	pair< set<char>::iterator, bool> wasInserted = uniqSymbs.insert(nextSymb);
	if (!wasInserted.second)
		{
		string estring;
		estring << "the SYMBOL " << nextSymb << " is already in the symbols list";
		throw NxsException(estring);
		}
	return true;
	}
	
void NxsDataType::CreateNewEquateEntry(
  char key,
  const string &s,
  bool polymorphicEquate,
  bool respectCase)
  	{
  	//std::cout << "adding equate " << key << " = " << s <<std::endl; 
  	if (!respectCase)
  		{
  		const char altCap = (islower(key) ? toupper(key) : (isupper(key) ? tolower(key) : '\0'));
  		if (altCap != '\0')
  			CreateNewEquateEntry(altCap, s, polymorphicEquate, true);
  		}
  	if (GetStateInBitCode(key) != NULL)
  		{
  		string errormsg;
  		errormsg << "The equate label " << key << " is already a valid symbol";
  		throw NxsException(errormsg);
  		}
  	if (equateCodes.find(key) != equateCodes.end())
  		{
  		string errormsg;
  		errormsg << "The equate label " << key << " is already an equate label";
  		throw NxsException(errormsg);
  		}
  		
  	unsigned argLen = (unsigned) s.length();
  	set<char> expansion;
  	for (unsigned k = 0; k < nCodeUnitsPerCharacter; ++k)
  		scratchSpace[k] = 0; 
  	for (unsigned i = 0; i < argLen; ++i)
  		{
  		if (IsNexusPunctuation(s[i]) || s[i] == '^')
  			{
  			if (s[i] != '+' && s[i] != '-')
  				{
  				string errormsg;
  				errormsg = "Equate value cannot be a punctuation mark or ^";
  				throw NxsException(errormsg);
  				}
  			}
  		const DataStorageType *n = GetStateInBitCode(s[i]);
  		if (n == NULL)
  			{
  			if (s[i] == missingSymbol)
  				n = missingCode;
  			else
  				{
				map<char, EquateValInfo>::const_iterator eIt = equateCodes.find(s[i]);
				if (eIt == equateCodes.end())
					{
					string errormsg;
					errormsg << "The equate value " << s[i] << " is not a known symbol (or equate label)";
					throw NxsException(errormsg);
					}
				n = eIt->second.first;
				if (eIt->second.second)
					polymorphicEquate = true;	//@ Is there a right way to do this? currently - if any of the expansion is polymorphic, consider the equate polymorphic
				}
  			}
  		for (unsigned m = 0; m < nCodeUnitsPerCharacter; ++m)
  			scratchSpace[m] |= n[m];
  		}
  	// now we put the coded stat in the allCodedStates array and store the pointer in 
  	//
  	symbolToCode[(int) key] = AppendCodedStateFromEquate(scratchSpace); //note that equates are case-sensitive so we don't call SetSymbolToCodePtr
  	equateCodes[(char) key] = make_pair<DataStorageType *, bool>(symbolToCode[(int) key], polymorphicEquate);
  	}
  
void NxsDataType::SwapInternalAliases(DataStorageType * begOfOldPtr, DataStorageType *begOfNewPtr, unsigned len)
  	{
	for(unsigned i = 0; i < len; i += nCodeUnitsPerCharacter)
		{
		for (EquateMap::iterator emIt = equateCodes.begin(); emIt != equateCodes.end(); ++emIt)
			{
			if (emIt->second.first == begOfOldPtr)
				emIt->second.first = begOfNewPtr;
			}
		replace(symbolToCode, symbolToCode + MAX_NSTATES_ALLOWED, begOfOldPtr, begOfNewPtr);
		if (scratchSpace == begOfOldPtr)
			scratchSpace = begOfNewPtr;
		if (missingCode == begOfOldPtr)
			missingCode = begOfNewPtr;
		begOfOldPtr	+= nCodeUnitsPerCharacter;
		begOfNewPtr	+= nCodeUnitsPerCharacter;
		}
  	NXS_ASSERT(Validate());
	}
	
bool NxsDataType::Validate() const
	{
	//std::cout << "symbols list = ";
	//for (SymbolsListType::const_iterator sym = symbolsList.begin(); sym != symbolsList.end(); ++sym)
	//	std::cout << *sym ;
	//std::cout << " should have length " << nStates << ". Reported to be "<< symbolsList.size() <<  std::endl;
	NXS_ASSERT(nStates + 1== symbolsList.size()); // symbolsList is Null terminated vector
	if (nStates == 0)
		return true;
	NXS_ASSERT(nCodeUnitsPerCharacter > 0);
	for (SymbolsListType::const_iterator sym = symbolsList.begin(); sym != symbolsList.end(); ++sym)
		{
		if (*sym)
			{
			//std::cout << "validating symbol"<< *sym << std::endl;
			NXS_ASSERT(symbolToCode[(unsigned)*sym]);
			if (!respectingCase)
				{
				if (islower(*sym))
					{
					NXS_ASSERT(symbolToCode[toupper(*sym)]);
					}
				else if (isupper(*sym))
					{
					NXS_ASSERT(symbolToCode[tolower(*sym)]);
					}
				}
			}
		}
	if (missingSymbol)
		{
		//std::cout << "missing symbol"<< missingSymbol << std::endl;
		NXS_ASSERT(!symbolToCode[(unsigned) missingSymbol]); // missing not in symbol to code lookup
		}
	if (gapSymbol)
		{
		//std::cout << "missing symbol"<< gapSymbol << std::endl;
		NXS_ASSERT(!symbolToCode[(unsigned)gapSymbol]); // gap not in symbolToCode lookup
		}
	for (EquateMap::const_iterator em = equateCodes.begin(); em != equateCodes.end(); ++em)
		{
		NXS_ASSERT(symbolToCode[(unsigned) em->first] == em->second.first);
		NXS_ASSERT(symbolToCode[(unsigned)em->first]);
		}
	
#	if CODED_STATES_AS_VEC
		NXS_ASSERT(codedStatesVec.size() >= nStates * nCodeUnitsPerCharacter);
#		if 0 //tucson2005: this code was previously inserted as a sanity check, but it wasn't compiling on 
			// windows so we have decided to just not check whether or not we are sane
			const int64_t begCodedStates = static_cast<const int64_t>(&codedStatesVec[0]);
			const int64_t endCodedStates = static_cast<const int64_t>(&codedStatesVec[codedStatesVec.size() -1]);
			if (missingCode)
				NXS_ASSERT(static_cast<int64_t>(missingCode) >= begCodedStates && static_cast<int64_t>(missingCode) <= endCodedStates);
			if (scratchSpace)
				NXS_ASSERT(static_cast<int64_t>(scratchSpace) >= begCodedStates && static_cast<int64_t>(scratchSpace) <= endCodedStates);
			for (unsigned i = 0; i < MAX_NSTATES_ALLOWED; ++i)
				{
				if (symbolToCode[i])
					NXS_ASSERT(static_cast<int64_t>(symbolToCode[i]) >= begCodedStates && static_cast<int64_t>(symbolToCode[i]) <= endCodedStates);
				}
			for (EquateMap::const_iterator em = equateCodes.begin(); em != equateCodes.end(); ++em)
				{
				NXS_ASSERT(static_cast<int64_t>(em->second.first) >= begCodedStates && static_cast<int64_t>(em->second.first) <= endCodedStates);
				}
#		endif

#	else
#		warn "NxsDataType::Validate() not fully written"
#	endif
	return true;
	}
/*--------------------------------------------------------------------------------------------------------------------------
|	Introducing a new equate causes havoc because of all of the aliasing of pointers in allCodedStates 
|	(pointed to by missingCode, scrathSpace, and the symbolToCode array.
*/
DataStorageType *NxsDataType::AppendCodedStateFromEquate(DataStorageType * newCode) //may require reallocating allCodedStates
	{
#	if CODED_STATES_AS_VEC
		std::vector<DataStorageType> newVecCodedStates = codedStatesVec;
		for (unsigned i = 0; i < nCodeUnitsPerCharacter;i++)
			 newVecCodedStates.push_back(newCode[i]);
		const unsigned oldSize = (unsigned)codedStatesVec.size();
		DataStorageType * orig = &codedStatesVec[0];
		codedStatesVec = newVecCodedStates;
		SwapInternalAliases(orig, &codedStatesVec[0], oldSize);
		DataStorageType *retPtr = &codedStatesVec[0] + oldSize;
		return retPtr;
#	else
		DataStorageType *newAllCodedStates = new DataStorageType[nCodeUnitsPerCharacter * (nCodedStatesStored + 1)];
		nxs_copy(allCodedStates, allCodedStates + (nCodeUnitsPerCharacter * nCodedStatesStored), newAllCodedStates);
		nxs_copy(newCode, newCode + nCodeUnitsPerCharacter * nCodedStatesStored, newAllCodedStates + (nCodeUnitsPerCharacter * nCodedStatesStored));
		DataStorageType *retPtr = newAllCodedStates + (nCodeUnitsPerCharacter * nCodedStatesStored);
		unsigned i = 0;
		SwapInternalAliases(allCodedStates, newAllCodedStates, nCodedStatesStored);
		++nCodedStatesStored;
		delete [] allCodedStates;
		allCodedStates = newAllCodedStates;
		return retPtr;
#	endif
	}

string NxsDataType::GetStatesAsNexusString(
  const DataStorageType *stateCode,
  bool polymorphic, 
  bool expandEquates) const
  	{
  	char c;
  	string temp = GetStatesAsNexusStringNoPunctuation(stateCode, polymorphic, &c, expandEquates);
  	if (c == 0)
  		return temp;
  	string retStr;
  	retStr.reserve(temp.size() + 2);
  	retStr << (c == 1 ? '(' : '{') << temp << (c == 1 ? ')' : '}');
  	return retStr;
  	}
  
string NxsDataType::GetStatesAsNexusStringNoPunctuation(
  const DataStorageType *stateCode,
  bool polymorphic,
  char *c,
  bool expandEquates) const
	{
	*c = 0;
	vector<unsigned> ords = GetOrdinationsOfStates(stateCode);
	//std::cout << "GetStatesAsNexusStringNoPunctuation call to GetOrdinationsOfStates with " << unsigned(stateCode[3]) << unsigned(stateCode[2]) << unsigned(stateCode[1]) << unsigned(stateCode[0]) << std::endl;
	//std::cout << "returned a vector of length ords = " << ords.size() << std::endl; 
	string s;
	if (ords.size() == 1)
		{
		NXS_ASSERT(ords[0] < symbolsList.size() - 1);	//last symbol is bogus '\0' for null terminated strings
		s = symbolsList[ords[0]];
		}
	else if (ords.empty() || ords.size() == nStates)
		{
		s = GetMissingChar();
		}
	else
		{
		if (!expandEquates)
			{
			for (EquateMap::const_iterator eIt = equateCodes.begin(); eIt != equateCodes.end(); ++eIt)
				{
				if (StatesComp(eIt->second.first, stateCode) == 0)
					{
					s = eIt->first;
					return s;
					}
				}
			}
		*c = (char) (polymorphic ? 1 : 2);
		//s << polymorphic ? '(' : '{';
		for(unsigned k = 0; k < ords.size(); ++k) 
			{
			//std::cout << " expanded to " << ords.size() << " states the ["<< k << "] = "<<unsigned(ords[k]) << std::endl;
			NXS_ASSERT(ords[k] < symbolsList.size() - 1);
			s << symbolsList[ords[k]];
			}
		//s << polymorphic ? ')' : '}';
		}
	return s;
	}
/*--------------------------------------------------------------------------------------------------------
|	Looks through the first nElements elements of the array c and returns the ordination of the first bit that is set to one.
|	-1 is means that none of the bits were on.
|	e.g.    an array 1 0 0 would return 0
|			an array 8 0 0 would return 3 
|			an array 0 1 0 would return 32 (assuming that  sizeof(DataStorageType) == 32) 
*/ 
unsigned IndexOfFirstOnBit(const DataStorageType *c, unsigned nElements)
	{
	
	DataStorageType key;
	unsigned retVal = 0;
	for (unsigned j = 0; j < nElements; ++j)
		{
		unsigned temp = c[j];
		for (unsigned i = 0; i < sizeof(DataStorageType); ++i, retVal += 8)
			{
			key = (DataStorageType) (temp & kAllDSTBits);
			if (key != 0)
				return (retVal + firstOnBit[key]);
			temp >>= 8;
			}
		}
	return UINT_MAX;
	}

					

// assumes the pattern will fit a symbol or equate
//
	
char NxsDataType::GetStateChar(
  const DataStorageType *s,
  bool	couldBeMissing) const
	{
	using std::equal;
	for (unsigned i = 0; i < nStates; ++i)
		{
#		if CODED_STATES_AS_VEC
			if (equal((&codedStatesVec[0]) + (i*nCodeUnitsPerCharacter), (&codedStatesVec[0]) + (1+i)*nCodeUnitsPerCharacter, s))
				return GetStateCharFromIndex(i);
#		else
			if (equal(allCodedStates + (i*nCodeUnitsPerCharacter), allCodedStates + (1+i)*nCodeUnitsPerCharacter, s))
				return GetStateCharFromIndex(i);
#		endif
		}
	const char missingChar = GetMissingChar();
	map<char, EquateValInfo>::const_iterator eIt = equateCodes.begin();
	for (; eIt != equateCodes.end(); ++eIt)
		{
		if ((couldBeMissing || eIt->first != missingChar) && equal(eIt->second.first, eIt->second.first + nCodeUnitsPerCharacter, s))
			return eIt->first;
		}
	if (couldBeMissing && equal(GetMissingBitCode(), GetMissingBitCode() + nCodeUnitsPerCharacter, s))
		return missingChar;
	throw XStateNotFound();
	}
	
void NxsDataType::AddSymbols(
  const string &s, 
  bool respectCase)
	{
	CheckSymbols(s, respectCase);
	string newSym;
	if (!symbolsList.empty())
		{
		SymbolsListType::const_iterator endIt = symbolsList.end();
  		--endIt;
  		NXS_ASSERT(*endIt == '\0');
  		for (vector<char>::iterator sLIt = symbolsList.begin(); sLIt != endIt; ++sLIt)
			newSym << *sLIt;
		}
	string oldSym = newSym;
	newSym << s;
	//	Store the old equates so that we can reconstitute them
	//
	typedef pair<string, bool> strEquateInfo;
	map<char, strEquateInfo> oldEquates;
	map<char, EquateValInfo>::iterator eIt = equateCodes.begin();
	for (; eIt != equateCodes.end(); ++eIt)
		{
		vector<unsigned> ords = GetOrdinationsOfStates(eIt->second.first);
		string eValStr;
		if (ords.size() == 1)
			{
			NXS_ASSERT(ords[0] < symbolsList.size() - 1); //-1 because last char in symbol list is bogus '\0'
			eValStr << symbolsList[ords[0]];
			}
		else if (ords.empty())
			eValStr << GetMissingChar();
		else if (ords.size() == nStates)
			eValStr = oldSym;
		else
			{
			for(unsigned k = 0; k < ords.size(); ++k) 
				{
				NXS_ASSERT(ords[k] < symbolsList.size() - 1); //-1 because last char in symbol list is bogus '\0'
				eValStr << symbolsList[ords[k]];
				}
			}
		oldEquates[eIt->first] = make_pair<string, bool>(eValStr, eIt->second.second);
		}
	char oldMissing = GetMissingChar();
	char oldMatch = GetMatchChar();
	char oldGap = GetMissingChar();
	
	ReplaceWithCheckedSymbols(newSym, respectCase);
	
	SetMissingChar(oldMissing);
	SetMatchChar(oldMatch);
	SetGapChar(oldGap);
	map<char, strEquateInfo>::iterator oeIt = oldEquates.begin();
	for (; oeIt != oldEquates.end(); ++oeIt)
		CreateNewEquateEntry(oeIt->first, oeIt->second.first, oeIt->second.second);
	}
		
void NxsDataType::Clear()
	{
	symbolsList.clear();
	symbolsList.push_back('\0');
	equateCodes.clear();
	nCodeUnitsPerCharacter = 1;
	nStates = 0;
#	if CODED_STATES_AS_VEC
		codedStatesVec.clear();
#	else
		nCodedStatesStored = 0;
		delete [] allCodedStates;
		allCodedStates = NULL;
#	endif
	//	return special chars to their defaults
	//
	missingSymbol = '?';	
	matchSymbol = gapSymbol = '\0';		// no defaults for these
	for (unsigned i = 0; i < MAX_NSTATES_ALLOWED; ++i)
		symbolToCode[i]  = NULL;
	scratchSpace = NULL; /*note scratchSpace is an alias to the end of contigCodedStates don't delete */
	missingCode = NULL;	/*note missingCode is an alias to the end of contigCodedStates don't delete */
	}
	
unsigned NxsDataType::CountStates(const DataStorageType *s) const
	{
	unsigned stateCount = 0;
	for (unsigned j = 0; j < nCodeUnitsPerCharacter; ++j)
		{
		unsigned temp = s[j];
		const DataStorageType singlekey = (DataStorageType) (temp & kAllDSTBits);
		stateCount += nOnesSet[singlekey];
		for (unsigned i = 1; i < sizeof(DataStorageType); ++i)
			{
			temp >>= 8;
			const unsigned char multikey = (DataStorageType) (temp & kAllDSTBits);
			stateCount += nOnesSet[multikey];
			}
		}
	return stateCount;	
	}
		

void NxsDataType::AddEquates(
  const string &rawEquate)
  	{
	string::const_iterator cIt = rawEquate.begin();
	while (cIt != rawEquate.end())
		{
		char keyCh = *cIt++;
	   	if (!IsLegalNexusChar(keyCh))
	   		ThrowIllegalNexusChar(keyCh, "an equate");
		if(*cIt != '=') 
			{
			string errormsg;
			errormsg << "Expecting '=' in EQUATE definition but found " << *cIt << " instead";
			throw NxsException(errormsg);
			}
		++cIt;
		if (cIt == rawEquate.end())
			throw NxsException("Expecting legal SYMBOLS between = and \" in the EQUATE definition");	
		bool poly = false;
		string vals;
		if (*cIt == '(')
	    	{
	    	poly = true;
	    	++cIt;
	    	while (*cIt != ')')
	    		{
	    		if (cIt == rawEquate.end())
	    			throw NxsException("Expecting ) before then end of the EQUATE definition");	
	    		vals << *cIt++;
	    		}
	    	if (vals.length() == 0)
	    		throw NxsException("Expecting legal SYMBOLS between the () in the EQUATE definition");	
	    	}
	    else if (*cIt == '{')
	    	{
	    	++cIt;
	    	while (*cIt != '}')
	    		{
	    		if (cIt == rawEquate.end())
	    			throw NxsException("Expecting } before then end of the EQUATE definition");	
	    		vals << *cIt++;
	    		}
	    	if (vals.length() == 0)
	    		throw NxsException("Expecting legal SYMBOLS between the {} in the EQUATE definition");	
	    	}
	    else
	    	vals << *cIt;
	    CreateNewEquateEntry(keyCh, vals, poly);
	    ++cIt;
		}
	}


string UserTypeDescription::GetMatrixDescription() const
	{
	string s;
	ncl::VerifyMatrixIsSquare(costMatrix);
	s << (unsigned)costMatrix.size() << "\n\t\t";
	NXS_ASSERT( costMatrix.size() == symbols.size());
	for (unsigned i = 0; i <  symbols.size(); ++i)
		s << symbols[i] << "\t";
	for (unsigned j = 0; j <  costMatrix.size(); ++j)
		{
		s << "\n\t[" << symbols[j] << "]\t";
		for (unsigned k = 0; k <  costMatrix.size(); ++k)
			{
			if (k == j)
				s << ".\t";
			else if (costMatrix[j][k] == DBL_MAX)
				s << "i\t";
			else
				s << costMatrix[j][k] << "\t";
			}
		}
				
	return s;
	}


