#ifndef NCL_NXSDATATYPE_H
#define NCL_NXSDATATYPE_H

#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/misc/utilities.hpp"
#include "pyphy/src/ncl/nxs_exception.hpp"
#include "pyphy/src/ncl/nxs_token.hpp"
#define CODED_STATES_AS_VEC 1

namespace ncl
{
typedef std::vector<VecDbl> 	DoubleMatrix;
/*--------------------------------------------------------------------------------------------------------------------------
|	Exception class thrown by VerifyMatrixIsSquare()
*/
class XNotSquareMatrix {};
void VerifyMatrixIsSquare(const DoubleMatrix &m);	// throws XNotSquareMatrix
bool CheckTriangeInequality(const DoubleMatrix &m, unsigned * const fromN, unsigned * const toN, unsigned * const intermed);
} //namespace ncl
#define MAX_NSTATES_ALLOWED CHAR_MAX
std::vector<unsigned> ConvertBitArrToVecUInt(const DataStorageType *, unsigned wordSize);
namespace ncl 
	{
	class NxsUnalignedMatrix;
	class NxsDiscreteMatrix;
	class NxsDiscreteMatrixBase;
	}
class NxsToken;
#if (MAX_NSTATES_ALLOWED > CHAR_MAX)
#	error do not define MAX_NSTATES_ALLOWED > CHAR_MAX (unless you change NxsDataType::symbolToCode)
#endif

/*--------------------------------------------------------------------------------------------------------------------------
|	Class used to translate between internal representations of of character states and the single character symbols used 
|	to display the state to the user.
|	Also maintains or provides state representations in bits through an array of DataStorageType's (a typedef for Unsigned 
|	char).
|	contains redundant information so that inlined functions can provide rapid translations between a state's index, 
|	bit representation, and symbol.
|
|	Possible extension: templates to control fundamental data storage type (as opposed to using the typedef DataStorageType)
|	for all NxsDataTypes.
|	specializations for types that require different numbers of DataStorageType to store their codes could speed up the 
|	simple types (by alleviated *'s and checking of nCodeUnitsPerCharacter, et).
*/
class NxsDataType 
	{	
	public:
		enum NxsDatatypesEnum
			{ 
				kStandard 	= 0, 
				kDNA		= 1, 
				kRNA		= 2, 
				kNucleotide = 3, 
				kProtein	= 4, 
				kContinuous = 5
			};

		class XStateNotFound {}; //thrown by GetStateChar if no legal state is found in the argument
		
		NxsDataType(const NxsDataType &r);
		NxsDataType(const NxsDatatypesEnum dataTypeCode, bool respectCase = false);
		~NxsDataType();
		NxsDataType &operator=(const NxsDataType &r);
		//	Accessors
		//
		NxsDatatypesEnum		GetDataTypeCategory() const
			{
			return dataTypeCategory;
			}
		char					GetGapChar() const;
		char					GetStateChar(const DataStorageType *, const bool couldBeMissingChar = true) const; //throws XStateNotFound
		char					GetMatchChar() const;
	   	const DataStorageType  *GetMissingBitCode() const;
	    char 					GetMissingChar() const;
		unsigned 		  		GetNumStates() const;
		unsigned 				GetNumElementsPerCode() const {return nCodeUnitsPerCharacter;}
		char					GetStateCharFromIndex(unsigned) const;
	    const DataStorageType  *GetStateInBitCode(char s) const;
	    const char			   *GetSymbols() const;
		bool					IsRespectingCase() const {return respectingCase;}
		bool					IsValidSymbol(char c) const;
	    std::map<char, std::string>	GetEquates() const;
	    
	    //Utilities
	    //
	    unsigned 				CountStates(const DataStorageType *) const;
		std::string				GetStatesAsNexusString(const DataStorageType *, bool polymorphic, bool expandEquates) const;
	    //std::string			GetStatesAsNexusString(const NxsDiscreteMatrixBase *dataMat, unsigned i, unsigned j) const;
	    std::vector<unsigned>	GetOrdinationsOfStates(const DataStorageType *) const;
		template<class T> void 	ReadNextState(T *dataMat, NxsToken& token, unsigned i, unsigned j) const;
		int						StatesComp(const DataStorageType *leftState,const DataStorageType *rightState) const;
		bool					FillWithIntersection(const DataStorageType *pattern, DataStorageType *inter, unsigned patLen) const;
		void					FillInCondLikeArray(const DataStorageType *pattern, double *arr) const;
		
		//Modifiers
		//
		void					AddSymbols(const std::string &s, bool);
		void					AddEquates(const std::string &e);
		void 					CreateNewEquateEntry(char key, const std::string &s, bool polymorphicEquate, bool respectCase = true);
		void					ReplaceSymbols(const std::string s, bool respectCase);
		void					SetGapChar(char );
		void					SetMatchChar(char );
		void					SetMissingChar(char );
		
	private:
		typedef DataStorageType  				  * DataStorageTypePtr;
		typedef std::pair<DataStorageType *, bool> 	EquateValInfo;
		typedef std::map<char, EquateValInfo> 		EquateMap;
		typedef std::vector<char> 					SymbolsListType;
		
		
		void 					AddStateChar(ncl::NxsDiscreteMatrix *dataMat, unsigned i, unsigned j, char ch, const NxsToken &token) const;
		void 					AddStateChar(ncl::NxsUnalignedMatrix *dataMat, unsigned i, unsigned j, char ch, const NxsToken &token) const;
		template<class T> void	AddStateCode(T *dataMat, unsigned i, unsigned j, char,const NxsToken &token) const;
		template<class T> void 	AddStateRange(T *dataMat, unsigned i, unsigned j, char fromC, char toC, const NxsToken &token) const;
		void 					Clear();
		bool 					CheckNewSymbol(char nextSymb, std::set<char> &uniqSymbs) const;
		void 					CheckSymbols(const std::string &s,  bool respectCase) const;
		void					Copy(const NxsDataType &r);
		std::string				GetEquateExpansion(const DataStorageType *) const;
		std::string				GetStatesAsNexusStringNoPunctuation(const DataStorageType *, bool polymorphic, char *, bool expandEquates) const;
		template<class T> void 	ReadNextMultipleStates(T *dataMat, unsigned i, unsigned j, NxsToken &token, bool parsingPoly) const;
		void 					ReplaceWithCheckedSymbols(const std::string s, bool respectCase);
		//void 					SetStateChar(ncl::NxsDiscreteMatrix *dataMat, unsigned i, unsigned j, char ch, const NxsToken &token) const;
		//void 					SetStateChar(ncl::NxsUnalignedMatrix *dataMat, unsigned i, unsigned j, char ch, const NxsToken &token) const;
		void 					SetSymbolToCodePtr(char c, DataStorageTypePtr p);
		DataStorageType		   *AppendCodedStateFromEquate(DataStorageType *); //may require reallocating allCodedStates
		void					SwapInternalAliases(DataStorageType *,DataStorageType *, unsigned);
		bool					Validate() const;
		unsigned 					nStates;		/* the number of possible states in this data type */ 
		unsigned					nCodeUnitsPerCharacter;	/* the length of DataStorageType arrays that store the state codes (this will be 1 + nStates/8) */
		DataStorageType			  * missingCode;	/* code with 1's for all of the states' bits (e.g. stores a pointer to a DataStorageType with the value 15 for DNA type) Alias to the end of contigCodedStates don't delete */
		char						missingSymbol;	/* character to be printed (or read) for missing data */ 
		char						matchSymbol;	/* character to be printed (or read) for data identical to the first taxon */
		char						gapSymbol;		/* character to be printed (or read) for inapplicable data */ 
		SymbolsListType				symbolsList;	/* the index to symbol mapping for the type  (length should be nStates) */
		DataStorageTypePtr			symbolToCode[MAX_NSTATES_ALLOWED]; /*the symbol to code * mapping.  Contains NULL for illegal state symbols */
		mutable DataStorageType	  * scratchSpace;	/* space for a state code used internally. Alias to the end of contigCodedStates don't delete */
		EquateMap					equateCodes;	/* mapping of an equate symbol to a std::pair<bit code, bool for whether the equate indicates polymorphism> */
#		if CODED_STATES_AS_VEC
			std::vector<DataStorageType>		codedStatesVec;
#		else
			DataStorageType			    *allCodedStates; /* bit codes for all of the characters in one array (in order the same order as the state's indices */
			unsigned					nCodedStatesStored;
#		endif
		bool						respectingCase;	/* true if 'a' is not the same state as 'A' */
		NxsDatatypesEnum			dataTypeCategory; 
	
		friend class ncl::NxsDiscreteMatrix;
		friend class ncl::NxsUnalignedMatrix;
	};	
	

typedef boost::shared_ptr<NxsDataType> DataTypeShPtr;

/*--------------------------------------------------------------------------------------------------------------------------
|	Struct-like class that stores a square costMatrix and a symbols list indicating the 
|	order of the symbols in the rows (and columns)
|	Note '\0' may be used in the symbols - this indicates a state that has no name (internal node
|	in a CSTree description).
|	Symbols not found in the NxsDataType's symbols list may be found as well (PAUP allows users to 
|	introduce symbols in NxsDataType descriptions).
*/
class UserTypeDescription
	{
	public:
		std::vector<char> symbols;
		ncl::DoubleMatrix costMatrix;
		
		std::string 		GetMatrixDescription() const;
	};

unsigned IndexOfFirstOnBit(const DataStorageType *c, unsigned nElements);

inline NxsDataType &NxsDataType::operator=(const NxsDataType &r) 
	{
	Copy(r);	//copy calls clear
	return *this;
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	
*/
inline bool NxsDataType::IsValidSymbol(char c) const
	{
	return (symbolToCode[(unsigned char) c] != NULL || (c != '\0' && ( c == matchSymbol || c == gapSymbol)));
	}
/*--------------------------------------------------------------------------------------------------------------------------
|	throws NxsExceptions if illegal or repeated symbols are passed in as `s'
*/
inline void NxsDataType::ReplaceSymbols(
  const std::string s, /*list of symbols in desired order (the order in this string determines the indexing of the states) */
  bool respectCase)
	{
	Clear();
	CheckSymbols(s, respectCase);
	ReplaceWithCheckedSymbols(s, respectCase);
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	Gets the symbol for the `ord'-th state (indexing starts at zero)
*/
inline char NxsDataType::GetStateCharFromIndex(
  unsigned ord) const
	{
	return symbolsList[ord];
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	CompFunction for two bit codes (returns 0 if the codes are equal, -1 is left is less, 1 if right is less)
*/
inline int NxsDataType::StatesComp(const DataStorageType *leftState, const DataStorageType *rightState) const
	{
	//	High index states appear later in the the code array, so we use "ReversePrecedence" array compare
	return ArrayCompareReversePrecedence<DataStorageType>(leftState, rightState, nCodeUnitsPerCharacter);
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the DataType's internal bit code pointer for the state `s'
*/
inline const DataStorageType *NxsDataType::GetStateInBitCode(char s) const
	{
	return symbolToCode[(unsigned char) s];
	}
	    
/*--------------------------------------------------------------------------------------------------------------------------
|	Returns the DataType's internal bit code pointer for the missing code
*/
inline const DataStorageType *NxsDataType::GetMissingBitCode() const 
	{
	return missingCode;
	}
	
inline unsigned NxsDataType::GetNumStates() const 
	{
	return nStates;
	}		

inline char NxsDataType::GetGapChar() const
	{
	return gapSymbol;
	}

inline char NxsDataType::GetMissingChar() const
	{
	return missingSymbol;
	}

inline char NxsDataType::GetMatchChar() const
	{
	return matchSymbol;
	}

inline void NxsDataType::SetGapChar( char g) 
	{
	gapSymbol = g;
	//SetSymbolToCodePtr(gapSymbol, NULL);
	}

inline void NxsDataType::SetMatchChar( char m) 
	{
	matchSymbol = m;
	}


/*--------------------------------------------------------------------------------------------------------------------------
|	Will return a bogus address if the datatype has no symbols (Num states = 0)
*/
inline const char* NxsDataType::GetSymbols() const
	{
	return &symbolsList[0];
	}

  	
/*--------------------------------------------------------------------------------------------------------------------------
|	Used to read groups of symbols surrounded by () (in the case of polymorphism) or {} (in the case of ambiguity)
*/
template<class T> void NxsDataType::ReadNextMultipleStates(
  T *dataMatrix,
  unsigned i, 
  unsigned j,
  NxsToken &token,
  bool parsingPoly) const
  	{
  	if (parsingPoly && j != UINT_MAX)
  		dataMatrix->SetPolymorphic(i, j);
  	char prev = 0;
	dataMatrix->ZeroState(i, j);
 	for(;;)
		{
		token.ReadSingleCharacter();
		if (token.GetTokenReference() == '~')
			{
			token.ReadSingleCharacter();
			AddStateRange<T>(dataMatrix, i, j, prev, token.GetTokenReference()[0], token);
			}
		else
			{
			prev = '\0';
			if (token.GetTokenReference() == '(')
				ReadNextMultipleStates<T>(dataMatrix, i, j, token, true);
			else if (token.GetTokenReference() == '{')
				ReadNextMultipleStates<T>(dataMatrix, i, j, token, false);
			else if ((token.GetTokenReference() == ')' && parsingPoly) || (token.GetTokenReference() == '}' && !parsingPoly))
				return;
			else if ((token.GetTokenReference() == '}' && parsingPoly) || (token.GetTokenReference() == ')' && !parsingPoly))
				{
				char breakChar = (parsingPoly ? ')' : '}');
				std::string errormsg;
				errormsg << "Improperly paired " << breakChar << " character in state label";
				throw NxsException(errormsg, token);
				}
			else
				{
				prev = token.GetTokenReference()[0];
				AddStateCode<T>(dataMatrix, i, j, token.GetTokenReference()[0], token);
				}
			}
		}
  	}

inline std::vector<unsigned> NxsDataType::GetOrdinationsOfStates(
  const DataStorageType *s) const
	{
	return ConvertBitArrToVecUInt(s, nCodeUnitsPerCharacter);	
	}


template<class T> inline void NxsDataType::AddStateCode(
  T *dataMatrix,
  unsigned i, 
  unsigned j, 
  char ch,
  const NxsToken &token) const
	{
	if (j == UINT_MAX)
		return;
	NXS_ASSERT(ch != '\0');
	AddStateChar(dataMatrix, i, j, ch, token);
	}

#endif
