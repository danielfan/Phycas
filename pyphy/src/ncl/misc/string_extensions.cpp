#include "phycas/force_include.h"
#include <cstdarg>
#include <boost/scoped_array.hpp>
#include "ncl/misc/string_extensions.hpp"
#include "ncl/nxs_exception.hpp"
using std::string;
using std::vector;
using std::va_list;
#	if defined(C_FUNCS_IN_STD_NAMESPACE)
		using std::strncpy;
		using std::strcpy;
		using std::strchr;
		using std::isgraph;
		using std::isalpha;
		using std::isupper;
		using std::toupper;
		using std::tolower;
		using std::size_t;
#	endif

std::string	Join(const VecString &v, const std::string & separator)
	{
	VecString::const_iterator vIt = v.begin();
	const VecString::const_iterator endIt = v.end();
	std::string retString;
	if (vIt != endIt)
		{
		retString = *vIt++;
		const bool hasSeparator = !separator.empty();
		for (; vIt != endIt; vIt++)
			{
			if (hasSeparator)
				retString << separator;
			retString << *vIt;
			}
		}
	return retString;
	}

bool IsALegalNexusLabelForObjectN(const std::string &s, unsigned n, bool *isANumber)
	{
	UInt u = 0;
	*isANumber = IsAnUnsigned(s, &u);
	if (*isANumber)
		return u == (n+1);
	return IsLegalNexusWord(s);
	}

void DblFormatter::RecreateFormatCode()
	{
	//@ perhaps we should be using %.g here?
	formatCode.clear();
	if (fieldWidth == UINT_MAX)
		{
		if (digitsAfterDecimal == UINT_MAX)
			formatCode = "%lf";
		else
			StrPrintF(formatCode, "%%.%dlf", digitsAfterDecimal);
		}
	else
		{
		if (digitsAfterDecimal == UINT_MAX)
			StrPrintF(formatCode, "%%%dlf", fieldWidth);
		else
			StrPrintF(formatCode, "%%%d.%dlf", fieldWidth, digitsAfterDecimal);
		}
	}


/*--------------------------------------------------------------------------------------------------------------------------
| Written to make it easy to initialize a vector of strings. Similar to the perl split function. Converts a string like
| this -- "A|bro|ken strin|g" -- to a vector of strings with four elements:  "A", "bro", "ken string", and "g".
*/
vector<string> SplitString(
  const string &str, char separator)     /* the string submitted for splitting */
	{
	vector<string> retVec;
	if (str.empty())
		return retVec;
	string ss;
	string::const_iterator p = str.begin();
	const string::const_iterator endIt = str.end();
	for (;;)
		{
		const bool done = (p == endIt);
		if (done || (*p == separator)) 
			{
			retVec.push_back(ss);
			ss.clear();
			if (done)
				return retVec;
			++p;
			if (p == str.end())
				return retVec;
			}
		ss << *p++;
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Converts the stored string to an unsigned int using the standard C function strtol, throwing NxsX_NotANumber if the 
|	conversion fails. Returns UINT_MAX if the number is too large to fit in an unsigned (or was a negative number).
*/
unsigned ConvertToUnsignedOrThrow(const string & s, ncl::IntegerConversion conv) X_SPEC_THROW(NxsX_NotANumber)
	{
	long l = ConvertToLongOrThrow(s, conv);
	if (l < 0)
		throw NxsX_NumberIsTooSmall();
	if (LONG_MAX > UINT_MAX && l > (long) UINT_MAX)
		throw NxsX_NumberIsTooLarge();
	return static_cast<unsigned> (l);
	}
/*--------------------------------------------------------------------------------------------------------------------------
|	Converts the stored string to an int using the standard C function strtol, throwing NxsX_NotANumber if the conversion 
|	fails. Returns INT_MAX if the number is too large to fit in an int or -INT_MAX if it is too small.
*/
int ConvertToIntOrThrow(const string & s, ncl::IntegerConversion conv) X_SPEC_THROW(NxsX_NotANumber)
	{
	long l = ConvertToLongOrThrow(s, conv);
	if (l >= INT_MAX)
		throw NxsX_NumberIsTooLarge();
	if (l <= -INT_MAX)
		throw NxsX_NumberIsTooSmall();
	return static_cast<int> (l);
	}



size_t CopyToCStr(const string & s, char *buffer, const size_t bufferLen)
	{
	assert(buffer != NULL);
	assert(bufferLen > 0);
	size_t len = s.length();
	if (len == 0)
		{
		buffer[0] = '\0';
		return 0;
		}
	if (len >= bufferLen)
		{
		strncpy(buffer, s.c_str(), bufferLen-1);
		buffer[bufferLen-1] = '\0';
		return  (bufferLen-1);
		}
	else
		{
		strcpy(buffer, s.c_str());
		return len;
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Returns a single-quoted version of the string.
*/ 
string GetSingleQuoted(const string & f)
	{
	string withQuotes;
	withQuotes.reserve(f.length() + 4);
	withQuotes = f;
	return ConvertToNexusSingleQuoted(withQuotes);
	}
	
unsigned replace_all_substr(std::string &target, const std::string & searchStr, const std::string & replaceStr)
	{
	unsigned nReplacements = 0;
	unsigned i = (unsigned)target.find(searchStr, 0);
	while (i != string::npos)
		{
		++nReplacements;
		target.replace(i, searchStr.length(), replaceStr);
		i = (unsigned)target.find(searchStr, i + (unsigned)replaceStr.length());
		}
	return i;
	}

string & ConvertToNexusSingleQuoted(string & f)
	{
	const string orig(f);
	f = '\'';
	for (std::string::const_iterator oIt = orig.begin(); oIt != orig.end(); ++oIt)
		{
		f << *oIt;
		if (*oIt == '\'')
			f << '\'';
		}
	f << '\'';
	return f;
	}

 
/*--------------------------------------------------------------------------------------------------------------------------
| Returns true if the string is a abbreviation (or complete copy) of the argument `s'.
*/
bool IsStdAbbreviation(
  const std::string &f, ///test string
  const string & s,   /* the string for which the stored string is potentially an abbreviation */
  bool 	respectCase)	     /* if true, comparison will be case-sensitive */
	{
	if (f.empty())
		return false;
		// s is the unabbreviated comparison string
	const unsigned slen = static_cast<unsigned long>(s.size());
	const unsigned flen = static_cast<unsigned long>(f.size());
	if (flen > slen)
		return false;

		// Examine each character in t and return false (meaning "not an abbreviation")
		// if at any point the corresponding character in s is different
	for (unsigned k = 0; k < flen; ++k)
		{
		if (respectCase)
			{
			if (f[k] != s[k])
				return false;
			}
		else if (toupper(f[k]) != toupper(s[k]))
			return false;
		}
	return true;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Returns true if `f' is a case-insensitive abbreviation (or complete copy) of `s' and the stored string 
| has all of the characters that are in the initial capitalized portion of `s'. For example if `s' is "KAPpa" then 
| "kappa", "kapp", or "kap" (with any capitalization pattern) will return true and all other strings will return false. 
| Always returns false if the stored string has length of zero.
*/
bool IsCapAbbreviation(
  const string & f,  /* the string to check */
  const string & s)   /* the abbreviation template to compare to */
	{
	if (f.empty())
		return false;

		// s is the unabbreviated comparison string
	const unsigned slen = static_cast<unsigned>(s.size());
	const unsigned flen = static_cast<unsigned>(f.size());
		// If the stored string is longer than s then it cannot be an abbreviation of s
	if (flen > slen)
		return false;

	unsigned k = 0;
	for (; k < slen; ++k) 
		{
		if (isupper(s[k]) || !isalpha(s[k]))
			{	// Still in the mandatory portion of the abbreviation
				// Note that non-alphabetic characters in the 
				// that are preceded by no lower-case letters are treated as mandatory
			if (k >= flen || toupper(f[k]) != s[k])
				return false;
			}
		else // Get here if we are no longer in the upper case portion of s 
			break;
		}

	// Check the lower case portion of s and any corresponding characters in t for mismatches
	// Even though the abbreviation is valid up to this point, it will become invalid if
	// any mismatches are found beyond the upper case portion of s
	//
	for (; k < flen; ++k)
		{
		if (toupper(f[k]) != toupper(s[k]))
			return false;
		}

	return true;
	}
/*-------------------------------------------------------------------------------------------------------------------------- 
| Returns true if the string needs to be surrounded by single-quotes to make it a single nexus token.
*/
NexusQuotingReq QuotesNeeded(const std::string & s) 
	{
	NexusQuotingReq nrq = kNoQuotesNeededForNexus;
	for (std::string::const_iterator sIt = s.begin(); sIt != s.end(); ++sIt)
		{
		if (!isgraph(*sIt))
			{
			if (*sIt != ' ')
				return kSingleQuotesNeededForNexus;
			nrq  = kUnderscoresSufficeForNexus;
			}
		else if (strchr("(){}\"-]/\\,;:=*`+<>", *sIt) != NULL)
			{
			// Get here if c is any NEXUS punctuation mark except left square bracket ([) or apostrophe (').
			// [ and ' never get returned as punctuation by NxsToken,
			// so we should never encounter them here. 
			//
			return (s.length() > 1 ? kSingleQuotesNeededForNexus : kNoQuotesNeededForNexus);
			}
		else if (strchr("\'[_", *sIt) != NULL)
			{
			// Get here if c is either an apostrophe or left square bracket. Quotes are needed if one of these
			// characters is all there is to this string
			//
			return kNoQuotesNeededForNexus;
			}
		}
	return nrq;
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Slow (catches exception if false)
*/
bool IsALong(const std::string &s, long * l)
	{
	if (l == NULL)
		{
		long nl;
		return IsALong(s, &nl);
		}
	try
		{
		*l = ConvertToLongOrThrow(s);
		}
	catch (...)
		{
		return false;
		}
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Slow (catches exception if false)
*/
bool IsAnUnsigned(const std::string &s, UInt * u) 
	{
	if (u == NULL)
		{
		UInt nu;
		return IsAnUnsigned(s, &nu);
		}
	try
		{
		*u = ConvertToUnsignedOrThrow(s);
		}
	catch (...)
		{
		return false;
		}
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
| 	Slow (catches exception if false)
*/
bool IsADouble(const std::string &s, double * d) 
	{
	if (d == NULL)
		{
		double nd;
		return IsADouble(s, &nd);
		}
	try
		{
		*d = ConvertToDoubleOrThrow(s);
		}
	catch (...)
		{
		return false;
		}
	return true;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Returns true if the stored string is a non-case-sensitive copy of the argument `s'. Note: will return true if both the
| stored string and `s' are empty strings.
*/
bool EqualsCaseInsensitive(
  const string &f,
  const string &s)   /* the comparison string */
	{
	unsigned k;
	unsigned slen = (unsigned)s.size();
	unsigned tlen = (unsigned)f.size();
	if (slen != tlen)
		return false;
	for (k = 0; k < tlen; ++k)
		{
		if ((char)toupper(f[k]) != (char)toupper(s[k]))
			return false;
		}
	return true;
	}



/*--------------------------------------------------------------------------------------------------------------------------
| Checks to see if the stored string begins with upper case letters and, if so, returns all of the contiguous capitalized
| prefix. If the stored string begins with lower case letters, an empty string is returned.
*/
string UpperCasePrefix(const string &s)
	{
	string x;
	unsigned i = 0;
	while (i < s.length() && isupper(s[i]))
		x << s[i++];
	return x;
	}


/*--------------------------------------------------------------------------------------------------------------------------
| Converts the stored string to a long using the standard C function strtol, throwing NxsX_NotANumber if the conversion 
| fails (or  NxsX_NumberIsTooLarge or NxsX_NumberIsTooSmall if the number is out of bounds)
*/
long ConvertToLongOrThrow(const string &s, ncl::IntegerConversion ) X_SPEC_THROW(NxsX_NotANumber)
	{
	long d = ConvertEntireStr<long>(s);
	if (d == LONG_MAX)
		throw NxsX_NotANumber();
	return d;
	}


/*--------------------------------------------------------------------------------------------------------------------------
| Converts the stored string to a double using the standard C function strtod, throwing NxsX_NotANumber if the conversion
| fails. throws NxsX_NumberIsTooLarge or NxsX_NumberIsTooSmall if the number is out of bounds
*/
double ConvertToDoubleOrThrow(const std::string & s) X_SPEC_THROW(NxsX_NotANumber)
	{
	double d = ConvertEntireStr<double>(s);
	if (d == DBL_MAX)
		throw NxsX_NotANumber();
	return d;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Transforms the vector of string objects by making them all lower case and then capitalizing the first portion of 
| them so that the capitalized portion is enough to uniquely specify each. Returns true if the strings are long enough 
| to uniquely specify each. Horrendously bad algorithm, but shouldn't be called often.
*/
bool SetToShortestAbbreviation(
  VecString  & strVec,		/* vector of string objects */
  bool			allowTooShort)  /* */
	{
	VecString upperCasePortion;
	unsigned long i;
	for (i = 0; i < strVec.size(); ++i)
		{
		// Change the next string to lower case
		//
		ToLower(strVec[i]);

		unsigned prefLen = 0;
		string pref;

		if (prefLen >= strVec[i].size())
			return false;
		pref << (char) toupper(strVec[i][prefLen++]);
		bool moreChars = true;

		// Keep adding letters from the current string until pref is unique.
		// Then add this pref to upperCasePortion (vector of previous prefs)
		//
		for (;moreChars;)
			{
			unsigned long prevInd = 0;
			for (; prevInd < upperCasePortion.size(); ++prevInd)
				{
				if (pref == upperCasePortion[prevInd])
					{
					//      Conflict  - both abbreviations need to grow
					//
					if (prefLen >= strVec[i].size())
						{
						if (allowTooShort)
							{
							if (prefLen < strVec[prevInd].size())
								upperCasePortion[prevInd] << (char) toupper(strVec[prevInd][prefLen]);
							moreChars = false;
							break;
							}
						else
							return false;
						}
					pref << (char) toupper(strVec[i][prefLen]);
					if (prefLen >= strVec[prevInd].size())
						{
						if (allowTooShort)
							{
							prevInd = 0;
							++prefLen;
							break;
							}
						else
							return false;
						}
					upperCasePortion[prevInd] << (char) toupper(strVec[prevInd][prefLen++]);
					prevInd = 0;
					break;
					}
				else
					{
					unsigned j;
					for (j = 0; j < prefLen; ++j)
						{
						if (pref[j] != upperCasePortion[prevInd][j])
							break;
						}
					if (j == prefLen)
						{
						//      pref agrees with the first part of another abbreviation, lengthen it.
						//
						if (prefLen >= strVec[i].size())
							{
							if (allowTooShort)
								{
								moreChars = false;
								break;
								}
							else
								return false;
							}
						pref << (char) toupper(strVec[i][prefLen++]);
						break;
						}
					}
				}
			if (prevInd == upperCasePortion.size() || !moreChars)
				{
				// Made it all the way through with no problems, add this 
				// prefix as command i's upper case portion
				//
				upperCasePortion.push_back(pref);
				break;
				}
			}
		}

	for (i = 0; i < strVec.size(); ++i)
		{
		for (unsigned j = 0; j < upperCasePortion[i].size(); ++j)
			strVec[i][j] = upperCasePortion[i][j];
		}

	return true;
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Returns a vector of string objects that match the entire `testStr'.
*/
VecString GetVecOfPossibleAbbrevMatches(
  const string  & testStr,	       /* string to match */
  const VecString & possMatches)   /* vector of possible matches */
	{
	VecString matches;
	for (unsigned i = 0; i < possMatches.size(); ++i)
		{
		if (Abbreviates(testStr, possMatches[i], ncl::kStringNoRespectCase))
			matches.push_back(possMatches[i]);
		}
	return matches;
	}

void FillVectorWithNumbers(VecString &v, unsigned fromN, unsigned endN)
	{
	v.reserve(v.size() + endN-fromN);
	string s;
	for (; fromN < endN; ++fromN)
		{
		s.clear();
		s << fromN;
		v.push_back(s);
		}
	}

/*--------------------------------------------------------------------------------------------------------------------------
| 	Returns true if the character is a space, tab, newline or < 7.
*/
inline bool IsNexusWhitespace(char ch)
	{
	return (ch == ' ' || ch == '\n' || ch == '\t' || ch < 7);
	}

const char *gIllegalChars = "()[]{}/\\,;:=*\'\"`<>^";
	
/*--------------------------------------------------------------------------------------------------------------------------
|	Function used when a single character (e.g. Gap character in the characters block) is being assigned by the user
|	and it cannot be whitespace or any one of the following ()[]{}/\\,;:=*\'\"`<>^
|	any of these characters will generate an NxsException.
*/
bool IsLegalNexusChar(
  char testC) 	/* character that the user has chosen */
	{
	return !(IsNexusWhitespace(testC) || strchr(gIllegalChars, testC) != NULL);
	}

void ThrowIllegalNexusChar(
  char 			/*testC*/, 	/* character that the user has chosen */
  const char  * name) /* name of the field that character is being assigned to (e.g. match, missing, gap) used in the NxsException if the character is illegal */ 
	{
	string s;
	s << name << " character cannot be whitespace (including _) or any of the following:  " << gIllegalChars;
	throw NxsException(s);	
	}

#if ! defined(NDEBUG) && defined(STOP_DEBUGGER_ON_TRIPPED_ASSERTS) //@pol -> mth defined(STOP...) ok? Was just STOP... before
	/*----------------------------------------------------------------------------------------------------------------------
	|	Called when a NXS_ASSERT fails.  Allows the debugger to stop and make the stack visible before the real assert is
	|	tripped
	*/
	bool	NxsAssertionFailed(const char *f, unsigned l)
		{
		int i = 0;
		++i;
			{
			std::ofstream failedAssertStream("LastFailedAssertion.txt");
			failedAssertStream << f << l << std::endl;
			failedAssertStream.close();
			}
		return i == 1;
		}
#endif
#if defined (BOOST_ENABLE_ASSERT_HANDLER) && 0 //Deprecated
	/*----------------------------------------------------------------------------------------------------------------------
	|	Defined whenever BOOST_ENABLE_ASSERT_HANDLER is defined.  This function is called when an assert trips in Boost 
	|	code.  It simply redirects the call to NxsAssertionFailed
	*/
	void boost::assertion_failed(char const * /*expr*/, char const * /*function*/, char const * file, long line)
		{
		NxsAssertionFailed(file, (unsigned) line);
		}
#endif


	// MakeStrPrintF and StrPrintF are almost identical (MakeStrPrintF would ideally be a wrapper around StrPrintF)
	// but the fact that they do multiple calls to va_start and va_end make it difficult to make one call the other
	//	The solution (for now at least) is to define the body of the function as a huge macro NCL_STR_PRINT_F_MACRO_HACK
	//
	// Comments for the NCL_STR_PRINT_F_MACRO_HACK
	// see  http://www.alphazed.co.uk/programming/vsnprintf.php for a discussion of problems with
	// the return value of vsnprintf.  Basically  we can trust that there was no truncation
	// iff -1 < nAdded < kInitialBufferSize - 1
	//
	//	Comment on the terminating for(;;)
	//  given that the last vsnprintf failed, we either double the buffer size or 
	//  increase the buffer size to the previous outsize (with padding for \0 and avoiding
	//  return code ambiguity.

#if defined(_MSC_VER)
#   define VSN_PRINT_F_NAME _vsnprintf
#elif defined (__MWERKS__)
#   define VSN_PRINT_F_NAME  std::vsnprintf
#else
#   define VSN_PRINT_F_NAME vsnprintf
#endif

#define NCL_STR_PRINT_F_MACRO_HACK  	\
	const int kInitialBufferSize = 512; \
	char buf[kInitialBufferSize]; \
	va_list argList; \
	va_start(argList, formatStr); \
	int nAdded = VSN_PRINT_F_NAME(buf, kInitialBufferSize, formatStr, argList); \
	va_end(argList); \
	if (nAdded > -1 && nAdded < kInitialBufferSize - 1) \
		{ \
		str.append(buf); \
		return str; \
		} \
	int nextSizeGuess = kInitialBufferSize; \
	int outsize = nAdded; \
	typedef boost::scoped_array<char> ScopedCharArrPtr; \
	for (;;)  \
		{  \
		nextSizeGuess = (outsize < nextSizeGuess + 1 ? nextSizeGuess * 2 : outsize + 2); \
		ScopedCharArrPtr dynBuffer(new char[nextSizeGuess]); \
		va_start(argList, formatStr); \
		outsize = VSN_PRINT_F_NAME(dynBuffer.get(), (unsigned long) nextSizeGuess, formatStr, argList);\
		va_end(argList); \
		if (outsize > -1 && outsize < nextSizeGuess - 1) \
			{ \
			str.append(dynBuffer.get()); \
			return str; \
			} \
		} 

std::string MakeStrPrintF(const char * formatStr, ...)
	{
	std::string str;
	NCL_STR_PRINT_F_MACRO_HACK
	}

/*--------------------------------------------------------------------------------------------------------------------------
| Appends a printf-style formatted string onto the end of this string and returns the number of characters added to the 
| string. For example, the following code would result in the string s being set to "ts-tv rate ratio = 4.56789":
|>
| double kappa = 4.56789;
| string s;
| StrPrintF(s, "ts-tv rate ratio = %.5f", kappa);
|>
*/
std::string & StrPrintF(string & str,
  const char * formatStr,	/* the printf-style format string */
  ...)					  /* other arguments referred to by the format string */
	{
	NCL_STR_PRINT_F_MACRO_HACK
	}
