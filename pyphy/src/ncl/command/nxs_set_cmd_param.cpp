#include "phycas/force_include.h"
#include "ncl/nxs_defs.hpp"
#include "ncl/command/nxs_set_cmd_param.hpp"
using std::pair;
using std::string;
unsigned GetPositiveInt(const string &nextTok, unsigned maxVal);

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
NxsIndexSetCmdOption::NxsIndexSetCmdOption(
  const string &n, 
  NxsIndexSet *manipVal, 
  const string &def, 
  MaxNumberSource maxP,
  LabelDecoder lab , 
  SetDecoder setP, 
  bool  shouldExpectName, 
  bool  persist, 
  string prefixWord,
  CmdPermissionLevel pLevel)
  	:SimpleCmdOptionInterface<NxsIndexSet>(n, manipVal, *manipVal, true, persist, pLevel),
  	labelProvider(lab),
  	setProvider(setP),
  	maxNumberProv(maxP),
  	expectingName(shouldExpectName),
  	allowNumbers(true),
  	pref(prefixWord),
  	defStr(def)
  	{
  	if (maxNumberProv)
  		maxIndex = maxNumberProv();
  	else
  		maxIndex = UINT_MAX;
  	}

VecString NxsIndexSetCmdOption::GetValidArgument()
	{
	RefreshMaxIndex();
	VecString v(2, "1");
	string s;
	s << maxIndex;
	v[1] = s;
	return v;
	}
/*----------------------------------------------------------------------------------------------------------------------
|	Calls the label to index provider to get the index associated with the name nextTok.  UINT_MAX is returned if the 
|	name isn't found.
*/
unsigned NxsIndexSetCmdOption::GetIndexFromLabel(
  const string & nextTok)
  	{
	if (labelProvider)
		{
		unsigned codedlabelInd = labelProvider(nextTok);
		if (codedlabelInd != UINT_MAX)
			return codedlabelInd;
		}
	return UINT_MAX;
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Tries to translate the string to an index.  First calls the label provider, if that doesn't work it uses
|	GetPositiveInt to cast the string to an index.
|	Assumes that RefreshMaxIndex has  been called
|	Can throw NxsX_MaxExceeded
*/
unsigned NxsIndexSetCmdOption::GetIndexFromToken(
  const string & nextTok)
	{
	unsigned i = GetIndexFromLabel(nextTok);
	if (i != UINT_MAX || !allowNumbers)
		return i;
	i = GetPositiveInt(nextTok, maxIndex + 1);
	if (i != UINT_MAX)
		--i;
	return  i;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads the vector notation for a set.  Assumes that the equals has NOT been read and maxIndex has been refreshed
*/
bool NxsIndexSetCmdOption::ReadVectorFormOfSet(
  NxsToken &token)
	{
	if (token.GetTokenReference() != '=')
		return FlagError(unrecognized, "=");
	currentValue->clear(); 
	
	for (unsigned ind = 0; ind <= maxIndex;ind++)
		{
		token.ReadSingleCharacter();
		if (token.GetTokenReference() == '1')
			currentValue->insert(ind);
		else if (token.GetTokenReference() != '0')
			return FlagError(unrecognized, "0 or 1");
		}
	++token;	
	return true;
	}  	

inline unsigned NxsIndexSetCmdOption::ReadStride(NxsToken &token, unsigned maxStride)
	{
	unsigned stride;
	try	
		{
		stride = GetPositiveInt(token.GetTokenReference(), maxStride);
		}
	catch (NxsX_MaxExceeded & )
		{
		stride = UINT_MAX;
		}
	return stride;
	}
/*----------------------------------------------------------------------------------------------------------------------
|	if nextTok holds an int > 0 and  <= maxVal, it will be returned as an int
|	otherwise UINT_MAX is returned if the value isn't an positive int
|	will throw NxsIndexSetCmdOption::NxsX_MaxExceeded if the string is > maxVal
*/
unsigned GetPositiveInt(const string &nextTok, unsigned maxVal)
	{
	UInt u;
	if (IsAnUnsigned(nextTok, &u) && u != UINT_MAX && u != 0)
		{
		if (u > maxVal)	
			throw NxsIndexSetCmdOption::NxsX_MaxExceeded(maxVal, u);
		return u;
		}
	return UINT_MAX;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reads the charset by calling either ReadVectorFormOfSet or ReadStandardFormOfSet.
|	If this NxsIndexSetCmdOption was told to expect a name, then it will read a name and the = sign before starting to 
|	read the set description.
*/
bool NxsIndexSetCmdOption::ReadValue(
 NxsToken &token,	/* the stream of tokens that are being read */
 bool equalsAlreadyRead) /* true if the equals sign has already been removed from the token stream (or is NOT expected) */
	{
	currentValue->clear();
	RefreshMaxIndex();
	string setName;
	if (expectingName)
		{
		if (!equalsAlreadyRead && !EatEqualsThenAdvance(token))
			return false;
		if (token.IsPunctuationToken())
			return FlagError(illegal_name, "it is a reserved punctuation character");
		else if (GetIndexFromLabel(token.GetTokenReference()) != UINT_MAX)
			return FlagError(illegal_name,"it is a label");
		unsigned x = GetPositiveInt(token.GetTokenReference(), INT_MAX);
		if (x != UINT_MAX && (x==0 || x-1 <= maxIndex))
			return FlagError(illegal_name,"it is a number that could also be a valid label");
		setName = token.GetTokenReference();
		++token;
		equalsAlreadyRead = false;
		}
	if (!equalsAlreadyRead)
		{
		if (token.GetTokenReference() == '(')
			{
			++token;
			
			bool vecMode = false;
			if (token.Equals("VECTOR"))
				{
				vecMode = true;
				++token;
				}
			else if (token.Equals("STANDARD"))
				++token;
			else
				return FlagError(unrecognized, "STANDARD or VECTOR");
			if (!EatWordThenAdvance(token, ")", ncl::kStringRespectCase))
				return false;
			if (vecMode)
				{
				if (!ReadVectorFormOfSet(token))
					return false;
				if (expectingName)
					currentValue->SetName(setName);
				return WasValidRead();
				}
			}
		if (!EatEqualsThenAdvance(token))
			return false;
		}
	if (!ReadStandardFormOfSet(token))
		return false;
	if (expectingName)
		currentValue->SetName(setName);
	return WasValidRead();
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Simply verifies that the set is not empty and that the indices are under the maximum value
*/
bool NxsIndexSetCmdOption::IsCurrentlyValid()
	{
	RefreshMaxIndex();
	if (currentValue->empty())
		return FlagError(query_for_error, MakeStrPrintF( "The %s is empty", GetDisplayType().c_str()));
	else if (currentValue->GetLast() > maxIndex )
		return FlagError(too_big, MakeStrPrintF("%ud", currentValue->GetLast()));
	return true;
	}

void NxsIndexSetCmdOption::ReturnValueToDefault()
	{
	if (defStr.empty())
		currentValue->clear();
	else
		{
		NxsToken token(defStr);
		bool readCorrectly = ReadValue(token, true);
		assert(readCorrectly);
		if (!readCorrectly)
			currentValue->clear();
		}
	}
		
/*----------------------------------------------------------------------------------------------------------------------
|	Reads the standard notation for a set.  Assumes that the equals HAS been read and maxIndex has been refreshed
*/
bool NxsIndexSetCmdOption::ReadStandardFormOfSet(
  NxsToken &token)
	{
	currentValue->clear();
	bool firstElement = true;
	try {
		for (;;)
			{
			bool setFound = false;
			if (setProvider != NULL)
				{
				const NxsIndexSet *setToAdd = setProvider(token.GetTokenReference());
				if (setToAdd != NULL)
					{
					setFound = true;
					++token;
					if (token.GetTokenReference() == '\\')
						{
						++token;
						unsigned stride = ReadStride(token, UINT_MAX);
						if (stride == UINT_MAX)
							return FlagError(illegal_modulus, setToAdd->GetName());
						unsigned charCount = 1;
						for (NxsIndexSet::const_iterator sIt = setToAdd->begin(); sIt != setToAdd->end(); ++sIt)
							{
							if ((charCount++ % stride) == 0 )
								currentValue->insert(*sIt);
							}
						++token;
						}
					else
						currentValue->insert(*setToAdd);
					}
				firstElement = false;
				}
			if (!setFound)
				{
				unsigned begLabelIndex = token.GetTokenReference() == '.' ? maxIndex : GetIndexFromToken(token.GetTokenReference());
				if (begLabelIndex != UINT_MAX)
					{
					++token;
					if (token.GetTokenReference() == '-')
						{
						++token;
						unsigned endLabelIndex = token.GetTokenReference() == '.' ? maxIndex : GetIndexFromToken(token.GetTokenReference());
						if (endLabelIndex == UINT_MAX || endLabelIndex <= begLabelIndex)
							return FlagError(illegal_range, MakeStrPrintF("%d - %d", 1 + begLabelIndex,  endLabelIndex + 1));
						++token;
						if (token.GetTokenReference() == '\\')
							{
							++token;
							unsigned stride = ReadStride(token, endLabelIndex - begLabelIndex);
							if (stride == UINT_MAX)
								return FlagError(illegal_modulus, MakeStrPrintF("%d", (int)stride));
							for (; begLabelIndex <= endLabelIndex; begLabelIndex+=stride)
								currentValue->insert(begLabelIndex);
							++token;
							}
						else
							currentValue->InsertRange(begLabelIndex, endLabelIndex);
						}
					else
						currentValue->insert(begLabelIndex);
					}
				else 	//	we've encountered something that isn't a legal set name, number or label. RETURN
					return (firstElement ? FlagError(unrecognized) : true);
				}
			firstElement = false;
			}
		}
	catch (NxsX_MaxExceeded & x)
		{
		string s;
		s << x.maxV;
		if (maxNumberProv)
			FlagError(too_big_labile, s);
		else
			FlagError(too_big, s);
		}
	return false;
	}

