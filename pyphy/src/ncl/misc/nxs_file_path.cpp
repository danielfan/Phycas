#include "phycas/force_include.h"
#include <iostream>
#include <fstream>
#include "ncl/nxs_defs.hpp"
#include "ncl/misc/nxs_file_path.hpp"
#include "ncl/misc/algorithm_extensions.hpp"
#include "ncl/nxs_exception.hpp"
#include "ncl/output/nxs_output.hpp"
#if defined (NCL_SUPPORT_DIR_LIST)
#   include <dirent.h>
#endif
#define SUPPORT_ABSOLUTE_PATHS

using std::list;
using std::ifstream;
using std::ofstream;
using std::string;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::strchr;
#endif
#if defined (NCL_SUPPORT_DIR_LIST)
	VecString NxsGetDirList(const string &dirName)
		{
		VecString retVec;
		if (NxsCheckForDir(dirName))
			{
			NxsInFilePath infP(dirName);
			const char *const natDirName = infP.GetNative().c_str();
			DIR * dir_p;
			dir_p = opendir(natDirName);
			struct dirent *dir_entry_p = readdir(dir_p);
			while(dir_entry_p !=  NULL)
				{
				retVec.push_back(string(dir_entry_p->d_name));
				dir_entry_p = readdir(dir_p);
				}
			closedir(dir_p);
			}
		return retVec;
		}
#endif	

NCL_TIME_T NxsFilePath::NxsGetModTime(const string &fileName)
	{
	NxsInFilePath infP(fileName);
	const char *const natFileName = infP.GetNative().c_str();
	struct NCL_STAT_STRUCT dStat;
	if (NCL_STAT(natFileName, &dStat) != 0)
		return NCL_TIME_T(0);
	return dStat.st_mtime;	/* time of last access */
	}
	
bool NxsMkDir(const string &dirName, NCL_MODE_T dirMode)
	{
	NxsInFilePath infP(dirName);
	const char *const natDirName = infP.GetNative().c_str();
	struct NCL_STAT_STRUCT dStat;
	return((NCL_STAT(natDirName, &dStat) != 0) && (NCL_MKDIR(natDirName, dirMode) == 0));
	}
	
bool NxsCheckForDir(const string &dirName)
	{
	NxsInFilePath infP(dirName);
	const char *const natDirName = infP.GetNative().c_str();
	struct NCL_STAT_STRUCT dStat;
	if (NCL_STAT(natDirName, &dStat) != 0)
		return false;
	return (NCL_S_ISDIR(dStat.st_mode));
	}
	
bool NxsCheckForOrMkDir(const string &dirName, NCL_MODE_T dirMode)
	{
	if (!NxsCheckForDir(dirName))
		return NxsMkDir(dirName, dirMode);
	return true;
	}

void NxsCheckForOrMkDirOrThrow(const string &dirName, NCL_MODE_T dirMode)
	{
	if (!NxsCheckForOrMkDir(dirName, dirMode))
		{
		string msg;
		msg << "Could not create the directory \"" << dirName << '\"';
		throw NxsException(msg);
		}
	}

void NxsMkDirOrThrow(const string &dirName, NCL_MODE_T dirMode)
	{
	if (!NxsMkDir(dirName, dirMode))
		{
		string msg;
		msg << "Could not create the directory \"" << dirName << '\"';
		throw NxsException(msg);
		}
	}

unsigned NxsFindIndexOfNextOrderedFile(
  const string &substitutionString, /*every * will in the substitutionString will be replaced with an index until a file is not found */
  unsigned firstIndexToCheck)
	{
	list<string> constPortion;
	split<char, string>(substitutionString, '*', &constPortion);
	if (constPortion.size() < 2)
		return 0;
	unsigned testVal = firstIndexToCheck;
	string testStr;
	for (;;)
		{
		list<string>::const_iterator cpIt = constPortion.begin();
		testStr = *cpIt++;
		for (; cpIt != constPortion.end(); ++cpIt)
			testStr << testVal << *cpIt;
		ifstream testStream;
		if (!NxsOpenInFile(testStr, testStream, false))
			return testVal;
		++testVal;
		}
	}

unsigned NxsFindIndexOfNextOrderedDir(
  const string &substitutionString, /*every * will in the substitutionString will be replaced with an index until a file is not found */
  unsigned firstIndexToCheck)
	{
	list<string> constPortion;
	split<char, string>(substitutionString, '*', &constPortion);
	if (constPortion.size() < 2)
		return 0;
	unsigned testVal = firstIndexToCheck;
	string testStr;
	for (;;)
		{
		list<string>::const_iterator cpIt = constPortion.begin();
		testStr = *cpIt++;
		for (; cpIt != constPortion.end(); ++cpIt)
			testStr << testVal << *cpIt;
		if (!NxsCheckForDir(testStr))
			return testVal;
		++testVal;
		}
	}


bool IsLegalFileOrDirChar(char c);
bool IsLegalFileOrDirWord(const string &s, const DirDescriptionFormat);
#if (NEW_PATH_STRING_IMPLEMENTATION)
	bool IsLegalUNIXFileOrDirWord(const string &s);
	bool IsLegalMacFileOrDirWord(const string &s);
	bool IsLegalDOSFileOrDirWord(const string &s);

	bool IsLegalUNIXFileOrDirWord(const string &s)
		{
		//@ THERE IS A LOT TO DO HERE (what are the restrictions on characters in UNIX file/dir names?)
		if (s.empty())	
			return false;
		string::const_iterator sIt = s.begin();
		if (*sIt == '/' || *sIt =='~')
			return false; //@ THIS is legal, but we are currently only supporting relative paths
		char prevChar = *sIt;
		++sIt;
		for (; sIt != s.end(); ++sIt)
			{
			if (*sIt == '/' && prevChar == '/')
				return false;
			prevChar = *sIt;
			}	
		return true;
		}
		
	bool IsLegalMacFileOrDirWord(const string &)
		{
		NXS_ASSERT(false);
		return false;
		}
		
	bool IsLegalDOSFileOrDirWord(const string &)
		{
		NXS_ASSERT(false);
		return false;
		}
		
#endif

bool NxsOpenOutFile(const string &fn, ofstream &oStream, OutModes o, bool canQuery)
	{
	NxsOutFilePath ofp(fn);
	if (!canQuery)
		ofp.SetQueryUserOnError(false);
	if (o == kReplaceOutMode)
		ofp.SetReplace();
	else if (o == kAppendOutMode)
		ofp.SetAppend();
	return ofp.Open(oStream);
	}
	
	
bool NxsOpenInFile(const string &fn, ifstream &iStream, bool canQuery)
	{
	NxsInFilePath ofp(fn);
	if (!canQuery)
		ofp.SetQueryUserOnError(false);
	return ofp.Open(iStream);
	}
	
const char internalDividerChar = '/';

//const char internalDividerChar = ':';
	
#if defined(NXS_USE_MAC_DIR_STYLE)
	const char nativeSeparator = ':';
#elif defined(NXS_USE_WINDOWS_DIR_STYLE)
	const char nativeSeparator = '\\';
#else
	const char nativeSeparator = '/';
#endif

//	this setting controls what will be accepted from the user.  
//	currently we're expecting the user to give directories in unix format
DirDescriptionFormat NxsFilePath::dirStringInterpretStyle = kUNIXDirStyle;

inline bool IsLegalFileOrDirChar(char c)
	{
	const char *illegal = " \t\n\\/:";
	return (strchr(illegal, c) == NULL);
	}
	
inline bool IsLegalFileOrDirWord(const string &s, const DirDescriptionFormat format) 
	{
#	if (NEW_PATH_STRING_IMPLEMENTATION)
		if (format == kUNIXDirStyle)
			return IsLegalUNIXFileOrDirWord(s);
		else if (format == kMacDirStyle)
			return IsLegalMacFileOrDirWord(s);
		return IsLegalDOSFileOrDirWord(s);
#	else
		return (!s.empty() &&  all(s.begin(), s.end(), IsLegalFileOrDirChar));
#	endif
	}

bool IsNull(char c);
bool IsDivider(char c);

inline bool IsNull(char c)
	{
	return c == '\0';
	}

inline bool IsDivider(char c)
	{
	return c == internalDividerChar;
	}
	
inline bool IsUpDirNotation(const string d)
	{
	return (d.length() == 2 && d[0] == '.' && d[1] == '.');
	}

void NxsFilePath::ShortenDirStack(list<string> &dirStack)
	{
	list<string>::iterator currIt = dirStack.begin();
#	if defined(SUPPORT_ABSOLUTE_PATHS) && !defined(NXS_USE_WINDOWS_DIR_STYLE) //POL, can you do the windows side of ABSOLUTE_PATHS
		if (currIt != dirStack.end() && currIt->length() == 0)
			++currIt;
#	endif
	while (currIt != dirStack.end())
		{
		if (currIt->length() < 1 || (currIt->length() == 1 && (*currIt)[0] == '.')) // remove empty elements and this (.) directories
			currIt = dirStack.erase(currIt);
		else 
			{
			if (IsUpDirNotation(*currIt))
				{
				if (currIt != dirStack.begin())
					{
					// shorten  "d/c/../"  to "d/"
					--currIt;
					if (IsUpDirNotation(*currIt))
						++currIt; // we can't erase ../../ 
					else
						{
						currIt = dirStack.erase(currIt);	//erase the previous directory name
						currIt = dirStack.erase(currIt);	// erase the .. directory name
						}
					}
				}
			++currIt;
			}
		}
	}
	
void NxsFilePath::TranslateDirStack(const list<string> &dirStack) const
	{
	list<string>::const_iterator dIt = dirStack.begin();
	const list<string>::const_iterator endIt = dirStack.end();
	
#	if defined(NXS_USE_UNIX_DIR_STYLE)
 		nativePath.clear();
 		if (dirStack.empty())
 			return;
 		if (dirStack.size() == 1 )
 			{
 			nativePath = (IsUpDirNotation(*dIt) ? "../" : *dIt);
 			return;
 			}
 		bool firstWord = true;
 		
 		if (dIt->empty())  // first character was /
 			{
 			// absolute path
 			firstWord = false;
 			}
 		bool lastWasUp = false;	
 		for (; dIt != endIt; ++dIt)
 			{
 			if (!firstWord)
 				nativePath << '/';
 			firstWord = false;
 			lastWasUp = IsUpDirNotation(*dIt);
 			if (lastWasUp)
 				nativePath << "..";
 			else
 				nativePath << *dIt;
 			}
 		if (lastWasUp)
 			nativePath  << '/';
#	elif defined(NXS_USE_MAC_DIR_STYLE)
 		nativePath.clear();
 		bool wordWritten = false;
 		if (dirStack.size() == 1 && !IsUpDirNotation(*dIt))
 			{
 			nativePath = *dIt;
 			return;
 			}
 		for (; dIt != endIt; ++dIt)
 			{
 			if (IsUpDirNotation(*dIt))
 				nativePath << ':';
 			else
 				{
 				wordWritten = true;
 				nativePath << ':' << *dIt;
 				}
 			}
 		if (!wordWritten)
 			nativePath  << ':';
#	elif defined(NXS_USE_WINDOWS_DIR_STYLE)
		//@POL Mark, submitting a path to _stat to see if it is a directory will fail if the path
		// ends in a backslash ('\\') character. What I'm not sure of is whether elsewhere you depend
		// on having a slash on the end (e.g. when combining with a file name)
		//
 		nativePath.clear();
 		for (; dIt != endIt; ++dIt)
 			{
 			if (IsUpDirNotation(*dIt))
 				nativePath += "..";
 			else
 				{
 				nativePath << '\\' << *dIt;
 				}
 			}
#	endif
	}

void NxsFilePath::CreateNative() const
	{
	nativePath = path;
	if (nativePath.empty())
		return;
#	if (NEW_PATH_STRING_IMPLEMENTATION)
		if (!nativePath.empty() && nativePath[0] != nativeSeparator)
			{
			list<string> dirStack;
			split<char, string>(nativePath, '/', &dirStack);
			ShortenDirStack(dirStack);
			TranslateDirStack(dirStack);
			}
#	else
		replace_if(nativePath.begin(), nativePath.end(), IsDivider, nativeSeparator);
#	endif
	}

bool NxsFilePath::IsDirectory() const
	{
	if (!this->Exists())
		return false;
	const char *const natDirName = this->GetNative().c_str();
	struct NCL_STAT_STRUCT dStat;
	if (NCL_STAT(natDirName, &dStat) != 0)
		{
		std::string msg;
		msg += "stat call failed for ";
		msg += natDirName;
		throw NxsException(msg);
		}
	bool retCode = NCL_S_ISDIR(dStat.st_mode);
	//std::cerr << "IsDirectory is returning "<< retCode << " for " << natDirName <<std::endl;
	return retCode;
	}
	
NxsFilePath::NxsFilePath(bool shouldBeDir)
	:sourceOfName(kProgrammer),
	isAbsolute(false),
	isDirty(true),
	isLegal(false),
	path(),
	pathAsEntered(),
	nativePath(),
	queryUserOnError(true)
	{
	assert(!shouldBeDir);	// file paths without filenames aren't supported yet
	}

NxsFilePath::NxsFilePath(const string &s, bool shouldBeDir )
	:sourceOfName(kProgrammer),
	isAbsolute(false),
	isDirty(true),
	isLegal(false),
	path(),
	pathAsEntered(),
	nativePath(),
	queryUserOnError(true)
	{
	assert(!shouldBeDir);	// file paths without filenames aren't supported yet
	if (ReadNewPathString(s))
		CreateNative();
	}
	
NxsFilePath::OpenOutputReturnCode NxsFilePath::OpenOutput(
  ofstream &outF, 
  bool replaceFile, 
  bool appendToFile) const
  	{
  	outF.close();
	outF.clear();
	if (Exists())
		{
		if (replaceFile)
			{
			if (appendToFile)
				return NxsFilePath::kBothReplaceAndAppend;
			outF.open(GetNative().c_str());
			if (outF.good())
				return NxsFilePath::kOpenAndReplacing;
			}
		else if (appendToFile)
			{
			outF.open(GetNative().c_str(), std::ios::app);
			if (outF.good())
				return NxsFilePath::kOpenAndAppending;
			}
		else
			return NxsFilePath::kNeitherReplaceOrAppend;
		}
	else if (IsLegalFileName())
		{
		outF.open(GetNative().c_str());
		if (outF.good())
			return NxsFilePath::kOpenedNewFile;
		}
	outF.close();
	outF.clear();
	return NxsFilePath::kCouldNotOpen;
  	}
  	
string NxsOutFilePath::ComposeReplaceAppendErrorMsg() const
	{
	assert ((GetReplace() && GetAppend()) || (!GetReplace() && !GetAppend()));
	
	string s;
	s << "The file " << GetFullName() << " already exists, and ";
	if (GetReplace() && GetAppend())
		s << "both replace and append are ";
	else
		s << "neither replace nor append is ";
	s << "set to \"true\"";
	return s;
	}
	
NxsFilePath::QueryResult NxsOutFilePath::QueryUserForReplaceAppend()
	{
	if (NxsOutput::GetUserQueryPtr() == NULL || !queryUserOnError)
		return kCouldNotQuery;
	userCorrection = kNoCorrection;
	unsigned userResponse = NxsOutput::GetUserQueryPtr()->UserChoice( "Error Opening File", ComposeReplaceAppendErrorMsg(),"Replace|Append|Cancel", 1, 1);
	if (userResponse == 2)
		return kUserCancelled;
	userCorrection = (userResponse == 0 ? kCorrectedToReplace : kCorrectedToAppend);
	return kProblemWasFixed;
  	}

NxsFilePath::QueryResult NxsOutFilePath::QueryUserForNewFile(const string &purposeOfFile)
	{
	NxsUserQuery * userQuery = NxsOutput::GetUserQueryPtr();
	if (userQuery == NULL || !queryUserOnError)
		return kCouldNotQuery;
	
	const string n = (purposeOfFile.empty() ? path : purposeOfFile);
	NxsOutFilePath p = userQuery->GetOutputFilePath(n);
	const std::string s = p.GetFullName(); //@temp should use p not s
	if (s.empty())
		return kUserCancelled;
	if (!ReadNewPathString(s))
		return QueryUserForNewFile(purposeOfFile);
	return kProblemWasFixed;
	}
	
NxsFilePath::QueryResult NxsInFilePath::QueryUserForNewFile(const string &purposeOfFile)
	{
	NxsUserQuery * userQuery = NxsOutput::GetUserQueryPtr();
	if (userQuery == NULL || !queryUserOnError)
		return kCouldNotQuery;
	
	const string n = (purposeOfFile.empty() ? path : purposeOfFile);
	NxsInFilePath p = userQuery->GetInputFilePath(n);
	const std::string s = p.GetFullName(); //@temp should use p not s
	if (s.empty())
		return kUserCancelled;
	if (!ReadNewPathString(s))
		return QueryUserForNewFile(purposeOfFile);
	return kProblemWasFixed;
	}
	  	
bool NxsOutFilePath::Open(ofstream &outF, string purposeOfFile)
	{
	NxsFilePath::OpenOutputReturnCode r = OpenOutput(outF, GetReplace(), GetAppend());
	wasOpened = isAppending = isReplacing = false;
	if (r == NxsFilePath::kCouldNotOpen || r == NxsFilePath::kBothReplaceAndAppend || r == NxsFilePath::kNeitherReplaceOrAppend)
		{
		if (r == NxsFilePath::kBothReplaceAndAppend || r == NxsFilePath::kNeitherReplaceOrAppend)
			{
			NxsFilePath::QueryResult res = QueryUserForReplaceAppend();
			if (res != kProblemWasFixed)
				return false;
			}
		else 
			{
			NxsFilePath::QueryResult fres = QueryUserForNewFile(purposeOfFile);
			if (fres != kProblemWasFixed)
				return false;
			}
		return Open(outF, purposeOfFile);
		}
	wasOpened = true;
	if (r == NxsFilePath::kOpenAndReplacing)
		isReplacing = true;
	if (r == NxsFilePath::kOpenAndAppending)
		isAppending = true;
	return true;
	}

bool NxsFilePath::Exists() const
	{
	ifstream testIn;
	bool retVal = OpenInput(testIn);
	testIn.close();
	return retVal;
	}

bool NxsFilePath::OpenInput(ifstream &inFStream) const
	{
	inFStream.close();
	inFStream.clear();
	if (!IsLegalFileName())
		return false;
	inFStream.open(GetNative().c_str());
	if (inFStream.good())
		return true;
	inFStream.close();
	inFStream.clear();
	return false;
	}
#if (NEW_PATH_STRING_IMPLEMENTATION)
	bool NxsFilePath::TranslateToInternalRep(string &)
		{
		if (dirStringInterpretStyle == kUNIXDirStyle)
			return true; //we are using the UNIX representation internally
		NXS_ASSERT(false);
		return false;
		}
#endif

bool NxsFilePath::ReadNewPathString(const string &s)
	{
	pathAsEntered = path = s;
	nativePath.clear();
	isDirty = true;
#	if (NEW_PATH_STRING_IMPLEMENTATION)
		if (!IsLegalFileName())
			return false;
		return TranslateToInternalRep(path);
#	else
		return IsLegalFileName();
#	endif
	}

bool NxsFilePath::IsLegalFileName() const
	{
	if (isDirty)
		isLegal  = IsLegalFileOrDirWord(GetFileName(), dirStringInterpretStyle);
	return isLegal;
	}

string NxsFilePath::GetFileName() const
	{
	string::size_type lastDivider = path.find_last_of(internalDividerChar);
	if (lastDivider == string::npos)
		return path;
	string fName;
	++lastDivider;
	string::size_type fLen = path.length() - lastDivider;
	fName.assign(path, lastDivider, fLen);
	return fName;
	}

bool NxsInFilePath::Open(ifstream &outF, string purposeOfFile)
	{
	if (OpenInput(outF))
		return  true;
	if (QueryUserForNewFile(purposeOfFile) != kProblemWasFixed)
		return false;
	return Open(outF, purposeOfFile);
	}
