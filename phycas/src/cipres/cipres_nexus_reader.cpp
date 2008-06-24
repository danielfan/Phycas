/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//#include "phycas/force_include.h"
#include <iostream>
#include <fstream>
#include <map>

#include "ncl/nxsexception.h"
#include "ncl/nxsdefs.h"
#include "ncl/nxsexception.h"
#include "ncl/nxscharactersblock.h"

#include "phycas/src/cipres/cipres_nexus_reader.hpp"
#include "phycas/src/cipres/CipresDataMatrixHelper.h"
#include "phycas/src/cipres/ConfigDependentHeaders.h"
#include "phycas/src/cipres/AllocateMatrix.hpp"
#include "phycas/src/cipres/CipresAssert.hpp"
#include "phycas/src/cipres/util_copy.hpp"
#if defined (NCL_SUPPORT_DIR_LIST)
#   include <dirent.h>
#endif
#define SUPPORT_ABSOLUTE_PATHS

using std::list;
using std::ifstream;
using std::ofstream;
using std::string;
using std::map;
using std::ifstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;

#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::strchr;
#endif

	

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
 				nativePath.append("/");
 			firstWord = false;
 			lastWasUp = IsUpDirNotation(*dIt);
 			if (lastWasUp)
 				nativePath.append("..");
 			else
 				nativePath.append(*dIt);
 			}
 		if (lastWasUp)
 			nativePath.append("/");
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
			nativePath.append(":");
 			if (!IsUpDirNotation(*dIt))
 				{
 				wordWritten = true;
 				nativePath.append(*dIt);
 				}
 			}
 		if (!wordWritten)
 			nativePath.append(":");
#	elif defined(NXS_USE_WINDOWS_DIR_STYLE)
		//@POL Mark, submitting a path to _stat to see if it is a directory will fail if the path
		// ends in a backslash ('\\') character. What I'm not sure of is whether elsewhere you depend
		// on having a slash on the end (e.g. when combining with a file name)
		//
 		nativePath.clear();
 		for (; dIt != endIt; ++dIt)
 			{
 			if (IsUpDirNotation(*dIt))
 				nativePath.append("..");
 			else
 				{
 				nativePath.append("\\");
 				nativePath.append(*dIt);
 				}
 			}
#	endif
	}

void NxsFilePath::CreateNative() const
	{
	nativePath = path;
	if (nativePath.empty())
		return;
	if (!nativePath.empty() && nativePath[0] != nativeSeparator)
		{
		list<string> dirStack;
		split<char, string>(nativePath, '/', &dirStack);
		ShortenDirStack(dirStack);
		TranslateDirStack(dirStack);
		}
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
	
NxsFilePath::NxsFilePath(const string &s, bool shouldBeDir )
	:sourceOfName(kProgrammer),
	isAbsolute(false),
	isDirty(true),
	path(),
	pathAsEntered(),
	nativePath(),
	queryUserOnError(true)
	{
	PHYCAS_ASSERT(!shouldBeDir);	// file paths without filenames aren't supported yet
	if (ReadNewPathString(s))
		CreateNative();
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
	inFStream.open(GetNative().c_str());
	if (inFStream.good())
		return true;
	inFStream.close();
	inFStream.clear();
	return false;
	}

bool NxsFilePath::ReadNewPathString(const string &s)
	{
	pathAsEntered = path = s;
	nativePath.clear();
	isDirty = true;
	return true;
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
	return OpenInput(outF);
	}


CipresNative::DiscreteMatrix * GetLastDiscreteMatrix(PhycasNexusReader & nexusReader, bool convertGapsToMissing = false)
	{
	return createNativeDiscreteMatrix(nexusReader, 0L, UINT_MAX, convertGapsToMissing);
	}
CipresNative::DiscreteMatrix * createNativeDiscreteMatrix(PhycasNexusReader & nexusReader, NxsTaxaBlock * taxaBlockPtr, unsigned int charBlockIndex, bool convertGapsToMissing)
	{
	NxsCharactersBlock * cb;
	if (charBlockIndex == UINT_MAX || taxaBlockPtr == NULL)
		cb = nexusReader.getActiveCharactersBlock();
	else
		cb = nexusReader.GetCharactersBlock(taxaBlockPtr, charBlockIndex);
	if (cb == 0L)
		return NULL;
	std::vector<const NxsDiscreteDatatypeMapper *> mappers = cb->GetAllDatatypeMappers();
	if (mappers.size() > 1)
		{
		std::string m("Characters block contains more than one datatype. Phycas cannot understand such matrices.");
		nexusReader.NexusWarn(m, NxsReader::SKIPPING_CONTENT_WARNING, 0, -1, -1);
		return NULL;
		}
	if (mappers.empty() || mappers[0] == NULL)
		{
		std::string m("Characters block does not contain a matrix with a valid mapping data structure. Please report this error to the developers of Phycas.");
		nexusReader.NexusWarn(m, NxsReader::SKIPPING_CONTENT_WARNING, 0, -1, -1);
		return NULL;
		}
	return new CipresNative::DiscreteMatrix(*cb, convertGapsToMissing);
	}


void PhycasNexusReader::NexusError(const std::string &msg, file_pos pos, unsigned line, unsigned col, CmdResult , NxsBlock* )
	{
	errorMsg = msg;
	filePositionOfError = pos;
	lineOfError = line;
	columnOfError = col;
	}
	
void PhycasNexusReader::ReadNxsFilePath(NxsInFilePath & filePath)
	{
	errorMsg.clear();
	const std::string fn = filePath.GetFullName();
	if (!filePath.Exists())
		{
		errorMsg = string("The file \'");
		errorMsg.append(fn);
		errorMsg.append("\' does not exist.");
		throw NxsException(errorMsg.c_str());
		}
	if (filePath.IsDirectory())
		{
		errorMsg  = string ("\'");
		errorMsg.append(fn);
		errorMsg.append("\' is a directory, not a file.");
		throw NxsException(errorMsg.c_str());
		}
	ifstream inStream;
	if (!filePath.Open(inStream))
		{
		errorMsg = string("The file \'");
		errorMsg.append(fn);
		errorMsg.append("\' could not be opened.");
		throw NxsException(errorMsg.c_str());
		}
	string msg;
	NxsToken token(inStream);	
	Execute(token);
	const unsigned nTaxaBlocks = GetNumTaxaBlocks();
	if (nTaxaBlocks == 0)
		{
		std::string m("No taxa blocks encountered in ");
		m.append(fn);
		m.append(" No information from the file was stored.");
		NexusWarn(m, NxsReader::SKIPPING_CONTENT_WARNING, 0, -1, -1);
		activeTaxaBlockIndex = UINT_MAX;
		activeCharactersBlockIndex = UINT_MAX;
		return;
		}
	else
		{
		activeTaxaBlockIndex = nTaxaBlocks - 1;
		if (nTaxaBlocks > 1)
			{
			std::string m("More than one taxa block encountered, only the last Taxa block is active");
			NexusWarn(m, NxsReader::SKIPPING_CONTENT_WARNING, 0, -1, -1);
			}
		trees.clear();
		}
	NxsTaxaBlock * activeTB = GetTaxaBlock(activeTaxaBlockIndex);
	PHYCAS_ASSERT(activeTB);

	const unsigned nCharBlocks = GetNumCharactersBlocks(activeTB);
	if (nCharBlocks == 0)
		activeCharactersBlockIndex = UINT_MAX;
	else
		{
		activeCharactersBlockIndex = nCharBlocks - 1;
		if (nCharBlocks > 1)
			{
			std::string m("More than one Characters block encountered, only the last Characters block is active");
			NexusWarn(m, NxsReader::SKIPPING_CONTENT_WARNING, 0, -1, -1);
			}
		}

	const unsigned nTreesBlocks = GetNumTreesBlocks(activeTB);
	for (unsigned tInd = 0; tInd < nTreesBlocks; ++tInd)
		{
		NxsTreesBlock * tb = GetTreesBlock(activeTB, tInd);
		std::vector<NxsFullTreeDescription> & currBlockTrees = tb->GetProcessedTrees();
		/* PELIGROSO -- we are stealing the trees from the first trees block */
		if (trees.empty())
			trees.swap(currBlockTrees);
		else
			std::copy(currBlockTrees.begin(), currBlockTrees.end(), std::back_inserter(trees));
		tb->Reset();
		}
	}

NxsTaxaBlock * PhycasNexusReader::getActiveTaxaBlock() const
	{
	if (activeTaxaBlockIndex == UINT_MAX)
		return NULL;
	return GetTaxaBlock(activeTaxaBlockIndex);
	}

NxsCharactersBlock * PhycasNexusReader::getActiveCharactersBlock() const
	{
	NxsTaxaBlock * activeTB = getActiveTaxaBlock();
	if (activeTB == 0L)
		return 0L;
	if (activeCharactersBlockIndex == UINT_MAX)
		return 0L;
	return GetCharactersBlock(activeTB, activeCharactersBlockIndex);
	}

unsigned PhycasNexusReader::GetNChar() const
	{
	NxsCharactersBlock * cb = getActiveCharactersBlock();
	if (cb == 0L)
		return 0;
	return cb->GetNumChar();
	}

const std::vector<NxsFullTreeDescription> & PhycasNexusReader::GetTrees() const
	{
	return trees;
	}


vector<string> PhycasNexusReader::GetTaxLabels() const
	{
	NxsTaxaBlock * activeTB = getActiveTaxaBlock();
	if (activeTB == 0L)
		return vector<string>();
	return activeTB->GetAllLabels();
	}
	


