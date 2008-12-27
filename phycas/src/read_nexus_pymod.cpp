/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
|																			  |
|  This program is free software; you can redistribute it and/or modify		  |
|  it under the terms of the GNU General Public License as published by		  |
|  the Free Software Foundation; either version 2 of the License, or		  |
|  (at your option) any later version.										  |
|																			  |
|  This program is distributed in the hope that it will be useful,			  |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of			  |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the			  |
|  GNU General Public License for more details.								  |
|																			  |
|  You should have received a copy of the GNU General Public License along	  |
|  with this program; if not, write to the Free Software Foundation, Inc.,	  |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.				  |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#include <iostream>
#include <fstream>

#include <boost/python.hpp>

#include "ncl/nxsexception.h"
#include "ncl/nxstreesblock.h"

#include "ncl/nxscxxdiscretematrix.h"

#if defined(USING_NUMARRAY)
// These next three lines required, otherwise get link error "unresolved external symbol _PyArrayHandle" 
// (at least in VC7)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#	include <boost/python/numeric.hpp>
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif

using namespace boost::python;

void translate(const NxsException & e)
	{
	// Use the Python 'C' API to set up an exception object
	PyErr_SetString(PyExc_Exception, e.what());
	}

#include "ncl/nxspublicblocks.h"
#include "phycas/src/phycas_nexus_reader.hpp"

//POL 
// I moved NxsFilePath declaration/definition to nxs_file_path.hpp and nxs_file_path.cpp, respectively
// I moved PhycasNexusReader declaration/definition to phycas_nexus_reader.hpp and phycas_nexus_reader.cpp, respectively
// 

#if 0 

#if defined(NXS_USE_MAC_DIR_STYLE)
	const char nativeSeparator = ':';
#elif defined(NXS_USE_WINDOWS_DIR_STYLE)
	const char nativeSeparator = '\\';
#elif defined(NXS_USE_UNIX_DIR_STYLE)
	const char nativeSeparator = '/';
#else
#	error must define a filesystem
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Similar to the Perl split function.  Copies the original container of T objects into the list outList as 
|	smaller containers broken at the specified value (which is omitted).
|	e.g. if the incoming string "usr/home/mine" and the splitAtVal was the char '/', then the resulting list would
|	have three strings "usr" "home" and "mine"
*/
template <typename T, class ORIG_CONTAINER>
void split(const ORIG_CONTAINER &origContainer, T splitAtVal, std::list<ORIG_CONTAINER> *outList)
	{
	typedef typename ORIG_CONTAINER::const_iterator OrigCIt;
	OrigCIt begIt = origContainer.begin();
	const OrigCIt endIt = origContainer.end();
	if (begIt == endIt)
		return;	//empty container sent
	OrigCIt copyToIt;
	do	{
		copyToIt = find(begIt, endIt, splitAtVal);
		if (begIt == copyToIt)	
			outList->push_back(ORIG_CONTAINER());
		else
			{
			outList->push_back(ORIG_CONTAINER(begIt,copyToIt));
			begIt = copyToIt;
			}
		++begIt;
		}
	while (copyToIt != endIt);
	}

#if defined(_MSC_VER)
#	define NCL_STAT				_stat
#	define NCL_STAT_STRUCT		_stat
#	define NCL_DIRMODE			0
#	define NCL_MODE_T			unsigned
#	define NCL_MKDIR(a,b)		_mkdir(a)
#	define NCL_S_ISDIR(a)		((_S_IFDIR & (a)) != 0)
#	define NCL_TIME_T			time_t
#	include <direct.h>
#else
#	define NCL_STAT				stat
#	define NCL_STAT_STRUCT		stat
#	define NCL_DIRMODE			(S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP| S_IROTH | S_IXOTH)
#	define NCL_MODE_T			mode_t
#	define NCL_MKDIR(a,b)		mkdir(a,b)
#	define NCL_S_ISDIR(a)		S_ISDIR(a)
#	define NCL_TIME_T			std::time_t
#endif


enum CmdResult
	{
	kCmdSucceeded,
	kCmdFailedSilent,
	kCmdFailedGenerateMessage
	};

enum DirDescriptionFormat  // names for the 3 systems of interpretting strings as paths
	{
	kUNIXDirStyle,
	kMacDirStyle,
	kDOSDirStyle
	};

class NxsFilePath
	{
	public:
		NxsFilePath(const std::string &, bool isDir = false);
		
		bool				Exists() const;
		std::string 		GetFileName() const;
		std::string 		GetFullName() const;
		const std::string & GetNative() const;
		NxsFilePath			GetParent() const;
		bool 				IsDirectory() const;
		bool				IsFile() const;
		void				SetQueryUserOnError(bool v) {queryUserOnError = v;}
		
		NxsFilePath & operator += (const NxsFilePath &);
		NxsFilePath & operator += (const std::string &);
	
		enum SourceOfFileName
			{
			kProgrammer = 0x01,		/* file path is not being used as command option */
			kDefaultOpt  = 0x02,		/* file path is controlled by a command option and the user has never entered anything */
			kUserSupplied = 0x03		/* file path was set by the user */
			};
		SourceOfFileName sourceOfName;
		bool				WasSetByUser() const {return (sourceOfName == kUserSupplied);}
		
		enum QueryResult
			{
			kCouldNotQuery,
			kProblemWasFixed,
			kUserCancelled
			};
		
		STATIC_DATA DirDescriptionFormat dirStringInterpretStyle; /* specifies the format that will be used to interpret strings as path descriptions */
	protected:
		enum OpenOutputReturnCode
			{
			kBothReplaceAndAppend,
			kOpenAndReplacing,
			kOpenAndAppending,
			kNeitherReplaceOrAppend,
			kOpenedNewFile,
			kCouldNotOpen
			};
		void			 	CreateNative() const;
		bool			 	OpenInput(std::ifstream &outF) const;
		bool				ReadNewPathString(const std::string &);
		void 				TranslateDirStack(const std::list<std::string> &dirStack) const;
		
		STATELESS_FUNC void			ShortenDirStack(std::list<std::string> &dirStack);
		
		bool				isAbsolute;
		bool				isDirty;
		std::string			path;
		std::string			pathAsEntered;
		mutable std::string	nativePath;
		bool				queryUserOnError;
		
	};


class NxsInFilePath : public NxsFilePath
	{
	public :
		NxsInFilePath() : NxsFilePath(false) {}
		NxsInFilePath(const std::string &s) : NxsFilePath(s, false) {}
		
		bool		Open(std::ifstream &outF, std::string purposeOfFile = std::string());
		
	};
bool IsUpDirNotation(const std::string d);

void NxsFilePath::ShortenDirStack(std::list<std::string> &dirStack)
	{
	std::list<std::string>::iterator currIt = dirStack.begin();
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
	
void NxsFilePath::TranslateDirStack(const std::list<std::string> &dirStack) const
	{
	std::list<std::string>::const_iterator dIt = dirStack.begin();
	const std::list<std::string>::const_iterator endIt = dirStack.end();
	
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
                if (dIt != dirStack.begin())
                    {
                    // Only add slash if in the middle of the path
 				    nativePath.append("\\");
                    }
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
		std::list<std::string> dirStack;
		split<char, std::string>(nativePath, '/', &dirStack);
		ShortenDirStack(dirStack);
		TranslateDirStack(dirStack);
		}
	}

inline const std::string &NxsFilePath::GetNative() const
	{
	if (isDirty)
		CreateNative();
	return nativePath;
	}
inline std::string NxsFilePath::GetFullName() const
	{
	return GetNative();
	}

bool NxsFilePath::OpenInput(std::ifstream &inFStream) const
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


inline bool IsUpDirNotation(const std::string d)
	{
	return (d.length() == 2 && d[0] == '.' && d[1] == '.');
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
	
NxsFilePath::NxsFilePath(const std::string &s, bool shouldBeDir )
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
	std::ifstream testIn;
	bool retVal = OpenInput(testIn);
	testIn.close();
	return retVal;
	}
	
bool NxsFilePath::ReadNewPathString(const std::string &s)
	{
	pathAsEntered = path = s;
	nativePath.clear();
	isDirty = true;
	return true;
	}

bool NxsInFilePath::Open(std::ifstream &outF, std::string purposeOfFile)
	{
	return OpenInput(outF);
	}

class PhycasNexusReader: public PublicNexusReader
	{
	public: 
		PhycasNexusReader(const int blocksToRead, NxsReader::WarningHandlingMode mode=NxsReader::WARNINGS_TO_STDERR) 
			:PublicNexusReader(blocksToRead, mode),
			activeTaxaBlockIndex(UINT_MAX),
			activeCharactersBlockIndex(UINT_MAX)
			{}
			
		/* Unlike the PublicNexusReader, the Phycas Nexus Reader owns the blocks, 
			and thus deletes them on destruction
		*/
		~PhycasNexusReader()
			{
			Clear();
			}
		void	ReadFilePath(const std::string &filePath)
			{
			NxsInFilePath fp(filePath);
			ReadNxsFilePath(fp);
			}
		std::vector<std::string> GetTaxLabels() const;
		const std::vector<NxsFullTreeDescription> & GetTrees() const;
		unsigned GetNChar() const;
		std::string GetErrorMessage()
			{
			return errorMsg;
			}
		virtual void Clear()
			{
			PublicNexusReader::DeleteBlocksFromFactories();
			}
		virtual void ClearUsedBlockList()
			{
			PublicNexusReader::ClearUsedBlockList();
			activeTaxaBlockIndex = UINT_MAX;
			activeCharactersBlockIndex = UINT_MAX;
			}
		NxsCharactersBlock * getActiveCharactersBlock() const;
		virtual void NexusError(const std::string &msg, file_pos pos, unsigned line, unsigned col, CmdResult , NxsBlock* );
	private:
		std::string errorMsg;
		unsigned activeTaxaBlockIndex;
		std::vector<NxsFullTreeDescription> trees;
		unsigned activeCharactersBlockIndex;
		file_pos filePositionOfError;
		unsigned lineOfError ;
		unsigned columnOfError ;


		void	ReadNxsFilePath(NxsInFilePath & filePath);
		NxsTaxaBlock * getActiveTaxaBlock() const;
	};

NxsCXXDiscreteMatrix * GetLastDiscreteMatrix(PhycasNexusReader & nexusReader, bool convertGapsToMissing);
NxsCXXDiscreteMatrix * createNativeDiscreteMatrix(PhycasNexusReader & nexusReader, NxsTaxaBlock * taxaBlockPtr, unsigned int charBlockIndex, bool convertGapsToMissing);

NxsCXXDiscreteMatrix * GetLastDiscreteMatrix(PhycasNexusReader & nexusReader, bool convertGapsToMissing = false)
	{
	return createNativeDiscreteMatrix(nexusReader, 0L, UINT_MAX, convertGapsToMissing);
	}

NxsCXXDiscreteMatrix * createNativeDiscreteMatrix(PhycasNexusReader & nexusReader, NxsTaxaBlock * taxaBlockPtr, unsigned int charBlockIndex, bool convertGapsToMissing)
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
	return new NxsCXXDiscreteMatrix(*cb, convertGapsToMissing);
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
		errorMsg = std::string("The file \'");
		errorMsg.append(fn);
		errorMsg.append("\' does not exist.");
		throw NxsException(errorMsg.c_str());
		}
	if (filePath.IsDirectory())
		{
		errorMsg  = std::string ("\'");
		errorMsg.append(fn);
		errorMsg.append("\' is a directory, not a file.");
		throw NxsException(errorMsg.c_str());
		}
	std::ifstream inStream;
	if (!filePath.Open(inStream))
		{
		errorMsg = std::string("The file \'");
		errorMsg.append(fn);
		errorMsg.append("\' could not be opened.");
		throw NxsException(errorMsg.c_str());
		}
	std::string msg;
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


std::vector<std::string> PhycasNexusReader::GetTaxLabels() const
	{
	NxsTaxaBlock * activeTB = getActiveTaxaBlock();
	if (activeTB == 0L)
		return std::vector<std::string>();
	return activeTB->GetAllLabels();
	}
	
#endif

BOOST_PYTHON_MODULE(_ReadNexus)
{
	class_<NxsFullTreeDescription>("FullTreeDescription", init<NxsFullTreeDescription>())
		.def("getName", &NxsFullTreeDescription::GetName, return_value_policy<copy_const_reference>())
		.def("getNewick", &NxsFullTreeDescription::GetNewick, return_value_policy<copy_const_reference>())
		.def("isRooted", &NxsFullTreeDescription::IsRooted)
		;

	class_<std::vector<NxsFullTreeDescription> >("VecFullTreeDescription", no_init)
		.def("__iter__",  iterator<std::vector<NxsFullTreeDescription> >())
		.def("size", &std::vector<NxsFullTreeDescription>::size)
		;
   
	class_<PhycasNexusReader>("NexusReaderBase", init<int>())
		.def("readFile", &PhycasNexusReader::ReadFilePath)
		.def("getNChar", &PhycasNexusReader::GetNChar)
		.def("getErrorMessage", &PhycasNexusReader::GetErrorMessage)
		.def("getTrees", &PhycasNexusReader::GetTrees, return_value_policy<copy_const_reference>())
		.def("getTaxLabels", &PhycasNexusReader::GetTaxLabels)
		.def("clear", &PhycasNexusReader::Clear)
		;
	
	def("getLastNexusDiscreteMatrix", GetLastDiscreteMatrix, return_value_policy<manage_new_object>());

	register_exception_translator<NxsException>(&translate);
}
