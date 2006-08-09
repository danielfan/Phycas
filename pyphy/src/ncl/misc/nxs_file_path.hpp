#if !defined (NCL_FILE_PATH_H)
#define NCL_FILE_PATH_H

//POL, I (MTH) tried to wrap unix directory listing code (in this file and the cpp file) that might not work on windows with
//  NCL_SUPPORT_DIR_LIST
//#define NCL_SUPPORT_DIR_LIST

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

// Many of the idea's for this class come from Boost's filesystem path:

//  boost/filesystem/path.hpp  -----------------------------------------------//

//  (C) Copyright Beman Dawes 2002. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.


//  See http://www.boost.org/libs/filesystem for documentation.

//----------------------------------------------------------------------------// 
#include "ncl/nxs_defs.hpp"
#include <fstream>
#include <list>
#include <sys/stat.h>

#define NEW_PATH_STRING_IMPLEMENTATION 1

//	Need to define a filesystem so we know whether to use / : or \ as directory separators
//		
#if !defined(NXS_USE_UNIX_DIR_STYLE) && !defined(NXS_USE_MAC_DIR_STYLE) && !defined(NXS_USE_WINDOWS_DIR_STYLE)
#	error must define a filesystem
#endif

// Note: const refers to the  path itself other fields are mutable.
enum DirDescriptionFormat  // names for the 3 systems of interpretting strings as paths
	{
	kUNIXDirStyle,
	kMacDirStyle,
	kDOSDirStyle
	};
enum OutModes
	{
	kUnspecifiedOutMode,
	kReplaceOutMode,
	kAppendOutMode
	};
	
class NxsFilePath
	{
	public:
		STATELESS_FUNC NCL_TIME_T NxsGetModTime(const std::string &dirName); //POL VC says time_t not in std
		
		NxsFilePath(bool isDir = false);
		NxsFilePath(const std::string &, bool isDir = false);
		
		bool				Empty() const;
		bool				Exists() const;
		std::string 		GetDirectory() const;
		std::string 		GetFileName() const;
		std::string 		GetFullName() const;
		const std::string & GetNative() const;
		NxsFilePath			GetParent() const;
		bool 				IsDirectory() const;
		bool				IsFile() const;
		bool				Remove() const;
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
		bool 			 	IsLegalFileName() const;
		OpenOutputReturnCode OpenOutput(std::ofstream &outF, bool replaceFile = false, bool appendToFile = false) const;
		bool			 	OpenInput(std::ifstream &outF) const;
		bool				ReadNewPathString(const std::string &);
		void 				TranslateDirStack(const std::list<std::string> &dirStack) const;
		
#		if (NEW_PATH_STRING_IMPLEMENTATION)
			STATIC_DATA_FUNC bool 	TranslateToInternalRep(std::string &); // assumes IsLegalFileName()
#		endif
		STATELESS_FUNC std::string 	RemoveLastDir(std::string *);
		STATELESS_FUNC void			ShortenDirStack(std::list<std::string> &dirStack);
		
		bool				isAbsolute;
		bool				isDirty;
		mutable bool		isLegal;
		std::string			path;
		std::string			pathAsEntered;
		mutable std::string	nativePath;
		bool				queryUserOnError;
		
		
		friend class NxsOutFileCmdOption;
		friend class NxsInFileCmdOption;
	};

class NxsOutFilePath : public NxsFilePath
	{
	public :
		NxsOutFilePath();
		NxsOutFilePath(const std::string &s);
		
		bool Open(std::ofstream &outF, std::string purposeOfFile = std::string());
			
	
		bool GetAppend() const;
		bool GetReplace() const;
		bool IsAppending() const;
		bool IsReplacing() const;
		bool WasOpened() const
			{
			return wasOpened;
			}
		
		void SetAppend(bool a = true) ;
		void SetReplace(bool r = true) ;
		
		bool	replaceFile, appendToFile;
		QueryResult QueryUserForNewFile(const std::string &);
		
	private:
		bool	isAppending, isReplacing, wasOpened;
		
		std::string ComposeReplaceAppendErrorMsg() const;
		
		NxsFilePath::QueryResult	QueryUserForReplaceAppend();
		
		enum 	ReplaceAppendUserCorrection
			{
			kNoCorrection,
			kCorrectedToReplace,
			kCorrectedToAppend
			};
		ReplaceAppendUserCorrection userCorrection;
		
		friend class NxsOutFileCmdOption;
		friend class NxsOutputCmdOption;
		
	};

bool NxsOpenOutFile(const std::string &fn, std::ofstream &, OutModes o, bool canQuery = true);
bool NxsOpenInFile(const std::string &fn, std::ifstream &, bool canQuery = true);
bool NxsMkDir(const std::string &DirName, NCL_MODE_T dirMode = NCL_DIRMODE);
void NxsMkDirOrThrow(const std::string &dirName, NCL_MODE_T dirMode = NCL_DIRMODE);
bool NxsCheckForDir(const std::string &dirName);
bool NxsCheckForOrMkDir(const std::string &dirName, NCL_MODE_T dirMode = NCL_DIRMODE);
void NxsCheckForOrMkDirOrThrow(const std::string &dirName, NCL_MODE_T dirMode = NCL_DIRMODE);
unsigned NxsFindIndexOfNextOrderedFile(const std::string &substitutionString, unsigned firstIndexToCheck = 0);
unsigned NxsFindIndexOfNextOrderedDir(const std::string &substitutionString, unsigned firstIndexToCheck = 0);
#if defined(NCL_SUPPORT_DIR_LIST)
	VecString NxsGetDirList(const std::string &dirName);
#endif
class NxsInFilePath : public NxsFilePath
	{
	public :
		NxsInFilePath() : NxsFilePath(false) {}
		NxsInFilePath(const std::string &s) : NxsFilePath(s, false) {}
		
		bool		Open(std::ifstream &outF, std::string purposeOfFile = std::string());
		QueryResult QueryUserForNewFile(const std::string &);
		
	};

inline NxsOutFilePath::NxsOutFilePath() 
  : NxsFilePath(false), 
  replaceFile(false), 
  appendToFile(false),
  isAppending(false), 
  isReplacing(false),
  wasOpened(false), 
  userCorrection(kNoCorrection) 
  	{}

inline NxsOutFilePath::NxsOutFilePath(const std::string &s) 
	:NxsFilePath(s, false), 
	replaceFile(false), 
	appendToFile(false), 
	isAppending(false), 
	isReplacing(false), 
  	wasOpened(false),
	userCorrection(kNoCorrection)
	{}


inline const std::string &NxsFilePath::GetNative() const
	{
	if (isDirty)
		CreateNative();
	return nativePath;
	}

inline NxsFilePath & NxsFilePath::operator+=(const std::string &s)
	{
	NxsFilePath p(s);
	return *this += p;
	}
	
inline std::string NxsFilePath::GetFullName() const
	{
	return GetNative();
	}
		
	
inline bool NxsFilePath::Empty() const
	{
	return path.empty();
	}

	
inline bool NxsFilePath::IsFile() const 
	{
	return Exists() && !IsDirectory();
	}

inline bool NxsOutFilePath::IsAppending() const
	{
	return isAppending;
	}

inline bool NxsOutFilePath::IsReplacing() const
	{
	return isReplacing;
	}

inline bool NxsOutFilePath::GetAppend() const
	{
	return (userCorrection == kNoCorrection ? appendToFile : (userCorrection == kCorrectedToAppend));
	}

inline bool NxsOutFilePath::GetReplace() const
	{
	return (userCorrection == kNoCorrection ? replaceFile : (userCorrection == kCorrectedToReplace));
	}

inline void NxsOutFilePath::SetAppend(bool a)
	{
	appendToFile = a;
	}

inline void NxsOutFilePath::SetReplace(bool r)
	{
	replaceFile = r;
	}

	
#endif

