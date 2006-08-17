#ifndef NCL_NXSDEFS_H
#define NCL_NXSDEFS_H

#include <string>

#if defined (__GNUC__) && (__GNUC__ > 2)
#	define NCL_USING_SSTREAM
#else
#	undef NCL_USING_SSTREAM
#endif
#if defined(NCL_USING_SSTREAM)
#	include <sstream>
#else
#	include <strstream>
#endif
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <cfloat>
#ifdef NDEBUG
#	define TESTING(ignore) error still debugging ignore
#else
#	define TESTING(ignore) 1
#	if defined (STOP_DEBUGGER_ON_TRIPPED_ASSERTS)
#		define BOOST_ENABLE_ASSERT_HANDLER
#	endif
#endif /* def NDEBUG */

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

//Keep memchk first
//
#if defined(MONITORING_ALLOCATION)
#	include "memchk.h"
#endif


#define NCL_NAME_AND_VERSION  "NCL version 2.1"
#define NCL_COPYRIGHT         "Copyright (c) 1999-2003 by Paul O. Lewis"
#define NCL_HOMEPAGEURL       "http://lewis.eeb.uconn.edu/ncl/"

// Maximum number of states that can be stored; the only limitation is that this
// number be less than the maximum size of an int (not likely to be a problem).
// A good number for this is 76, which is 96 (the number of distinct symbols
// able to be input from a standard keyboard) less 20 (the number of symbols
// symbols disallowed by the NEXUS standard for use as state symbols)
//
#define NCL_MAX_STATES         76

#if defined(__MWERKS__) || defined(__DECCXX) || defined(_MSC_VER)
	typedef long		file_pos;
#else
	typedef std::streampos	file_pos;
#endif

#undef	SUPPORT_OLD_NCL_NAMES
typedef	unsigned int					UInt;
typedef std::vector<bool>				VecBool;
typedef std::vector<char>				VecChar;
typedef std::vector<UInt>				VecUInt;
typedef std::vector<int>				VecInt;
typedef std::vector<double>				VecDbl;
typedef std::pair<int, int>				PairInt;
typedef std::pair<UInt, UInt>			PairUInt;
typedef std::pair<double, double>		PairDbl;
typedef std::vector<std::string>		VecString;
typedef std::vector<VecString>			VecVecString;
typedef std::vector<VecUInt>			VecVecUInt;
typedef VecVecUInt::const_iterator		VecVecUInt_ConstIt;
typedef std::set<unsigned>				SetUInt;
typedef std::map<unsigned, std::string>		LabelMap;
typedef std::map<unsigned, VecString>		StateLabelMap;
		
// The following typedefs are simply for maintaining compatibility with existing code.
// The names on the right are deprecated and should not be used.
//
typedef	VecBool							BoolVect;
typedef VecUInt							UVec;
typedef PairDbl							DblPair;
typedef SetUInt							IntSet;
typedef std::map<std::string, IntSet>	IntSetMap;
typedef VecVecString					AllelesVect;
typedef VecString						LabelList;
typedef VecString						StrVec;
typedef VecString						vecStr;
typedef std::map<unsigned, StrVec>		LabelListBag;
typedef std::map<std::string, std::string>	AssocList;
typedef std::string						NxsString;

class NxsIndexSet;
typedef std::map<std::string, NxsIndexSet> 			PartitionDescription;
typedef std::map<std::string, PartitionDescription> 	NamedPartitions;
typedef unsigned char 							DataStorageType; /* binary encodings of character data are arrays of DataStorageType */

//	The standard is not clear on what to do if '' is encountered.  These 
//	are (almost certainly) user errors.  defining NXS_THROW_IF_EMPTY_TOKENS will
//	throw an exception if these are encountered.  This simplifies code by letting us assume that
// 	NxsToken::GetNextToken will always contain a token (except at the beginning/end of files)
//	 
//
#define NXS_THROW_IF_EMPTY_TOKENS

//	deriving from boost::noncopyable means that an object's copy constructor and equals operator cannot be used.
//	Despite the fact that boost::noncopyable is trivial and inlined, inheriting from it does affect the executable size
//	 (at least under MWerks and Apple's g++).  NON_COPYABLE and  AND_NON_COPYABLE allow one to derive from boost::noncopyable
//	in debug mode (this should be good enough, because it helps catch compile time errors).
//	use "class XXX NON_COPYABLE" if the class is not derived from any other base and
//	"class XXX : public YYY AND_NON_COPYABLE" if there are other base classes
//
#include <boost/noncopyable.hpp>
#if defined (NDEBUG)
#	define NON_COPYABLE 
#	define AND_NON_COPYABLE 
#	if defined(IGNORE_NXS_ASSERT)
#		define NXS_ASSERT(test)
#	else
#		define NXS_ASSERT(test) assert(test)
#	endif
#else
#	define NON_COPYABLE : public boost::noncopyable
#	define AND_NON_COPYABLE , public boost::noncopyable
#	if defined(IGNORE_NXS_ASSERT)
#		define NXS_ASSERT(test)
#	elif defined(STOP_DEBUGGER_ON_TRIPPED_ASSERTS) //@pol -> mth defined(STOP...) ok? Was just STOP... before
		bool	NxsAssertionFailed(const char *, unsigned); 
#		define NXS_ASSERT(test) ((test) ? ((void) 0) : (NxsAssertionFailed(__FILE__, __LINE__) ? assert(test) : assert(test)))
#	else
#		define NXS_ASSERT(test) (assert(test))
#	endif
#endif
#include "ncl/misc/compile_assert.hpp"

#define NCL_COMMAND_MAXLEN 512
#undef NCL_USING_BLOCK_SUPPLIERS

#define NEW_NXS_CMD_MGR 1

#if NEW_NXS_CMD_MGR
	//NEW vs OLD block and reader are inexorably linked because the old reader relies on the next block pointer that is only present in the old block
#	define NEW_NXS_BLOCK_AND_READER 1
#else
#	define NEW_NXS_BLOCK_AND_READER 0
#endif

#define NEW_NXS_TOKEN 1

#define OLD_NXS_CMD_MGR !NEW_NXS_CMD_MGR
#define OLD_NXS_TOKEN !NEW_NXS_TOKEN
#define OLD_NXS_BLOCK_AND_READER !NEW_NXS_BLOCK_AND_READER

//	if NCL_NXS_THROW_UNDEFINED NxsX_UndefinedException will be thrown as the result
//	of programmer mis-use or the library (in cases that only arise at run-time)
//	
#define NCL_NXS_THROW_UNDEFINED

#if NEW_NXS_TOKEN
#	if OLD_NXS_BLOCK_AND_READER
#		define SUPPORT_NXS_TOKEN_LABILE_FLAGS
#	else
#		undef SUPPORT_NXS_TOKEN_LABILE_FLAGS
#	endif
#	define GET_NEXT_TOKEN ReadToken
#else
#	define GET_NEXT_TOKEN GetNextToken
#endif
#undef DECLARING_ALL_EXCEPTIONS

#if defined DECLARING_ALL_EXCEPTIONS
#	define X_SPEC_THROW(x) throw(x)
#else
#	define X_SPEC_THROW(x) 
#endif

enum CmdResult
	{
	kCmdSucceeded,
	kCmdFailedSilent,
	kCmdFailedGenerateMessage
	};
	
enum CmdPermissionLevel 
	{
	kCmdPermBasicUser, 
	kCmdPermAdvancedUser, 
	kCmdPermRegisteredUsers, 
	kCmdPermBetaTesters, 
	kCmdPermDeveloper
	};
		

#define IMPLEMENTING_NEW_CHARS_BLOCK 1
#define IMPLEMENTING_OLD_CHARS_BLOCK !IMPLEMENTING_NEW_TREES_BLOCK

#define SUPPORT_TREE_OUTPUT

#if !defined(NXS_USE_UNIX_DIR_STYLE) && !defined(NXS_USE_MAC_DIR_STYLE) && !defined(NXS_USE_WINDOWS_DIR_STYLE)
#	define NXS_USE_UNIX_DIR_STYLE
#endif
#if !defined(NXS_IO_SUPPORT_INTERACTIVE_MODE) && ! defined(NXS_IO_NON_INTERACTIVE_MODE)
#	define NXS_IO_SUPPORT_INTERACTIVE_MODE
#endif

#include "ncl/output/output_fwd_decl.hpp"

#if 0
//@POL 26-Oct-2005 Mark, I had to disable this when compiling pyphy modules because it trips on
// some Microsoft internal headers (e.g. xlocmon) that are invoked when <boost/format>
// is #included
template<typename T>
std::string & operator+=(std::string &, const T & d); //temp
template<typename T>
std::string & operator+=(std::string & s, const T & d) //temp
	{
	// use of a potentially dangerous string += operator 
	// will match this function and trip the following compile time assert
	//  If this happens you should use the << operator instead of +=
	// previously we supported string += int in ncl/misc/string_extensions.hpp
	// this could cause problems if that header was not included early enough 
	// because compilers would match to the std namespace's string += char operator
	// (without warning in the case of g++) resulting in appending the ascii interpretation 
	// of the integer instead of the number in string form.
	//	To avoid such problems we can append with the << operator (in ncl/misc/string_extensions.hpp)
	// this templated function detects any remnants of the old += style (or catches us if we relapse).
	COMPILE_TIME_ASSERT(0, BadStringOperator);
	return s;
	}
#endif

#endif
