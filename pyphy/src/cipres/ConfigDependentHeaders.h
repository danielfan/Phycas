#if ! defined(CONFIG_DEPENDENT_HEADERS_H)
#define CONFIG_DEPENDENT_HEADERS_H

/**
 * This file uses the macros defined in config.h (presumably created 
 * by running configure) to include some fundamental headers.
 * (TL - is this cipres_config.h, not config.h ?)
 * 
 * This has been tested on Mac and Windows (Visual Studio .NET 2003).
 * Additional accommodations may be needed for some flavors of *nix.
 */

#if defined (HAVE_CONFIG_H)
#	if !defined (CIPRES_CONFIG_H)
#		include "CipresCommlib/cipres_config.h"
#		undef PACKAGE
#		undef PACKAGE_BUGREPORT
#		undef PACKAGE_NAME 
#		undef PACKAGE_STRING
#		undef PACKAGE_TARNAME
#		undef PACKAGE_VERSION
#	endif
#else
	// omniORB is our default ORB, if CIPRES_USING_OMNI_ORB is not 0
#	if !defined(CIPRES_USING_OMNI_ORB)
#       define CIPRES_USING_OMNI_ORB 1
#	endif
#endif //defined (HAVE_CONFIG_H)

#if !defined(CIPRES_USING_OMNI_ORB)
#	define CIPRES_USING_OMNI_ORB 0
#endif
#if ! CIPRES_USING_OMNI_ORB
#	error "only omniORB is currently supported"
#endif

#if defined(_MSC_VER)
#	undef	HAVE_COMPILE_TIME_DISPATCH
#	if CIPRES_USING_OMNI_ORB
		// These VC7 warnings were being generated by OmniORB
		//  warning C4311: 'type cast' : pointer truncation from 'char *' to '_CORBA_Sequence_String::ptr_arith_t'
		//  warning C4312: 'type cast' : conversion from '_CORBA_Sequence_WString::ptr_arith_t' to '_CORBA_WChar *' of greater size
		//  warning C4267: 'initializing' : conversion from 'size_t' to '_CORBA_ULong', possible loss of data
#		pragma warning( disable : 4311 4312 4267)
#	endif
#else
#	define HAVE_COMPILE_TIME_DISPATCH
#endif

	/* For typedefs like uint8_t */
#if HAVE_INTTYPES_H
#	include <inttypes.h>
#elif HAVE_STDINT_H
#	include <stdint.h>
#elif defined(_MSC_VER) && _MSC_VER >= 1200
#	include <basetsd.h>
	typedef   INT8 int8_t;
	typedef  UINT8 uint8_t;
	typedef  INT64 int64_t;
	typedef UINT64 uint64_t;
#elif defined(_MSC_VER)
	typedef signed char int8_t;
	typedef unsigned char uint8_t;
	typedef long long int64_t;
	typedef unsigned long long uint64_t;
	//typedef unsigned int uint32_t;
	//typedef unsigned short uint16_t;
#elif defined(_WIN32)
#	include <stdint.h>
#endif

	/* For size_t */
#if defined(HAVE_STDDEF_H)
#	include <stddef.h>
#elif defined(_MSC_VER)
//#	include <cstddef>
//#include <cstdlib>
#endif 


/* what follows is some regrettable macro hackery for including header files while
	created by different IDL compilers (ACE and OmniORB).  
	ACE wants to produce 2 different headers (S and C), while omniORB just produces one,
	so we have to use different include statements.
*/
#if CIPRES_USING_OMNI_ORB
#	define CIPRES_USING_IORTABLE 0
#	if defined(_MSC_VER)
		// avoiding Microsoft VC7 warning C4290: C++ exception specification ignored except to indicate a function is not __declspec(nothrow)
#		define ACE_THROW_SPEC(X)
#	else
#		define ACE_THROW_SPEC(X) throw X
#	endif
#else 
#	error "only CIPRES_USING_OMNI_ORB is supported"
#endif 


#endif
