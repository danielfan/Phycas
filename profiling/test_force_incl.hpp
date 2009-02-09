#define PHYCAS_ASSERT(x)
#define X_SPEC_THROW(x)
#define STATIC_DATA   static 
#define STATIC_CONST  static
#define STATELESS_FUNC static 
#define POLPY_NEWWAY  1
#define HAVE_STDDEF_H 1
#define FILE_AND_LINE
#define MONITORING_ALLOCATION 1
#define INSTANTIATE_MEMCHK 1

#if defined (__APPLE__)
#	define NXS_USE_MAC_DIR_STYLE
#	define HAVE_STDINT_H 1
#elif defined(_WIN32)
#	undef HAVE_STDINT_H
#	define NXS_USE_WINDOWS_DIR_STYLE
#	if defined(_MSC_VER)
#		if _MSC_VER >= 1500
#			include <cstdio>
#			if !defined(vsnprintf)
#				define vsnprintf _vsnprintf_s
#			endif
#			define sprintf sprintf_s
#		endif
#		if _MSC_VER >= 1200
#			include <basetsd.h>
			typedef   INT8 int8_t;
			typedef  UINT8 uint8_t;
			typedef  INT64 int64_t;
			typedef UINT64 uint64_t;
#		else
			typedef signed char int8_t;
			typedef unsigned char uint8_t;
			typedef long long int64_t;
			typedef unsigned long long uint64_t;
#		endif
#	endif
#else
#	define NXS_USE_UNIX_DIR_STYLE
#endif

#include <sys/stat.h>

#include <vector>
typedef std::vector<double> VecDbl;
