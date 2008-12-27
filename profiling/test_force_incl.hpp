#define PHYCAS_ASSERT(x)
#define X_SPEC_THROW(x)
#define STATIC_DATA   static 
#define STATIC_CONST  static
#define STATELESS_FUNC static 
#define POLPY_NEWWAY  1
#define HAVE_STDINT_H 1
#define HAVE_STDDEF_H 1
#define FILE_AND_LINE
#define MONITORING_ALLOCATION 1
#define INSTANTIATE_MEMCHK 1

#if defined (__APPLE__)
#	define NXS_USE_MAC_DIR_STYLE
#elif defined(_WIN32)
#	define NXS_USE_WINDOWS_DIR_STYLE
#else
#	define NXS_USE_UNIX_DIR_STYLE
#endif

#include <sys/stat.h>

#include <vector>
typedef std::vector<double> VecDbl;
