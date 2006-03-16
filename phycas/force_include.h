#if !defined (FORCE_INCLUDE_HEADERS_H)
#define FORCE_INCLUDE_HEADERS_H
// mth, added because I can't find the force include header option for bjam 
#if defined(FORCE_INCLUDE_CONSOLE_CONFIG)
#	include "phycas/phycas_console_config.h"
#elif defined(FORCE_INCLUDE_SOCKET_CONFIG)
#	include "phycas/phycas_config.h"
#elif defined(FORCE_INCLUDE_MTH_SOCKET_CONFIG)
#	include "phycas/mth_dev_phycas_config.h"
#elif defined(FORCE_INCLUDE_MTH_CONSOLE_CONFIG)
#	include "phycas/mth_dev_phycas_config.h"
#elif defined(FORCE_INCLUDE_POL_SOCKET_CONFIG)
#	include "phycas/pol_dev_phycas_socket_config.h"
#elif defined(FORCE_INCLUDE_POL_CONSOLE_CONFIG)
#	include "phycas/pol_dev_phycas_console_config.h"
#elif defined (FORCE_INCLUDE_EXTENDING_PYTHON_CONFIG)
#	include "phycas/pol_pyphy_config.h"
#elif ! defined(PHYCAS_CONFIG_H)
#	error "Must define FORCE_INCLUDE_(MTH_|POL_)?(CONSOLE|SOCKET)_CONFIG or FORCE_INCLUDE_EXTENDING_PYTHON_CONFIG use your compiler's force include option to include the correct header"
#endif

#endif

