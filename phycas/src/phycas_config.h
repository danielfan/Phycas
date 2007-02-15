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

#ifndef PHYCAS_CONFIG_H
#define PHYCAS_CONFIG_H
  
//
// below here from old phypy_config.h
//
#define PYTHON_ONLY	// for code that only makes sense when exported to Python (most code can be used to construct a pure C++ program)
#define POL_PHYCAS	// for POL-specific code related to PHYCAS (this one needs to go)
//#define POL_PYPHY	// for POL-specific code related to boost python (this one needs to go too)
#define NO_IDL_TYPES
#define POLPY_NEWWAY  1
#define POLPY_OLDWAY  !(POLPY_NEWWAY)
//#define USING_NUMARRAY
//#define INTERFACE_WITH_CIPRES
//
// above here from old phypy_config.h
//

#if defined(PYTHON_ONLY)
#	include <Python.h>
#endif

#if defined(NDEBUG)
#	define FILE_AND_LINE
#	define PHYCAS_ASSERT(expr)
#elif defined(_MSC_VER) // Microsoft Visual Studio C++ compiler issues error C1055 "out of keys" unless __LINE__ is not used in asserts or phycas::Uniform constructor
#	define FILE_AND_LINE
#	include <boost/assert.hpp>
#	define PHYCAS_ASSERT(expr)  if (!(expr)) boost::assertion_failed((const char *)#expr, (const char *)__FUNCTION__, (const char *)"" /* __FILE__ */, 0L /* __LINE__*/)
#	define IGNORE_NXS_ASSERT	// used in nxs_defs.hpp to cause NXS_ASSERT to expand to nothing
#else
#	define FILE_AND_LINE  __FILE__,__LINE__
#	include <cassert>
#	define PHYCAS_ASSERT(expr)  assert(expr)
#endif

#define NXS_SUPPRESS_OUTPUT

#if defined(__MWERKS__) && defined(__APPLE__)
#	include <MSLCarbonPrefix.h>
#endif

/* to make cipres services supplied by phycas reentrant, MTH is on a crusade against statics 
The following defines will help me wade through harmless, class-level functions and 
  and static data and functions that hinder thread-safety.
*/
 /* STATIC_DATA - for use with static data members and local variables that are declared static */
#define STATIC_DATA static 
 /* STATIC_DATA_FUNC -  for use to declare static functions that rely on static data */
#define STATIC_DATA_FUNC static 
 /* STATELESS_FUNC - for functions that are part of a class's namespace, but are static AND do not 
   access in static data */

#define STATELESS_FUNC static 
  /* STATIC_SINGLETON -  used for singletons that are declared as statics: */
#define STATIC_SINGLETON static
  /* STATIC_CONST_ACCESSOR - used to declare functions that refer to static data that is constant
  these functions will be thread-safe despite using static data */
#define STATIC_CONST_ACCESSOR static
  /* STATIC_CONST declares data as const and static */
#define STATIC_CONST static

/*currently we only support a console or suppressing output */
#define CONSOLE_PHOREST 
#if defined (CONSOLE_PHYCAS_NO_SOCKET)
#   if !defined(NCL_USE_STD_OUTPUT)
#		define NCL_USE_NXS_CONSOLE_OUTPUT
#   endif
#   undef NCL_SOCKET_IO
#   define NCL_PRINT_COMMAND_STATE_TO_HIDDEN 0
#elif defined(NXS_SUPPRESS_OUTPUT)
//@POL 27-Oct-2005 Mark, I get lots of errors when I try to use this code
#   undef NCL_SOCKET_IO
#	define NCL_USER_SUPPLIED_OUTPUT
#	define NCL_SUPPRESS_OUTPUT
#   define NCL_PRINT_COMMAND_STATE_TO_HIDDEN 0
#	define NCL_USER_OUTPUT_FWD_DECLARATIONS "phycas/src/ncl/output/nxs_suppress_fwd_decl.hpp"
#	define NCL_USER_OUTPUT_HEADER "phycas/src/ncl/output/nxs_suppress_output_stream.hpp"
#	define NCL_USER_OUTPUT_MGR_HEADER "phycas/src/ncl/output/nxs_suppress_output_mgr.hpp"
#else
#	define NCL_SOCKET_IO
#	define NCL_USER_SUPPLIED_OUTPUT
#   define NCL_PRINT_COMMAND_STATE_TO_HIDDEN 1
#	define NCL_USER_OUTPUT_FWD_DECLARATIONS "phycas/src/ncl/output/nxs_xml_socket_fwd_decl.hpp"
#	define NCL_USER_OUTPUT_HEADER "phycas/src/ncl/output/nxs_xml_socket_output_stream.hpp"
#	define NCL_USER_OUTPUT_MGR_HEADER "phycas/src/ncl/output/nxs_xml_socket_output_mgr.hpp"
#endif //defined (CONSOLE_PHYCAS_NO_SOCKET)

#if defined (CORBA_CLIENT_PHYCAS) || defined(CORBA_SERVER_PHYCAS)
#	define CORBA_PHYCAS
#endif

	// if this is just a CIPRES module, then we don't read PHYCAS blocks
#if !defined (CIPRES_USING_PHYCAS_CODE_BASE)
#	define READ_PHYCAS_BLOCK
#endif

#if defined (__APPLE__)
#	define MAC_PHOREST
#	define MACH_O_PHOREST
#elif defined(_WIN32)
#	define WIN_PHOREST
#	define _WIN32_WINNT   0x0501	// assuming target OS is Windows XP or later, see phycas/misc/stopwatch.hpp
#else
#	define LINUX_PHOREST	//POL added
#endif

	// try to read _run.nex
#define DEBUGGING_FROM_RUN_NEX
#define STOP_DEBUGGER_ON_TRIPPED_ASSERTS 1

#if defined(__MWERKS__)
#	define C_FUNCS_IN_STD_NAMESPACE
#	define HAVE_PRAGMA_UNUSED
#	if defined(MAC_PHOREST)
#		define BOOST_DISABLE_THREADS
#	endif
#	if __ide_target("d_phycas_mac_console") || __ide_target("d_PhycasWinConsole") || __ide_target("d_phycas_mac_socket")|| __ide_target("mth_dev_mac_console")|| __ide_target("dummy")
#		undef NDEBUG 
#	else
#		define NDEBUG
#	endif
#	if __ide_target("d_phycas_mac_memchk")
#		define MONITORING_ALLOCATION	// Special target for instantiating memory checker (it really slows things down, and so it shouldn't be on for all debugging
#	endif
#endif	/* defined(__MWERKS__) */

#if defined (WIN_PHOREST)
	//@@ this misspecification of the directory style is a nxsfilepath bug that needs to be fixed
	//#	define NXS_USE_WINDOWS_DIR_STYLE
#	undef NXS_USE_WINDOWS_DIR_STYLE
#elif defined (LINUX_PHOREST) 
#	define NXS_USE_UNIX_DIR_STYLE
#elif defined (MAC_PHOREST)
#	if defined (MACH_O_PHOREST)
#		define NXS_USE_UNIX_DIR_STYLE
#	else
#		define NXS_USE_MAC_DIR_STYLE
#	endif
#else
#	error need to define a platform
#endif

#include "phycas/src/ncl/singleton.hpp"

class Phylodiversity;
class DistributionManager;
class PhoSampler;
class PhoMCMC;
class PWD;
class LongOperationManager;

#if defined(USE_LOKI_SINGLETON)
	typedef Loki::SingletonHolder<DistributionManager, Loki::CreateUsingNew, Loki::SingletonWithLongevity> DistributionManagerSingletonHolder;
	typedef Loki::CreateUsingNew<DistributionManager> DistributionManagerCreator;
	typedef Loki::SingletonHolder<PhoSampler, Loki::CreateUsingNew, Loki::SingletonWithLongevity> SamplerSingletonHolder;
	typedef Loki::CreateUsingNew<PhoSampler> PhoSamplerCreator;
	typedef Loki::SingletonHolder<PhoMCMC, Loki::CreateUsingNew, Loki::SingletonWithLongevity> MCMCSingletonHolder;
	typedef Loki::CreateUsingNew<PhoMCMC> PhoMCMCCreator;
	typedef Loki::SingletonHolder<PWD, Loki::CreateUsingNew, Loki::SingletonWithLongevity> PWDSingletonHolder;
	typedef Loki::CreateUsingNew<PWD> PWDCreator;
	typedef Loki::SingletonHolder<LongOperationManager, Loki::CreateUsingNew, Loki::SingletonWithLongevity> LongOperationManagerSingletonHolder;
	typedef Loki::CreateUsingNew<LongOperationManager> LongOperationManagerCreator;
	typedef Loki::SingletonHolder<Phylodiversity, Loki::CreateUsingNew, Loki::SingletonWithLongevity> PhylodiversitySingletonHolder;
	typedef Loki::CreateUsingNew<Phylodiversity> LongOperationManagerCreator;
#else
	typedef ncl::SingletonHolder<DistributionManager> DistributionManagerSingletonHolder;
	typedef ncl::SingletonHolder<DistributionManager> DistributionManagerCreator;
	typedef ncl::SingletonHolder<PhoSampler> SamplerSingletonHolder;
	typedef ncl::SingletonHolder<PhoSampler> PhoSamplerCreator;
	typedef ncl::SingletonHolder<PhoMCMC> MCMCSingletonHolder;
	typedef ncl::SingletonHolder<PhoMCMC> PhoMCMCCreator;
	typedef ncl::SingletonHolder<PWD> PWDSingletonHolder;
	typedef ncl::SingletonHolder<PWD> PWDCreator;
	typedef ncl::SingletonHolder<LongOperationManager> LongOperationManagerSingletonHolder;
	typedef ncl::SingletonHolder<LongOperationManager> LongOperationManagerCreator;
	typedef ncl::SingletonHolder<Phylodiversity> PhylodiversitySingletonHolder;
	typedef ncl::SingletonHolder<Phylodiversity> PhylodiversityCreator;

#endif

const unsigned int kNxsDistributionManagerLongevity = 40;
const unsigned int kNxsTreesManagerLongevity = 40;
const unsigned int kNxsCharactersManagerLongevity = 40;
const unsigned int kNxsTaxaManagerLongevity = 50;
const unsigned int kNxsOutputManagerLongevity = 100;
#include "phycas/src/phycas.h"

#endif	/* __PHOREST_CONFIG_H */
