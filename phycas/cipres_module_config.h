#if !defined CIPRES_MODULE_CONFIG_H
#define CIPRES_MODULE_CONFIG_H

#if defined(_MSC_VER)
//tucson2005: this needs to come before #include "CipresCommlib/ConfigDependentHeaders.h"
// because of need to define CIPRES_USING_OMNI_ORB
// because of need to define ACE_THROW_SPEC
#	include "./windows_cipres_config.h"
#endif

#include "CipresCommlib/ConfigDependentHeaders.h"

#if defined (HAVE_CONFIG_H)
#	undef PACKAGE
#	undef PACKAGE_BUGREPORT
#	undef PACKAGE_NAME 
#	undef PACKAGE_STRING
#	undef PACKAGE_TARNAME
#	undef PACKAGE_VERSION
#elif !defined(_MSC_VER) //tucson2005: windows config file included above
#	error Need to run configure
#endif

#if defined(WRONG_CONFIG_H_INCLUDED)
#	error Include paths should be set up to include the config.h created when configuring CipresCommLib, not the file created from PHYCAS_ROOT/cipres_services !
#endif

#define CIPRES_USING_PHYCAS_CODE_BASE
#define NXS_IO_NON_INTERACTIVE_MODE
#define NXS_SUPPRESS_OUTPUT

#define SCM_MODULE
#define POL_TEMP
#define CORBA_SERVER_PHYCAS
#include "./phycas_config.h"

#endif
