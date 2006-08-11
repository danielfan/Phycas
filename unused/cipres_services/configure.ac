# This file is processed by autoconf to create a Makefile for the C++ portion of the CIPRES comm lib.
# This file was written by hand (by mth using PAUP's configure.ac and http://autotoolset.sourceforge.net/tutorial.html as guides).

AC_INIT(phycas_cipres_services, 0.1)
# Safety-check for srcdir
AC_CONFIG_SRCDIR(../ncl/nxs_exception.cpp)

if  test -n "$PHYCAS_ROOT" ; then
	:
else
	AC_MSG_ERROR([PHYCAS_ROOT must be set to the parent of the phycas and cipres_services directory])
fi
if  test -n "$NCL_ROOT" ; then
	:
else
	AC_MSG_WARN([NCL_ROOT is not set; it should be the parent of the ncl directory.  Assuming that it is the same as PHYCAS_ROOT...])
	NCL_ROOT="$PHYCAS_ROOT"
fi

if  test -n "$BOOST_ROOT" ; then
	:
else
	AC_MSG_ERROR([BOOST_ROOT must be set to the parent of the IDL directory])
fi

# the Makefile that we use imports the  from CIPRES/framework/C++/config/Makefile that is produced by 
# running configure for the CIPRES library
# So we need to know the path to the build dir (not necessarily the source or install dir) for CIPRES
# The precedence rules for the possible values for this setting are:
#	1	use the CIPRES_LIB_BUILD_DIR= arg to configure
#	2	use the CIPRES_LIB_BUILD_DIR environmental var
#	3	add /C++ to the CIPRES_LIB_ROOT  environmental var
#	4	quit with an error
#
AC_ARG_VAR(CIPRES_LIB_BUILD_DIR, [path to parent of CipresIDL, CipresCommlib, and config directories of CIPRES.  This directories should be configured and built for this platform.])

if test -n "$CIPRES_LIB_BUILD_DIR" ; then
	:
else
	if test -n "$CIPRES_LIB_ROOT" ; then
		AC_MSG_WARN([CIPRES_LIB_BUILD_DIR is not set and was not specified using an arg to configure; Assuming that it is the same as CIPRES_LIB_ROOT/C++...])
		CIPRES_LIB_BUILD_DIR="$CIPRES_LIB_ROOT/C++"
	else
		AC_MSG_ERROR([CIPRES_LIB_BUILD_DIR must be used to specify the path to the built version of C++ directory that creates the cipres lib])
	fi
fi

# Get build-platform and target info
AC_CANONICAL_BUILD
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE

#--------------------------------------------------------------------------------------------------#
#  Preliminaries                                                                                   #
#--------------------------------------------------------------------------------------------------#

# Output configuration macros to config.h using template from config.h.in
#	bogus destination for config.h, because we are currently just using the 
#	autoconf process for Makefile creation (we use CIPRES_LIB_ROOT/C++/config.h)
AM_CONFIG_HEADER(bogus_config.h)


# compilers macros used by autoheader
AC_PROG_CC
AC_PROG_CXX

# required because we are building a library 
AC_PROG_RANLIB


AC_PROG_INSTALL

#--------------------------------------------------------------------------------------------------#
#  Headers                                                                                         #
#--------------------------------------------------------------------------------------------------#

AC_DEFINE(WRONG_CONFIG_H_INCLUDED, 1, make sure we are including the config.h in CIPRES_LIB_ROOT/C++)


#--------------------------------------------------------------------------------------------------#
#  Set CXXFLAGS                                                                                    #
#--------------------------------------------------------------------------------------------------#

# Initialize CXXFLAGS to an empty string to prevent it from defaulting to "-g -O2"
CXXFLAGS=""

# Check for debugging mode.
AC_ARG_ENABLE(debugging,
	AC_HELP_STRING([--enable-debugging],[build for debugging]), , [enable_debugging=no])
if test "$enable_debugging" = yes; then
	AC_MSG_NOTICE([
	   *** NOTE: debugging is enabled; optimization is suppressed!
	   ])
fi

# Initialize optimization flag in case it doesn't get set below.
CXXFLAGS_OPTIM_SPEED="-O"

# "-g" may not work with some compilers, but end users shouldn't be
if test "$enable_debugging" = yes; then
	CXXFLAGS_OPTIM_SPEED="-O0"
	CXXFLAGS="-g"
fi

if test "$CC" = "icc"; then
#	Intel C compiler for Linux
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-O3"
		CXXFLAGS_OPTIM_SIZE="-O2"
	fi
elif test "$CC" = "ccc"; then
#	Compaq C compiler for Linux
	if test "x$arch" = "x"; then
		arch="host"
	fi
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-fast -inline speed -arch $arch"
		CXXFLAGS_OPTIM_SIZE="-fast -inline size -unroll 1 -arch $arch"
	fi
elif test "$CC" = "xlc"; then
#	IBM XL C compiler
	if test "x$arch" = "x"; then
		arch="auto"
	fi
	if test "$enable_debugging" = no; then
		CXXFLAGS_OPTIM_SPEED="-O3 -qarch=$arch"
		CXXFLAGS_OPTIM_SIZE="-O3 -qarch=$arch"
	fi
elif test "x$GCC" = "xyes" ; then
	if test "$enable_debugging" = yes; then
#		Suppress warnings about possibly uninitialized variables but show everything else (used for
#   	development, but these warnings should also not trip for release builds)
		CXXFLAGS_WARNINGS="-Wall -Wno-uninitialized -Wimplicit -Wmissing-prototypes"
	else
#		Just suppress warnings about possibly uninitialized variables	
		CXXFLAGS_WARNINGS="-Wno-uninitialized"
		CXXFLAGS_OPTIM_SPEED="-O3 -ffast-math"
		CXXFLAGS_OPTIM_SIZE="-O2 -ffast-math"
	fi
	case "$build_os" in
	 	darwin*) CXXFLAGS_WARNINGS="$CXXFLAGS_WARNINGS -Wno-long-double";;
	 	*);;
	esac	
fi
if test "x$CXXFLAGS_OPTIM_SIZE" = "x"; then		#FIXME: just test for empty? here and elsewhere
	CXXFLAGS_OPTIM_SIZE=$CXXFLAGS_OPTIM_SPEED
fi


# libraries and include directories  needed to compile services that link to libCipres and libCiphycas
CIPHYCAS_DIR="\$(PHYCAS_ROOT)/cipres_services/ciphycas"
CIPHYCAS_LIBFLAG="Ciphycas"
CIPHYCAS_PREFIX_FILE="\$(PHYCAS_ROOT)/phycas/cipres_module_config.h"


CPPFLAGS="-I$PHYCAS_ROOT -I$BOOST_ROOT -I$NCL_ROOT"
# Variables to substitute in Makefile:
AC_SUBST(CXXFLAGS_OPTIM_SPEED)
AC_SUBST(CXXFLAGS_OPTIM_SIZE)
AC_SUBST(CXXFLAGS_WARNINGS)

AC_SUBST(CIPHYCAS_DIR)
AC_SUBST(CIPHYCAS_LIBFLAG)
AC_SUBST(CIPHYCAS_PREFIX_FILE)
AC_SUBST(PHYCAS_ROOT)
AC_SUBST(NCL_ROOT)
AC_SUBST(CIPRES_LIB_BUILD_DIR)

AC_CONFIG_FILES([Makefile
		ciphycas/Makefile
		ciphycas/RealMakefile
		read_nexus/Makefile
		read_nexus/RealMakefile
		tree_improve_client/Makefile
		tree_improve_client/RealMakefile
		tree_merge/Makefile
		tree_merge/RealMakefile ])
AC_OUTPUT
