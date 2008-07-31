#!/bin/sh
################################################################################################
# This file sourced whenever phycas is run through the Phycas.app 
# The Phycas.app executes ~/.phycas/run_phycas.sh which sources this file before invoking python.
# 
# The default behavior is for this file to source the file "$HOME/.phycas/phycas_gui_env.sh" 
#  which is overwritten every time a Phycas.app window is opened so that it always reflects 
#  the current location of Phycas.app 
# 
# If you would like to augment the default behavior of Phycas you can:
#     - define (and export) the environmental variable PYTHONINTERPRETER 
#       to point to the full path to the python executable that you would like to use. 
#     - add the commnd:
# export USES_I_PYTHON=1
#       to make the Phycas.app bundle take you into the IPython environment.
# 
# To use  different version of the phycas python package with the Phycas.app bundle
# you can comment out the line
#    source "$env_settings_path"
# below, and make sure that  PHYCAS_PYTHONPATH and PHYCAS_DYLD_LIBRARY_PATH are
#   set correctly for the version of phycas that you would like to use
################################################################################################
env_settings_path="$HOME/.phycas/phycas_gui_env.sh"
if test -f "$env_settings_path"
then
  source "$env_settings_path" || exit 1
fi

if ! test -z "$PHYCAS_PYTHONPATH"
then
	if test -z "$PYTHONPATH"
	then
		PYTHONPATH="$PHYCAS_PYTHONPATH"
	else
		PYTHONPATH="$PYTHONPATH:$PHYCAS_PYTHONPATH"
	fi
fi	
if ! test -z "$PHYCAS_DYLD_LIBRARY_PATH"
then
	if test -z "$DYLD_LIBRARY_PATH"
	then
		DYLD_LIBRARY_PATH="$PHYCAS_DYLD_LIBRARY_PATH"
	else
		DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$PHYCAS_DYLD_LIBRARY_PATH"
	fi
fi	
if ! test -z "$PHYCAS_PATH"
then
	if test -z "$PATH"
	then
		PATH="$PHYCAS_PATH"
	else
		PATH="$PATH:$PHYCAS_PATH"
	fi
fi	
export PYTHONPATH
export DYLD_LIBRARY_PATH
export PATH

export USES_I_PYTHON=1
