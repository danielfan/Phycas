#!/bin/sh
mac_os_dir=`dirname "$0"`
contents_os_dir=`dirname "$mac_os_dir"`
resources_dir="$contents_os_dir/Resources"
PHYCAS_GUI_RESOURCES_DIR="$resources_dir"
export PHYCAS_GUI_RESOURCES_DIR
phycas_dir="$resources_dir/phycas"
if test -d "$phycas_dir"
then
	if test -z "$PYTHONPATH"
	then
		PYTHONPATH="$resources_dir"
	else
		PYTHONPATH="$PYTHONPATH:$resources_dir"
	fi
	if test -z "$DYLD_LIBRARY_PATH"
	then
		DYLD_LIBRARY_PATH="$resources_dir:$phycas_dir/Conversions"
	else
		DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:$resources_dir:$phycas_dir/Conversions"
	fi
	PATH="$mac_os_dir:$PATH"
else
	echo "$0: $phycas_dir does not exist!"
	exit 1
fi
env_settings_path="$HOME/.phycas/phycas_gui_env.sh"
if ! test -d "$HOME/.phycas"
then
	mkdir "$HOME/.phycas"
fi
if test -d "$HOME/.phycas"
then
	echo "#!/bin/sh" > "$env_settings_path"
	echo "PYTHONPATH=\"$PYTHONPATH\"" >> "$env_settings_path"
	echo "DYLD_LIBRARY_PATH=\"$DYLD_LIBRARY_PATH\"" >> "$env_settings_path"
	echo "PATH=\"$PATH\"" >> "$env_settings_path"
	echo "export PYTHONPATH" >> "$env_settings_path"
	echo "export DYLD_LIBRARY_PATH" >> "$env_settings_path"
	echo "export PATH" >> "$env_settings_path"
	active_env_path="$HOME/.phycas/active_phycas_env.sh"
	if ! test -f "$active_env_path"
	then
		echo "#!/bin/sh" > "$active_env_path"
		echo "################################################################################################" >> "$active_env_path"
		echo "# This file sourced whenever phycas is run through the Phycas.app " >> "$active_env_path"
		echo "# The Phycas.app executes ~/.phycas/run_phycas.sh which sources this file before invoking python." >> "$active_env_path"
		echo "# " >> "$active_env_path"
		echo "# The default behavior is for this file to source the file \"\$HOME/.phycas/phycas_gui_env.sh\" " >> "$active_env_path"
		echo "#  which is overwritten every time a Phycas.app window is opened so that it always reflects " >> "$active_env_path"
		echo "#  the current location of Phycas.app " >> "$active_env_path"
		echo "# " >> "$active_env_path"
		echo "# If you would like to augment the default behavior of Phycas you can:" >> "$active_env_path"
		echo "#     - define (and export) the environmental variable PYTHONINTERPRETER " >> "$active_env_path"
		echo "#       to point to the full path to the python executable that you would like to use. " >> "$active_env_path"
		echo "#     - add the commnd:" >> "$active_env_path"
		echo "# export USES_I_PYTHON=1" >> "$active_env_path"
		echo "#       to make the Phycas.app bundle take you into the IPython environment." >> "$active_env_path"
		echo "# " >> "$active_env_path"
		echo "# To use  different version of the phycas python package with the Phycas.app bundle" >> "$active_env_path"
		echo "# you can comment out the line" >> "$active_env_path"
		echo "#    source \"\$env_settings_path\"" >> "$active_env_path"
		echo "# below, and make sure that  PYTHONPATH and DYLD_LIBRARY_PATH are set correctly" >> "$active_env_path"
		echo "# for the version of phycas that you would like to use" >> "$active_env_path"
		echo "################################################################################################" >> "$active_env_path"
		echo >> "$active_env_path"
		echo "env_settings_path=\"\$HOME/.phycas/phycas_gui_env.sh\"" >> "$active_env_path"
		echo "if test -f \"\$env_settings_path\"" >> "$active_env_path"
		echo "then" >> "$active_env_path"
		echo "  source \"\$env_settings_path\" || exit 1" >> "$active_env_path"
		echo "fi" >> "$active_env_path"
	fi	
	run_cmd_path="$HOME/.phycas/run_phycas.sh"
	if ! test -f "$run_cmd_path"
	then
		echo "#!/bin/sh" > "$run_cmd_path"
		echo "active_env_path=\"\$HOME/.phycas/active_phycas_env.sh\"" >> "$run_cmd_path"
		echo "if test -f \"\$active_env_path\"" >> "$run_cmd_path"
		echo "then" >> "$run_cmd_path"
		echo "  source \"\$active_env_path\" || exit 1" >> "$run_cmd_path"
		echo "fi" >> "$run_cmd_path"
		echo "if test -z \"\$PYTHONINTERPRETER\"" >> "$run_cmd_path"
		echo "then" >> "$run_cmd_path"
		echo "   PYTHONINTERPRETER=python" >> "$run_cmd_path"
		echo "fi" >> "$run_cmd_path"
		echo "if test \$# -eq 0" >> "$run_cmd_path"
		echo "then" >> "$run_cmd_path"
		echo "    if ! test \"\$USES_I_PYTHON\" = \"1\"" >> "$run_cmd_path"
		echo "    then" >> "$run_cmd_path"
		echo "        \$PYTHONINTERPRETER -i -c \"from phycas import *\"" >> "$run_cmd_path"
		echo "    else" >> "$run_cmd_path"
		echo "        \$PYTHONINTERPRETER -c \"import IPython ; IPython.Shell.start().mainloop()\" -c \"from phycas import *\"" >> "$run_cmd_path"
		echo "    fi" >> "$run_cmd_path"
		echo "else" >> "$run_cmd_path"
		echo "    \$PYTHONINTERPRETER -i \$@" >> "$run_cmd_path"
		echo "fi" >> "$run_cmd_path"
		chmod +x "$run_cmd_path"
	fi
	if ! test -f "$run_cmd_path"
	then
		echo "The mandatory file $run_cmd_path does not exist, and cannot be created (perhaps there is a directory in the way or you do not have permissions)."
		exit 1
	fi
	if test $# -eq 0
	then
		"$run_cmd_path"
	else
		bd=`dirname $1`
		cd "$bd" || exit 1
		"$run_cmd_path" $@
	fi
else
	echo "The mandatory directory $HOME/.phycas does not exist, and cannot be created (perhaps there is a file in the way or you do not have permissions to write to the directory $HOME)."
	exit 1
fi

