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
	echo "export PHYCAS_PYTHONPATH=\"$PYTHONPATH\"" >> "$env_settings_path"
	echo "export PHYCAS_DYLD_LIBRARY_PATH=\"$DYLD_LIBRARY_PATH\"" >> "$env_settings_path"
	echo "export PHYCAS_PATH=\"$PATH\"" >> "$env_settings_path"
	active_env_path="$HOME/.phycas/active_phycas_env.sh"
	if ! test -f "$active_env_path"
	then
		cp "$resources_dir/active_phycas_env.sh" "$active_env_path"
	fi	
	startup_py_path="$HOME/.phycas/startup.py"
	if ! test -f "$startup_py_path"
	then
		if test -f "$resources_dir/startup.py"
		then
			cp "$resources_dir/startup.py" "$startup_py_path"
		fi
	fi	
	run_cmd_path="$HOME/.phycas/run_phycas.sh"
	if ! test -f "$run_cmd_path"
	then
		cp "$resources_dir/run_phycas.sh" "$run_cmd_path"
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

