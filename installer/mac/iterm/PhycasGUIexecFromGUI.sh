#!/bin/sh
mac_os_dir=`dirname "$0"`
contents_os_dir=`dirname "$mac_os_dir"`
resources_dir="$contents_os_dir/Resources"
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
	PATH="$PATH:$mac_os_dir"
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
	
	run_cmd_path="$HOME/.phycas/run_phycas.sh"
	if ! test -f "$run_cmd_path"
	then
		echo "#!/bin/sh" > "$run_cmd_path"
		echo "env_settings_path=\$HOME/.phycas/phycas_gui_env.sh" >> "$run_cmd_path"
		echo "if test -f \$env_settings_path" >> "$run_cmd_path"
		echo "then" >> "$run_cmd_path"
		echo "  source \$env_settings_path || exit 1" >> "$run_cmd_path"
		echo "fi" >> "$run_cmd_path"
		echo "if test -z \"\$PYTHONINTERPRETER\"" >> "$run_cmd_path"
		echo "then" >> "$run_cmd_path"
		echo "  if test -z \"\$USES_I_PYTHON\"" >> "$run_cmd_path"
		echo "  then" >> "$run_cmd_path"
		echo "   PYTHONINTERPRETER=python" >> "$run_cmd_path"
		echo "  else" >> "$run_cmd_path"
		echo "    PYTHONINTERPRETER=ipython" >> "$run_cmd_path"
		echo "  fi" >> "$run_cmd_path"
		echo "fi" >> "$run_cmd_path"
		echo "\$PYTHONINTERPRETER -i -c \"from phycas import *\"" >> "$run_cmd_path"
		chmod +x "$run_cmd_path"
	fi
	if ! test -f "$run_cmd_path"
	then
		echo "The mandatory file $run_cmd_path does not exist, and cannot be created (perhaps there is a directory in the way or you do not have permissions)."
		exit 1
	fi
	$run_cmd_path
else
	echo "The mandatory directory $HOME/.phycas does not exist, and cannot be created (perhaps there is a file in the way or you do not have permissions to write to the directory $HOME)."
	exit 1
fi

