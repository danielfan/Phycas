#!/bin/sh
active_env_path="$HOME/.phycas/active_phycas_env.sh"
if test -f "$active_env_path"
then
  source "$active_env_path" || exit 1
fi
if test -z "$PYTHONINTERPRETER"
then
   PYTHONINTERPRETER=python
fi
if test -z "$USES_I_PYTHON"
then
	echo
	echo "Testing for IPython..."
	echo
	echo "You can decrease the amount of time that phycas takes to launch by adding the line:"
	echo "   export USES_I_PYTHON=1"
	echo "or"
	echo "   export USES_I_PYTHON=0"
	echo "to the end of the file ~/.phycas/active_phycas_env.sh"
	$PYTHONINTERPRETER -c "import IPython" >/dev/null 2>&1
	if test $? -eq 0
	then
		USES_I_PYTHON=1
	else
		USES_I_PYTHON=0
	fi
fi
if test $# -eq 0
then
	if test "$USES_I_PYTHON" = "1"
	then
		$PYTHONINTERPRETER -c "from IPython.Shell import IPShell ; IPShell(['-i', '-c','from phycas import *']).mainloop()"
	else
		$PYTHONINTERPRETER -i -c "from phycas import *"
	fi
else
	if test $# -eq 1
	then
		if test "$USES_I_PYTHON" = "1"
		then
			$PYTHONINTERPRETER -c "from IPython.Shell import IPShell ; IPShell(['-i', \"$1\"]).mainloop()"
		else
			$PYTHONINTERPRETER -i $1
		fi
	else
		$PYTHONINTERPRETER -i $@
	fi
fi
