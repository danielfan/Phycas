#!/bin/sh"
active_env_path="$HOME/.phycas/active_phycas_env.sh"
if test -f "$active_env_path"
then
  source "$active_env_path" || exit 1
fi
if test -z "$PYTHONINTERPRETER"
then
   PYTHONINTERPRETER=python
fi
if test $# -eq 0
then
    if ! test "$USES_I_PYTHON" = "1"
    then
        $PYTHONINTERPRETER -i -c "from phycas import *"
    else
        $PYTHONINTERPRETER -c "import IPython ; IPython.Shell.start().mainloop()" -c "from phycas import *"
    fi
else
    $PYTHONINTERPRETER -i $@
fi
