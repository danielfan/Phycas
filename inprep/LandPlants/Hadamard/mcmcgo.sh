#!/bin/bash

export PYTHONPATH="$HOME/Documents/Projects/pdevgit/Phycas"
export DYLD_LIBRARY_PATH="$PYTHONPATH/phycas/Conversions"

python hadamardMCMC.py

