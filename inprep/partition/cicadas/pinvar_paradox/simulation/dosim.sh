#!/bin/bash

./seq-gen -l100  -n1 -s100 -mHKY -t2 -z730107 -x paupblock.txt -on < tree_with_names.txt > small.nex
./seq-gen -l5000 -n1 -s1   -mHKY -t2 -z135799 -x paupblock.txt -on < tree_with_names.txt > large.nex
cat large.nex small.nex > ../analysis/concatenated.nex
