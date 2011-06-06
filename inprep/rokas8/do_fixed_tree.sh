#!/bin/sh
if test -z $1
then
    echo "Expecting a tree number as an argument" 
    exit 1
fi
n="$1"
sd=$(expr $n % 100)
if ! test -d "sd${sd}"
then
    mkdir "sd${sd}" || exit 2
fi
if ! test -d "sd${sd}/tree${n}"
then
    mkdir "sd${sd}/tree${n}" || exit 2
fi
cd "sd${sd}/tree${n}"
tree_str=$(cat ../../alltrees.nex | grep "PAUP_${n} " | awk '{print $5}' | sed -e 's/;//')
if test -z ${tree_str}
then
    echo "Could not find tree PAUP_${n} in alltrees.nex"
    exit 3
fi
cat ../../rokas8.py | sed -e "s/rokas8T1000C.nex/..\/..\/rokas8T1000C.nex/" | sed -e "s/TREE_DESCRIPTION_HERE/\'${tree_str}\'/" > rokas8.py
python rokas8.py $1
