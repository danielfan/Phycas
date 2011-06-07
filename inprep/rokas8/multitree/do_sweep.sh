#!/bin/sh
cv="$1"
if test -z $cv
then
    echo "Expecting a minimum coefficient of variation to apply to the parameters in the reference distribution"
    exit 1
fi
cap="$2"
if test -z $cap
then
    echo "Expecting a maximum split frequency to apply to every internal edge in the reference distribution"
    exit 1
fi
set -x
for s in 1 0.999 0.99 0.98 0.97 0.96 0.95 0.9 0.8 0.7 0.6 0.5 0.25 0.1
do
    dir="mincv${cv}_scale${s}_cap${cap}"
    if ! test -d "$dir"
    then
        mkdir "${dir}" || exit 
    fi
    fn="${dir}/refdist.txt"
    head -n1 ../pilot/refdist.txt | python 1010511/alter_edge_lengths.py -i --scale=${s} --cap=${cap} | awk '{print $1}' > "${fn}" || exit
    cat ../pilot/refdist_params.txt >> "${fn}" || exit
    # master_edgelen is not used, but currently required - @TODO
    echo 'master_edgelen = Gamma(1.0, 1.0)' >> "${fn}" || exit
    cat ../pilot/sump.refdist.txt >> "${fn}" || exit
done
set +x
echo "Assured that parameters had a coefficient of variation of at least ${cv}"
echo "and capped the upper limit on a split probability at ${cap}"
