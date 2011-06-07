#!/bin/sh
set -x
for s in 1 0.999 0.99 0.98 0.97 0.96 0.95 0.9 0.8 0.7 0.6 0.5 0.25 0.1
do
    if ! test -d scale${s}
    then
        mkdir scale${s}
    fi
    head -n1 ../pilot/refdist.txt | python 1010511/alter_edge_lengths.py -i --scale=${s} | awk '{print $1}' > scale${s}/refdist.txt
    cat ../pilot/refdist_params.txt >> scale${s}/refdist.txt
    echo 'master_edgelen = Gamma(1.0, 1.0)' >> scale${s}/refdist.txt
    cat ../pilot/sump.refdist.txt >> scale${s}/refdist.txt
done

