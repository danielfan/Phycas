This directory contains results from several different multi-tree steppingstone
    sampling runs.  The focal tree in the reference distribution in these
    runs differs by setting the internal branch lengths to be x*scaler where
    x is the probability estimated in the refdist by the pilot run.
    
These are produced by:
    1. git clone git://gist.github.com/1010511.git
    
    2. running commands like:

head -n1 ../pilot/refdist.txt | python gist-1010511/alter_edge_lengths.py -i --scale=0.5 | awk '{print $1}' > scale0.5/refdist.txt
cat ../pilot/refdist_params.txt >> scale0.5/refdist.txt


    
    
