This analysis estimated the log marginal likelihood separately for each of the 105 possible trees for the 6 taxon Rokas data set. This allows the overall log marginal likelihood to be computed for comparision with the multi-tree-topology steppingstone approach. The Jukes-Cantor model was used to minimize the number of variables.

The newway.py script in this folder runs a brief multi-tree stepping-stone analysis on the first 1000 sites of the Rokas et al. 2003 data using the JC model. It should get a value close to -4182.55668757, which was obtained usiong the brute force approach (see hardway folder).

Manifest
--------
_readme.txt 
--> this file

newway.py 
--> Phycas script that does the analysis

newway.sh 
--> shell script used to run newway.py on kwyjibo (Mark Holder's cluster) - not really necessary since this runs on a laptop in a couple of minutes

rokas6first1000.nex 
--> the data file (first 1000 sites of the Rokas et al. 2003 yeast data set)

sumt_ref_dist.txt 
--> the reference distribution file (see stdway folder for details)
