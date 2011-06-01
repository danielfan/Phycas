This analysis estimated the log marginal likelihood separately for each of the 105 possible trees for the 6 taxon Rokas data set. This allows the overall log marginal likelihood to be computed for comparision with the multi-tree-topology steppingstone approach. The Jukes-Cantor model was used to minimize the number of variables.

On Mark's cluseter (kwyjibo), ran this command:

qsub -t 1-105 brute6taxon.sh

After the 105 runs finished, used this command

tail -n 1 ss.sump.*

to spit out last line of each of the output files (last line gives the log marginal likelihood). Applied this BBEdit regex to convert to format expected by summarize.py:

Find: ==> ss.sump.(\d+).txt <==\r\s(-\d+.\d+) Generalized Stepping Stone method\r

Replace: \1  \2

Manifest
--------
_readme.txt --> this file
rokas6first1000.nex --> first 1000 sites of the Rokas et al. (2003) yeast data set
brute6taxon.sh --> shell script designed for submission to kwyjibo's SGE
brute6taxon.py --> the Phycas script to do the analyses
data_first1000jc.txt --> tree number and log marginal likelihood
summarize.py --> computes overall marginal likelihood from data_first1000jc.txt