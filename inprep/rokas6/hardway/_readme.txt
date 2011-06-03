This analysis estimated the log marginal likelihood separately for each of the 105 possible trees for the 6 taxon Rokas data set. This allows the overall log marginal likelihood to be computed for comparision with the multi-tree-topology steppingstone approach. The Jukes-Cantor model was used to minimize the number of variables.

On Mark's cluster (kwyjibo), ran this command:

qsub -t 1-105 brute6taxon.sh

After the 105 runs finished, used this command

tail -n 1 ss.sump.*

to spit out last line of each of the output files (last line gives the log marginal likelihood). Applied this BBEdit regex to convert to format expected by summarize.py:

Find: ==> ss.sump.(\d+).txt <==\r\s(-\d+.\d+) Generalized Stepping Stone method\r

Replace: \1  \2

Saved the 105 lines thus produced to the file data_first1000jc.txt.

Finally, ran summarize.py to read in data_first1000jc.txt. Here is the estimate (the majority of the file, showing the log marginal likelihood and posterior probability of each of the 105 trees, has been omitted):

Overall log marginal likelihood (based on 105 distinct tree topologies):
  log(marginal likelihood) = -4182.55668757

Manifest
--------
_readme.txt 
--> this file

rokas6first1000.nex 
--> first 1000 sites of the Rokas et al. (2003) yeast data set

brute6taxon.sh 
--> shell script designed for submission to kwyjibo's SGE

brute6taxon.py 
--> the Phycas script to do the analyses

data_first1000jc.txt 
--> tree number and log marginal likelihood

summarize.py 
--> computes overall marginal likelihood from data_first1000jc.txt

output.96.txt 
--> output from the analysis of tree 96, which has the highest posterior probability (needed for the reference distribution details, which are copied below)

Reference distribution details (from tree 96):
  edgelen_1003 = Gamma(53.94841, 0.00116)
  edgelen_1 = Gamma(25.43606, 0.00130)
  edgelen_1002 = Gamma(30.21331, 0.00118)
  edgelen_2 = Gamma(77.51640, 0.00118)
  edgelen_1001 = Gamma(25.52571, 0.00170)
  edgelen_3 = Gamma(33.92934, 0.00153)
  edgelen_1000 = Gamma(17.87053, 0.00211)
  edgelen_4 = Gamma(63.23979, 0.00150)
  edgelen_5 = Gamma(173.60092, 0.00164)
