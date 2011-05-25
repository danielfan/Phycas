# This file was originally used to test SAMC, but now there is SAMCTest for that

import sys,os
from phycas import *

p = Phycas()
try:
    d = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p.data_file_name = getPhycasTestData("green.nex")
except NameError:
    p.data_file_name = os.path.join("phycas", "Tests","Data", "green.nex")
p.default_model = 'jc'

# Create a pairwise distance matrix. Here is the matrix PAUP* computes for comparison:
#
#                        1        2        3        4        5        6        7        8        9
#  1 Chara               -
#  2 Conocephalum  0.17472        -
#  3 Bazzania      0.20250  0.14700        -
#  4 Sphagnum      0.21371  0.14325  0.16794        -
#  5 Osmunda       0.23451  0.18157  0.17084  0.16890        -
#  6 Picea         0.23557  0.19147  0.17961  0.16025  0.17375        -
#  7 Iris          0.23451  0.18650  0.19147  0.17863  0.18059  0.14325        -
#  8 Asplenium     0.26795  0.23451  0.22300  0.22092  0.18552  0.20352  0.21885        -
#  9 Nicotiana     0.23874  0.19346  0.20149  0.17570  0.18256  0.15265  0.10590  0.22508        -
# 10 Avena         0.23346  0.19446  0.21166  0.17765  0.18948  0.15360  0.11574  0.21473  0.13674
dmatrix = p.pairwiseDistanceMatrix()    # this function does not yet exist

# Build a stepwise-addition tree using the specified addition sequence.
# Here is PAUP's stepwise addition ME tree:
#
# Heuristic search completed
#    Total number of rearrangements tried = 0
#    Score of best tree(s) found = 0.83775
#    Number of trees retained = 1
#    Time used = 0.01 sec
#
#   Translate
#       1 Chara,
#       2 Conocephalum,
#       3 Bazzania,
#       4 Sphagnum,
#       5 Osmunda,
#       6 Picea,
#       7 Iris,
#       8 Asplenium,
#       9 Nicotiana,
#       10 Avena
#       ;
# tree PAUP_1 = [&U] (1:0.111657,2:0.063067,(3:0.076958,(4:0.071965,((5:0.072360,8:0.113155):0.013395,(6:0.069457,((7:0.049397,9:0.056506):0.008506,10:0.064784):0.017809):0.014771):0.011382):0.012152):0.010430);

addseq = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
print p.stepwiseAdditionMETree(dmatrix, addseq)    # this function does not yet exist
