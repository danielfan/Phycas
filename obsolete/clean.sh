#!/bin/sh

rm -f ./phypy/*.pyc 

echo "Removing so files from ./phypy/Conversions"
rm ./phypy/Conversions/_Conversions.so*

echo "Removing so files from ./phypy/DataMatrix"
rm ./phypy/DataMatrix/_DataMatrixBase.so*

echo "Removing so files from ./phypy/Likelihood"
rm ./phypy/Likelihood/_LikelihoodBase.so*

echo "Removing so files from ./phypy/Phylogeny"
rm ./phypy/Phylogeny/_Phylogeny.so*

echo "Removing so files from ./phypy/ProbDist"
rm ./phypy/ProbDist/_ProbDist.so*

echo "Removing so files from ./phypy/ReadNexus"
rm ./phypy/ReadNexus/_ReadNexus.so*

echo "Done."
