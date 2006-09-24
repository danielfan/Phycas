#!/bin/sh

function abort()
{
echo "Shell script aborted because of failed example"
exit 1
}

TESTDIR=`pwd`
rm -f $TESTDIR/runall_diffs.txt
cd ../Examples

echo "****************************" >> $TESTDIR/runall_diffs.txt
echo "*** Running ExplorePrior ***" >> $TESTDIR/runall_diffs.txt
echo "****************************" >> $TESTDIR/runall_diffs.txt
echo " " >> $TESTDIR/runall_diffs.txt
cd ExplorePrior
python ExplorePrior.py
diff nodata.nex.p reference_output/nodata.nex.p >> $TESTDIR/runall_diffs.txt 
if [ $? -ne 0 ] 
then
  abort
fi
diff nodata.nex.t reference_output/nodata.nex.t >> $TESTDIR/runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..

echo " " >> $TESTDIR/runall_diffs.txt
echo "***************************" >> $TESTDIR/runall_diffs.txt
echo "*** Running FixedParams ***" >> $TESTDIR/runall_diffs.txt
echo "***************************" >> $TESTDIR/runall_diffs.txt
echo " " >> $TESTDIR/runall_diffs.txt
cd FixedParams
python FixedParams.py
diff params.p reference_output/params.p >> $TESTDIR/runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
diff trees.t reference_output/trees.t >> $TESTDIR/runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..

echo " " >> $TESTDIR/runall_diffs.txt
echo "****************************" >> $TESTDIR/runall_diffs.txt
echo "*** Running GelfandGhosh ***" >> $TESTDIR/runall_diffs.txt
echo "****************************" >> $TESTDIR/runall_diffs.txt
echo " " >> $TESTDIR/runall_diffs.txt
cd GelfandGhosh
python GelfandGhosh.py
diff ggout.txt reference_output/ggout.txt >> $TESTDIR/runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..

echo " " >> $TESTDIR/runall_diffs.txt
echo "******************************" >> $TESTDIR/runall_diffs.txt
echo "*** Running LikelihoodTest ***" >> $TESTDIR/runall_diffs.txt
echo "******************************" >> $TESTDIR/runall_diffs.txt
echo " " >> $TESTDIR/runall_diffs.txt
cd LikelihoodTest
python LikelihoodTest.py
diff simulated.nex reference_output/simulated.nex >> $TESTDIR/runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
diff check.nex reference_output/check.nex >> $TESTDIR/runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..

#echo " " >> $TESTDIR/runall_diffs.txt
#echo "**********************" >> $TESTDIR/runall_diffs.txt
#echo "*** Running Phycas ***" >> $TESTDIR/runall_diffs.txt
#echo "**********************" >> $TESTDIR/runall_diffs.txt
#echo " " >> $TESTDIR/runall_diffs.txt
#cd Phycas
#python Phycas.py
#diff green.nex.p reference_output/green.nex.p >> $TESTDIR/runall_diffs.txt 
#if [ $? -ne 0 ] 
#then
#  abort
#fi
#diff green.nex.t reference_output/green.nex.t >> $TESTDIR/runall_diffs.txt
#if [ $? -ne 0 ] 
#then
#  abort
#fi
#cd ..

echo " " >> $TESTDIR/runall_diffs.txt
echo "**************************" >> $TESTDIR/runall_diffs.txt
echo "*** Running Polytomies ***" >> $TESTDIR/runall_diffs.txt
echo "**************************" >> $TESTDIR/runall_diffs.txt
echo " " >> $TESTDIR/runall_diffs.txt
cd Polytomies
python Polytomies.py
diff analHKY.nex.p reference_output/analHKY.nex.p >> $TESTDIR/runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
diff analHKY.nex.t reference_output/analHKY.nex.t >> $TESTDIR/runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
diff simHKY.nex reference_output/simHKY.nex >> $TESTDIR/runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..

echo " " >> $TESTDIR/runall_diffs.txt
echo "*************************" >> $TESTDIR/runall_diffs.txt
echo "*** Running Simulator ***" >> $TESTDIR/runall_diffs.txt
echo "*************************" >> $TESTDIR/runall_diffs.txt
echo " " >> $TESTDIR/runall_diffs.txt
cd Simulator
python Simulator.py
diff simulated.nex reference_output/simulated.nex >> $TESTDIR/runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..
