#!/bin/sh

function abort()
{
echo "Shell script aborted because of failed example"
exit 1
}

rm -f runall_diffs.txt

echo "****************************" >> runall_diffs.txt
echo "*** Running ExplorePrior ***" >> runall_diffs.txt
echo "****************************" >> runall_diffs.txt
echo " " >> runall_diffs.txt
cd ExplorePrior
python ExplorePrior.py
diff nodata.nex.p reference_output/nodata.nex.p >> ../runall_diffs.txt 
if [ $? -ne 0 ] 
then
  abort
fi
diff nodata.nex.t reference_output/nodata.nex.t >> ../runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..

echo " " >> runall_diffs.txt
echo "***************************" >> runall_diffs.txt
echo "*** Running FixedParams ***" >> runall_diffs.txt
echo "***************************" >> runall_diffs.txt
echo " " >> runall_diffs.txt
cd FixedParams
python FixedParams.py
diff params.p reference_output/params.p >> ../runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
diff trees.t reference_output/trees.t >> ../runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..

echo " " >> runall_diffs.txt
echo "****************************" >> runall_diffs.txt
echo "*** Running GelfandGhosh ***" >> runall_diffs.txt
echo "****************************" >> runall_diffs.txt
echo " " >> runall_diffs.txt
cd GelfandGhosh
python GelfandGhosh.py
diff ggout.txt reference_output/ggout.txt >> ../runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..

echo " " >> runall_diffs.txt
echo "******************************" >> runall_diffs.txt
echo "*** Running LikelihoodTest ***" >> runall_diffs.txt
echo "******************************" >> runall_diffs.txt
echo " " >> runall_diffs.txt
cd LikelihoodTest
python LikelihoodTest.py
diff simulated.nex reference_output/simulated.nex >> ../runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
diff check.nex reference_output/check.nex >> ../runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..

#echo " " >> runall_diffs.txt
#echo "**********************" >> runall_diffs.txt
#echo "*** Running Phycas ***" >> runall_diffs.txt
#echo "**********************" >> runall_diffs.txt
#echo " " >> runall_diffs.txt
#cd Phycas
#python Phycas.py
#diff green.nex.p reference_output/green.nex.p >> ../runall_diffs.txt 
#if [ $? -ne 0 ] 
#then
#  abort
#fi
#diff green.nex.t reference_output/green.nex.t >> ../runall_diffs.txt
#if [ $? -ne 0 ] 
#then
#  abort
#fi
#cd ..

echo " " >> runall_diffs.txt
echo "**************************" >> runall_diffs.txt
echo "*** Running Polytomies ***" >> runall_diffs.txt
echo "**************************" >> runall_diffs.txt
echo " " >> runall_diffs.txt
cd Polytomies
python Polytomies.py
diff analHKY.nex.p reference_output/analHKY.nex.p >> ../runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
diff analHKY.nex.t reference_output/analHKY.nex.t >> ../runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
diff simHKY.nex reference_output/simHKY.nex >> ../runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..

echo " " >> runall_diffs.txt
echo "*************************" >> runall_diffs.txt
echo "*** Running Simulator ***" >> runall_diffs.txt
echo "*************************" >> runall_diffs.txt
echo " " >> runall_diffs.txt
cd Simulator
python Simulator.py
diff simulated.nex reference_output/simulated.nex >> ../runall_diffs.txt
if [ $? -ne 0 ] 
then
  abort
fi
cd ..

