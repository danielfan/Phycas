#! /bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -q general.q
export PYTHONPATH="/share/apps/lib/python2.4/site-packages"
export LD_LIBRARY_PATH="/opt/gridengine/lib/lx26-amd64:/opt/gridengine/lib/lx26-amd64:/share/apps/lib/python2.4/site-packages/phycas/Conversions"
python brute6taxon.py $SGE_TASK_ID

