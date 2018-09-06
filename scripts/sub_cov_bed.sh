#!/bin/bash
# megahit
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G

SCRIPT_DIR=$1;shift
COV=$1;shift
OUTDIR=$1;shift

TAB=$(sed 's/.*\///' <<<$COV|sed 's/\..*//').tab

# change to session temp folder
cd $TMP

python $SCRIPT_DIR/cov_bed.py $COV $TAB

cp $TAB $OUTDIR/.
