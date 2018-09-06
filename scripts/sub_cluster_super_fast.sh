#!/bin/bash
# megahit
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=0G

INDIR=$1;shift
OUTDIR=$1;shift
ID=$1;shift

FILE=$(ls $INDIR|sed -n "$SGE_TASK_ID p")

cd $TMP

usearch -cluster_fast $INDIR/$FILE -id $ID -uc $FILE.uc

mkdir -p $OUTDIR

cp *.uc $OUTDIR/.