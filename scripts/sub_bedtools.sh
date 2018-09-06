#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G


DIR=$1;shift

#GFF=$1;shift
SF=$1;shift
BAM=$1;shift
OUTDIR=$1;shift

GFF=$(sed -n -e "$SGE_TASK_ID p" $DIR/L${SF}_splitfiles.txt)
F=$(sed 's/.*\///' <<<$BAM|sed 's/\..*//')
G=$(sed 's/.*\///' <<<$GFF)

cd $TMP

bedtools coverage -b $BAM -a $GFF -counts >$F.$G.cov

cp $F.$G.cov $OUTDIR/. 