#!/bin/bash
# megahit
#$ -S /bin/bash
#$ -cwd

OUTDIR=$1;shift
PREFIX=$1;shift

cd $TMP

cap3 $@ >/dev/null

#cp * $OUTDIR/$PREFIX/.