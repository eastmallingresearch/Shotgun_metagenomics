#!/bin/bash
# megahit
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G

PROC=$1;shift
OUTDIR=$1;shift
PREFIX=$1;shift

# change to session temp folder
cd $TMP

samtools -O BAM -o $PREFIX -T $PREFIX -@ $PROC $@

cp $PREFIX.bam $OUTDIR/.
