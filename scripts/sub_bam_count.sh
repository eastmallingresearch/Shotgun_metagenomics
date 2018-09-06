#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G
#$ -pe smp 2

SCRIPT_DIR=$1;shift
BAM=$1;shift
GFF=$1;shift
OUT=$1;shift

F=$(sed 's/.*\///' <<<$BAM|sed 's/\..*//')

cd $TMP

samtools view $BAM|$SCRIPT_DIR/bam_scaffold_count.pl $GFF $@ > $F.cov

mkdir -p $OUT

cp * $OUT/.