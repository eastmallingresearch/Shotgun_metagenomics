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

megahit \
 -o $TMP/output \
 -t $PROC \
 --kmin-1pass \
 --out-prefix $PREFIX \
 $@ 

mkdir -p $OUTDIR/$PREFIX

cd $TMP/output

pigz -9 -p 12 -r *

cp -r * $OUTDIR/$PREFIX/