#!/bin/bash
# spades
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G

PROC=$1;shift
OUTDIR=$1; shift
FORWARD=$1; shift
REVERSE=$1; shift
PREFIX=$1;shift

# change to session temp folder
cd $TMP

metaspades.py \
 --meta \
 --only-assembler \
 -o $TMP \
 -1 $FORWARD \
 -2 $REVERSE \
 -t $PROC \
 $@ 

mkdir -p $OUTDIR/$PREFIX

pigz -9 -p 12 -r *

cp * -r  $OUTDIR/$PREFIX/