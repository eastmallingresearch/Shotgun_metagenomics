#!/bin/bash
# megahit
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 24

OUTDIR=$1;shift

# change to session temp folder
cd $TMP

cp -L $@ .
gunzip *
for f in *;do e=$(echo $e,$f);done
e=$(echo $e|sed 's/,//')

megahit \
 -o $TMP/output \
 -t 23 \
 --kmin-1pass \
 --out-prefix everything \
 --k-min=27 --k-step 10 --k-max 127 \
 -r $e

mkdir -p $OUTDIR/everything

cd $TMP/output

pigz -9 -p 12 -r *

echo cp -r * $OUTDIR/everything/

