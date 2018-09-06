#!/bin/bash
# Mnimio
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=2G

OUTDIR=$1; shift

cd $TMP

Minimo $@

cp *.fa $OUTDIR/.