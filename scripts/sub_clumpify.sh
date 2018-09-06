#!/bin/bash
# clumpify
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=2G
#$ -pe smp 12

SCRIPT_DIR=$1; shift
OUTDIR=$1; shift

# change to session temp folder
cd $TMP
zcat $@ > temp.fq

NAME=$(sed 's/.*\///' <<<$1)

clumpify.sh \
 in=temp.fq \
 out=$NAME.fq.gz \

cp *.fq.gz $OUTDIR/.

