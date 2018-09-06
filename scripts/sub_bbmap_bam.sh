#!/bin/bash
# megahit
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G

PROC=$1; shift
OUTDIR=$1; shift
PREFIX=$1; shift
REF=$1; shift
FORWARD=$1; shift
REVERSE=${1:-NOTHING}; shift

# change to session temp folder
cd $TMP

F=$(sed 's/.*\///' <<<$FORWARD)

bbmap.sh \
 in1=$FORWARD \
 in2=$REVERSE \
 out=$F.bam \
 t=$PROC $@ -Xmx31g
#| samtools sort -O bam -o ${PREFIX}.${F}.bam -T TEMP_${PREFIX} -@ $PROC -m 2G

mkdir -p $OUTDIR

cp * $OUTDIR/.

