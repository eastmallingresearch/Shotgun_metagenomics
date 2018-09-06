#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G

PROC=$1;shift
LOC=$1;shift
DIR=$1;shift
SF=$1;shift
DATABASE=$1;shift

DATA=$(sed -n -e "$SGE_TASK_ID p" $DIR/L${SF}_splitfiles.txt)
RANDNAME=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n1)

cd $TMP

echo -e "name\treference\nX\t$DATA" >${RANDNAME}_metafile.txt

functionalAnnotation.py -m ${RANDNAME}_metafile.txt -db $DATABASE -n $PROC -o $RANDNAME $@

cp $RANDNAME $DIR/ -r