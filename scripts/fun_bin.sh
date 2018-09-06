#!/bin/bash

#########################################################
#														#
# hmm search of database using functionalAnnotation.py	#
#														#
#########################################################

PROC=$1;shift
LOC=$1;shift	
FILES=$1;shift
DATABASE=$1;shift

SCRIPT_DIR=$(readlink -f ${0%/*})

cd $LOC
dir=`mktemp -d -p $LOC`
for f in $FILES
do
#	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  <$f > ${f}2
#	sed -i -e '1d' ${f}2
#	sed -i -e 's/ .*//' ${f}2	
	split -l 20000 $f -a 4 -d $dir/${f}.
	cd $dir
	find $PWD -name "$f*" >L${f}_splitfiles.txt
	TASKS=$(wc -l L${f}_splitfiles.txt|awk -F" " '{print $1}')
    qsub -pe smp $PROC -t 1-$TASKS:1 -tc 50 -h $SCRIPT_DIR/sub_fun_bin.sh $PROC $LOC $dir $f $DATABASE $@
	cd ..
done	

