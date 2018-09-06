#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G

read -r -d '' HELP << EOM
#########################################################################
#                                                                       #
#   Wrapper script for calculating coverage                             #
#                                                                       #
#   usage: coverage.sh -p <program> [options]                           #
#                                                                       #
#   -p <metaspades|megahit>                                             #
#                                                                       #
#   assemble.sh -p metaspades Forward Reverse Output PREFIX             #
#   assemble.sh -p megahit Output PREFIX -1 Forward -2 Reverse          #
#                                                                       #
#########################################################################
EOM

function print_help {
	echo;echo "$HELP" >&1;echo;
	exit 0
}

if [ $# -eq 2 ];
then
   print_help
fi

OPTIND=1

while getopts ":hsp:" options; do
	case "$options" in
	s)
	  SCRIPT_DIR=$OPTARG
	  ;;
	p)
	  program=$OPTARG
	  break
	  ;;
	h)  
	  print_help
	  exit 0
	  ;;
	esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if [[ -z "$SCRIPT_DIR" ]]; then 
	SCRIPT_DIR=$(readlink -f ${0%/*})
fi


case $program in

bam_count.sh|bam_count)
	REGSERVER=$1;shift
	qsub -l h=$REGSERVER $SCRIPT_DIR/sub_bam_count.sh $SCRIPT_DIR $@
	exit 0
;;
bedtools.sh|bedtools)	
	REGSERVER=$1;shift
	LOC=$1;shift
	GFF=$1;shift
	dir=`mktemp -d -p $LOC`
	cd $LOC
	split -l 100000 $GFF -a 4 -d $dir/${GFF}.
	cd -
	cd $dir
	find $PWD -name "$GFF*" >L${GFF}_splitfiles.txt
	TASKS=$(wc -l L${GFF}_splitfiles.txt|awk -F" " '{print $1}')
	qsub -l h=$REGSERVER -t 1-$TASKS:1 -tc 25  $SCRIPT_DIR/sub_bedtools.sh $dir $GFF $@
	#$SCRIPT_DIR/sub_bedtools.sh $dir $GFF $@
	cd -
	exit 0
;;
*)
	echo "Invalid assembly program: $program" >&2
	exit 1
esac
