#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G

read -r -d '' HELP << EOM
#########################################################################
#																		#
#	Wrapper script for partioning (clumping) reads						#
#																		#
#	usage: partition.sh -p <program> [options]							#
#																		#
#	-p <clumpify> 													    #
#																		#
#	partition.sh -p clumpify Outdir merged_reads unmerged_reads 	 	#
#	 																	#
#	 																	#
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

clumpify|clumify|clump)
	qsub $SCRIPT_DIR/sub_clumpify.sh $SCRIPT_DIR $@
	exit 0
;;
*)
	echo "Invalid merge program: $program" >&2
	exit 1
esac
