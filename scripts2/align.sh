#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G

read -r -d '' HELP << EOM
#########################################################################
#                                                                       #
#   Wrapper script for assmbling data                                   #
#                                                                       #
#   usage: align.sh -p <program> [options]                           #
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

bbmap|BBMap)	
	PROC=$1;shift
	REGSERVER=$1;shift
	qsub -pe smp $PROC -l h=$REGSERVER  $SCRIPT_DIR/sub_bbmap_bam.sh $PROC $@
	exit 0
;;
*)
	echo "Invalid assembly program: $program" >&2
	exit 1
esac
