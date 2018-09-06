#!/bin/bash
#$ -S /bin/bash
#$ -cwd

SCRIPT_DIR=$(readlink -f ${0%/*})

#==========================================================
#	Set Help (add message between EOM blocks
#==========================================================	
read -r -d '' HELP << EOM
#############################################################
#                                                           #
#	Metagenomics pipeline	for Illumina data               #
#                                                           #
#	usage: PIPELINE.sh -c <program> [options]               #
#                                                           #
#############################################################

 -c <program>	Program can be any of the defined programs
 -h		display this help and exit	
EOM


function print_help {
	echo;echo "$HELP" >&1;echo;
	exit 1
}

if [ $# -eq 0 ];
then
   print_help
fi

#==========================================================
#	Set command line switch options
#==========================================================

OPTIND=1 

while getopts ":hs:c:" options; do
	case "$options" in
	s)
	    SCRIPT_DIR=$OPTARG
	    ;;
	c)  
 	    program=$OPTARG
	    break
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    exit 1
      	    ;;
	h)  
	    print_help
	    exit 0
 	    ;;
	?) 
	    echo "Invalid option: -$OPTARG" >&2
	    echo "Call PIPELINE with -h switch to display usage instructions"
	    exit 1
	    ;;
	esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

#==========================================================
#	Set (-c) programs
#==========================================================

case $program in

wibble)	
	echo $SCRIPT_DIR 
	echo "Wibbling..."
	exit 1
;;
split|splitfq|splitfq.sh)
	$SCRIPT_DIR/splitfq.sh $@
	exit 0
;;
trim|trim.sh)
	$SCRIPT_DIR/trim.sh $@
	exit 0
;;
clean|clean.sh)
	$SCRIPT_DIR/clean.sh $SCRIPT_DIR $@
	exit 0
;;
filter|filter.sh)
	$SCRIPT_DIR/filter.sh $@
	exit 0
;;
concat|concatonate|cat|concat.sh)
	qsub $SCRIPT_DIR/submit_concat.sh $@
	exit 0
;;
normalise|normalise.sh)
	$SCRIPT_DIR/normalise.sh $@
	exit 0
;;
dereplicate|dereplicate.sh)
	$SCRIPT_DIR/dereplicate.sh $@
	exit 0
;;
correct)
	qsub $SCRIPT_DIR/submit_correct.sh $SCRIPT_DIR $@
	exit 0
;;
interleave|inter)
	qsub $SCRIPT_DIR/submit_interleave.sh $@
	exit 0
;;
assemble|assemble.sh)
	$SCRIPT_DIR/assemble.sh $@
	exit 0
;;
post|post.sh)
	qsub $SCRIPT_DIR/submit_post.sh $SCRIPT_DIR $@
	exit 0
;;
merge|merge.sh)
	$SCRIPT_DIR/merge.sh $@
	exit 0
;;
partition|partition.sh)
	$SCRIPT_DIR/partition.sh $@
	exit 0
;;
align|align.sh)
	$SCRIPT_DIR/align.sh $@
	exit 0
;;
coverage|coverage.sh)
	$SCRIPT_DIR/coverage.sh $@
	exit 0
;;
cov_bed)
	SERVER=$1;shift
	qsub -l h=$SERVER $SCRIPT_DIR/sub_cov_bed.sh $SCRIPT_DIR $@
	exit 0
;;
cluster_super_fast)
  SERVER=$1;shift
  CORES=$1;shift
  FILELOC=$1
  TASKS=$(ls $FILELOC/*.fasta|wc -l)
  qsub -l h=$SERVER -t 1-$TASKS:1 -tc $CORES $SCRIPT_DIR/sub_cluster_super_fast.sh $@
  exit 0
;;
megafilt|MEGAFILT)
    qsub -l h=balcklace01,h=blacklace11 $SCRIPT_DIR/mega_duk.sh $SCRIPT_DIR $@
   #$SCRIPT_DIR/mega_duk.sh $SCRIPT_DIR $@
	exit 0
;;
mega|MEGA)
	
	STRAW_DN=$1; shift
	SAMPLE=$1; shift
	READS=$1; shift

	JOBNAME=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n1)

	qsub -N MEGA_SPLIT_${JOBNAME}_1 $SCRIPT_DIR/submit_splitfq.sh $STRAW_DN/../raw/${SAMPLE}_1.fq.gz $READS $STRAW_DN/split
	  
 	qsub -N MEGA_SPLIT_${JOBNAME}_2 $SCRIPT_DIR/submit_splitfq.sh $STRAW_DN/../raw/${SAMPLE}_2.fq.gz $READS $STRAW_DN/split

 	qsub -hold_jid MEGA_SPLIT_${JOBNAME}_1 -hold_jid MEGA_SPLIT_${JOBNAME}_2 $SCRIPT_DIR/MEGA_BEHOLD.sh $SCRIPT_DIR $STRAW_DN $SAMPLE $JOBNAME

	exit 0
;;
TEST)
	echo "test program run with options:" $@
	exit 0
;;

*)
	echo "Invalid program: $program" >&2
	exit 1
esac
