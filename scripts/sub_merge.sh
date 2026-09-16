#!/bin/bash -l
#SBATCH -J merge
#sBATCH --mem=10000
#SBATCH -o merge"%j".out

#SCRIPT_DIR=$1; shift
OUTDIR=$1; shift
FORWARD=$1; shift
REVERSE=$1; shift
OUTPREFIX=$1; shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

F=$(sed 's/.*\///' <<<$FORWARD)
R=$(sed 's/.*\///' <<<$REVERSE)

~/pipelines/common/functional-analysis/carnelian/util/fileMerge \
 . \
$OUTPREFIX \
 .fq.gz \
 $FORWARD \
 $REVERSE 

mkdir -p $OUTDIR

cp * $OUTDIR/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
