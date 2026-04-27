#!/bin/bash -l
#SBATCH -J carnelian
#sBATCH --mem=10000
#SBATCH -o carnelian"%j".out

#SCRIPT_DIR=$1; shift
OUTDIR=$1; shift
FA_DIR=$1; shift
MODEL_DIR=$1; shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

conda activate carnelian

~/pipelines/common/functional-analysis/carnelian/carnelian.py annotate \
$@ \
$FA_DIR \
$MODEL_DIR \
. \
FragGeneScan

mkdir -p $OUTDIR

cp -r * $OUTDIR/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
