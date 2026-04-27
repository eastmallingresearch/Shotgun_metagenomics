#!/bin/bash -l
#SBATCH -J metabat
#sBATCH --mem=40000
#SBATCH -o metabat_"%j".out

A=$1;shift
OUT=$1;shift
#BAMS=$1;shift


# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

runMetaBat.sh  --unbinned -m 1500 -x 0 --minCVSum 0.5 $A $@

mkdir -p $OUT

cp -r * $OUT/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
