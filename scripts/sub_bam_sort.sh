#!/bin/bash -l
#SBATCH -J bam_sort
#sBATCH --mem=25000
#SBATCH -o bam_sort_"%j".out

PROC=$1;shift
OUTDIR=$1;shift
PREFIX=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# execute sorting 
samtools sort -O BAM -o $PREFIX -T $PREFIX -@ $PROC $@

cp $PREFIX $OUTDIR/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
