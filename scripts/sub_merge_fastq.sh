#!/bin/bash -l
#SBATCH -J mergefastq
#sBATCH --mem=500G
#SBATCH -o merge_fastq_"%j".out

S=$1;shift
R=$1;shift
OUT=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# micromamba run -n binchicken
echo Sample $S Read $R

cat ${S}*_L*_${R}.* > cleaned.fastq.gz

echo "catted"

mkdir -p $OUT

X=$(basename $S)

cp cleaned.fastq.gz $OUT/${X}_${R}.cleaned.fastq.gz

cd ..

rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
