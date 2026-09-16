#!/bin/bash -l
#SBATCH -J metabat
#sBATCH --mem=40000
#SBATCH -o metabat_"%j".out

ASSEMBLY=$1;shift
DEPTH=$1;shift
SAMPLE=$1;shift
OUT=$1;shift
MINL=$1;shift
CORES=$1;shift
SEED="${1:-$RANDOM}"
#BAMS=$1;shift


# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

metabat2 -i $ASSEMBLY -a $DEPTH -m 2500 -t $CORES -o out.bin --seed $SEED --unbinned

mkdir -p $OUT/$SAMPLE

cp -r * $OUT/$SAMPLE/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
