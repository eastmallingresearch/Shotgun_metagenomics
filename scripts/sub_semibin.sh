#!/bin/bash -l
#SBATCH -J semibin
#sBATCH --mem=60000
#SBATCH -o semibin_"%j".out

ASSEMBLY=$1;shift
SAMPLE=$1;shift
BAMS=$1;shift
OUT=$1;shift
CORES=$1;shift
#BAMS=$1;shift


# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}


PIXI_PROJECT=$APPS/semibin2

pixi run --manifest-path $PIXI_PROJECT/pixi.toml \
SemiBin2 single_easy_bin \
  -i $ASSEMBLY \
  -b $BAMS \
  -o semibin2.bin \
  --sequencing-type short_read \
  --engine gpu \
  -t $CORES
 
mkdir -p $OUT/$SAMPLE

cp -r * $OUT/$SAMPLE/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
