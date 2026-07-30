#!/bin/bash -l
#SBATCH -J dRep
#sBATCH --mem=40000
#SBATCH -o dRep_"%j".out

BINDIR=$1;shift
GINFO=$1;shift
OUT=$1;shift
CORES=$1;shift

echo $CORES

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}


echo "Running dRep"

dRep dereplicate drep_output \
  -g $BINDIR/*.fa \
  --genomeInfo $GINFO \
  --S_algorithm fastANI \
  -pa 0.9 \
  -sa 0.95 \
  -nc 0.30 \
  -comp 50 \
  -con 10 \
  -p $CORES


mkdir -p $OUT

cp -r * $OUT/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r	
