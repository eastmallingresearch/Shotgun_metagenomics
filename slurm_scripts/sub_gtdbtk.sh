#!/bin/bash -l
#SBATCH -J gtdbtk
#sBATCH --mem=40000
#SBATCH -o gtdbtk_"%j".out

BINDIR=$1;shift
OUT=$1;shift
CORES=$1;shift

echo $CORES

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}


echo "Running gtdbtk"

gtdbtk classify_wf \
  --genome_dir $BINDIR/ \
  --out_dir . \
  -x fa \
  --cpus $CORES

mkdir -p $OUT

cp -r * $OUT/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r	
