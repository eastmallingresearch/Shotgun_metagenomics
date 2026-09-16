#!/bin/bash -l
#SBATCH -J GUNC
#sBATCH --mem=40000
#SBATCH -o GUNC_"%j".out

BINDIR=$1;shift
GUNCDB=$1;shift
OUT=$1;shift
CORES=$1;shift


# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}


echo "Running GUNC"

gunc run \
  --input_dir $BINDIR/ \
  --file_suffix .fa \
  --db_file $GUNCDB \
  --threads $CORES \
  --out_dir .


mkdir -p $OUT

cp -r * $OUT/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
