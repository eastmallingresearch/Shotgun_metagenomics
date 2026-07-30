#!/bin/bash -l
#SBATCH -J checkm2
#sBATCH --mem=40000
#SBATCH -o checkm2_"%j".out

BINS=$1;shift
DB=$1;shift
OUT=$1;shift
SAMPLE=$1;shift
CORES=$1;shift


# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# make MAGScot output directory

mkdir checkm2

echo "Running checkm2"

pixi run --manifest-path $APPS/checkm2/pixi.toml \
  checkm2 predict --threads $CORES --input $BINS/ --output-directory checkm2 --database_path $DB -x fa
  #/path/to/shared/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd
   
mkdir -p $OUT/

cp -r checkm2 $OUT/$SAMPLE/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
