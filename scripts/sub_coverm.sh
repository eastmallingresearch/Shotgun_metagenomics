#!/bin/bash -l
#SBATCH -J coverm
#sBATCH --mem=40000
#SBATCH -o coverm_"%j".out

BAMDIR=$1;shift
DREPDIR=$1;shift
OUT=$1;shift
CORES=$1;shift

echo $CORES

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}


echo "Running dRep"

coverm genome \
  --bam-files $BAMDIR/*.sorted.bam \
  --genome-fasta-directory $DREPDIR/drep_output/dereplicated_genomes/ \
  -x fa \
  -m count \
  --min-covered-fraction 0 \
  -t $CORES \
  -o abundance.tsv

#  -m count relative_abundance trimmed_mean covered_fraction \
#  --min-covered-fraction 10 \
#  -t $CORES \
#  -o abundance.tsv

mkdir -p $OUT

cp -r * $OUT/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r	
