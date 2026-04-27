#!/bin/bash -l
#SBATCH -J kaiju_protein
#sBATCH --mem=80G
#SBATCH -o kaiju_protein_count_"%j".out

INFILE=$1;shift
OUTSUFFIX=$1;shift
OUT=$1;shift

SCRIPT_DIR=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

SCRIPT_DIR/kaiju_out_to_prot_id.pl $INFILE > prot.counts

mkdir -p $OUT

echo "copying protein counts"
cp prot.counts $OUT/${OUTSUFFIX}.prot.counts

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
