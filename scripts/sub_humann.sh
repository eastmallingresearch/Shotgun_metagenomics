#!/bin/bash -l
#SBATCH -J kaiju
#sBATCH --mem=60000
#SBATCH -o humann_"%j".out


FORWARD=$1;shift
REVERSE=$1;shift
ID=$1;shift
OUTFOLDER=$1;shift
SCRIPTDIR=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

echo "concatinating sequence files"
zcat -f -- $FORWARD $REVERSE > $ID.fastq

conda activate humann3

echo "Running humann pipeline"
humann -i $ID.fastq -o . $@

mkdir -p $OUTFOLDER

rm $ID.fastq 

echo "Copying data"
cp  -r *  $OUTFOLDER/. 

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
