#!/bin/bash -l
#SBATCH -J humann3
#sBATCH --mem=60000
#SBATCH -o humann_"%j".out

INFILE=($(sed -n -e "$SLURM_ARRAY_TASK_ID p" $1));shift

FORWARD="${INFILE[0]}"
REVERSE="${INFILE[1]}"
ID="${INFILE[2]}"
OUTFOLDER="${INFILE[3]}"

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

echo "concatinating sequence files"
zcat -f -- $FORWARD $REVERSE > $ID.fq

# sleep 10s

conda activate humann3

echo "Running humann pipeline this will take a while..."
humann -i $ID.fq -o . $@

mkdir -p $OUTFOLDER

rm $ID.fq 

echo "Copying data"
cp  -r *  $OUTFOLDER/. 

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
