#!/bin/bash -l
#SBATCH -J whokaryote
#sBATCH --mem=40000
#SBATCH -o whokaryote_"%j".out

ASSEMBLY=$1;shift
OUT=$1;shift
PRODIGAL=$1;shift
SAMPLE=$1;shift
MINSIZE=$1;shift
CORES=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}


echo "Running whokaryote"
conda activate whokaryote
whokaryote.py --contigs $ASSEMBLY --outdir $SAMPLE --prodigal_file $PRODIGAL  --minsize $MINSIZE --threads $CORES --f

mkdir -p $OUT

cp -r $SAMPLE $OUT/.

echo "Finished whokaryote"

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
