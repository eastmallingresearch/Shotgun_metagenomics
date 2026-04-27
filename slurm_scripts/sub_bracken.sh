#!/bin/bash -l
#SBATCH -J bracken
#sBATCH --mem=20000
#SBATCH -o braken_"%j".out

KRAKEN_DB=$1;shift
KRAKEN_REPORT=$1;shift
READ_LENGTH=$1;shift
OUTPUT_PREFIX=$1;shift
OUTDIR=$1;shift

# Variables
RANKS=("D" "P" "C" "O" "F" "G" "S") # Desired ranks: Domain, Phylum, etc.


# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}


# Loop through taxonomic ranks
for RANK in "${RANKS[@]}"; do
    bracken -d $KRAKEN_DB \
            -i $KRAKEN_REPORT \
            -o ${OUTPUT_PREFIX}_${RANK}.txt \
            -r $READ_LENGTH \
            -l $RANK
    echo "Bracken completed for rank $RANK"
done

mkdir -p $OUTDIR

cp *.* $OUTDIR/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
