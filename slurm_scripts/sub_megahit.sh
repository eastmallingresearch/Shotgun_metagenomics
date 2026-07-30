#!/bin/bash -l
#SBATCH -J megahit
#sBATCH --mem=500G
#SBATCH -o megahit_"%j".out

R1=$1;shift
R2=$1;shift
CORES=$1;shift
OUT=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# micromamba run -n binchicken
echo "Running megahit hit with files $R1"
megahit \
 -1 $R1 \
 -2 $R2 \
 --presets meta-large \
 -m 0.85 \
 -t $CORES \
 -o moobles

echo "making output directory $OUT" 
mkdir -p $OUT

cp moobles -r  $OUT/.

cd ..

rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
