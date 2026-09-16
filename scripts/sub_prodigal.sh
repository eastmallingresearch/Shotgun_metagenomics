#!/bin/bash -l
#SBATCH -J prodigal
#sBATCH --mem=40000
#SBATCH -o prodigal_"%j".out

ASSEMBLY=$1;shift
OUT=$1;shift
SAMPLE=$1;shift
CORES=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# make prodigal working directory
mkdir tmp_workfolder

echo "Running prodigal with $CORES workers"

cat $ASSEMBLY|parallel -j $CORES --block 999k --recstart '>' --pipe \
  prodigal -p meta -a tmp_workfolder/sample.{#}.faa -d tmp_workfolder/sample.{#}.ffn -f gff -o tmp_workfolder/sample.{#}.gff
  
cat tmp_workfolder/sample.*.faa > $SAMPLE.prodigal.faa
cat tmp_workfolder/sample.*.ffn > $SAMPLE.prodigal.ffn
cat tmp_workfolder/sample.*.gff > $SAMPLE.prodigal.gff

rm -r tmp_workfolder  

mkdir -p $OUT

cp $SAMPLE.* $OUT/.

echo "completed successsfully"

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
