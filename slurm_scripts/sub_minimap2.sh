#!/bin/bash -l
#SBATCH -J minimap
#sBATCH --mem=500G
#SBATCH -o minimap_"%j".out
REF=$1;shift
OUT=$1;shift
R1=$1;shift
R2=$1;shift
S=$1;shift
CORES=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}


echo "Running minmap hit for sample $S with index $REF"
minimap2 -ax sr --split-prefix=tmp_${S} \
  -t $CORES \
  $REF \
  $R1 \
  $R2 \
  |samtools view -u - | samtools sort -@ 8 -m 1G -o $S.sorted.bam
  
echo "Indexing output"  
samtools index $S.sorted.bam

echo "making output directory $OUT" 
mkdir -p $OUT

cp $S.*  $OUT/.

cd ..

rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
