#!/bin/bash -l
#SBATCH -J bam_count
#sBATCH --mem=40000
#SBATCH -o bam_sort_"%j".out

SCRIPT_DIR=$1;shift
BAM=$1;shift
GFF=$1;shift
OUT=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

F=$(sed 's/.*\///' <<<$BAM|sed 's/\..*//')

echo $F

samtools view $BAM|$SCRIPT_DIR/bam_scaffold_count.pl $GFF $@ > $F.cov

mkdir -p $OUT

cp * $OUT/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
