#!/bin/bash -l
#SBATCH -J bbmerge
#sBATCH --mem=85000
#SBATCH -o bbmerge"%j".out

#SCRIPT_DIR=$1; shift
OUTDIR=$1; shift
FORWARD=$1; shift
REVERSE=$1; shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

F=$(sed 's/.*\///' <<<$FORWARD)
R=$(sed 's/.*\///' <<<$REVERSE)

~/pipelines/common/bbtools/bbmerge.sh \
 -Xmx80g \
 in1=$FORWARD \
 in2=$REVERSE \
 out=$F.merged.fq.gz \
 outu1=$F.unmerged.fq.gz \
 outu2=$R.unmerged.fq.gz \
 $@

mkdir -p $OUTDIR

cp *.fq.gz $OUTDIR/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
