#!/bin/bash -l
#SBATCH -J kaiju_table
#sBATCH --mem=12000
#SBATCH -o kaiju_table__"%j".out

NODES=$1;shift
NAMES=$1;shift
OUTFILE=$1;shift
OUT=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}


kaiju2table -t $NODES -n $NAMES -r species -l superkingdom,phylum,class,order,family,genus,species -o $OUTFILE $@

mkdir -p $OUT

cp $OUTFILE $OUT/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
