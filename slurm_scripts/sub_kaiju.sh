#!/bin/bash -l
#SBATCH -J kaiju
#sBATCH --mem=200000
#SBATCH -o kaiju_"%j".out

NODES=$1;shift
NAMES=$1;shift
DB=$1;shift
#BINS=$1;shift
OUTFILE=$1;shift
OUT=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}


kaiju -t $NODES -f $DB -o k.out $@ #-z 20 -v -i $forward -j $reverse 

kaiju-addTaxonNames -t $NODES -n $NAMES -r superkingdom,phylum,class,order,family,genus,species -i k.out -o $OUTFILE

mkdir -p $OUT

cp $OUTFILE $OUT/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
