#!/bin/bash -l
#SBATCH -J kaiju_table
#sBATCH --mem=80G
#SBATCH -o kaiju_correct_counts_"%j".out

INFILE=$1;shift
OUTSUFFIX=$1;shift
OUT=$1;shift

SCRIPT_DIR=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

awk -F"\t" '{print gsub(/,/,",",$5) "\t" $5}' < $INFILE > S.new_counts

echo "kaiju file successfully read"

awk -F"\t" '{if($1>0){print $2}}' S.new_counts|$SCRIPT_DIR/kaiju_taxon_to_hash.pl > S.pcounts 

echo "P hash created"

 
awk -F"\t" '{if($1==1){print $2}}' S.new_counts|$SCRIPT_DIR/kaiju_taxon_to_hash.pl > S.ncounts

echo "N hash created"

Rscript $SCRIPT_DIR/prop_reads.R $OUTSUFFIX S.pcounts S.ncounts S.new_counts 


mkdir -p $OUT

echo "copying corrected counts"
cp corrected.counts $OUT/${OUTSUFFIX}.corrected.counts

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
