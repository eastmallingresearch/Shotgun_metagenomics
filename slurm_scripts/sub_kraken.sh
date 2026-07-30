#!/bin/bash -l
#SBATCH -J kraken
#SBATCH -o kraken_"%j".out

DBPATH=$1;shift
FORWARD=$1;shift
REVERSE=$1;shift
OUTPUT=$1;shift
REPORT=$1;shift
OUTDIR=$1;shift
CORES=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}


kraken2 --db $DBPATH \
        --threads $CORES \
        --paired \
		--gzip-compressed \
        --output $OUTPUT \
        --report $REPORT \
        $FORWARD $REVERSE
#--classified-out classified_reads#.fq \
#--unclassified-out unclassified_reads#.fq \

mkdir -p $OUTDIR

cp $OUTPUT $OUTDIR/.
cp $REPORT $OUTDIR/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
