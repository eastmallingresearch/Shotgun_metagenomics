#!/bin/bash -l
#SBATCH -J megafilt
#sBATCH --mem=40000
#SBATCH -o megafilt_"%j".out

TRUSEQ=$1; shift
PHIXREF=$1; shift
RIBOKMERS=$1; shift
OUTDIR=$1; shift
FORWARD=$1; shift
REVERSE=$1; shift
RNA=${1:-true}; shift
TRIML=( ktrim=l k=23 mink=11 hdist=1 tpe tbo t=10 )
TRIMR=( ktrim=r k=23 mink=11 hdist=1 tpe tbo t=10 )
PHIX=( k=31 hdist=1 t=10 )
RRNA=( k=31 t=10 )

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

F=$(sed 's/.*\///' <<<$FORWARD)
R=$(sed 's/.*\///' <<<$REVERSE)

# remove forward adapters
~/pipelines/common/bbtools/bbduk.sh threads=10 in1=$FORWARD in2=$REVERSE out1=O1F out2=O1R ref=$TRUSEQ ${TRIML[@]}

# remove reverse adapters
~/pipelines/common/bbtools/bbduk.sh threads=10 in1=O1F in2=O1R out1=O2F out2=O2R ref=$TRUSEQ ${TRIMR[@]}

# remove phix
~/pipelines/common/bbtools/bbduk.sh threads=10 in1=O2F in2=O2R out1=$F.filtered.fq.gz out2=$R.filtered.fq.gz ref=$PHIXREF ${PHIX[@]}

# remove rRNA
if [ "$RNA" = true ]; then
  ~/pipelines/common/bbtools/bbduk.sh threads=10 in1=$F.filtered.fq.gz in2=$R.filtered.fq.gz out1=$F.O3F.fq.gz out2=$R.O3R.fq.gz outm1=$F.rRNA.fq.gz outm2=$R.rRNA.fq.gz ref=$RIBOKMERS ${RRNA[@]} stats=$F.stats.txt
  sleep 8m
  rm $F.filtered.fq.gz
  rm $R.filtered.fq.gz
  mv $F.O3F.fq.gz $F.filtered.fq.gz
  mv $R.O3R.fq.gz $R.filtered.fq.gz
fi

mkdir -p $OUTDIR/stats

cp *.fq.gz $OUTDIR/.
cp $F.stats.txt $OUTDIR/stats/.


cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
