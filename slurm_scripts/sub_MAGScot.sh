#!/bin/bash -l
#SBATCH -J MAGScot
#sBATCH --mem=40000
#SBATCH -o MAGScot_"%j".out

ASSEMBLY=$1;shift
MAGDIR=$1;shift
TSVIN=$1;shift
OUT=$1;shift
SAMPLE=$1;shift
CORES=$1;shift
PRODIGAL="${1:-0}";shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

if [ "$PRODIGAL" = "0" ]; then
  echo "Running prodigal with $CORES workers"
  # make prodigal working directory
  mkdir tmp_workfolder
  cat $ASSEMBLY|parallel -j $CORES --block 999k --recstart '>' --pipe \
  prodigal -p meta -a tmp_workfolder/sample.{#}.faa -d tmp_workfolder/sample.{#}.ffn -o tmpfile
  
  cat tmp_workfolder/sample.*.faa > sample.prodigal.faa
  cat tmp_workfolder/sample.*.ffn > sample.prodigal.ffn
  rm -r tmp_workfolder tmpfile  
 
else
  echo "Copying prodigal files"
  cp $PRODIGAL sample.prodigal.faa
fi

echo "Running hmmer"

# HMM search against the GTDB r207 marker sets shipped in the MAGScoT repo's hmm/ folder
hmmsearch -o sample.hmm.tigr.out --tblout sample.hmm.tigr.hit.out --noali --notextw --cut_nc --cpu $CORES \
  $MAGDIR/hmm/gtdbtk_rel207_tigrfam.hmm sample.prodigal.faa
hmmsearch -o sample.hmm.pfam.out --tblout sample.hmm.pfam.hit.out --noali --notextw --cut_nc --cpu $CORES \
  $MAGDIR/hmm/gtdbtk_rel207_Pfam-A.hmm sample.prodigal.faa
  

grep -v "^#" sample.hmm.tigr.hit.out | awk '{print $1"\t"$3"\t"$5}' > sample.tigr
grep -v "^#" sample.hmm.pfam.hit.out | awk '{print $1"\t"$4"\t"$5}' > sample.pfam
cat sample.pfam sample.tigr > sample.hmm

# debug copying hmm files  
#mkdir -p $OUT/MAGScot
#cp sample.hmm  $OUT/MAGScot/$SAMPLE.hmm


echo "Running MAGScot"

# make MAGScot output directory
mkdir magout

Rscript $MAGDIR/MAGScoT.R \
  -i $TSVIN \
  --hmm sample.hmm \
  -o magout/magscot \
  -p default \
  --min_markers 25 \
  --min_sharing 0.8
  
echo "Creating new bins"
$MAGDIR/magscot_bin_extractor.py magout/magscot.refined.contig_to_bin.out $ASSEMBLY bins  

mkdir -p $OUT/$SAMPLE

cp -r * $OUT/$SAMPLE/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
