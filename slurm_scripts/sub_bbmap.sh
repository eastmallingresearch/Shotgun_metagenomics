#!/bin/bash -l
#SBATCH -J bbmap
#sBATCH --mem=40000
#SBATCH -o bbmap"%j".out

#SCRIPT_DIR=$1; shift
REF=$1; shift
OUTDIR=$1; shift
FORWARD=$1; shift
REVERSE=${1:-NOTHING}; shift
#MEMORY=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}
F=$(sed 's/.*\///' <<<$FORWARD)
R=$(sed 's/.*\///' <<<$REVERSE)

~/pipelines/common/bbtools/bbmap.sh \
 threads=20 \
 in1=$FORWARD \
 in2=$REVERSE \
 outu1=$F.cleaned.fq.gz \
 outu2=$R.cleaned.fq.gz \
 outm1=$F.unclean.fq.gz \
 outm2=$R.uncleaned.fq.gz \
 path=$REF/ \
 $@

mkdir -p $OUTDIR/unclean

cp *.cleaned.fq.gz $OUTDIR/.
cp *.unclean.fq.gz $OUTDIR/unclean/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r


#sbatch --mem=40000 -p medium -c 20 /home/gdeakin/projects/BlueC_coir/metagenomics_pipeline/scripts/slurm/sub_bbmap.sh /home/gdeakin/projects/BlueC_coir/metagenomics_pipeline/common/resources/contaminants/bbmap_human_strawberry_coconut /home/gdeakin/projects/BlueC_coir/data/metagenomics/cleaned /home/gdeakin/projects/BlueC_coir/data/metagenomics/filtered/K_V6_EKDN230032044-1A_HJJTFDSX7_L4_1.fq.gz.filtered.fq.gz /home/gdeakin/projects/BlueC_coir/data/metagenomics/filtered/K_V6_EKDN230032044-1A_HJJTFDSX7_L4_2.fq.gz.filtered.fq.gz 50g minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 t=8 

#-Xmx 