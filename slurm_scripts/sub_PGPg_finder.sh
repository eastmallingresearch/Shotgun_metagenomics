#!/bin/bash -l
#SBATCH -J PGPg
#SBATCH -o PGPg_"%j".out

FORWARD=$1;shift
REVERSE=$1;shift
OUTDIR=$1;shift
DB=$1;shift
CORES=$1;shift

# make session temp directory
mkdir $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

# change to session temp folder
cd $TMPDIR/${SLURM_JOB_USER}_${SLURM_JOBID}

mkdir OUTPUT
mkdir INDIR

# diamond can't handle paired end reads, need to merge first

F=$(sed 's/.*\///' <<<$FORWARD)

echo $F

~/pipelines/common/bbtools/bbmerge.sh \
 -Xmx40g \
 in1=$FORWARD \
 in2=$REVERSE \
 out=$F.merged.fa

# DIAMOND 

# Default values
min_identity=30
min_query_cover=30
diamond_extra=""
min_score=""
evalue="1e-5"
diamond_mode=""

# Parse options
#while true; do
#    case "$1" in
#        --piden) min_identity=$2; shift 2 ;;
#        --qcov) min_query_cover=$2; shift 2 ;;
#        --extra) diamond_extra=$2; shift 2 ;;
#        --bitscore) min_score=$2; shift 2 ;;
#        --evalue) evalue=$2; shift 2 ;;
#        --dmode) diamond_mode=$2; shift 2 ;;
#        -h|--help) display_help; exit 0 ;;
#        --) shift; break ;;
#        *) display_help; exit 1 ;;
#    esac
#done


# Check if DIAMOND database exists
#script_dir=$(dirname "$(dirname "$(readlink -f "$0")")")
#diamond_db= "$APPS/PGPg_finder/database/metagenome.dmnd"

if [ ! -f "$DB" ]; then
    echo "Error: DIAMOND database not found at $DB."
    exit 1
else
    echo "DIAMOND database found at $DB"
fi

# Create a file to store gene counts
gene_counts_file=OUTPUT/gene_counts.txt
echo -e "Sample\tID\tCount" > "$gene_counts_file"
echo "Created gene counts file: $gene_counts_file"

# Run DIAMOND
echo "Running DIAMOND for PLaBAse alignment for ${sample}"
diamond blastx -d "$DB" \
  -q $F.merged.fa \
  -o OUTPUT/${F}_diamond.txt \
  -k 1 \
  -p $CORES \
  -e "$evalue" \
  --id "$min_identity" \
  --query-cover "$min_query_cover" \
  $( [ -n "$min_score" ] && echo "--min-score $min_score" ) \
  $diamond_mode \
  $diamond_extra


#python ~/apps/PGPg_finder/PGPg_finder.py -w $PROGRAM \
#        -i INDIR \
#        -o OUTPUT \
#        -t $THREADS

mkdir -p $OUTDIR

cp OUTPUT/* $OUTDIR/.

cd ..
rm ${SLURM_JOB_USER}_${SLURM_JOBID} -r
