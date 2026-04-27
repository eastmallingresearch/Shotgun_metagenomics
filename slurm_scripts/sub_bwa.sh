#!/bin/bash
#SBATCH --job-name=bwa_human
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --partition=medium
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

REF=$1; shift
OUTDIR=$1; shift
FORWARD=$1; shift
REVERSE=${1:-NOTHING}; shift

THREADS=${SLURM_CPUS_PER_TASK:-1}

# Load modules if needed on your cluster
# module purge || true
# module load bwa-mem2 samtools || true

# Use node-local scratch if available
TMP_WORK="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/bwa_human_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMP_WORK"
cd "$TMP_WORK"

F=$(basename "$FORWARD")
R=$(basename "$REVERSE")

SAMPLE_R1="${F%.fq.gz}"
SAMPLE_R2="${R%.fq.gz}"

echo "Reference: $REF"
echo "Output dir: $OUTDIR"
echo "Forward: $FORWARD"
echo "Reverse: $REVERSE"
echo "Threads: $THREADS"
echo "Working dir: $TMP_WORK"

# Align to human reference
bwa-mem2 mem -t "$THREADS" "$REF" "$FORWARD" "$REVERSE" \
  | samtools view -@ 4 -b -o sample.vs_human.bam -

# Keep pairs where both mates are unmapped
samtools view -@ 4 -u -f 12 -F 2304 sample.vs_human.bam > sample.nonhuman_pairs.bam

# Name-collate and convert back to FASTQ
samtools collate -u -O sample.nonhuman_pairs.bam \
  | samtools fastq -@ 4 \
      -1 "${SAMPLE_R1}.nonhuman.fq.gz" \
      -2 "${SAMPLE_R2}.nonhuman.fq.gz" \
      -0 /dev/null \
      -s /dev/null \
      -n -

mkdir -p "$OUTDIR"/human_mapped
mkdir -p "$OUTDIR"/logs

cp "${SAMPLE_R1}.nonhuman.fq.gz" "$OUTDIR"/
cp "${SAMPLE_R2}.nonhuman.fq.gz" "$OUTDIR"/
cp sample.vs_human.bam "$OUTDIR"/human_mapped/"${SAMPLE_R1}.vs_human.bam"

echo "Done."
