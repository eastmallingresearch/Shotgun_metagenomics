#!/usr/bin/env bash
#SBATCH --job-name=uniref90_sum
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/uniref90_sum_%A_%a.out
#SBATCH --error=logs/uniref90_sum_%A_%a.err

set -euo pipefail

FILE_LIST="$1"
MAPS_DB="$2"
OUTDIR="${3:-uniref90_sums}"

mkdir -p "$OUTDIR" logs

INFILE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")"

if [[ -z "$INFILE" ]]; then
    echo "No input file for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

BASE="$(basename "$INFILE")"
OUTFILE="${OUTDIR}/${BASE}.UniRef90.sum.tsv"

SCRATCH_ROOT="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
WORKDIR="$(mktemp -d "${SCRATCH_ROOT}/uniref90_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_XXXXXX")"

trap 'rm -rf "$WORKDIR"' EXIT

WORKDB="${WORKDIR}/work.sqlite"
FLAT_TSV="${WORKDIR}/flat.tsv"
SQL_OUT="${WORKDIR}/out.tsv"

# Make a simple local symlink so the SQLite ATTACH path is clean.
ln -s "$(realpath "$MAPS_DB")" "${WORKDIR}/maps.db"

# Normalize the flat file to tab-delimited two-column input.
# Assumes: column 1 = UniProtID, column 2 = numeric value, no header.
awk 'BEGIN { OFS="\t" } NF >= 2 { print $1, $2 }' "$INFILE" > "$FLAT_TSV"

# If your files have a header, use this instead:
# awk 'BEGIN { OFS="\t" } NR > 1 && NF >= 2 { print $1, $2 }' "$INFILE" > "$FLAT_TSV"

sqlite3 "$WORKDB" <<SQL
.bail on

CREATE TABLE flat (
    UniProtID TEXT NOT NULL,
    value REAL NOT NULL
);

.mode tabs
.import ${FLAT_TSV} flat

CREATE INDEX flat_UniProtID_idx ON flat(UniProtID);

ATTACH DATABASE 'file:${WORKDIR}/maps.db?mode=ro' AS maps;

.headers on
.mode tabs
.once ${SQL_OUT}

WITH flat_summed AS (
    SELECT
        UniProtID,
        SUM(value) AS value_sum
    FROM flat
    GROUP BY UniProtID
)
SELECT
    m.UniRef90,
    SUM(f.value_sum) AS value_sum
FROM flat_summed AS f
JOIN maps.idmapping_updated AS m
    ON m.UniProtID = f.UniProtID
WHERE m.UniRef90 IS NOT NULL
  AND m.UniRef90 <> ''
GROUP BY m.UniRef90
ORDER BY m.UniRef90;
SQL

mv "$SQL_OUT" "$OUTFILE"

echo "Wrote $OUTFILE"
