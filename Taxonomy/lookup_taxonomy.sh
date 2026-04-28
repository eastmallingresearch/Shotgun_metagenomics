#!/usr/bin/env bash
set -euo pipefail

DB="${1:-ncbi_taxonomy.sqlite}"
TAXIDS="${2:-taxids.txt}"
OUT="${3:-taxonomy_results.tsv}"

tmp_taxids="$(mktemp)"
awk 'NF && $1 !~ /^#/ { print $1 }' "$TAXIDS" > "$tmp_taxids"

sqlite3 "$DB" <<SQL > "$OUT"
.headers on
.mode tabs

DROP TABLE IF EXISTS temp_queries;
CREATE TEMP TABLE temp_queries (
  tax_id INTEGER
);

.import $tmp_taxids temp_queries

WITH RECURSIVE lineage(query_tax_id, tax_id, parent_tax_id, rank, depth) AS (
  SELECT
    q.tax_id,
    n.tax_id,
    n.parent_tax_id,
    n.rank,
    0
  FROM temp_queries q
  LEFT JOIN nodes n ON n.tax_id = q.tax_id

  UNION ALL

  SELECT
    l.query_tax_id,
    n.tax_id,
    n.parent_tax_id,
    n.rank,
    l.depth + 1
  FROM lineage l
  JOIN nodes n ON n.tax_id = l.parent_tax_id
  WHERE l.parent_tax_id != l.tax_id
    AND l.depth < 200
),

named_lineage AS (
  SELECT
    l.query_tax_id,
    l.tax_id,
    l.parent_tax_id,
    l.rank,
    l.depth,
    COALESCE(nm.name, 'taxid:' || l.tax_id) AS name
  FROM lineage l
  LEFT JOIN names nm ON nm.tax_id = l.tax_id
)

SELECT
  q.tax_id AS input_taxid,

  COALESCE((SELECT name FROM names WHERE tax_id = q.tax_id), 'NA') AS taxon_name,
  COALESCE((SELECT rank FROM nodes WHERE tax_id = q.tax_id), 'NA') AS taxon_rank,

  COALESCE((SELECT name FROM named_lineage WHERE query_tax_id = q.tax_id AND rank = 'superkingdom' ORDER BY depth LIMIT 1), 'NA') AS superkingdom,
  COALESCE((SELECT name FROM named_lineage WHERE query_tax_id = q.tax_id AND rank = 'kingdom'      ORDER BY depth LIMIT 1), 'NA') AS kingdom,
  COALESCE((SELECT name FROM named_lineage WHERE query_tax_id = q.tax_id AND rank = 'phylum'       ORDER BY depth LIMIT 1), 'NA') AS phylum,
  COALESCE((SELECT name FROM named_lineage WHERE query_tax_id = q.tax_id AND rank = 'class'        ORDER BY depth LIMIT 1), 'NA') AS class,
  COALESCE((SELECT name FROM named_lineage WHERE query_tax_id = q.tax_id AND rank = 'order'        ORDER BY depth LIMIT 1), 'NA') AS tax_order,
  COALESCE((SELECT name FROM named_lineage WHERE query_tax_id = q.tax_id AND rank = 'family'       ORDER BY depth LIMIT 1), 'NA') AS family,
  COALESCE((SELECT name FROM named_lineage WHERE query_tax_id = q.tax_id AND rank = 'genus'        ORDER BY depth LIMIT 1), 'NA') AS genus,
  COALESCE((SELECT name FROM named_lineage WHERE query_tax_id = q.tax_id AND rank = 'species'      ORDER BY depth LIMIT 1), 'NA') AS species,

  COALESCE(
    (SELECT group_concat(part, '; ')
     FROM (
       SELECT rank || ':' || name AS part
       FROM named_lineage
       WHERE query_tax_id = q.tax_id
       ORDER BY depth DESC
     )),
    'NA'
  ) AS full_taxonomy
  
FROM temp_queries q;
SQL

rm -f "$tmp_taxids"

echo "Wrote: $OUT"
