#!/usr/bin/env bash
# creates a sqllite database for taxonomy lookup.
# requires a copy of names.dmp and nodes.dmp from NCBI

set -euo pipefail

DB="ncbi_taxonomy.sqlite"

rm -f "$DB" names_sci.tsv nodes.tsv

# Scientific names from names.dmp
awk '
BEGIN {
  FS = "[[:space:]]*\\|[[:space:]]*"
  OFS = "\t"
}
function clean(x) {
  gsub(/^[ \t]+|[ \t|]+$/, "", x)
  return x
}
{
  tax_id = clean($1)
  name = clean($2)
  name_class = clean($4)

  if (name_class == "scientific name") {
    print tax_id, name
  }
}
' names.dmp > names_sci.tsv

# Taxonomy tree from nodes.dmp
awk '
BEGIN {
  FS = "[[:space:]]*\\|[[:space:]]*"
  OFS = "\t"
}
function clean(x) {
  gsub(/^[ \t]+|[ \t|]+$/, "", x)
  return x
}
{
  print clean($1), clean($2), clean($3)
}
' nodes.dmp > nodes.tsv

sqlite3 "$DB" <<'SQL'
PRAGMA journal_mode = OFF;
PRAGMA synchronous = OFF;
PRAGMA temp_store = MEMORY;

CREATE TABLE names (
  tax_id INTEGER PRIMARY KEY,
  name TEXT
);

CREATE TABLE nodes (
  tax_id INTEGER PRIMARY KEY,
  parent_tax_id INTEGER,
  rank TEXT
);

.mode tabs
.import names_sci.tsv names
.import nodes.tsv nodes

CREATE INDEX idx_nodes_parent ON nodes(parent_tax_id);
CREATE INDEX idx_nodes_rank ON nodes(rank);
CREATE INDEX idx_names_name ON names(name);

ANALYZE;
SQL

rm -f names_sci.tsv nodes.tsv

echo "Built reusable taxonomy database: $DB"
