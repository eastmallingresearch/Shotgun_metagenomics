#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# download_metagenome_functional_dbs.sh
#
# Bulk downloads for a soil/metagenome functional annotation/enrichment database.
#
# Suggested use:
#
#   chmod +x download_metagenome_functional_dbs.sh
#   ./download_metagenome_functional_dbs.sh /path/to/annotation_dbs
#
# Optional switches:
#
#   RUN_HUGE=1 ./download_metagenome_functional_dbs.sh db
#       Also download very large sequence/match files.
#
#   RUN_KEGG=1 ./download_metagenome_functional_dbs.sh db
#       Download a small set of KEGG REST mapping dumps.
#       KEGG REST is academic-use only and should be rate-limited.
#
#   RUN_RECURSIVE=1 ./download_metagenome_functional_dbs.sh db
#       Recursively mirror selected FTP directories, e.g. AMRFinderPlus latest.
#
#   RUN_CARD=1 ./download_metagenome_functional_dbs.sh db
#       Try CARD direct download endpoint.
#
#   RUN_TCDB_PAGE=1 ./download_metagenome_functional_dbs.sh db
#       Save TCDB download page for manual/dynamic links.
#
###############################################################################

DBROOT="${1:-metagenome_annotation_dbs}"

RUN_HUGE="${RUN_HUGE:-0}"
RUN_KEGG="${RUN_KEGG:-0}"
RUN_RECURSIVE="${RUN_RECURSIVE:-0}"
RUN_CARD="${RUN_CARD:-0}"
RUN_TCDB_PAGE="${RUN_TCDB_PAGE:-1}"

# Conservative wget settings:
# -c resumes partial downloads.
# --tries/--waitretry tolerate flaky FTP/HTTP mirrors.
# --no-verbose keeps logs readable.
WGET="wget -c --tries=10 --timeout=120 --waitretry=20 --retry-connrefused --no-verbose"

mkdir -p "$DBROOT"

log() {
    printf '\n[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}

get() {
    # get URL OUTDIR
    url="$1"
    outdir="$2"
    mkdir -p "$outdir"
    log "Downloading: $url"
    (
        cd "$outdir"
        $WGET "$url"
    )
}

try_get() {
    # try_get URL OUTDIR
    # Non-fatal download helper for URLs that occasionally move or are generated.
    url="$1"
    outdir="$2"
    mkdir -p "$outdir"
    log "Trying: $url"
    (
        cd "$outdir"
        $WGET "$url"
    ) || {
        log "WARNING: failed or unavailable: $url"
        return 0
    }
}

get_as() {
    # get_as URL OUTFILE
    url="$1"
    outfile="$2"
    mkdir -p "$(dirname "$outfile")"
    log "Downloading: $url -> $outfile"
    $WGET -O "$outfile" "$url"
}

###############################################################################
# 1. UniProt / UniRef / UniParc backbone
###############################################################################

log "UniProt ID mapping backbone"

UNIPROT_IDMAP_BASE="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping"

get "$UNIPROT_IDMAP_BASE/README"                   "$DBROOT/uniprot/idmapping"
get "$UNIPROT_IDMAP_BASE/idmapping.dat.gz"         "$DBROOT/uniprot/idmapping"
get "$UNIPROT_IDMAP_BASE/idmapping_selected.tab.gz" "$DBROOT/uniprot/idmapping"

# Optional full UniProtKB records.
# Useful if you later want to parse DR lines, EC/Rhea annotations, reviewed status,
# keywords, comments, etc. TrEMBL is very large.
if [ "$RUN_HUGE" = "1" ]; then
    log "UniProtKB full flat files; Swiss-Prot is manageable, TrEMBL is huge"
    UNIPROT_KB_BASE="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"
    get "$UNIPROT_KB_BASE/uniprot_sprot.dat.gz"  "$DBROOT/uniprot/knowledgebase"
    get "$UNIPROT_KB_BASE/uniprot_trembl.dat.gz" "$DBROOT/uniprot/knowledgebase"
fi

# Optional UniRef90 XML if you need cluster membership beyond the mapping table.
if [ "$RUN_HUGE" = "1" ]; then
    log "UniRef90 XML; huge, but contains cluster membership"
    UNIREF90_BASE="https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90"
    get "$UNIREF90_BASE/uniref90.xml.gz" "$DBROOT/uniprot/uniref90"
fi

###############################################################################
# 2. Gene Ontology: ontology + UniProt-GOA + external2go maps
###############################################################################

log "Gene Ontology and GOA"

# GO ontology for enrichment parent propagation / slims / DAG-aware analysis.
get "https://current.geneontology.org/ontology/go-basic.obo" "$DBROOT/go/ontology"
get "https://current.geneontology.org/ontology/go.json"      "$DBROOT/go/ontology"

# Full upstream UniProt-GOA GAF.
# Very useful if you want a cleaner long-form UniProtKB -> GO table than the
# semicolon-packed GO column in idmapping_selected.tab.gz.
get "https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz" "$DBROOT/go/goa"

# GOA external2go mappings. These are small and useful for crosswalking
# InterPro, EC, Rhea, UniProt keywords/locations, etc. to GO.
GO_EXT_BASE="https://ftp.ebi.ac.uk/pub/databases/GO/goa/external2go"

for f in \
    interpro2go \
    hamap2go \
    uniprotkb_kw2go \
    uniprotkb_sl2go \
    pfam2go
do
    try_get "$GO_EXT_BASE/$f" "$DBROOT/go/external2go"
done

###############################################################################
# 3. InterPro: integrated protein family/domain layer
###############################################################################

log "InterPro protein mappings and metadata"

INTERPRO_BASE="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release"

get "$INTERPRO_BASE/protein2ipr.dat.gz"        "$DBROOT/interpro"
get "$INTERPRO_BASE/entry.list"                "$DBROOT/interpro"
get "$INTERPRO_BASE/interpro2go"               "$DBROOT/interpro"
get "$INTERPRO_BASE/ParentChildTreeFile.txt"   "$DBROOT/interpro"
get "$INTERPRO_BASE/interpro.xml.gz"           "$DBROOT/interpro"

# match_complete has region-level match data and can be large.
if [ "$RUN_HUGE" = "1" ]; then
    get "$INTERPRO_BASE/match_complete.dat.gz" "$DBROOT/interpro"
fi

###############################################################################
# 4. Pfam: native HMM/domain resources
###############################################################################

log "Pfam HMMs and metadata"

PFAM_BASE="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"

get "$PFAM_BASE/Pfam-A.hmm.gz"        "$DBROOT/pfam"
get "$PFAM_BASE/Pfam-A.hmm.dat.gz"    "$DBROOT/pfam"
get "$PFAM_BASE/Pfam-A.clans.tsv.gz"  "$DBROOT/pfam"
get "$PFAM_BASE/Pfam-A.regions.tsv.gz" "$DBROOT/pfam"
get "$PFAM_BASE/relnotes.txt"         "$DBROOT/pfam"
get "$PFAM_BASE/md5_checksums"        "$DBROOT/pfam"

# Optional huge Pfam files.
if [ "$RUN_HUGE" = "1" ]; then
    get "$PFAM_BASE/Pfam-A.full.gz"  "$DBROOT/pfam"
    get "$PFAM_BASE/Pfam-A.fasta.gz" "$DBROOT/pfam"
    get "$PFAM_BASE/pfamseq.gz"      "$DBROOT/pfam"
fi

###############################################################################
# 5. Rhea + ENZYME: reaction and EC layer
###############################################################################

log "Rhea reaction mappings and ENZYME EC metadata"

RHEA_BASE="https://ftp.expasy.org/databases/rhea/tsv"

# rhea-tsv.tar.gz contains most non-UniProt TSV mappings.
get "$RHEA_BASE/rhea-tsv.tar.gz" "$DBROOT/rhea"

# UniProt mappings are separate and important for UniProtKB -> reaction joins.
get "$RHEA_BASE/rhea2uniprot_sprot.tsv"     "$DBROOT/rhea"
get "$RHEA_BASE/rhea2uniprot_trembl.tsv.gz" "$DBROOT/rhea"

# A few key TSVs explicitly, useful even if you do not unpack the tarball yet.
for f in \
    rhea2ec.tsv \
    rhea2go.tsv \
    rhea2kegg_reaction.tsv \
    rhea2metacyc.tsv \
    rhea2reactome.tsv \
    rhea-directions.tsv \
    rhea-relationships.tsv
do
    try_get "$RHEA_BASE/$f" "$DBROOT/rhea"
done

# EC number metadata.
get "https://ftp.expasy.org/databases/enzyme/enzyme.dat" "$DBROOT/enzyme"
try_get "https://ftp.expasy.org/databases/enzyme/enzclass.txt" "$DBROOT/enzyme"

# Optional ChEBI ontology/metadata; useful for interpreting Rhea participants.
CHEBI_BASE="https://ftp.ebi.ac.uk/pub/databases/chebi"
try_get "$CHEBI_BASE/ontology/chebi.obo" "$DBROOT/chebi"
try_get "$CHEBI_BASE/flat_files/compounds.tsv.gz" "$DBROOT/chebi"
try_get "$CHEBI_BASE/flat_files/database_accession.tsv.gz" "$DBROOT/chebi"
try_get "$CHEBI_BASE/flat_files/relation.tsv.gz" "$DBROOT/chebi"

###############################################################################
# 6. eggNOG: orthology + predicted functional summaries
###############################################################################

log "eggNOG v7 files"

# The eggNOG v7 web download page exposes these file names. Direct file URLs
# have varied between releases, so this section first tries common direct paths.
# If they fail, it saves the downloads page for manual link extraction.
EGGNOG_DIR="$DBROOT/eggnog"
mkdir -p "$EGGNOG_DIR"

#https://eggnogdb.org/public/eggnog7/e7.og_info_kegg_go.tsv.gz

# Try the likely direct URLs.
for f in \
    e7.og_info_kegg_go.tsv.gz \
    e7.protein_families.tsv.gz \
    e7.clu2ogs.tsv.gz \
    e7.taxid_info.tsv.gz
do
    try_get "https://eggnogdb.org/downloads/$f" "$EGGNOG_DIR"
done

# Save the download page as a fallback for extracting actual hrefs.
try_get "https://eggnogdb.org/downloads/" "$EGGNOG_DIR"

# Optional huge eggNOG sequence/tree files.
if [ "$RUN_HUGE" = "1" ]; then
    for f in \
        e7.proteins.fa.gz \
        e7.og_fasta_sequences.tar \
        e7.trees.tsv.gz
    do
        try_get "https://eggnogdb.org/downloads/$f" "$EGGNOG_DIR"
    done
fi

###############################################################################
# 7. OrthoDB: orthologous groups and cross-references
###############################################################################

log "OrthoDB v12 files"

ORTHODB_BASE="https://data.orthodb.org/v12/download/odb_data_dump"

for f in \
    README.txt \
    odb12v2_species.tab.gz \
    odb12v2_levels.tab.gz \
    odb12v2_level2species.tab.gz \
    odb12v2_OGs.tab.gz \
    odb12v2_OG_xrefs.tab.gz \
    odb12v2_OG_pairs.tab.gz
do
    get "$ORTHODB_BASE/$f" "$DBROOT/orthodb/v12"
done

# Very large OrthoDB protein/gene mapping files.
if [ "$RUN_HUGE" = "1" ]; then
    for f in \
        odb12v2_genes.tab.gz \
        odb12v2_gene_xrefs.tab.gz \
        odb12v2_OG2genes.tab.gz \
        odb12v2_aa_fasta.gz \
        odb12v2_og_aa_fasta.gz
    do
        get "$ORTHODB_BASE/$f" "$DBROOT/orthodb/v12"
    done
fi

###############################################################################
# 8. dbCAN / CAZy-like CAZyme annotation resources
###############################################################################

log "dbCAN CAZyme database files"

DBCAN_BASE="https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total"

# Core files described by run_dbCAN docs.
try_get "$DBCAN_BASE/CAZy.dmnd"      "$DBROOT/dbcan"
try_get "$DBCAN_BASE/dbCAN.hmm"      "$DBROOT/dbcan"
try_get "$DBCAN_BASE/dbCAN_sub.hmm"  "$DBROOT/dbcan"

# The official run_dbCAN command is often the most reliable way to get the full
# moving database set:
#
#   run_dbcan database --db_dir "$DBROOT/dbcan" --aws_s3
#
# This script keeps to wget, but the command above is preferred once run_dbCAN is
# installed.

###############################################################################
# 9. TCDB: transporter classification
###############################################################################

log "TCDB"

# TCDB uses web download links/buttons for several mapping files. Save the page;
# use the links on that page for current files. This avoids hard-coding dynamic
# script URLs that change.
if [ "$RUN_TCDB_PAGE" = "1" ]; then
    try_get "https://saierlab.biosci.ucsd.edu/download/" "$DBROOT/tcdb"
    try_get "https://tcdb.org/" "$DBROOT/tcdb"
fi

###############################################################################
# 10. MEROPS: peptidases/proteases
###############################################################################

log "MEROPS peptidase/protease files"
MEROPS_BASE="ftp.ebi.ac.uk/pub/databases/merops/current_release"

# These filenames are listed on the MEROPS download page. Some servers redirect;
# try_get makes this non-fatal.
for f in \
    dnld_list.txt \
    pepunit.lib \
    protease.lib \
    merops_scan.lib \
    meropsrefs.txt
do
    try_get "$MEROPS_BASE/$f" "$DBROOT/merops"
done

# Save download page too, because release tarball names change, e.g.
# meropsweb124.tar.gz in release 12.4.
try_get "https://www.ebi.ac.uk/merops/download_list.shtml" "$DBROOT/merops"

###############################################################################
# 11. NCBI PGAP HMMs / NCBIfam / TIGRFAM / PRK-style models
###############################################################################

log "NCBI PGAP HMMs / protein family models"

# Current HMM directory can contain release tarballs, seed alignments, metadata,
# and subdirectories. Recursive mirror is optional to avoid accidental huge pulls.
if [ "$RUN_RECURSIVE" = "1" ]; then
    mkdir -p "$DBROOT/ncbi_hmm/current"
    (
        cd "$DBROOT/ncbi_hmm/current"
        wget -c -r -np -nH --cut-dirs=1 -R "index.html*" \
            --tries=10 --timeout=120 --waitretry=20 --retry-connrefused \
            "https://ftp.ncbi.nlm.nih.gov/hmm/current/"
    )
else
    # Common high-level files; non-fatal because filenames can change by release.
    for f in \
        README \
        release.txt \
        hmm_PGAP.HMM.tgz \
        hmm_PGAP.SEED.tgz \
        for-interpro.tsv \
        hmm_PGAP.LIB
    do
        try_get "https://ftp.ncbi.nlm.nih.gov/hmm/current/$f" "$DBROOT/ncbi_hmm/current"
    done
fi

###############################################################################
# 12. AMRFinderPlus database: AMR/stress/virulence specialist layer
###############################################################################

log "NCBI AMRFinderPlus database"

AMRFINDER_BASE="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"

if [ "$RUN_RECURSIVE" = "1" ]; then
    mkdir -p "$DBROOT/amrfinderplus/latest"
    (
        cd "$DBROOT/amrfinderplus/latest"
        wget -c -r -np -nH --cut-dirs=5 -R "index.html*" \
            --tries=10 --timeout=120 --waitretry=20 --retry-connrefused \
            "$AMRFINDER_BASE/"
    )
else
    for f in \
        ReferenceGeneCatalog.txt \
        ReferenceGeneHierarchy.txt \
        AMRProt.fa \
        AMRProt-mutation.fa \
        AMR_CDS.fa \
        AMR_DNA.fa \
        fam.tsv \
        changes.txt \
        version.txt
    do
        try_get "$AMRFINDER_BASE/$f" "$DBROOT/amrfinderplus/latest"
    done
fi

###############################################################################
# 13. CARD: optional AMR layer
###############################################################################

if [ "$RUN_CARD" = "1" ]; then
    log "CARD latest data"
    # CARD's latest/data endpoint usually returns the current data archive.
    # Content-Disposition may provide the versioned filename.
    mkdir -p "$DBROOT/card"
    (
        cd "$DBROOT/card"
        wget -c --content-disposition --tries=10 --timeout=120 \
            "https://card.mcmaster.ca/latest/data" || true
    )
fi

###############################################################################
# 14. CDD: optional NCBI conserved domain metadata
###############################################################################

log "NCBI CDD metadata"

CDD_BASE="https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd"

for f in \
    README \
    cddid.tbl.gz \
    cddid_all.tbl.gz \
    cddannot.dat.gz \
    cddannot_generic.dat.gz \
    cddmasters.fa.gz \
    bitscore_specific.txt \
    family_superfamily_links
do
    try_get "$CDD_BASE/$f" "$DBROOT/ncbi_cdd"
done

# Optional huge CDD model archives.
if [ "$RUN_HUGE" = "1" ]; then
    for f in \
        cdd.tar.gz \
        fasta.tar.gz \
        acd.tar.gz
    do
        try_get "$CDD_BASE/$f" "$DBROOT/ncbi_cdd"
    done
fi

###############################################################################
# 15. KEGG REST mappings: optional and academic-use only
###############################################################################

if [ "$RUN_KEGG" = "1" ]; then
    log "KEGG REST mappings; academic-use only, rate-limited"

    mkdir -p "$DBROOT/kegg"

    # KEGG asks users of rest.kegg.jp to limit requests to <=3 calls/sec.
    # This section is a handful of calls only, with sleeps.
    get_as "https://rest.kegg.jp/list/ko"       "$DBROOT/kegg/list_ko.tsv"
    sleep 1
    get_as "https://rest.kegg.jp/list/pathway"  "$DBROOT/kegg/list_pathway.tsv"
    sleep 1
    get_as "https://rest.kegg.jp/list/module"   "$DBROOT/kegg/list_module.tsv"
    sleep 1
    get_as "https://rest.kegg.jp/list/reaction" "$DBROOT/kegg/list_reaction.tsv"
    sleep 1
    get_as "https://rest.kegg.jp/list/enzyme"   "$DBROOT/kegg/list_enzyme.tsv"
    sleep 1

    get_as "https://rest.kegg.jp/link/pathway/ko" "$DBROOT/kegg/ko_to_pathway.tsv"
    sleep 1
    get_as "https://rest.kegg.jp/link/module/ko"  "$DBROOT/kegg/ko_to_module.tsv"
    sleep 1
    get_as "https://rest.kegg.jp/link/brite/ko"   "$DBROOT/kegg/ko_to_brite.tsv"
    sleep 1
    get_as "https://rest.kegg.jp/link/reaction/ko" "$DBROOT/kegg/ko_to_reaction.tsv"
    sleep 1
    get_as "https://rest.kegg.jp/link/enzyme/ko" "$DBROOT/kegg/ko_to_enzyme.tsv"
    sleep 1
fi

###############################################################################
# 16. BioCyc / MetaCyc note
###############################################################################

cat > "$DBROOT/biocyc_metacyc_NOTE.txt" <<'EOF'
BioCyc / MetaCyc:
  Bulk flat files are generally downloaded after accepting the BioCyc/MetaCyc
  license terms from the BioCyc download site. 

Useful files after download commonly include:
  pathways.dat
  reactions.dat
  enzrxns.dat
  proteins.dat
  genes.dat
  compounds.dat
  classes.dat

Load them as a separate licensed data layer rather than attempting to fetch them
blindly with wget.
EOF

###############################################################################
# 17. Manifest
###############################################################################

log "Writing manifest"

find "$DBROOT" -type f -printf '%p\t%s\n' \
    | sort \
    > "$DBROOT/download_manifest.path_size.tsv"

cat > "$DBROOT/README.downloads.txt" <<EOF
Downloaded on: $(date -Is)

Root:
  $DBROOT

Useful next steps:
  1. Verify large downloads:
       find "$DBROOT" -name "*.gz" -size -10k -print
       find "$DBROOT" -name "*.html*" -print

  2. Build SQLite from:
       uniprot/idmapping/idmapping.dat.gz
       go/goa/goa_uniprot_all.gaf.gz
       interpro/protein2ipr.dat.gz
       rhea/rhea2uniprot_*.tsv*
       eggnog/e7.*.tsv.gz
       orthodb/v12/*.tab.gz

  3. Specialist local annotation databases:
       pfam/Pfam-A.hmm.gz
       dbcan/dbCAN.hmm, dbCAN_sub.hmm, CAZy.dmnd
       ncbi_hmm/current/*
       amrfinderplus/latest/*
       merops/*

  4. KEGG:
       Only downloaded when RUN_KEGG=1.
       Respect KEGG academic-use and rate-limit terms.

  5. BioCyc/MetaCyc:
       See biocyc_metacyc_NOTE.txt.

Manifest:
  download_manifest.path_size.tsv
EOF

log "Done"
