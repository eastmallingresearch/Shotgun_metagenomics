# Functional analysis and taxonomy 

The are now quite a few options in this space, Megan +DIAMOND, MetaLAFFA (uses Diamond), SqueezeMeta, Humann3.
  
However none of them seem more suited than using the output from Kaiju (which is orders of magnitude faster than DIAMOND - there's MMseqs2 which might be a good replacement).
  

Kaiju outputs a list of GenBank/Refseq protein IDs. These can be mapped to other systems to unlock various functional annotations


## Mapping files

### Core identifier / GO layer
| Layer                    | Main files / locations                                                                                                                                                                                                      | Notes                                                                                                                                                                                             |
| ------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| UniProt ID spine         | `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz`                                                                                                                    | Your main 3-column backbone: `UniProtKB-AC, ID_type, ID`.                                                                                                                                         |
| UniProt selected mapping | `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz`                                                                                                           | Wide version with `RefSeq`, `EMBL-CDS`, `GO`, `UniRef90`, etc.                                                                                                                                    |
| UniProtKB full records   | `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz``````https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz` | Use only if you want to parse richer `DR`, `CC`, EC, Rhea, keywords, comments, reviewed status.                                                                                                   |
| GO annotations           | `https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz`                                                                                                                                                 | Full UniProt-GOA annotation file; good for UniProtKB accession → GO. GOA says these files include manual and computational annotations and are released roughly every four weeks. ([EMBL-EBI][1]) |
| GO ontology              | `https://current.geneontology.org/ontology/go-basic.obo``````https://current.geneontology.org/ontology/go.json`                                                                                                             | Needed for parent-term propagation, slim mapping, enrichment rollups.                                                                                                                             |
| GO external mappings     | `https://ftp.ebi.ac.uk/pub/databases/GO/goa/external2go/`                                                                                                                                                                   | Includes mappings like `interpro2go`, `ec2go`, `rhea2go`, etc.; GOA describes these as mappings from InterPro, EC, HAMAP, UniProt keywords, locations, etc. to GO. ([EMBL-EBI][2])                |

[1]: https://www.ebi.ac.uk/GOA/?utm_source=chatgpt.com "GOA | European Bioinformatics Institute"
[2]: https://www.ebi.ac.uk/GOA/downloads.html?utm_source=chatgpt.com "Downloads | European Bioinformatics Institute"


### Protein family / domain layer
| Layer                          | Main files / locations                                                                                                                                                                                                                                | Notes                                                                                                                                                                                                                          |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| InterPro protein mappings      | `https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz`                                                                                                                                                                     | UniProt protein → InterPro entry. This is probably your best primary domain/family mapping. InterPro bulk downloads include precomputed InterPro data and prior releases via FTP. ([interpro-documentation.readthedocs.io][1]) |
| InterPro match-level detail    | `https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/match_complete.dat.gz`                                                                                                                                                                  | Use if you need coordinates, member database signatures, domain locations.                                                                                                                                                     |
| InterPro metadata              | `https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list``````https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro.xml.gz``````https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/ParentChildTreeFile.txt` | Entry names, types, hierarchy.                                                                                                                                                                                                 |
| InterPro → GO                  | `https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go`                                                                                                                                                                            | Also mirrored through GOA `external2go`.                                                                                                                                                                                       |
| Pfam native HMMs               | `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz``````https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz`                                                                                            | Use if you want to run HMMER yourself. Pfam’s current release directory contains current flat files, and Pfam mappings are also available through InterPro/Pfam FTP structure. ([pfam-docs.readthedocs.io][2])                 |
| CDD / COG / PRK / NCBI domains | `https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/`                                                                                                                                                                                                          | Optional. CDD includes NCBI-curated domains plus imported Pfam, SMART, COG, PRK and TIGRFAM models. ([NCBI][3])                                                                                                                |

[1]: https://interpro-documentation.readthedocs.io/en/latest/download.html?utm_source=chatgpt.com "How to download InterPro data? — InterPro Documentation"
[2]: https://pfam-docs.readthedocs.io/en/latest/ftp-site.html?utm_source=chatgpt.com "FTP Site — Pfam Documentation"
[3]: https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml?utm_source=chatgpt.com "Conserved Domains Database (CDD) and Resources"

### Enzyme / reaction / pathway layer
| Layer                            | Main files / locations                                                                                                                                                                                  | Notes                                                                                                                                                                                      |
| -------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Rhea reaction mappings           | `https://ftp.expasy.org/databases/rhea/tsv/`                                                                                                                                                            | Key files below. Rhea provides TSV cross-references to EC, UniProtKB, GO, KEGG reaction, MetaCyc, Reactome, etc. ([rhea-db.org][1])                                                        |
| Rhea → UniProt                   | `https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_sprot.tsv``````https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_trembl.tsv.gz`                                                            | UniProtKB accession → precise reaction annotation via Rhea.                                                                                                                                |
| Rhea → EC                        | `https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv`                                                                                                                                                 | Better than EC alone, because Rhea can represent reactions without EC numbers.                                                                                                             |
| Rhea → GO                        | `https://ftp.expasy.org/databases/rhea/tsv/rhea2go.tsv`                                                                                                                                                 | Useful for MF propagation.                                                                                                                                                                 |
| Rhea → KEGG / MetaCyc / Reactome | `https://ftp.expasy.org/databases/rhea/tsv/rhea2kegg_reaction.tsv``````https://ftp.expasy.org/databases/rhea/tsv/rhea2metacyc.tsv``````https://ftp.expasy.org/databases/rhea/tsv/rhea2reactome.tsv`     | Good pathway bridge.                                                                                                                                                                       |
| Rhea relationships               | `https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv``````https://ftp.expasy.org/databases/rhea/tsv/rhea-relationships.tsv`                                                                   | Needed to collapse reaction directions / related reaction IDs.                                                                                                                             |
| ENZYME EC metadata               | `https://ftp.expasy.org/databases/enzyme/enzyme.dat`                                                                                                                                                    | EC number → name, catalytic activity, Swiss-Prot pointers. ExPASy says ENZYME is the EC nomenclature database and provides this single downloadable file. ([enzyme.expasy.org][2])         |
| KEGG                             | `https://rest.kegg.jp/link/pathway/ko``````https://rest.kegg.jp/link/module/ko``````https://rest.kegg.jp/link/brite/ko``````https://rest.kegg.jp/link/reaction/ko``````https://rest.kegg.jp/link/ec/ko` | KEGG has REST endpoints, but the public REST service is academic-use only and rate-limited to 3 calls/sec; for true local bulk use, KEGG FTP requires licensing. ([genome.jp][3])          |
| BioCyc / MetaCyc                 | `https://biocyc.org/download.shtml`                                                                                                                                                                     | Useful files after download: `pathways.dat`, `reactions.dat`, `enzrxns.dat`, `proteins.dat`, `genes.dat`, `pathways.col`, `enzymes.col`. Pathway Tools documents these flat-file formats.  |

[1]: https://www.rhea-db.org/help/download?utm_source=chatgpt.com "Rhea help results"
[2]: https://enzyme.expasy.org/?utm_source=chatgpt.com "Expasy - ENZYME"
[3]: https://www.genome.jp/kegg/rest/keggapi.html?utm_source=chatgpt.com "KEGG API Manual"

### Orthology / metagenome-friendly function layer
| Layer   | Main files / locations                   | Notes                                                                                                                                                                                                                                |
| ------- | ---------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| eggNOG  | `https://eggnogdb.org/downloads/`        | Useful files include `e7.og_info_kegg_go.tsv.gz`, `e7.protein_families.tsv.gz`, and `e7.proteins.fa.gz`. eggNOG’s download page describes OG info files with proteins, KEGG KO terms, symbols and GO slim terms. ([eggnogdb.org][1]) |
| OrthoDB | `https://data.orthodb.org/v12/download/` | Useful files: `odb12v1_OG2genes.tab.gz`, `odb12v1_OG_xrefs.tab.gz`, `odb12v1_gene_xrefs.tab.gz`, `odb12v1_OGs.tab.gz`. The OrthoDB v12 download page lists these bulk mappings.                                                      |
| OMA     | `https://omabrowser.org/oma/current/`    | Optional. OMA bulk downloads include orthology relationships, OMA groups/HOGs, sequences, annotations, and mappings to UniProt/RefSeq/EntrezGene. ([omabrowser.org][2])                                                              |

[1]: https://eggnogdb.org/downloads/?utm_source=chatgpt.com "eggNOG"
[2]: https://omabrowser.org/oma/uses/?utm_source=chatgpt.com "Access the OMA Data"


### Soil-relevant specialist layers
| Layer                                | Main files / locations                                                                                                            | Notes                                                                                                                                                                                                                                                                          |
| ------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| dbCAN / CAZy-style CAZyme annotation | Automatic: `run_dbcan database --db_dir db --aws_s3` Manual base: `https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/` | dbCAN provides `CAZy.dmnd`, `dbCAN.hmm`, `dbCAN_sub.hmm`, substrate mappings, PUL files, and related databases. The docs recommend the `run_dbcan database` command for installation. ([Project name not set][1])                                                              |
| TCDB transporters                    | `https://tcdb.org/` and `https://saierlab.biosci.ucsd.edu/download/`                                                              | TCDB has download links for all protein sequences, annotated substrates, PDB mappings, and RefSeq/UniProt accession to TCID mappings. ([saierlab.biosci.ucsd.edu][2])                                                                                                          |
| MEROPS proteases                     | `https://www.ebi.ac.uk/merops/download_list.shtml`                                                                                | Useful files: `dnld_list.txt`, `pepunit.lib`, `protease.lib`, `merops_scan.lib`, and the SQL dump. MEROPS explicitly provides accession lists and FASTA libraries for peptidases/inhibitors. ([EMBL-EBI][3])                                                                   |
| NCBIfam / TIGRFAM / PGAP HMMs        | `https://ftp.ncbi.nlm.nih.gov/hmm/`                                                                                               | NCBI says current TIGRFAM models are distributed as part of current PGAP HMMs, alongside NCBIFAMs and PRK-derived models. ([NCBI][4])                                                                                                                                          |
| AMRFinderPlus                        | `https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/`                                   | Key files include `ReferenceGeneCatalog.txt`, `ReferenceGeneHierarchy.txt`, `AMRProt.fa`, and internal hierarchy/HMM files. The NCBI AMRFinderPlus database docs recommend `ReferenceGeneCatalog.txt` for database content. ([GitHub Wiki][5])                                 |
| CARD                                 | `https://card.mcmaster.ca/latest/data`                                                                                            | CARD provides reference DNA/protein sequences, detection models, ARO ontology and RGI support. For pipelines, commonly used files include `card.json`, protein/nucleotide FASTA model files, and ARO metadata after extraction. ([McMaster Antibiotic Resistance Database][6]) |

[1]: https://run-dbcan.readthedocs.io/en/latest/getting_started/database_description.html?utm_source=chatgpt.com "Database Description — dbcan"
[2]: https://saierlab.biosci.ucsd.edu/download/?utm_source=chatgpt.com "Download | Saier Lab"
[3]: https://www.ebi.ac.uk/merops/download_list.shtml?utm_source=chatgpt.com "MEROPS - the Peptidase Database"
[4]: https://www.ncbi.nlm.nih.gov/genome/annotation_prok/tigrfams/?utm_source=chatgpt.com "TIGRFAMs at NCBI"
[5]: https://github-wiki-see.page/m/ncbi/amr/wiki/AMRFinderPlus-database?utm_source=chatgpt.com "AMRFinderPlus database - ncbi/amr GitHub Wiki"
[6]: https://card.mcmaster.ca/?utm_source=chatgpt.com "The Comprehensive Antibiotic Resistance Database"



The uniprot mapping file contains the majority of useful accessions

```shell
# UNIPROT selected mapping file
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
# columns: UniProtKB-AC,UniProtKB-ID,GeneID,RefSeq,GI,PDB,GO,UniRef100,UniRef90,UniRef50,UniParc,PIR,NCBI-taxon,MIM,UniGene,PubMed,EMBL,EMBL-CDS,Ensembl,Ensembl_TRS,Ensembl_PRO,Additional PubMed
# Uniprot complete mpping file
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
```
Possible mappings to download:
idmapping.dat.gz
goa_uniprot_all.gaf.gz
go-basic.obo
protein2ipr.dat.gz
interpro2go
rhea2uniprot_sprot.tsv
rhea2uniprot_trembl.tsv.gz
rhea2ec.tsv
rhea2go.tsv
rhea2kegg_reaction.tsv
enzyme.dat
eggNOG OG info / protein family files
dbCAN databases
TCDB mappings
MEROPS accession / FASTA libraries
NCBI PGAP HMMs
AMRFinderPlus database


Useful final annotations:
EMBL-CDS
UniRef90 / UniRef50
GO
eggNOG
KEGG
NCBI_TaxID
RefSeq
UniParc

Sequence / clustering:
  UniRef90
  UniRef50
  UniParc
  UniProtKB accession
  NCBI_TaxID

General function:
  GO
  InterPro
  Pfam
  eggNOG / COG
  KEGG KO

Enzyme / reaction / pathway:
  EC
  Rhea
  KEGG reaction / module / pathway
  BioCyc / MetaCyc
  UniPathway

Soil-relevant specialist function:
  CAZy / dbCAN
  TCDB
  MEROPS
  TIGRFAMs / NCBIfam
  CARD / AMRFinderPlus, optional but useful

Context / interpretation:
  NCBI_TaxID
  source database accession type
  reviewed/unreviewed UniProt status, if you can add it

To extract from id mapping:
UniProtKB-ID
UniParc
UniRef100
UniRef90
UniRef50
RefSeq
RefSeq_NT
EMBL
EMBL-CDS
GO
GeneID
Gene_Name
Gene_OrderedLocusName
Gene_ORFName
Gene_Synonym
NCBI_TaxID
eggNOG
KEGG
BioCyc
UniPathway
Reactome
PlantReactome
TCDB
MEROPS
PHI-base
STRING
OMA
OrthoDB
HOGENOM
GeneTree
PATRIC
VEuPathDB
EnsemblGenome
EnsemblGenome_PRO
EnsemblGenome_TRS

### Mapping files

```shell
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
    ec2go \
    rhea2go \
    metacyc2go \
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
try_get "$CHEBI_BASE/Flat_file_tab_delimited/compounds.tsv.gz" "$DBROOT/chebi"
try_get "$CHEBI_BASE/Flat_file_tab_delimited/database_accession.tsv.gz" "$DBROOT/chebi"
try_get "$CHEBI_BASE/Flat_file_tab_delimited/relation.tsv.gz" "$DBROOT/chebi"

###############################################################################
# 6. eggNOG: orthology + predicted functional summaries
###############################################################################

log "eggNOG v7 files"

# The eggNOG v7 web download page exposes these file names. Direct file URLs
# have varied between releases, so this section first tries common direct paths.
# If they fail, it saves the downloads page for manual link extraction.
EGGNOG_DIR="$DBROOT/eggnog"
mkdir -p "$EGGNOG_DIR"

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

ORTHODB_BASE="https://data.orthodb.org/v12/download"

for f in \
    README.txt \
    odb12v1_species.tab.gz \
    odb12v1_levels.tab.gz \
    odb12v1_level2species.tab.gz \
    odb12v1_OGs.tab.gz \
    odb12v1_OG_xrefs.tab.gz \
    odb12v1_OG_pairs.tab.gz
do
    get "$ORTHODB_BASE/$f" "$DBROOT/orthodb/v12"
done

# Very large OrthoDB protein/gene mapping files.
if [ "$RUN_HUGE" = "1" ]; then
    for f in \
        odb12v1_genes.tab.gz \
        odb12v1_gene_xrefs.tab.gz \
        odb12v1_OG2genes.tab.gz \
        odb12v1_aa_fasta.gz \
        odb12v1_og_aa_fasta.gz
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

MEROPS_BASE="https://www.ebi.ac.uk/merops/download"

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
        hmm_PGAP.HMM.tsv \
        hmm_PGAP.LIB.tgz
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

```

## Kaiju functional analysis

### [Taxonomy](../master/Taxonomy/README.md)

Run the kaiju taxonomy scripts first 

### Functional analysis 

The output data contains a huge number of proteins, almost all which can not be distinguished at the functional level - they may have different functions, but the details are not yet available. Due to this massive duplication, it makes sense to shrink the data down to unique protein functions.  

Below is a short bit of R code to do this. Not all steps are strictly neccessary (the DESeq stuff can all be dropped if it's not going to be used)

```R
## Load required libraries
library(DESeq2)
library(data.table)
library(tidyverse)

# uncomment below 2 lines to install metafuncs package
# library(devtools)
# install_github("eastmallingresearch/Metabarcoding_pipeline/scripts")
library(metafuncs)


# Custom functions
dt_to_df <-
function (DT, row_names = 1) 
{
  DF <- as.data.frame(DT)
  row.names(DF) <- DF[, row_names]
  DF <- DF[, -row_names]
  DF
}

# host <- ifelse(Sys.info()[['sysname']]=="Windows","/","~")

#=============================================================================== -->
#       Load data
#=============================================================================== -->

# Load data

countData   <- fread("prot.counts.txt")
colData   <- fread("colData.txt") # this is just a meta data file, could just use countData row names
accession <- fread("prot.acc.prots.named.out")

#=============================================================================== -->
#       Create DESeq object to calculate initial size factors
#=============================================================================== -->

# row_names column of countData object
row_names <- 1 

# creates a dds object and also subsets and orders colData by colnames of countData
dds <- DESeqDataSetFromMatrix(dt_to_df(countData,row_names), dt_to_df(colData)[names(countData)[-row_names],], ~1)

# get size factors
sf <- sizeFactors(estimateSizeFactors(dds))


#=============================================================================== -->
##       Filter reads
#============================================================================ -->

# memory efficient row deletion
delete <- function(DT, del.idxs) {           
  keep.idxs <- setdiff(DT[, .I], del.idxs);  # select row indexes to delete
  cols = names(DT);
  DT.subset <- data.table(DT[[1]][keep.idxs]); # this is the subsetted table
  setnames(DT.subset, cols[1]);
  for (col in cols[2:length(cols)]) {
    DT.subset[, (col) := DT[[col]][keep.idxs]];
    DT[, (col) := NULL];  # delete
  }
  return(DT.subset);
}


# simple filter to remove anything with a count less than 10
# best test this works first before mashing the full table..
#test <- countData[1:9,]
#rowSums(test[,-1])
#test <- delete(test,which(rowSums(test[,-1])<10)) 
#rowSums(test[,-1])

##### why is this here before combining accessions????
#remove <- which(rowSums(countData[,-1])<10)
#cat("removing",length(remove),"out of",nrow(countData),"total proteins\n")
#countData <- delete(countData,remove) # this will still take a while...
#countData <- countData[-1,]

# Combine accessions - remove blanks
blank_accession <- accession[V2=="",]
accession <- accession[V2!="",]
accession[,full:=gsub(";.*","",V2)]

# quick test
CD <- head(countData,4000)
test <- head(accession,4000)
CD <- test[CD,on=c("V1"="ProtID")]
cols=names(CD)[-1:-3]
CD1 <- CD[, lapply(.SD, sum, na.rm=TRUE), by=full,.SDcols=cols ]
# end test

# remove some issues in the accessions, e.g. remove the partial names from hypothetical proteins
# this may not be complete - update as required
accession[,full:=gsub("\001[a-zA-Z0-9]*\\.. "," ",full)]
accession[,full:=gsub("^(hypothetical protein) ([a-zA-Z0-9_\\.]*)$","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)( [a-zA-Z0-9_\\.\\-]*)$","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)( [a-zA-Z0-9_\\.\\-]* \\(plasmid\\)$)","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)( \\(plasmid\\)$)","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)(, partial \\(plasmid\\)$)","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)(, partial$)","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)( [a-zA-Z0-9_\\.\\-]*, partial$)","\\1",full)]

# merge countData and accessions
countData <- accession[countData,on=c("V1"="ProtID")] # this will take a bit (not so long about 3mins)
cols=names(countData)[-1:-3]

# combine counts on full name - these are proteis whoch can not be differentiated by function
countData <- countData[, lapply(.SD, sum, na.rm=TRUE), by=full,.SDcols=cols ] # this should be reasonably fast

# first row is the stuff which is unnamed, but I guess could still be kept..
# update using ref symatics
countData[is.na(full), full := "Unknown function"]

# remove hypothetical proteins - or not (I think I'll keep them in for now)
#idx <- which(countData$full=="hypothetical protein")
#countData <- delete(countData,idx)

# write data
sf <- sizeFactors(dds) 

# reorder columns - not necessary, but will save having to call again
colData <- colData[names(countData)[-1],on="Sample_name"]

# add size factors to metadata
colData[,sizefactors:=sf]

fwrite(countData,"countData.small.txt",sep="\t")
fwrite(colData,"colData.sf.txt",sep="\t",row.names = T)
fwrite(accession,"accessions.small.txt",sep="\t")
sfD <- data.table(sf=sf)
fwrite(sfD,"sizeFactors.txt",row.names = T)
```

There are a lot of resources available which ara related to UniProtID and IPR ID.   
The pipeline produces a UniProtID for each hit in the metagenome. Below are some scripts for working with UniProtIDs

#### Extract Uniprot ID from countData
```R
library(data.table)
counts <- fread("countData.txt") # from above - the full count table
ID <- counts[,1,drop=F] # gets the UniProtID
ID[,ProtID:=gsub("\\..*","",ProtID)] # removes the protein version (this is not necessary)
fwrite(ID,"UniProtIDs.txt") # write to file
```

#### Uniprot databases

There are several useful databases curretly stored at:
[ebi link](https://www.ebi.ac.uk/interpro/download/)
The protein to ipr database will be used here  


There's also the UniProt mapping file
[UniMap](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz)
  
```shell
wget https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz

wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
zcat idmapping.dat.gz |grep "RefSeq" >ncbi_uniprot.txt # extract ncbi IDs

```

#### Merge Uniprot database with IDs
```R
library(data.table)
ID <- fread("UniProtIDs.txt") # load IDS

map <- fread("ncbi_uniprot.txt",header=F) # load mappings

setnames(map,c("ProtID","ref","RefSeqID"))
map[,RefSeqID:=gsub("\\..*","",RefSeqID)] # remove .[0-9] freom refseqIDs
setkey(map,ProtID,RefSeqID) # add keys to map

uniprot <- fread("protein2ipr.dat") # load uniprot database
setnames(uniprot,c("ProtID","IPRID","Description","FAMID","start","end")) # add col names to uniprot data

# Best to add keys to these databases as they're big
setkey(ID,"ProtID)
setkey(uniprot,"ProtID") # this will take a bit to load into memory
# merge 
DT <- uniprot[ID,on=c("ProtID")] # left join ID on uniprot - this is fast with keying

# Extract unampped IDs (some of these will be ncbi IDs)
unmapped <- DT[is.na(IPRID),1,drop=F]
setkey(unmapped,"ProtID") # add a key
mapped <- map[unmapped] # left join unmapped on map (should match RefSeqIDs in the ID file [and link them to Uniprot ID])
uni_mapped <- mapped[complete.cases(mapped),]
setkey(uni_mapped,"ProtID","RefSeqID")

# lots of IPR IDs pointing at the same accession (different start and end points) - dups can be removes

DT[,c("start","end"):=NULL] # remove start and end columns
setkey(DT, ProtID, IPRID) # set a key
DT <- unique(DT,by = key(DT)) # retain unique 

# map extracts onto uniprot
DT2 <- uniprot[uni_mapped,on=c("ProtID")] # left join
DT2[,c("start","end","ref"):=NULL] # remove start and end columns
setkey(DT2, ProtID, IPRID) 
DT2 <- unique(DT2,by = key(DT2)) # retain unique 

DT <- DT[complete.cases(DT),] # remove missing values from DT
DT[,RefSeqID:=NA]

DT <- rbindlist(DT,DT2)

# write file
fwrite(DT,"count_ipr.txt",sep="\t")

```

#### Modify accessions file

```r
# add maps in accessions
accessions <- fread("accessions.txt")

accession[,full:=gsub(";.*","",V2)] # add full and remove altnames

accession[,full:=gsub("\001[a-zA-Z0-9]*\\.. "," ",full)]
accession[,full:=gsub("^(hypothetical protein) ([a-zA-Z0-9_\\.]*)$","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)( [a-zA-Z0-9_\\.\\-]*)$","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)( [a-zA-Z0-9_\\.\\-]* \\(plasmid\\)$)","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)( \\(plasmid\\)$)","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)(, partial \\(plasmid\\)$)","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)(, partial$)","\\1",full)]
accession[,full:=gsub("^(hypothetical protein)( [a-zA-Z0-9_\\.\\-]*, partial$)","\\1",full)]

accession[,full:=gsub("RecName: Full=","",full)]
# set full to all lower case
accession[,full:=tolower(full)]

accession[,V1:=gsub("\\..*","",V1)]

# now add in the UniProt and NCBI IDs from above
setkey(accession,V1)

map <- fread("ncbi_uniprot.txt",header=F) # load mappings
setnames(map,c("ProtID","ref","RefSeqID")) # set the names
map[,RefSeqID:=gsub("\\..*","",RefSeqID)] # remove .[0-9] freom refseqIDs
setkey(map,ProtID,REfSeqID) # add keys to map

NCBIs <- map[accession,on=c("RefSeqID"="V1")] # left join map on the accessions
#UNIs  <- map[accession,on=c("ProtID"="V1")]
NCBIs <- unique(NCBIs,by="RefSeqID") # there are a few duplicate UniProt names matching to NCBI IDs
#UNIs<-unique(UNIs, by="ProtID")

NCBIs[is.na(ProtID), ProtID := RefSeqID] # update NCBIs column, substitute NA Uniprot values with RefSeqID, this is correct it just seems to be the wrong way round
sum(accession$V1!=NCBIs$RefSeqID) # should = 0, if it the two tables are in the same order and we can copy the ProtIds over

accession[,ProtID:=NCBIs$ProtID] # add new ProtID column
setcolorder(accession,c("ProtID","full"))
fwrite(accession[,1:2],"uniprot_accessions.txt",sep="\t") # this file should be able to map back any descriptions from countData.small to an accession number

```

## END KAIJU


## Humann3

There's also humann3 which can do functional analysis - requires full alignment so will take a while to run.
The default alignment method uses Bowtie2 (which is no longer recommended for alignment) to create BAM files.
I'll install the default method for testing if it's any use before implementing a more appropriate alignment method.

NOTE: I'm not convinced that Humann gives us anything extra compared to Kaiju, certainly not for the complex samples, e.g. soil, that we deal with. I have implemented steps to extract protein annotations from Kaiju which is now in the above. This will give us exactly the same data as from Humann, but I'll need to implement some of their pathway analysis stuff.  

The pathway analysis is the most useful part of Humann3, it's just a shame it produces so few useful results in it's current implementation.  
I've yet to fully integrate it into the Kaiju pathway, but have made some progress.

### Installation

Couple of methods provided, I'd skip the conda version as it can result in glib issues

The default installations will need patching to the latest levels of the software.

Using conda environment

```shell

conda create --name humann3 -y python=3.7
conda activate humann3
```

This will probably fail on the HPC due to incompatible glib versions

```shell
conda config --add channels biobakery
conda install -y humann -c biobakery
```

Alternative is to install using pip excluding the binary enable build from source
Will still have to install metaphlan - try conda build

Beware there's a bug in MetaPhlan3 which renders it unusable (hardcoded link to a non existant dropbox location.. hum, auto download of "stuff" - good way to get users to trust you.)

Installing MetaPhlan 4.x should fix this issue


```shell
pip install humann --no-binary :all:
pip install humann --upgrade

# this will install an alternative build of metaphlan with correct glib version
mamba install -c conda-forge -c bioconda metaphlan
```

### Test build

Requires ~ 60G of memory

```shell
humann_test

humann_databases --download chocophlan DEMO humann_dbs
humann_databases --download uniref DEMO_diamond humann_dbs

mkdir tests
cd tests
wget https://github.com/biobakery/humann/raw/master/examples/demo.fastq.gz

humann -i demo.fastq.qz -o sample_results
```

The test run takes a long time to run due to the mapping step by Bowtie - we're talking hours here for a test fastq of 21,000 reads.

### Update databases
```shell
# pangenome
humann_databases --download chocophlan full [/PATH/TO/DATABASES] --update-config yes

# proteins
humann_databases --download uniref uniref90_diamond [/PATH/TO/DATABASES] --update-config yes

# annotation
humann_databases --download utility_mapping full [/PATH/TO/DATABASES] --update-config yes
```

### Running
```shell
humann -i sample_reads.fastq -o sample_results
```

#### Paired end reads
According to the authors of Humann3 the best method for dealing  with paired end data is simply to concatenate reather than merge the reads.

```shell
for FR in $PROJECT_FOLDER/data/$RUN/cleaned/*_1.fq.gz.filtered.fq.gz.cleaned.fq.gz; do
 RR=$(sed 's/_1.fq.gz/_2.fq.gz/' <<< $FR)
 S=$(sed 's/\(.*\/\)\(.*_1\)\(\..*\)/\2/' <<< $FR)
sbatch --mem=60000 -p long -c 20 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_humann.sh \
 $FR \
 $RR \
 ${S} \
 $PROJECT_FOLDER/data/$RUN/humann/ \
 $@
done
```

There are several output files from Humann, the most 


## Carnelia

### Installing Carnelian

You may need to edit the Makefile before running the final make. Remove the -static flag.

Requires python 2.7

```shell
conda create --name carnelian python=2.7
conda activate carnelian
```

```shell
conda install -c conda-forge vowpalwabbit
pip install -U scikit-learn
pip install pandas
pip install biopython==1.76

git clone https://github.com/snz20/carnelian
cd carnelian/util/ext
tar -zxf gdl-1.1.tar.gz
cd gdl-1.1
sh autogen.sh
./configure --prefix=$PWD/GDL/ #/home/xxx/local
make && make install
cd ../..
make
```
You may need to modify the Makefile in util folder (remove static) and possiblt add extra libraries.
Finally modify the library path to include the path where the new libraries have been installed.

### Install FragGeneScan
```shell
conda install -c bioconda fraggenescan
```

### Get pre-built model (vowpalwrabit 8.11)
```shell
# move to root of carnelian folder
cd ..
# make data directory
mkdir data
cd data
# prebuilt models
wget http://bergerlab-downloads.csail.mit.edu/carnelian/EC-2010-DB-model.tar.gz &
wget http://bergerlab-downloads.csail.mit.edu/carnelian/cog-model.tar.gz &

tar -zxf EC-2010-DB-model.tar.gz
tar -zxf cog-model.tar.gz
cd ..
```

## build model database
```shell
# Only necessary for new models (requires fasta file) - the above are prebuilt 
python2 ./carnelian.py train -k 8 --num_hash 4 -l 30 -c 5 data/EC-2010-DB models
```

### run annotation

carnelian searches for fasta files in a directory (files must have .fasta file type)

```shell
# convert fq 2 fa
for FR in $PROJECT_FOLDER/data/merged/*; do
  S=$(sed 's/\(.*\/\)\(.*_1\)\(\..*\)/\2/' <<< $FR)
  sbatch --mem=10000 -p short -c 1 $PROJECT_FOLDER/metagenomics_pipeline/common/scripts/sub_fq_2_fa.sh \
  $PROJECT_FOLDER/data/fasta \
  $FR \
  $RR \
  $S -n -z
done

# move files into unique folders
cd $PROJECT_FOLDER/data/fasta
for f in *; do
 S=$(sed 's/_.*//' <<<$f)
 mkdir $S
 mv $f $S/.
done

# set model location
MODEL=$PROJECT_FOLDER/metagenomics_pipeline/common/functional-analysis/carnelian/models/EC-2010-DB


# run annotation pipeline
# syntax
# python2 carnelian.py annotate -k 8 -n 20 <FASTA_DIRECTORY> <MODEL_DIRECTORY> <DIRECTORY?> FragGeneScan 

for DIR in $PROJECT_FOLDER/data/fasta/*/; do
  S=$(sed 's/fasta/carnelian/' <<< $DIR)
  sbatch --mem=80000 -p long -c 20 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_carnelian_annotate.sh \
  $S \
  $DIR \
  $MODEL \
  -k 8 -n 20
done
```

### Sum the annotations into count files
Carnelian has it's own method for doing this. It produces raw reads and one's normalised for length of annotation - the normalised reads are useful (maybe) for comparison between annotations.

It's pretty trivial to do the summing outside Carnelian.

```
#carnelian abundance labels_dir abundance_matrix_dir sampleinfo_file data/EC-2010-DB/ec_lengths.tsv
cd $PROJECT_FOLDER/data/carnelian
find . -type f -name *.label|xargs -I% mv % .
printf %b '#!/usr/bin/perl -s -w\nmy %annot_hash;\nwhile(<STDIN>) {\n  chomp;\n  $annot_hash{$_}++;\n}\nforeach (keys %annot_hash) {\n  print "$_\t$annot_hash{$_}\n";\n}\n' >counts.pl
chmod 755 counts.pl
for F in *.label; do
  S=$(sed 's/label/counts/' <<< $F)
  ./counts.pl < $F > $S &
done
```

Then combine counts into a countData object
```R
library(data.table)
library(tidyverse)

myfiles <- list.files(".",".*counts$",full.names=F,recursive=F)

# load count files
qq <- lapply(myfiles,fread) 

# apply names to appropriate list columns (enables easy joining of all count tables)
invisible(lapply(seq_along(qq),function(i) setnames(qq[[i]],"V2",gsub("_.*","",myfiles[[i]]))))

# merge count tables (full join)
countData <- Reduce(function(...) {merge(..., all = TRUE)}, qq)

# NA to 0
countData <- countData[,lapply(.SD, function(x) {x[is.na(x)] <- "0" ; x})]

# set the names
setnames(countData,"V1","EC")

# write the final table
fwrite(countData,"countData",sep="\t",quote=F,col.names=T)
```


## Kraken

```shell
# refseq
wget -c https://refdb.s3.climb.ac.uk/kraken2-microbial/hash.k2d &
wget -c https://refdb.s3.climb.ac.uk/kraken2-microbial/opts.k2d &
wget -c https://refdb.s3.climb.ac.uk/kraken2-microbial/taxo.k2d &

# full 
wget -c https://refdb.s3.climb.ac.uk/maxikraken2_1903_140GB/hash.k2d &
wget -c https://refdb.s3.climb.ac.uk/maxikraken2_1903_140GB/opts.k2d &
wget -c https://refdb.s3.climb.ac.uk/maxikraken2_1903_140GB/taxo.k2d &

# build data base
kraken-build --download-taxonomy --db kraken
kraken-build --add-to-library nr.gz --db kraken

```

## Sort bam files
```shell
for BAM in $PROJECT_FOLDER/data/assembled/aligned/megahit/*.bam; do
 PREFIX=$(echo $f|sed -e 's/\..*//')
  sbatch --mem-per-cpu 4500M -c 10 \
 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_bam_sort.sh \
 10 $PROJECT_FOLDER/data/sorted $PREFIX $BAM
done
```
## Run metabat

I need to scriptify this at some stage

```shell
# get list of bam filed for each assembly
BAM=$(for f in $PROJECT_FOLDER/data/assembled/aligned/sorted/$P1*; do echo $f; done|tr  '\n' ' ')

# Run metabat
sbatch --mem=40000 -p medium $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_metabat.sh \
$PROJECT_FOLDER/data/assembled/megahit/$PREFIX/$PREFIX.contigs.fa.gz \
$PROJECT_FOLDER/data/taxonomy_binning/${PREFIX}_BINS $BAMS
```

## Taxonomy assignment


Kaiju needs plenty of memory to load the nr database 
The sub_kaiju script is hard coded to use 20 processors for classification

```shell
sbatch --mem=120000 -p medium -c 20 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_kaiju.sh \
$PROJECT_FOLDER/data/kaiju/nodes.dmp \
$PROJECT_FOLDER/data/kaiju/names.dmp \
$PROJECT_FOLDER/data/kaiju/nr_euk/kaiju_db_nr_euk.fmi \
${PREFIX}.kaiju.out \
$PROJECT_FOLDER/data/taxonomy/$PREFIX \
-e 100 -m 100 -E 0.0000000001 -z 20 -v \
-i $PROJECT_FOLDER/data/taxonomy_binning/${PREFIX}_BINS/${PREFIX}.bins.fa

# drop  various columns from the output
awk -F"\t" '{print $1,$2,$4,$6,$8}' OFS="\t" < ATTINGHAM.kaiju.out >ATTINGHAM.k2.out
```

Taxon names can be added to the kaiju output using kaiju-addTaxonNames
```shell
kaiju-addTaxonNames -t $PROJECT_FOLDER/data/kaiju/nodes.dmp -n $PROJECT_FOLDER/data/kaiju/names.dmp -r superkingdom,phylum,class,order,family,genus,species -i ${PREFIX}.kaiju.out -o ${PREFIX}.names.out
```

### Add protein names to the bins

I've done this by using sqlite - it may not be the best method but it is fairly speedy and easy to implement. 
There is a slight problem that a sqlite query can have a maximum of 9999 'or' statements.

#### Set-up database
```shell
cd $PROJECT_FOLDER/data/kaiju/nr_euk
zgrep ">.*?\[" -oP nr.gz |sed 's/..$//'|sed 's/>//'|sed 's/MULTIGENE: //'|sed 's/ /|/' >nr.names
sqlite3 nr.db "CREATE TABLE nr(acc TEXT PRIMARY KEY, desc TEXT)"
sqlite3 -separator "|" nr.db ".import nr.names nr" 2>/dev/null
```

#### Extract protein names

This is a bit rubbish - better to use a temp table

```shell
# create multiple sql scripts with a maximum of 9999 terms
cd $PROJECT_FOLDER/data/taxonomy/$PREFIX
awk -F"\t" '{print $6}' OFS="," ${PREFIX}.kaiju.out|sed 's/.$//'|awk -F"," '{ for(i = 1; i <= NF; i++) { print "acc=\x27"$i"\x27 OR"; } }'|sed '$s/OR//'|split -l 9999

# loop through and execute the list of sql scripts
for f in x*; do 
 sed -i -e '$s/OR//' $f
 sed -i -e '1s/acc/SELECT * FROM nr WHERE acc/' $f
 sqlite3 $PROJECT_FOLDER/data/kaiju/nr_euk/nr.db <$f >> ${PREFIX}.prots.out
done
```
Temp table method

Download accessions2taxid for proteins and create a sqlite table from it:
ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz


```shell
for FILE in DR_13_EKDN230032028-1A_HKM3WDSX7_L4_1.pcounts ; do 
  echo "create temp table prots(acc text primary key,count text);" >query.sql
  echo ".mode csv" >>query.sql
  echo ".separator ;" >>query.sql
  #echo ".import 'test.txt' prots" >>query.sql
  echo ".import '$FILE' prots" >>query.sql
  echo "SELECT prots.acc,desc,count FROM prots LEFT JOIN nr ON prots.acc = nr.acc;" >>query.sql
  sqlite3 /data/data2/scratch2/deakig/nr_euk/nr.db <query.sql > ../$FILE.prots.named.out;
done


```



### Merge bin taxonomy and protein names

```R
library(tidyverse)
library(data.table)

PREFIX=xxx

dat <- fread(paste0(PREFIX,".names.out"),fill=T,sep="\t")
dat[,c("bin","contig"):=tstrsplit(V2,".fa.",fixed=T)]
dat[,acc:=sub(",.*","",V6)]
dat[,c("kingdom","phylum","class","order","family","genus","species"):=tstrsplit(V8,";",fixed=T)]

prot <- fread(paste0(PREFIX,".prots.out"),header=F)
setnames(prot,c("acc","protein"))
prot <- unique(prot)
dat <- prot[dat,on=c("acc==acc")]

dat[,(5:10):=NULL]

setnames(dat,c("V1","V2"),c("assigned","fullname"))

fwrite(dat,paste0(PREFIX,".taxoprot.txt"),quote=F,row.names=F,col.names=T)

```

## Count bin hits

### Generate gff file for all bins
```shell
cd $PROJECT_FOLDER/data/taxonomy_binning/${PREFIX}_BINS
$PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/awk_bin_to_gff.sh  < $PREFIX.bins.fa > $PREFIX.gff &
```

### count overlapping features

```shell
for BAM in $PROJECT_FOLDER/data/sorted/$P1*; do
  sbatch --mem 40000 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_bam_count.sh \
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm \
  $BAM \
  $PROJECT_FOLDER/data/taxonomy/$PREFIX/${PREFIX}.gff \
  $PROJECT_FOLDER/data/taxonomy/$PREFIX/map
done
```

#### Remove unused fields from cov output files 

```shell
cd $PROJECT_FOLDER/data/taxonomy/$PREFIX/map
for F in *.cov; do
  O=$(sed 's/_.*_L/_L/' <<<$F|sed 's/_1\.cov/.tab/')
  awk -F"\t" '{
   sub("ID=","",$(NF-1));
   sub(/fa\..*/,"fa",$(NF-1));
   print $1,$(NF-1),$NF 
  }' OFS="\t" $F > $O &
done 
```

### Merge count data

```R
library(data.table)

# location of files to load (assuming run after previuos step without leaving directory)
tmpdir <- "." # paste0(args[1],"/")

# load count files
qq <- lapply(list.files(tmpdir ,"*.tab",full.names=T),function(x) fread(x,sep="\t"))

# get the sample names  
names <- sub("_1\\.tab","",list.files(tmpdir ,"*.tab",full.names=F,recursive=F))

# kaiju taxonomy
names<-sub("([A-Z0-9]*)(_N.*_L)(.)(.*)","\\1_\\3",list.files(tmpdir ,"A.*",full.names=F,recursive=F))

# aggregate by domain
qa <- lapply(qq,function(DT) DT[,sum(V3),by = V2])

# apply names to appropriate list columns (enables easy joining of all count tables)
qa <- lapply(seq(1:length(qa)),function(i) {X<-qa[[i]];colnames(X)[2] <- names[i];return(X)})

# merge count tables (full join)
countData <- Reduce(function(...) {merge(..., all = TRUE)}, qa)

# rename first column
setnames(countData,"V2","Bin")

# NA to 0
countData <- countData[,lapply(.SD, function(x) {x[is.na(x)] <- "0" ; x})]

# write table
fwrite(countData,"countData",sep="\t",quote=F,row.names=F,col.names=T)

# without aggregation
qa <- qq

# merge contig and bin names
qa <- lapply(seq(1:length(qa)),function(i) {X<-qa[[i]];X[,sub_bin:=paste(V2,V1,sep=".")];X[,c("V1","V2"):=NULL];return(X)})

#kaiju
qa <- lapply(seq(1:length(qa)),function(i) {X<-qa[[i]];X[,sub_bin:=paste0(taxon_name,taxon_id)];X[,c("file","percent","taxon_name","taxon_id"):=NULL];return(X)})

# It's possible some of the sub bin names are duplicated...
qa <- lapply(qq,function(DT) DT[,sum(V3),by = sub_bin])

# apply names to appropriate list columns (enables easy joining of all count tables)
qa <- lapply(seq(1:length(qa)),function(i) {X<-qa[[i]];colnames(X)[2] <- names[i];return(X)})

# merge count tables (full join)
countData <- Reduce(function(...) {merge(..., all = TRUE)}, qa)

# NA to 0
countData <- countData[,lapply(.SD, function(x) {x[is.na(x)] <- "0" ; x})]

# write table
fwrite(countData,"sub_bin.countData",sep="\t",quote=F,row.names=F,col.names=T)

```

