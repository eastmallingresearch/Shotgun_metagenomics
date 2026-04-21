I have two pipeline kraken and kaiju for taxonomic profiling.  
Both are useful for soil metagenomics, and it is useful to run both

# Kraken pipeline
The CropDiversit HPC maintains several version of kraken databases  
Will check how to use these

## Database options

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

# prebuilt options
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20260226.tar.gz

```

## Run Kraken

```shell
for FR in $PROJECT_FOLDER/data/cleaned/*_1*.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  S=$(sed 's/\(.*\/\)\(.*_1\)\(\..*\)/\2/' <<< $FR)
  sbatch --mem=250000 -p himem -c 20 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_kraken.sh \
  $PROJECT_FOLDER/PATHTOKRAKENDB \
  $FR \
  $RR \
  ${S}.kraken.out \
  ${S}.kraken.report.out \
  $PROJECT_FOLDER/data/taxonomy/kraken/ \
  20
done
```
## Estimate abundances

Braken is used to estimate abundances from kraken reports...

However kraken seem to already do this  - and the output is more useful.   
The below is the standard bracken pipeline. It's fast, but doesn't give the full taxonomy, only the level bracken was run at (e.g. only gives species, not the full taxonomy).   

```shell
for KR in $PROJECT_FOLDER/data/taxonomy/kraken/*.report.out; do
  S=$(sed 's/\(.*\/\)\(.*_1\)\(\..*\)/\2/' <<< $KR)
  sbatch --mem=20000 -p short $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_bracken.sh \
  $PROJECT_FOLDER/PATHTOKRAKENDB \
  $KR \
  200 \ # read length
  ${S} \
  $PROJECT_FOLDER/data/taxonomy/bracken/ 
done

```

The output from kraken does include the full and the recacalculated (brakenised) read counts.

The below R script will wrangle the kraken output into count and taxonomy files
```R
#CODE TO COMBINE countData and taxData into single files (one for each)
library(data.table)

# Set the path
PATH <- "."

file_names <- list.files(PATH, "species.out", full.names = TRUE)

# Taxonomic ranks in hierarchical order
ranks <- c("R","R1","R2","K","P","C","O","F","G","S")

# Helper function:
# update a rank, then clear all lower ranks
clear_lower <- function(curr, rank_name, rank_levels) {
  i <- match(rank_name, rank_levels)
  
  if (!is.na(i) && i < length(rank_levels)) {
    for (j in (i + 1):length(rank_levels)) {
      curr[[rank_levels[j]]] <- NA_character_
    }
  }
  
  curr
}

# Function to process one kraken/bracken report file
process_one_file <- function(myfile) {
  
  # read tab-separated file
  DT <- fread(myfile, sep = "\t", header = FALSE, strip.white = TRUE, quote = "")
  
  # keep only taxonomy ranks of interest
  x <- DT[V4 %in% ranks]
  
  # current taxonomy path
  current <- setNames(as.list(rep(NA_character_, length(ranks))), ranks)
  
  # store species rows
  out <- vector("list", 0L)
  
  for (i in seq_len(nrow(x))) {
    rk <- x$V4[i]
    tx <- x$V6[i]
    
    current[[rk]] <- tx
    current <- clear_lower(current, rk, ranks)
    
    if (rk == "S") {
      row_out <- c(
        as.list(x[i]),
        current
      )
      out[[length(out) + 1L]] <- row_out
    }
  }
  
  # return empty objects safely if no species rows found
  if (length(out) == 0L) {
    sample_name <- gsub("_1\\..*", "", basename(myfile))
    
    countData <- data.table(TaxID = character(), tmp = numeric())
    setnames(countData, "tmp", sample_name)
    
    taxData <- data.table(
      TaxID = character(),
      R = character(), R1 = character(), R2 = character(),
      K = character(), P = character(), C = character(),
      O = character(), F = character(), G = character(), S = character()
    )
    
    return(list(countData = countData, taxData = taxData))
  }
  
  # bind all species rows
  result <- rbindlist(out, fill = TRUE)
  
  # remove unused columns
  result[, c("V1", "V3", "V4", "V6") := NULL]
  
  # rename columns:
  # V2 = count column
  # V5 = TaxID
  sample_name <- gsub("_1\\..*", "", basename(myfile))
  setnames(result, c("V2", "V5"), c(sample_name, "TaxID"))
  
  # split into count and taxonomy
  countData <- result[, .(TaxID, get(sample_name))]
  setnames(countData, c("TaxID", sample_name))
  
  taxData <- unique(result[, .(TaxID, R, R1, R2, K, P, C, O, F, G, S)])
  
  list(countData = countData, taxData = taxData)
}

# Process all files
res_list <- lapply(file_names, process_one_file)

# Extract count and taxonomy tables
count_list <- lapply(res_list, `[[`, "countData")
tax_list   <- lapply(res_list, `[[`, "taxData")

# Combine all count tables by TaxID
countData <- Reduce(
  function(x, y) merge(x, y, by = "TaxID", all = TRUE),
  count_list
)

# Replace missing counts with zero
for (j in seq_along(countData)) {
  if (names(countData)[j] != "TaxID") {
    set(countData, which(is.na(countData[[j]])), j, 0)
  }
}

# Combine taxonomy tables into one, remove duplicates
taxData <- unique(rbindlist(tax_list, fill = TRUE), by = "TaxID")

# Optional: keep taxData in same order as countData
taxData <- taxData[match(countData$TaxID, taxData$TaxID)]

# write out tables
fwrite(taxData,"../kraken.taxData",sep="\t")
fwrite(countData,"../kraken.countData",sep="\t")

```



# Kaiju pipeline
Kaiju needs a database - for taxonomy probably best to stick to the nr_euk database. Other, or custom options are available 


## Download and setup database  

```shell
kaiju-makedb -s nr_euk

```

## Run kaiju
```shell
for FR in $PROJECT_FOLDER/data/cleaned/*_1*.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  S=$(sed 's/\(.*\/\)\(.*_1\)\(\..*\)/\2/' <<< $FR)
  sbatch --mem=120000 -p medium -c 20 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_kaiju.sh \
  $PROJECT_FOLDER/data/kaiju/nodes.dmp \
  $PROJECT_FOLDER/data/kaiju/names.dmp \
  $PROJECT_FOLDER/data/kaiju/nr_euk/kaiju_db_nr_euk.fmi \
  ${S}.kaiju.out \
  $PROJECT_FOLDER/data/kaiju_taxonomy/ \
  -z 20 -v \
  -i $FR \
  -j $RR
done
```
