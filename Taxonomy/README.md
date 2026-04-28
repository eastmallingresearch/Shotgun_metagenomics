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
  sbatch --mem=250000 -p himem -c 20 $PROJECT_FOLDER/metagenomics_pipeline/scripts/sub_kraken.sh \
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
  sbatch --mem=20000 -p short $PROJECT_FOLDER/metagenomics_pipeline/scripts/sub_bracken.sh \
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
  sbatch --mem=200000 -p himem -c 20 $PROJECT_FOLDER/metagenomics_pipeline/scripts/sub_kaiju.sh \
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

#### Species counts

The below works but it is better to use the script for corrected counts further down

Use kaiju tools to create a table of counts at the species rank - this can then be manipulated in R  
The kaiju2table program is fast.

```shell
for K in $PROJECT_FOLDER/data/kaiju_taxonomy/${P1}*.out; do
S=$(sed 's/\(.*\/\)\(.*_1\)\(\..*\)/\2/' <<< $K)
sbatch --mem=12000 -p short -c 1 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_kaiju_table.sh \
 $PROJECT_FOLDER/data/kaiju/nodes.dmp \
 $PROJECT_FOLDER/data/kaiju/names.dmp \
 ${S}.kaiju.counts \
 $PROJECT_FOLDER/data/kaiju_results/ \
 $K  
done

# kaiju2table -t $PROJECT_FOLDER/data/kaiju/nodes.dmp -n $PROJECT_FOLDER/data/kaiju/names.dmp -r species -l superkingdom,phylum,class,order,family,genus,species -o ${K}.counts $K &
```

#### Correct counts 

Something goes here...

Plan was to implement something like DiTASic, but DiTASic is not going to work with Kaiju. I may be able to adapt the model they use to apply to a protein database - won't be easy though. In fact it is full of problems, giving up on the idea for now

Will assign multimapping reads based on the proportion of uniqueliy mapped reads per taxon.

```shell
# correct counts
for K in $PROJECT_FOLDER/data/kaiju_taxonomy/${P1}*.out; do
S=$(sed 's/\(.*\/\)\(.*_1\)\(\..*\)/\2/' <<< $K)
sbatch --mem=20G -p medium -c 1 $PROJECT_FOLDER/metagenomics_pipeline/scripts/sub_kaiju_correct_counts.sh \
 $K \
 $S \
 $PROJECT_FOLDER/data/kaiju_results \
 $PROJECT_FOLDER/metagenomics_pipeline/scripts
done
```

#### Produce counts and taxonomy

```R
library(data.table)
library(tidyverse)

file_suffix <- gsub("\\..*","",list.files(".",".*corrected.counts$",full.names=F,recursive=F))

# load count files
qq <- lapply(file_suffix,function(i) fread(paste0(i,".corrected.counts"))) 

# apply names to appropriate list columns (enables easy joining of all count tables)
invisible(lapply(seq_along(qq),function(i) setnames(qq[[i]],"tot",file_suffix[i])))
invisible(lapply(qq,function(DT)DT[,c("prop","V2"):=NULL]))
# merge count tables (full join)
countData <- Reduce(function(...) {merge(..., all = TRUE)}, qq)

# NA to 0
countData <- countData[,lapply(.SD, function(x) {x[is.na(x)] <- "0" ; x})]
# 
# # add OTU column
# countData[,OTU:=paste0("OTU",1:nrow(countData))]
# 
# count_cols <- names(countData)[-1]
# countData[,(count_cols):=lapply(.SD,as.numeric),.SDcols=count_cols]
setnames(countData,"V1","taxon_id")
```


Using sql is one way of doing this 
Use names.dmp and nodes.dmp from kaiju download (or directly from NCBI - they're taxonomy file) to create sqlite database

create-sqlite.sh can create a taxonomy database from names.dmp and nodes.dmp

```shell
create_sqlite_taxonomy.sh
```

The name and nodes files need to be in the same folder as the sqlite script.  

The names and nodes files have an odd separator field \t|\t possibly. It's best to remove tabs from the files before making the databases

```shell
sed -i 's/\t//g' names.dmp > names_tabless.dmp
sed -i 's/\t//g' nodes.dmp > nodes_tabless.dmp

```

Imports the files into two tables named as above (the script will throw up a lot of info messages as both cmp files contain a lot of irrelevent columns - there's no way to suppress these messages in sqlite [as far as I know])


```shell
./create-sqlite.sh taxonomy.db
```

Below is a sql script to query the taxonomy database  - but it's possibly easier to do everything in R rather than via sqlite directly

```sql
-- recursive query to get all parents of id
WITH RECURSIVE
  taxonomy(i) AS (
    VALUES(2927976)
    UNION
    SELECT parent_taxID FROM nodes,taxonomy
    WHERE nodes.taxID = taxonomy.i
  )
  SELECT nodes.rank,nodes.taxID,name FROM nodes
  INNER JOIN names ON 
  nodes.taxID = names.taxID
  WHERE nodes.taxID IN taxonomy AND names.name_class="scientific name";

-- no rank|1|root
-- superkingdom|2|Bacteria
-- genus|44675|Geothrix
-- phylum|57723|Acidobacteriota
-- no rank|131567|cellular organisms
-- class|533205|Holophagae
-- order|574975|Holophagales
-- family|574976|Holophagaceae
-- no rank|2647902|unclassified Geothrix
-- species|2927976|Geothrix sp. Red802
```

It's a lot easier to do the querying via R than sqlite (well I think so)
The below runs the sql query via sqldf - the whole lot could probably be replaced with some, megaverbose(TM) dplyr code.

```r
library(tidyverse)
library(data.table)
library(sqldf)

countData <- fread("countData")

query <- function(i){
  paste0("WITH RECURSIVE
  taxonomy(i) AS (
    VALUES(",i,")
    UNION
    SELECT parent_taxID FROM nodes,taxonomy
    WHERE nodes.taxID = taxonomy.i
  )
  SELECT nodes.rank,name FROM nodes
  INNER JOIN names ON 
  nodes.taxID = names.taxID
  WHERE nodes.taxID IN taxonomy 
  AND names.name_class='scientific name' AND 
    (rank='species' OR
     rank='genus' OR
     rank='family' OR
     rank='order' OR
     rank='class' OR
     rank='phylum' OR
     rank='kingdom' OR
     rank='superkingdom')
   ")
}
fetch <- function(i){
  X<- sqldf(query(i),connection=con)
  X<-setNames(data.frame(t(X[,-1])),X[,1])
  setDT(X)
  
  cols <- c(taxon_id=i,
			superkingdom = NA_real_, 
            kingdom = NA_real_, 
            phylum = NA_real_,
            class = NA_real_,
            order = NA_real_,
            family = NA_real_,
            genus = NA_real_,
            species = NA_real_)
  
  X<-add_column(X, !!!cols[setdiff(names(cols), names(X))])
  setcolorder(X,c("taxon_id","superkingdom","kingdom","phylum","class","order","family","genus","species"))
  X
}
con <- DBI::dbConnect(RSQLite::SQLite(),"taxonomy.db",flags=SQLITE_RO)

taxData <- rbindlist(apply(countData[,1],1,fetch))
taxData[,taxon_id:=countData$taxon_id]
setcolorder(taxData,"taxon_id")
fwrite(taxData,"taxData",sep="\t")

# dplyr version
# X <- tbl(con,"names") # produces some weird tibble like structure which is not subsetable - tbl(con,"names")[,1] - results in an error
# presumably there's some very verbose dplyr syntax (maybe even as verbose as sql...) to query these tbl classes
# o.k the tbl constructs a query, need to use collect(query) to return results as a table
# Y <-tbl(con,"nodes")
# Z <- inner_join(X,Y,y="taxID")
# show_query(Z)
# nms <- collect(X)
# nds <- collect(Y)
# setDT(nms)
# setDT(nds)

```

### Extract protein accessions from Kaiju ouptut 
```shell
for f in *.kaiju.out; do
  S=$(echo $f|awk -F"\/" '{print $NF}')
  grep ^C $f|awk -F"\t" '{split($6,a,",");print a[1]}' >> $S.protein_accessions.txt
done
```

Some code goes here to concatenate the counts and accessions
