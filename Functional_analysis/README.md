# Functional analysis

This is all in need of updating.  
I have a working pipeline using Kaiju, but there are other methods available which can be of use.  
However most of them are, so far, not worth implementing for complex microbiomes (i.e. anything that is not model organism related).  
This will change as more data is added to public databases.  

Carnelian (no longer developed???) and Humann3 are worth exploring. Humann3 does work, but the couple of times I've used it, it didn't produce any useful results - it also relies on full alignment (alignment files can be deleted after run), and is just generally slow. Has potential to be useful.

## Kaiju pipeline

### Taxonomy 

This needs to be done first - I'll move and link this at some stage

Kaiju will need setting up first - details below

First step is to set-up an nr database

(The database is stored in $PROJECT_FOLDER/data/kaiju for the following scripts)
```shell
kaiju-makedb -s nr_euk 
```

``

Run kaiju against paired end reads
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

#### Species counts

The below works but it is better to use the script for corrected counts further down

Best bet is to use kaiju tools to create a table of counts at the species rank - this can then be manipulated in R  
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
sbatch --mem=80G -p short -c 1 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_kaiju_correct_counts.sh \
 $K \
 $S \
 $PROJECT_FOLDER/data/kaiju_results \
 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm
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
#!/bin/bash

# create-sqlite.sh 
# create-sqlite.sh [NAME.db]

dbfile=$1
shift

sqlite3 <<EOT
.open $dbfile
CREATE TABLE names (
	taxID INT, 
	name VARCHAR(300), 
	unique_name VARCHAR(300), 
	name_class VARCHAR(300)
);
	
CREATE TABLE nodes (
	taxID INT, 
	parent_taxID INT, 
	rank VARCHAR(300)
);

.shell echo Importing names
.separator '|'
.import names_tabless.dmp names
.shell echo Indexing names.
CREATE INDEX name_idx ON names(taxID);
CREATE INDEX name_class ON names (name,name_class);

.shell echo Importing nodes
.separator '|'
.import names_tabless.dmp nodes
.shell echo Indexing nodes.
CREATE INDEX nodes_idx ON nodes(taxID);

EOT
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

### Functional analysis 

The output data contains a huge number of proteins, almost all which can not be distinguished at the functional level - they may have different functions, but the details are not yet available. Due to this massive duplication, it makes sense to shrink the data down to unique protein functions.  

Below is a short bit of R code to do this. Not all steps are stricktly neccessary (the DESeq stuff can all be dropped if it's not going to be used)

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
```shell
wget https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz
```

#### Merge Uniprot database with IDs
```R
library(data.table)
ID <- fread("UniProtIDs.txt") # load IDS
uniprot <- fread("protein2ipr.dat") # load uniprot database
setnames(uniprot,c("ProtID","IPRID","Description","FAMID","start","end")) # add col names to uniprot data

# Best to add keys to these databases as they're big
setkey(ID,"ProtID)
setkey(uniprot,"ProtID") # this will take a bit to load into memory

# merge 
DT <- uniprot[ID,on=c("ProtID")] # left join ID on uniprot - this is fast with keying

# lots of IPR IDs pointing at the same accession (different start and end points) - dups can be removes

DT[,c("start","end"):=NULL] # remove start and end columns
setkey(DT, ProtID, IPRID) # set a key
DT <- unique(DT,by = key(DT)) # retain unique 

# write file
fwrite(DT,"count_ipr.txt",sep="\t")
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

