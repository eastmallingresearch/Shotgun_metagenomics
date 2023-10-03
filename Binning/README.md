## Functional binning
This pipeline is based on the HirBin (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3686-6) pipeline  
I've replaced most of the code as it was not capable of dealing with soil metagenomics in a reasonable time (still only uses Perl and R - both pretty slow). Instead some of the processes use a lot of memory.

HirBin identifies functional domains (Pfam/Tigrfam etc.), then adds an additional step to cluster the bins into a set of sub-bins. The identification of the domains is the only bit left of the HirBin pipeline I haven't needed to rewrite (uses HMMER to do the actual id of the domains).

The pipeline described below is an example taken from the oak decline project

I'm likely also to implement an assembly free functional binning pipeline - Metakallisto or Carnelian (this was actually the original plan, but I got sidetracked into doing assemblies).  


```shell
# set some variables
PROJECT_FOLDER=~/projects/metagenomics
PREFIX=MYSAMPLES
P1=${PREFIX:0:1}
ln -s ~/pipeline/metagenomics $PROJECT_FOLDER/metagenomics_pipeline
```

### Annotation
Annotating an assembly uses functionalAnnotaion.py (HirBin)  

```shell
functionalAnnotation.py -m METADATA_FILE -db DATABASE_FILE -e EVALUE_CUTOFF -n N -p MAX_ACCEPTABLE_OVERLAP
```

fun_bin.sh uses functionAnnotation.py to find domians, but splits assembly into 20,000 droid chunks for running on cluster.   

![#f03c15](https://placehold.it/15/f03c15/000000?text=+) Assembly can be gz compressed or uncompressed
```shell
$PROJECT_FOLDER/metagenomics_pipeline/scripts/fun_bin.sh \
 1 $PROJECT_FOLDER/data/assembled/$PREFIX \
 $PREFIX.contigs.fa \
 ~/pipelines/common/resources/pfam/Pfam-A.hmm \
 -e 1e-03
```

#### Concatenate annotation output
The output from fun_bin.sh will be in a large number of sub directories inside a tmp.* directory in the (above) $PROJECT_FOLDER/data/assembled/megahit/$PREFIX directory

The below scripts concatenate the output

```shell
find -type f -name X.gff|head -n1|xargs -I% head -n1 % >$PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.gff
find -type f -name X.gff|xargs -I% grep -v "##" % >>$PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.gff
find -type f -name X.pep|xargs -I% cat % >$PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.pep
find -type f -name X.hmmout|xargs -I% grep -v "#" % | \
awk -F" " '($21~/^[0-9]+$/) && ($20~/^[0-9]+$/) {print $4,$1,$20,$21,$3,$7}' OFS="\t"| \
$PROJECT_FOLDER/metagenomics_pipeline/scripts/subbin_domain_extractor.pl \
> $PROJECT_FOLDER/data/binning/$PREFIX/$PREFIX.domains
```

#### Get gff with different MAX_ACCEPTABLE_OVERLAP - example below will produce gff with all overlapping features
```shell
# this is a script from HirBin
con_coor.py -p 1 -o $PREFIX.2.gff -d $PREFIX.pep -m $PREFIX.hmmout
```

#### Sub bins

##### Extract functional domains as aa strings
arg[4] = number of chunks to read the data in (larger has lower memory footprint, but will be slower)  
arg[5] = T/F for parallel processing (currently hard coded to 10 cores) - this will slurp up a lot of memory
```shell
Rscript subbin_fasta_extractor.R $PREFIX.domains $PREFIX.pep $PREFIX_hirbin_output 100 T
```

##### Cluster extracted domains
```shell
PIPELINE.sh - c cluster_super_fast SERVERS CORES INPUT_DIR OUTPUT_DIR CLUS_ID
```
```shell
$PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c cluster_super_fast \
  blacklace[01][0-9].blacklace 100 \
  $PROJECT_FOLDER/data/binning/$PREFIX/${PREFIX}_clustering/forClustering \
  $PROJECT_FOLDER/data/binning/$PREFIX/${PREFIX}_clustering/clust0.7 \
  0.7 
```

##### concatenate clustering output
reduced.txt (and the \*.uc files) is a map of domains (column 2) to sub bins (column 1)
```shell
cat $PROJECT_FOLDER/data/binning/${PREFIX}_clustering/clust0.7/*.uc > $PROJECT_FOLDER/data/binning/$PREFIX/reduced.txt
```

### Mapping
Mapping with BBMAP - anything that outputs a bam file can be substituted
align reads to assembly - will need to index first
```shell
bbmap.sh ref=$PREFIX.contigs.fa.gz usemodulo=t #k=11
```
Then map
```shell
for FR in $PROJECT_FOLDER/data/fastq/$P1*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c align -p bbmap \
  16 blacklace[01][0-9].blacklace \
  $PROJECT_FOLDER/data/assembled/aligned/megahit \
  $PREFIX \
  $PROJECT_FOLDER/data/assembled/megahit/$PREFIX/${PREFIX}.contigs.fa.gz \
  $FR \
  $RR \
  maxindel=100 \
  unpigz=t \
  touppercase=t \
  path=$PROJECT_FOLDER/data/assembled/megahit/$PREFIX/ \
  usemodulo=T \ 
  -Xmx31g
done
```

### Count coverage
bam_scaffold_count.pl will output a cov file (bedtools output file). It's vastly faster than bedtools for counting overlapping features (with no minimum overlap). If minimum overall > 1 is required use bedtools.
NOTE: by default this counts unique F + R and paired overlap (*i.e.* if both overlap it gets a score of 1), bedtools counts F + R  (i.e. if both overlap it gets a score of 2). Set BEDTOOLS to true to count as per bedtools

```shell
bam_scaffold_count.pl GFF_FILE [BEDTOOLS] < SAM_LINE
samtools view bam_file|~/pipelines/metagenomics/scripts/bam_scaffold_count.pl $PREFIX.gff > bam_counts.cov
samtools view bam_file|~/pipelines/metagenomics/scripts/bam_scaffold_count.pl $PREFIX.gff true> bedtools_counts.cov
```
Get counts - pipeline
```shell
for BAM in $PROJECT_FOLDER/data/assembled/aligned/megahit/$P1*.bam; do
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c coverage -p bam_count \
  blacklace[01][0-9].blacklace \
  $BAM \
  $PROJECT_FOLDER/data/$PREFIX/${PREFIX}.gff \
  $PROJECT_FOLDER/data/assembled/counts/megahit \
done
```

#### Bin counts (deprecated - sub bin counter will do this as well)
Get bin counts from the cov format domain counts and output as count matrix
```shell
Rscript $PROJECT_FOLDER/metagenomics_pipeline/scripts/cov_count.R cov_file_location "$PREFIX.countData"
```

#### Sub bin counts
Create tab files from cov files (originally for consistency with HirBin, but now to reduce processing in R)  
I could get bam_scaffold_count.pl to output in this format, or just bung the script below into the bam_count pipeline to produce both...
```shell
for F in *.cov; do
  O=$(sed 's/_.*_L/_L/' <<<$F|sed 's/_1\.cov/.tab/')
  awk -F"\t" '{sub("ID=","",$(NF-1));OUT=$1"_"$(NF-1)"_"$4"_"$5;print OUT,$(NF-1),$4,$5,$NF}' OFS="\t" $F > $O &
done 
```

Parse the tab files to create countData matrix
```shell
Rscript subbin_parser.R reduced.txt tab_file_location $PREFIX
```

## Asembly free functional binning
Carnelian might work

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

## Humann3

There's also humann3 which can do functional analysis - requires full alignment so will take a while to run.
The default alignment method uses Bowtie2 (which is no longer recommended for alignment) to create BAM files.
I'll install the default method for testing if it's any use before implementing a more approriate alignment method.

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


# Taxonomy binning

Taxonomy binning uses a mashup of various pipelines. I did try and implement Anvio, but it is vastly too slow (and memory hungry) for the size of data involved in soil metegenomics.  

The binning is done by Metabat which requires sorted bam files (as produced above).   
![#f03c15](https://placehold.it/15/f03c15/000000?text=+) Note - cluster jobs are written for slurm not grid engine

Kaiju by itself maybe of use for binning raw (or filtered) reads - use as below

## Kaiju binning and taxonomy

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

### Species counts

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

### Correct counts with DiTASiC
Something goes here...
It does - DiTASic is not going to work with Kaiju. But I may be able to adapt the model they use to apply to a protein database - won't be easy though. In fact it is full of problems, giving up on the idea for now

Will assign multimapping reads based on the proportion of uniqueliy mapped reads per taxon.

This has all been scripted now...

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


#### IGNORE - THIS IS SCRIPTED IN THE ABOVE

Details of how the counts are corrected

Will need counts for all taxon entries, and which are multi hits.
```shell
for K in $PROJECT_FOLDER/data/kaiju_taxonomy/${P1}*.out; do
  S=$(sed 's/\(.*\/\)\(.*_1\)\(\..*\)/\2/' <<< $K)
  awk -F"\t" '{print gsub(/,/,",",$5) "\t" $5}' < $K > ${S}.new_counts  
done
```

A quick perl script to add all taxons to a hash

This is in 

```perl
#!/usr/bin/perl -s -w
use List::Util 'sum';
my %taxon_hash; 
while(<STDIN>) {
  chomp;
  my @array=split /,/,$_;
  foreach(@array) {
    $taxon_hash{$_}++;
  }
}

my $taxon_sum = sum values %taxon_hash;

foreach (keys %taxon_hash) {
  print "$_\t$taxon_hash{$_}\n";
}
```
add them together
```shell
for f in *.new_counts; do
  S=$(sed 's/\..*//' <<< $f)
  awk -F"\t" '{if($1>0){print $2}}' $f|./test.pl > ${S}.pcounts & 
done
```
get the totals without multi-hits as well
```shell
for f in *.new_counts; do
  S=$(sed 's/\..*//' <<< $f)
  awk -F"\t" '{if($1==1){print $2}}' $f|./test.pl > ${S}.ncounts & 
done
```
Then add proportional values for each multi-hit to the totals (in R for easy multiprocessing - eats memory though)
```R
library(data.table)
#   library(parallel)
library(tidyverse)

del.rows <- function(DT, del.idxs) {           # pls note 'del.idxs' vs. 'keep.idxs'
  keep.idxs <- setdiff(DT[, .I], del.idxs);  # select row indexes to keep
  cols = names(DT);
  DT.subset <- data.table(DT[[1]][keep.idxs]); # this is the subsetted table
  setnames(DT.subset, cols[1]);
  for (col in cols[2:length(cols)]) {
    DT.subset[, (col) := DT[[col]][keep.idxs]];
    DT[, (col) := NULL];  # delete
  }
  return(DT.subset);
}

file_suffix <- gsub("\\..*","",list.files(".",".*pcounts$",full.names=F,recursive=F))

freadFun <- function(i) fread(cmd=paste0("cat header.txt ",i,".new_counts"),fill=T,header=T)

pcounts    <- lapply(file_suffix,function(i) fread(paste0(i,".pcounts"))) 
ncounts    <- lapply(file_suffix,function(i) fread(paste0(i,".ncounts"))) 
# new_counts <- mclapply(file_suffix,freadFun,mc.cores=10)
fwrite(as.data.table(t(as.data.table(paste0("V",1:23)))),"header.txt",sep=",")
new_counts <- lapply(file_suffix,freadFun)

new_counts <- lapply(new_counts,function(DT) del.rows(DT,c(1,which(DT[,1]<2))))
del_cols <- lapply(new_counts,function(DT)names(DT)[c(1,length(names(DT)))])
invisible(lapply(seq_along(new_counts),function(i){new_counts[[i]][,(del_cols[[i]]):=NULL];new_counts[[i]][,seq:=(1:nrow(new_counts[[i]]))]}))
lapply(pcounts,function(DT) DT[,V1:=as.character(V1)])
lapply(ncounts,function(DT) DT[,V1:=as.character(V1)])

multi_hits <- lapply(new_counts,melt,id.vars="seq")
invisible(lapply(multi_hits,function(DT){DT[,variable:=NULL];DT[,value:=as.factor(value)];setnames(DT,"value","V1")}))
multi_hits <- lapply(multi_hits,function(DT) del.rows(DT,which(DT[,2]=="")))

final_counts <- lapply(seq_along(multi_hits),function(i) {
  X <- pcounts[[i]][multi_hits[[i]],on="V1"]
  X[,prop:=V2/sum(V2),by=seq]
  X<- X[,lapply(.SD,sum),by=V1,.SDcols="prop"]
  X <- merge(X, ncounts[[i]], all=TRUE)
  X[,tot:=rowSums(.SD, na.rm = TRUE), .SDcols =c("V2","prop")]
  #ncounts[[i]][X,tot:=V2+prop,on="V1"]
})

lapply(seq_along(final_counts),function(i) fwrite(final_counts[[i]],paste0(file_suffix[[i]],".corrected_counts"),sep="\t"))

```

END IGNORE

### Produce counts and taxonomy

This needs editing, the taxonomy is not correct. It uses the output from Kaiju to assign the taxonomy - this is not complete as it excludes various taxa. Possibly only assignes to species levels, rather than higher taxonomic ranks. Will need to assign taxonomy directly from the names.dmp and nodes.dmp files.


RUN THIS PART

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

THE BIT BELOW IS NOT WORKING CORRCTLY

```r
# load taxonomy data (and false counts)
qq <- lapply(file_suffix,function(i) fread(paste0(i,".kaiju.counts"))) 
invisible(lapply(qq,function(DT)DT[,c("file","percent","reads"):=NULL]))

# the below does not find all the taxonomy included in the output.
taxData <- Reduce(function(...) {merge(..., all = TRUE)}, qq)
taxData[,taxon_id:=as.character(taxon_id)]
fwrite(taxData,"taxData.txt",sep=";",quote=F,row.names=F,col.names=F)
fread("taxData.txt",fill=T)
taxData[,V9:=NULL]
setnames(taxData,c("taxon_id","kingdom","phylum","class","order","family","genus","species"))
fwrite(taxData,"taxData.txt",sep=";",quote=F)

fwrite(countData,paste0("countData"),sep="\t",quote=F,row.names=F,col.names=T)

```

END NOT WORKING

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
The below runs the sql query via sqldf - the whole lot could probably be replaced with aome dplyr code.

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
