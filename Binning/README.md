## Functional binning
This pipeline is based on the HirBin (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3686-6) pipeline  
I've replaced most of the code as it is not capable of dealing with soil metagenomics in a reasonable time. Some of the processes use a lot of memory

HirBin identifies functional domains (Pfam/Tigrfam etc.), then adds an additional step to cluster the bins into a set of sub-bins. The identification of the domains is the only bit left of the HirBin pipeline I haven't needed to rewrite (uses HMMER to do the actual id of the domains).

The pipeline described below is an example taken from the oak decline project
```shell
# set some variables
PROJECT_FOLDER=~/projects/Oak_decline/metagenomics
PREFIX=BIGWOOD
P1=${PREFIX:0:1}
ln -s ~/pipeline/metagenomics $PROJECT_FOLDER/metagenomics_pipeline
```

### Annotation
Annotating an assembly uses functionalAnnotaion.py (HirBin),

```shell
functionalAnnotation.py -m METADATA_FILE -db DATABASE_FILE -e EVALUE_CUTOFF -n N -p MAX_ACCEPTABLE_OVERLAP
```

fun_bin.sh uses functionAnnotation.py to find domians, but splits assembly file into 20,000 droid chunks for running on cluster.  

```shell
$PROJECT_FOLDER/metagenomics_pipeline/scripts/fun_bin.sh \
 1 $PROJECT_FOLDER/data/assembled/megahit/$PREFIX \
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
bam_scaffold_count.pl will output a cov file (bedtools output file). Vastly faster than bedtools for counting overlapping features (with no minimum overlap). If minimum overall > 1 is required use bedtools.
NOTE: this counts unique F + R and paired overlap (*i.e.* if both overlap it gets a score of 1), bedtools counts F + R  (i.e. if both overlap it gets a score of 2)  
bam_scaffold_count.pl can be run in two modes, bin counting or domain hit counting mode. Add cov (or anything) as a second flag for domain hit counting - this is necessary to get sub bin counts.

```shell
samtools view bam_file|~/pipelines/metagenomics/scripts/bam_scaffold_count.pl $PREFIX.gff > bam_counts.txt
samtools view bam_file|~/pipelines/metagenomics/scripts/bam_scaffold_count.pl $PREFIX.gff cov> bam_file.cov
```
Get counts - pipeline
```shell
for BAM in $PROJECT_FOLDER/data/assembled/aligned/megahit/$P1*.bam; do
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c coverage -p bam_count \
  blacklace[01][0-9].blacklace \
  $BAM \
  $PROJECT_FOLDER/data/$PREFIX/${PREFIX}.gff \
  $PROJECT_FOLDER/data/assembled/counts/megahit \
  cov
done
```

#### Bin counts
```shell
Rscript $PROJECT_FOLDER/metagenomics_pipeline/scripts/cov_count.R "." "$P1.*\\.cov" "$PREFIX.countData"
```

#### Sub bin counts
Create tab files from cov files (originally for consistency with HirBin, but now to reduce processing in R, could be run in parallel)
```shell
for F in $P1*.cov; do
  O=$(sed 's/_.*_L/_L/' <<<$F|sed 's/_1\.cov/.tab/')
  awk -F"\t" '{sub("ID=","",$(NF-1));OUT=$1"_"$(NF-1)"_"$4"_"$5;print OUT,$(NF-1),$4,$5,$NF}' OFS="\t" $F > $O
done 
```

Parse the tab files to create countData matrix
```shell
Rscript subbin_parser.R reduced.txt *.tab $PREFIX.countData.sub_bins
```

