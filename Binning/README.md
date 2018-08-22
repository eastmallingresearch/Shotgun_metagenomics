## functional binning with HirBin
 I've had to hack some of the HirBin scripts (specifically clusterbinstosubbins.py) as it doesn't work in current format
 also it uses usearch for clustering, while this is good the free 32bit version will almost certainly run out of memory for any sort of 
 soil metagenome assembly. I will probably change this to use vsearch which will give output the same as usearch

 HirBin does an hmm alignment of an assembly to a protein domain database 
 this will take a looong time unless the assembly and preferbly the hmm database is divided into chunks

 it has three/four steps

### Annotate
annotate uses functionalAnnotaion.py, but splits input file into 20,000 droid chunks for running on cluster (25 concurrent jobs)
functionalAnnotation.py -m METADATA_FILE -db DATABASE_FILE -e EVALUE_CUTOFF -n N -p MAX_ACCEPTABLE_OVERLAP

```shell
$PROJECT_FOLDER/metagenomics_pipeline/scripts/fun_bin.sh \
 1 $PROJECT_FOLDER/data/assembled/megahit/$PREFIX \
 $PREFIX.contigs.fa \
 ~/pipelines/common/resources/pfam/Pfam-A.hmm \
 -e 1e-03
```
#### concatenate annotate output
```shell
find -type f -name X.gff|head -n1|xargs -I% head -n1 % >$PREFIX.gff
find -type f -name X.gff|xargs -I% grep -v "##" % >>$PREFIX.gff
find -type f -name X.pep|xargs -I% cat % >$PREFIX.pep
find -type f -name X.hmmout|head -n1|xargs -I% head -n3 % >$PREFIX.hmmout   
find -type f -name X.hmmout|xargs -I% grep -v "#" % >>$PREFIX.hmmout
find -type f -name X.hmmout|head -n1|xargs -I% tail -n10 % >>$PREFIX.hmmout

grep -v "#" $PREFIX.hmmout|awk -F" " '($21~/^[0-9]+$/) && ($20~/^[0-9]+$/) {print $4,$1,$20,$21,"+",$7}' OFS="\t" > $PREFIX.hmm.cut
awk -F"\t" '{print $1}' $PREFIX.hmm.cut|sort|uniq > $PREFIX.domains # this is better method as some domians may not be present in gff due to filtering
# cut -f9 $PREFIX.gff|sort|uniq|sed 's/ID=//'|tail -n +2 > $PREFIX.2.domains # the tail bit gets rid of the first line of output an d the grep removes errors in the output
```
#### Get gff with different MAX_ACCEPTABLE_OVERLAP - example below will produce gff with all overlapping features
```shell
con_coor.py -p 1 -o $PREFIX.2.gff -d $PREFIX.pep -m $PREFIX.hmmout
```
## mapping
mapping is not implemented very well in HirBin, will do this seperately with bbmap
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
  path=$PROJECT_FOLDER/data/assembled/megahit/$PREFIX/ 
  usemodulo=T 
done
```

## Count coverage
bedtools code is inefficient at getting over-lapping counts (if min overlap is set to 1)
I've written something in perl which is way less memory hungry and takes about a millionth of the time to run
output is not a cov file but just counts per domain - not certain the sub-binning is worth while (could modify bam_count to return a cov/tab file to implement this step)
takes about ten minutes on a single core to run, could easily get it to produce a cov file
bam_scaffold_count.pl will output a cov file rather than counts per domain, it is memory hungry ~10G for large (>2G) gff files  
e.g.
```shell
samtools view bam_file|~/pipelines/metagenomics/scripts/bam_scaffold_count.pl $PREFIX.gff > bam_counts.txt
samtools view bam_file|~/pipelines/metagenomics/scripts/bam_scaffold_count.pl $PREFIX.gff cov> bam_file.cov
```
Get counts
```shell
for BAM in $PROJECT_FOLDER/data/assembled/aligned/megahit/$P1*.bam; do
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c coverage -p bam_count \
  blacklace[01][0-9].blacklace \
  $BAM \
  $PROJECT_FOLDER/data/assembled/megahit/$PREFIX/${PREFIX}.gff \
  $PROJECT_FOLDER/data/assembled/counts/megahit \
  cov
done
```

Convert per sample count to count matrix
```shell
Rscript $PROJECT_FOLDER/metagenomics_pipeline/scripts/cov_count.R "." "$P1.*\\.cov" "$PREFIX.countData"
```
## Sub binning 
I've hacked around with a few of the HirBin settings - for speed mostly and for consistency (or the lack of) in domain names
will require a tab (converted cov) file from bam_scaffold_count.pl  
e.g. 
```shell
awk -F"\t" '{sub("ID=","|",$(NF-1));OUT=$1$(NF-1)":"$4":"$5":"$7;print OUT,$NF}' OFS="\t" x.cov > x.tab
```
All samples
```shell
for F in $P1*.cov; do
  O=$(sed 's/_.*_L/_L/' <<<$F|sed 's/_1\.cov/.tab/')
  awk -F"\t" '{sub("ID=","|",$(NF-1));OUT=$1$(NF-1)":"$4":"$5":"$7;print OUT,$NF}' OFS="\t" $F > $O
done 
```

#### create the required metadata file
```shell
echo -e \
"Name\tGroup\tReference\tAnnotation\tCounts\Domain\n"\
"$PREFIX\tSTATUS\t$PREFIX.pep\t$PREFIX.hmm.cut\tEMPTY\t$PREFIX.domains" > $PREFIX.metadata.txt
```

### Clustering and counting
The extractSequences module of clusterBinsToSubbins is inpractical for large pep files - the R script subbin_fasta_extractor.R should be used in it's palace
```shell
Rscript subbin_fasta_extractor.R $PREFIX.domains $PREFIX.pep $PREFIX_hirbin_output
# clusterBinsToSubbins.py -m metadata.txt -id 0.7 --onlyClustering -f -o  $PREFIX_hirbin_output # this will create the sub bins - use subbin_fasta_extractor in preference to this 
clusterBinsToSubbins.py -m metadata.txt -id 0.7 --reClustering --onlyClustering -f -o  $PREFIX_hirbin_output# clustering without sub bin extraction (no parsing)
clusterBinsToSubbins.py -m metadata.txt -id 0.95 --reClustering -f -o  $PREFIX_hirbin_output# recluster at a different identity plus parsing
clusterBinsToSubbins.py -m metadata.txt -id 0.7 --onlyParsing -f -o  $PREFIX_hirbin_output# this will make count files for $PREFIX.tab to the bins and sub bins 
```