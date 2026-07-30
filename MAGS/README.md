# Common pipeline

## Set up environment

```shell
PROJECT_FOLDER=~/projects/MYPROJECT
ln -s ~/pipelines/metagenomics $PROJECT_FOLDER/metagenomics_pipeline

mkdir $PROJECT_FOLDER/data
mkdir $PROJECT_FOLDER/data/fastq
mkdir $PROJECT_FOLDER/data/filtered
mkdir $PROJECT_FOLDER/data/cluster
mkdir $PROJECT_FOLDER/data/assembled

```

## Adapter removal and contaminant filtering
```shell
for FR in $PROJECT_FOLDER/data/fastq/*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  sbatch --mem=20000 -p medium -c 10 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/mega_duk.sh \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/adapters/truseq.fa \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/contaminants/phix_174.fa \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/contaminants/ribokmers.fa.gz \
  $PROJECT_FOLDER/data/filtered \
  $FR \
  $RR \
  false
done 
```

## Human contaminant removal
```shell
for FR in "$PROJECT_FOLDER"/data/filtered/*_1*.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< "$FR")
  sbatch \
    --mem=40000 \
    -p medium \
    -c 20 \
    "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/slurm/sub_bwa.sh \
    "$PROJECT_FOLDER"/metagenomics_pipeline/common/resources/contaminants/bwa_human_ref/hs1.fa \
    "$PROJECT_FOLDER"/data/cleaned \
    "$FR" \
    "$RR"
done
```
## Assembly

```shell
mkdir $PROJECT_FOLDER/data/cluster/assemble
cd $PROJECT_FOLDER/data/cluster/assemble

sbatch \
  --mem=500G \
  -p himem \
  -c 64 \
  "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_megahit.sh \
    $R1 \
	$R2 \
    64 \
	"$PROJECT_FOLDER"/data/assembled
```

# MAG (prokaryote) pipeline

This pipeline was built to produce a combined set of MAGs from soil derived from six assemblies of six apple orchards.

```shell
ORCHARD=CAS
```

The assembly is likely to produce (especially for complex soils) a huge number of short <1Kb contigs. 
For producing MAGs it is advisable to filter out shorter contigs - a suggested minimum size (enforced by later steps) is 2.5Kb

## Filtering to min 2.5Kb contig length
```shell
mkdir -P $PROJECT_FOLDER/data/MAGs/filtered_assembly
seqkit seq -m 2500 $PROJECT_FOLDER/data/assembled/$ORCHRD/final.contigs.fa > $PROJECT_FOLDER/data/MAGs/filtered_assembly/$ORCHARD/orchard.2.5kb.fa
```

The pipeline can only handle prokaryotic contigs (this is a limit of many of the chosen-best in class-tools. 
TODO create eukaryote pipeline

## Split contigs into eukaryote/prokaryote

### Run prodigal
```shell
for d in "$PROJECT_FOLDER"/data/MAGs/filtered_assembly/*; do
  ASSEMBLY=$d/orchard.2.5kb.fa
  ORCHARD=$(basename $d)
  CORES=24
  sbatch --mem=60G -p long -c $CORES  "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_prodigal.sh \
    $ASSEMBLY \
    "$PROJECT_FOLDER"/data/MAGs/prodigal \
    $ORCHARD \
    $CORES
done
 ```

### Run whokaryote

whokaryote classifies each contig into prokaryote, eukaryote or unknown, and outputs as fasta file for each 

```shell
for d in "$PROJECT_FOLDER"/data/MAGs/filtered_assembly/*; do
  ASSEMBLY=$d/orchard.2.5kb.fa
  ORCHARD=$(basename $d)
  CORES=24
  sbatch --mem=60G -p long -c $CORES "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_whokaryote.sh \
    $ASSEMBLY \
    "$PROJECT_FOLDER"/data/MAGs/whokaryote \
	"$PROJECT_FOLDER"/data/MAGs/prodigal/$ORCHARD.prodigal.gff \
    $ORCHARD \
	2500 \
	$CORES
done
```

   
## Prokaryote pipeline

### Align

Aligning using minimap - will need to index each of the prokaryote contig files

#### Indexing

```shell
minimap2 -x sr -d "$PROJECT_FOLDER"/data/MAGs/whokaryote/$ORCHARD/prokaryotes.mmi "$PROJECT_FOLDER"/data/MAGs/whokaryote/$ORCHARD/prokaryotes.fasta
```

#### Run alignment

```shell
for R1 in "$PROJECT_FOLDER"/data/cleaned/*_1.cleaned.fastq.gz; do
 R2=$(sed 's/_1/_2/g' <<<$R1)
 S=$(basename $R1|sed 's/_1.*//')
 ORCHARD=$(sed 's/_.*//' <<<$S)
 ORCHARD=$(sed 's/CS/CD/' <<<$ORCHARD)
 sbatch \
  --mem=60G \
  -p long \
  -c 24 \
  "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_minimap2.sh \
   $PROJECT_FOLDER/data/MAGs/whokaryote/${ORCHARD}/prokaryotes.mmi \
   "$PROJECT_FOLDER"/data/MAGs/aligned \
   $R1 \
   $R2 \
   $S \
   24
done 
```

### Make depth FILES

Need a single concatenated file for all samples in an assembly

```shell
INP=$(find $PROJECT_FOLDER/data/MAGs $FILES -type f -name "${ORCHARD}*.bam")
jgi_summarize_bam_contig_depths --outputDepth $PROJECT_FOLDER/data/MAGs/depths/$ORCHARD.depth.txt $INP
```

### BINNING

#### metabat

```shell
for F in "$PROJECT_FOLDER"/data/MAGs/depths/*.txt; do
 ORCHARD=$(basename $F|sed 's/\..*//')
 sbatch \
  --mem=36G \
  -p long \
  -c 16 \
  "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_metabat.sh \
   $PROJECT_FOLDER/data/MAGs/whokaryote/${ORCHARD}/prokaryotes.fasta \
   $F \
   $ORCHARD \
   $PROJECT_FOLDER/data/MAGs/binning/metabat \
   2500 \
   16 \
   1288
done 
```

#### Semibin
```shell
for F in "$PROJECT_FOLDER"/data/depths/*.txt; do
  SAMPLE=$(basename $F .depth.txt)
  echo sbatch \
    --mem=80G \
    -p gpu \
    -c 64 \
    --gpus=1 \
    "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_semibin.sh \
    $PROJECT_FOLDER/data/whokaryote/$SAMPLE/prokaryotes.fasta \
    $SAMPLE \
    $PROJECT_FOLDER/data/aligned/${SAMPLE}'\*.sorted.bam' \
    $PROJECT_FOLDER/data/MAGs/binning/semibin \
    64
done 
```


#### MAGScot prokaryote bin refiner

Most bin refiners discard fungal genomes - hence why we split earlier
TODO: find refiner which can handle eukaryote bins


##### Generate metabat tsv files

```shell
for d in $PROJECT_FOLDER/data/MAGs/binning/metabat/*; do
  for f in $d/*.fa; do
    bin=$(basename "$f" .fa)
    grep ">" "$f" | sed 's/^>//' | awk -v b="$bin" '{print b"\t"$1"\tmetabat2"}'
  done > $d.metabat2.contigs_to_bin.tsv
done
```

##### semibin tsv files

```shell
for d in $PROJECT_FOLDER/data/MAGs/binning/semibin/*; do
  for f in $d/semibin2.bin/output_bins/*.fa.gz; do
    bin=$(basename "$f" .fa.gz)
    zcat "$f" | grep ">" | sed 's/^>//' | awk -v b="$bin" '{print b"\t"$1"\tsemibin2"}'
  done > $d.semibin2.contigs_to_bin.tsv
done
```

##### concatenate tsv files per orchard

```shell
cd $PROJECT_FOLDER/data/MAGs/binning
for d in $PROJECT_FOLDER/data/MAGs/binning/semibin/*tsv; do
  s=$(sed 's/semibin/metabat/g' <<<$d) 
  o=$(basename "$d" .semibin2.contigs_to_bin.tsv)
  cat $d $s > $o.contigs_to_bin.tsv
done
```

#### Run MAGScot

```shell
for d in "$PROJECT_FOLDER"/data/MAGs/whokaryote/*; do
  ASSEMBLY=$d/prokaryotes.fasta
  ORCHARD=$(basename $d)
  CORES=24
  sbatch --mem=60G -p long -c $CORES "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_MAGScot.sh \
    $ASSEMBLY \
    ~/apps/MAGScoT/ \
	"$PROJECT_FOLDER"/data/MAGs/${ORCHARD}.contigs_to_bin.tsv \
	"$PROJECT_FOLDER"/data/MAGs/MAGScot \
    $ORCHARD \
	$CORES
done
```

#### Bin checking and filtering

All of these steps are useful, but in a complex soil metagenome they may not work as intended.  

##### Checkm2
```shell
for d in "$PROJECT_FOLDER"/data/MAGs/MAGScot/*; do
  BINS=$d/bins
  ORCHARD=$(basename $d)
  CORES=24
  sbatch --mem=60G -p long -c $CORES "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_checkm2.sh \
    $BINS \
    ~/apps/checkm2/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
	"$PROJECT_FOLDER"/data/MAGs/MAGScot \
	$ORCHARD \
	$CORES
done
```

Extract combined bin set
Need to decide on what bins to retain based contamination and completeness  
Standard is duel filter, discard <50 complete and > 30% contamination  some simple R code could do this   
But, this needs checking depending on the complexity of the the biome

##### Rename MAGScot bins
```shell
for d in "$PROJECT_FOLDER"/data/MAGs/MAGScot/*; do
  ORCHARD=$(basename $d)
  cd $ORCHARD/bins
  rename "s/^/${ORCHARD}_/" *.fa
  cd ../..
 done
```

#####  Rename checkm2 files

Also copies the per orchard files to the checkm2_files folder

```shell
for d in "$PROJECT_FOLDER"/data/MAGs/MAGScot/*; do
  ORCHARD=$(basename $d)
  cd $ORCHARD/checkm2
	rename "s/^/${ORCHARD}_/" *.tsv
  cd ../..
 done

mkdir "$PROJECT_FOLDER"/data/MAGs/MAGScot/checkm2_files
  
cd "$PROJECT_FOLDER"/data/MAGs/MAGScot/checkm2_files
find "$PROJECT_FOLDER"/data/MAGs/MAGScot/ -type f -name *quality_report.tsv|xargs -I % cp % .

sed -i -e 's/^[^N]/CAS_m/' ${ORCHARD}_quality_report.tsv

```


##### filter and keep good bins

```shell
mkdir -P "$PROJECT_FOLDER"/data/MAGs/MAGScot/combined/good
mkdir "$PROJECT_FOLDER"/data/MAGs/MAGScot/combined/dodgy

cd "$PROJECT_FOLDER"/data/MAGs/MAGScot
R
```

```R
library(data.table)

# get checkm files
files <- list.files(
  path = "./checkm2_files",
  pattern = ".*quality_report\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

# load files
qq <- lapply(files, fread)

# name the list by  basename
names(qq) <- basename(files)

# add an "orchard" column to each data.table
lapply(names(qq),function(i) qq[[i]][,orchard:=gsub("_.*","",i)])


# combine checkm2 files
dat <- rbindlist(qq)

# select out dodgy bins1
dodge <- dat[Contamination>10,]

# select good bins
good <-dat[Contamination <=10,]

# how many bins in each category
c(All=nrow(dat),Good=nrow(good),Dodgy=nrow(dodge))

# copy files 
file.copy(paste0(good$orchard,"/bins/",good$Names,".fa"),"combined/good/")
file.copy(paste0(dodge$orchard,"/bins/",dodge$Names,".fa"),"combined/good/")

# write out checkm2 combined good results
fwrite(good,"combined_filtered_bins.txt")
```

##### GUNC filtering- checking (this is not so useful perhaps in complex soil metagenomes)

```shell
gunc download_db $APPS/gunc_db --db progenomes_2.1
```

```shell
CORES=24
sbatch --mem=60G -p long -c $CORES "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_GUNC.sh \
    "$PROJECT_FOLDER"/data/MAGs/MAGScot/combined/good \
    ~/apps/GUNC/gunc_db_progenomes2.1.dmnd \
	"$PROJECT_FOLDER"/data/MAGScot/combined/GUNC \
	$CORES
```

#### Dereplication with dRep	


##### dRep pre-run steps

```shell

cd "$PROJECT_FOLDER"/data/MAGs/MAGScot/checkm2_files

echo "genome,completeness,contamination" > drep_genomeInfo.csv
tail -n +2 *.tsv | awk -F'\t' '{print $1".fa,"$2","$3}' >> drep_genomeInfo.csv
```

##### Run dRep
```shell
CORES=24
sbatch --mem=60G -p long -c $CORES "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_dRep.sh \
  "$PROJECT_FOLDER"/data/MAGs/combined/good \
  "$PROJECT_FOLDER"/data/MAGs/checkm2_files/drep_genomeInfo.csv \
  "$PROJECT_FOLDER"/data/MAGs/dRep \
  $CORES
```

The final set of MAGs are in dRep/drep_output/dereplicated_genomes.  
 
#### Abundance

##### Build minimap index	

First rename all the contigs per orchard to ensure they're unique (essential) and traceable (useful - maybe), otherwise some of the downstream steps fail. Then concatenate them 


```shell
# check unique
cd "$PROJECT_FOLDER"/data/MAGs/dRep/drep_output/dereplicated_genomes/
for f in *.fa; do
  ORCHARD=$(basename $f|sed 's/_.*//')
  sed -i -e "s/^>/>${ORCHARD}_/" $f
done

# cat MAGs
cat "$PROJECT_FOLDER"/data/MAGs/dRep/drep_output/dereplicated_genomes/*.fa > "$PROJECT_FOLDER"/data/MAGs/dRep/genome_catalog.fa

#useful to archive the MAGs, might as well do it now
pigz -c "$PROJECT_FOLDER"/data/MAGs/dRep/genome_catalog.fa > "$PROJECT_FOLDER/data/MAGs/mags.tar.gz
```

Run indexing

```shell
minimap2 -x sr -I 10G -d "$PROJECT_FOLDER"/data/MAGs/dRep/genome_catalog.mmi "$PROJECT_FOLDER"/data/MAGs/dRep/genome_catalog.fa
```

##### Map reads back to new index
```shell
for R1 in "$PROJECT_FOLDER"/data/cleaned/*_1.cleaned.fastq.gz; do
 R2=$(sed 's/_1/_2/g' <<<$R1)
 S=$(basename $R1|sed 's/_1.*//')
 ORCHARD=$(sed 's/_.*//' <<<$S)
 ORCHARD=$(sed 's/CS/CD/' <<<$ORCHARD)
 sbatch \
  --mem=60G \
  -p long \
  -c 24 \
  "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_minimap2.sh \
   "$PROJECT_FOLDER"/data/MAGs/dRep/genome_catalog.mmi \
   "$PROJECT_FOLDER"/data/MAGs/mag_aligned \
   $R1 \
   $R2 \
   $S \
   24
done 
```

##### Calculate coverage

```shell
CORES=24
sbatch --mem=60G -p himem -c $CORES "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_coverm.sh \
  "$PROJECT_FOLDER"/data/MAGs/mag_aligned \
  "$PROJECT_FOLDER"/data/MAGs/dRep \
  "$PROJECT_FOLDER"/data/MAGs/abundance \
  $CORES
```  
  
####Taxonomy

```shell
CORES=24
sbatch --mem=150G -p himem -c $CORES "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_gtdbtk.sh \
  "$PROJECT_FOLDER"/data/MAGs/drep_output/dereplicated_genomes \
  "$PROJECT_FOLDER"/data/MAGs/taxonomy \
  $CORES
```


