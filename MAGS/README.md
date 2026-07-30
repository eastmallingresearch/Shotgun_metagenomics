# Set up environment

```shell
PROJECT_FOLDER=~/projects/MYPROJECT
ln -s ~/pipelines/metagenomics $PROJECT_FOLDER/metagenomics_pipeline

mkdir $PROJECT_FOLDER/data
mkdir $PROJECT_FOLDER/data/fastq
mkdir $PROJECT_FOLDER/data/filtered
mkdir $PROJECT_FOLDER/data/cluster

```

# Adapter removal and contaminant filtering
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

# Human contaminant removal
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
#
```shell
for F in "$PROJECT_FOLDER"/data/split_cleaned/*L2_1*.fq.gz; do
 S=$(echo $F|sed 's/_EKDN.*//')
 sbatch --mem=100M -p short "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_merge_fastq.sh $S 1 "$PROJECT_FOLDER"/data/cleaned
 sbatch --mem=100M -p short "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_merge_fastq.sh $S 2 "$PROJECT_FOLDER"/data/cleaned
done 
```

# Assembly

```shell
mkdir $PROJECT_FOLDER/data/runfiles/assemble
cd $PROJECT_FOLDER/data/runfiles/assemble
 
FILES=$PROJECT_FOLDER/data/cleaned
R1=$(find "$FILES" -type f -name '*_1.fq.gz.filtered.nonhuman.fq.gz' | paste -sd, - )
R2=$(sed 's/_1/_2/g' <<<$R1)

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
	 
### #Suggest filtering at min 2.5Kb contig length
```shell
seqkit seq -m 2500 $PROJECT_FOLDER/data/assembled/ASB/final.contigs.fa > $PROJECT_FOLDER/data/assembled/assemble.2.5kb.fa
```

# Possibly split into eukaryote/prokaryote
# first step is to run prodigal
```shell
for d in "$PROJECT_FOLDER"/data/assembled/*; do
  ASSEMBLY=$d/orchard.2.5kb.fa
  ORCHARD=$(basename $d)
  CORES=24
  sbatch --mem=60G -p long -c $CORES  "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_prodigal.sh \
    $ASSEMBLY \
    "$PROJECT_FOLDER"/data/prodigal \
    $ORCHARD \
    $CORES
done
 ```

# second step is to run whokaryote
```shell
for d in "$PROJECT_FOLDER"/data/assembled/*; do
  ASSEMBLY=$d/orchard.2.5kb.fa
  ORCHARD=$(basename $d)
  CORES=24
  sbatch --mem=60G -p long -c $CORES "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_whokaryote.sh \
    $ASSEMBLY \
    "$PROJECT_FOLDER"/data/whokaryote \
	"$PROJECT_FOLDER"/data/prodigal/$ORCHARD.prodigal.gff \
    $ORCHARD \
	2500 \
	$CORES
done
```
   
# ~4% eukaryote reads

#### Prokaryote workflow

# Align

# Using minimap - will need to index those orchrd files first
```shell
minimap2 -x sr -d prokaryotes.mmi prokaryotes.fasta
```

```shell
for R1 in "$PROJECT_FOLDER"/data/cleaned/C[SD]*_1.cleaned.fastq.gz; do
 R2=$(sed 's/_1/_2/g' <<<$R1)
 S=$(basename $R1|sed 's/_1.*//')
 ORCHARD=$(sed 's/_.*//' <<<$S)
 ORCHARD=$(sed 's/CS/CD/' <<<$ORCHARD)
 sbatch \
  --mem=60G \
  -p long \
  -c 24 \
  "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_minimap2.sh \
   $PROJECT_FOLDER/data/whokaryote/${ORCHARD}/prokaryotes.mmi \
   "$PROJECT_FOLDER"/data/aligned \
   $R1 \
   $R2 \
   $S \
   24
done 
```

# Make depth FILES
```shell
INP=$(find $FILES -type f -name "ASP*.bam")
jgi_summarize_bam_contig_depths --outputDepth ASP.depth.txt $INP

INP=$(find $FILES -type f -name "CAS*.bam")
jgi_summarize_bam_contig_depths --outputDepth CAS.depth.txt $INP

INP=$(find $FILES -type f -name "C[DS]*.bam")
jgi_summarize_bam_contig_depths --outputDepth CD.depth.txt $INP 

INP=$(find $FILES -type f -name "HER*.bam")
jgi_summarize_bam_contig_depths --outputDepth HER.depth.txt $INP

INP=$(find $FILES -type f -name "MOR*.bam")
jgi_summarize_bam_contig_depths --outputDepth MOR.depth.txt $INP

INP=$(find $FILES -type f -name "ASB*.bam")
jgi_summarize_bam_contig_depths --outputDepth ASB.depth.txt $INP
```

# BINNING

## metabat

```shell
for F in "$PROJECT_FOLDER"/data/depths/*.txt; do
 ORCHARD=$(basename $F|sed 's/\..*//')
 sbatch \
  --mem=36G \
  -p long \
  -c 16 \
  "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_metabat.sh \
   $PROJECT_FOLDER/data/whokaryote/${ORCHARD}/prokaryotes.fasta \
   $F \
   $ORCHARD \
   $PROJECT_FOLDER/data/binning/metabat \
   2500 \
   16 \
   1288
done 
```

## Semibin
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
    $PROJECT_FOLDER/data/binning/semibin \
    64
done 
```

## BIN REFINEMENT
Most bin refiners discard fungal genomes - probably best to filter bins first into pro/eukaryote

### MAGScot prokaryote refiners
generate tsv files

```shell
for d in $PROJECT_FOLDER/data/binning/metabat/*; do
  for f in $d/*.fa; do
    bin=$(basename "$f" .fa)
    grep ">" "$f" | sed 's/^>//' | awk -v b="$bin" '{print b"\t"$1"\tmetabat2"}'
  done > $d.metabat2.contigs_to_bin.tsv
done
```

```shell
for d in $PROJECT_FOLDER/data/binning/semibin/*; do
  for f in $d/semibin2.bin/output_bins/*.fa.gz; do
    bin=$(basename "$f" .fa.gz)
    zcat "$f" | grep ">" | sed 's/^>//' | awk -v b="$bin" '{print b"\t"$1"\tsemibin2"}'
  done > $d.semibin2.contigs_to_bin.tsv
done
```

```shell
cd $PROJECT_FOLDER/data/binning
```

```shell
for d in $PROJECT_FOLDER/data/binning/semibin/*tsv; do
  s=$(sed 's/semibin/metabat/g' <<<$d) 
  o=$(basename "$d" .semibin2.contigs_to_bin.tsv)
  cat $d $s > $o.contigs_to_bin.tsv
done
```

```shell
pixi run --manifest-path ~/envs/magscot/pixi.toml \
  Rscript $MAGScoT_folder/MAGScoT.R \
  -i orchard1.contigs_to_bin.tsv \
  --hmm orchard1.hmm \
  -o orchard1_magscot \
  -p bac120+ar53 \
  --min_markers 25 \
  --min_sharing 0.8
```

```shell
for d in "$PROJECT_FOLDER"/data/whokaryote/*; do
  ASSEMBLY=$d/prokaryotes.fasta
  ORCHARD=$(basename $d)
  CORES=24
  sbatch --mem=60G -p long -c $CORES "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_MAGScot.sh \
    $ASSEMBLY \
    ~/apps/MAGScoT/ \
	"$PROJECT_FOLDER"/data/binning/${ORCHARD}.contigs_to_bin.tsv \
	"$PROJECT_FOLDER"/data/binning/MAGScot \
    $ORCHARD \
	$CORES
done
```

### Bin checking (checkm2)

```shell
for d in "$PROJECT_FOLDER"/data/binning/MAGScot/*; do
  BINS=$d/bins
  ORCHARD=$(basename $d)
  CORES=24
  sbatch --mem=60G -p long -c $CORES "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_checkm2.sh \
    $BINS \
    ~/apps/checkm2/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
	"$PROJECT_FOLDER"/data/binning/MAGScot \
	$ORCHARD \
	$CORES
done
```

extract combined bin set
need to decide on what bins to retain based contamination and completeness
duel filter, discard <50 complete and > 30% contamination
some simple R code is probably best

### recover files (may need renaming..)
```shell
for d in "$PROJECT_FOLDER"/data/binning/MAGScot/*; do
  ORCHARD=$(basename $d)
  cd $ORCHARD/bins
  rename "s/^/${ORCHARD}_/" *.fa
  cd ../..
 done
```

#### filter and keep dodgy, but checkable (in R)

```R
dodge <- dat[Contamination>10,]
# MAG set
dat <-dat[Contamination <=10,]
dim(dat)

mkdir "$PROJECT_FOLDER"/data/binning/MAGScot/combined
mkdir "$PROJECT_FOLDER"/data/binning/MAGScot/dodgy
```

### GUNC
```shell
gunc download_db $APPS/gunc_db --db progenomes_2.1
```

```shell
CORES=24
sbatch --mem=60G -p long -c $CORES "$PROJECT_FOLDER"/metagenomics_pipeline/scripts/sub_GUNC.sh \
    "$PROJECT_FOLDER"/data/binning/combined/good \
    ~/apps/GUNC/gunc_db_progenomes2.1.dmnd \
```
	"$PROJECT_FOLDER"/data/binning/combined/GUNC \
	$CORES
