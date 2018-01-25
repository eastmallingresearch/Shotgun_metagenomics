# Shotgun_metagenomics

## Preprocessing

```shell
# add project folders
$PROJECT_FOLDER=~/projects/myproject
ln -s ~/pipelines/metagenomics $PROJECT_FOLDER/metagenomics_pipeline

mkdir $PROJECT_FOLDER/data
mkdir $PROJECT_FOLDER/data/fastq
mkdir $PROJECT_FOLDER/data/trimmed
mkdir $PROJECT_FOLDER/data/filtered
mkdir $PROJECT_FOLDER/data/normalised
mkdir $PROJECT_FOLDER/data/merged
```


### Adapter trimming (trimmomatic)
Edit the *FR* and the sed to the required format
```shell
for FR in $PROJECT_FOLDER/data/fastq/*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c trim \
  $FR \
  $RR \
  $PROJECT_FOLDER/data/trimmed \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/adapters/truseq.fa \
  4
done
```
### Synthetic construct/contaminant filtering (BBduc/bowtie)
```shell
for FR in $PROJECT_FOLDER/data/trimmed/*_1.fq.gz.trimmed.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c filter -p bbduk \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/contaminants/phix_174.fa \
  $PROJECT_FOLDER/data/cleaned \
  $FR \
  $RR 
done
```

### Human contaminant removal (BBMap)

### Normalization (BBNorm)

### paired read merge (BBMerge/USEARCH)

### Quality check (BLASTn)
