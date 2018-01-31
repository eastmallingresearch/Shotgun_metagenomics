# Shotgun_metagenomics

## Preprocessing
The workflow should include at the least adapter trimming and filtering for phix/contamination. Normalisation, error correction and merging are dependent on the data and/or the assebley pipeline. Trimmomatic can also trim for quality if the data is of poor quality.

The bbtools preprocessing pipeline has a number of good options for many of these tasks.

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


### Adapter trimming (trimmomatic - extra command line options, e.g. quality trimming, will be passed to trimmomatic)
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
### Synthetic construct/contaminant filtering (BBduc)
```shell
for FR in $PROJECT_FOLDER/data/trimmed/*_1.fq.gz.trimmed.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c filter -p bbduk \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/contaminants/phix_174.fa \
  $PROJECT_FOLDER/data/filtered \
  $FR \
  $RR \
  k=31 \
  hdist=1
done
```

### Human contaminant removal (BBMap)
```shell
for FR in $PROJECT_FOLDER/data/filtered/*_1.fq.gz.trimmed.fq.gz.filtered.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c filter -p bbmap \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/contaminants/bbmap_human \
  $PROJECT_FOLDER/data/cleaned \
  $FR \
  $RR \
  minid=0.95 \
  maxindel=3 \
  bwr=0.16 \
  bw=12 \
  quickmatch \
  fast \
  minhits=2 \
  t=8
done
```
### Normalization and error correction (BBNorm)
```shell
for FR in $PROJECT_FOLDER/data/cleaned/*_1.fq.gz.trimmed.fq.gz.filtered.fq.gz.cleaned.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c normalise -p bbnorm \
  $PROJECT_FOLDER/data/corrected \
  $FR \
  $RR  \
  target=100 \
  min=5 \
  ecc=t \
  passes=1 \
  bits=16 prefilter
done
```
### Paired read merge (BBMerge)
```shell
for FR in $PROJECT_FOLDER/data/corrected/*_1.fq.gz.trimmed.fq.gz.filtered.fq.gz.cleaned.fq.gz.corrected.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c merge -p bbmerge-auto \
  $PROJECT_FOLDER/data/merged \
  $FR \
  $RR  \
  rem k=62 \
  extend2=50 \
  t=12 \
  vstrict
done
```

# output can be catted to form a 
```shell

```


## Assembly
metaspades and megahit are two decent options

### metaspades
Metaspades can only run on paired reads (no option to use single and/or merged pairs, or multiple libraries)
```shell
for FR in $PROJECT_FOLDER/data/corrected/*_1.fq.gz.trimmed.fq.gz.filtered.fq.gz.cleaned.fq.gz.corrected.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p metaspades \
  $PROJECT_FOLDER/data/assembled \
  $FR \
  $RR  \
  $PREFIX \
  -k 21,33,55,77
```


