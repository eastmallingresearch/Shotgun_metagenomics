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
### Adapter removal/quality filtering and contaminant filtering
BBTools has good options for doing all of this. 

I've merged adapter removal, phix filtering and rRNA removal into a single operation using BBDuk (though it has to run multiple times, rather than a single passthrough). To modify any settings will require editing the mega_duk.sh script. Alternatively the three operations can be run seperately using bbduk (PIPELINE.sc -c bbduk). 

rRNA removal is probably best left for metatranscriptomic data. The final false in the below script prevents this step from running

#### Adapter/phix/rRNA removal
Runs all three of the options in "Filtering full options" shown at bottom
```shell
for FR in $PROJECT_FOLDER/data/trimmed/*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c MEGAFILT \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/adapters/truseq.fa \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/contaminants/phix_174.fa \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/contaminants/ribokmers.fa.gz \
  $PROJECT_FOLDER/data/filtered \
  $FR \
  $RR \
  false
done  
```
bbduk command line arguments used:  
adapter removal forward; ktrim=l k=23 mink=11 hdist=1 tpe tbo t=10
adapter removal reverse; ktrim=r k=23 mink=11 hdist=1 tpe tbo t=10
phix filtering; k=31 hdist=1 t=4
rRNA filtering; k=31 t=4 

#### Human contaminant removal (BBMap)
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
#### Normalization and error correction (BBNorm)
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
#### Paired read merge (BBMerge)
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

#### Partition reads (Clumpify)
```shell
for UR in $PROJECT_FOLDER/data/merged/*_1.fq.gz.trimmed.fq.gz.filtered.fq.gz.cleaned.fq.gz.corrected.fq.gz.unmerged.fq.gz; do
  MR=$(sed 's/merged/unmerged/' <<< $UR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c partition -p clumpify \
  $PROJECT_FOLDER/data/clumped \
  $MR \
  $UR 
done
```

#### rename files (could have implemented in the jobs above)
```shell
find $PROJECT_FOLDER/data -type f -name *.fq.gz|rename 's/(.*_[12]).*(\.[a-zA-Z]+\.fq\.gz$)/$1$2/'
```
