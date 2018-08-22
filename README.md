# Shotgun_metagenomics

## Preprocessing

## Assembly
metaspades and megahit are two decent options

#### metaspades
Metaspades can only run on paired reads (no option to use single and/or merged pairs, or multiple libraries)
```shell
for FR in $PROJECT_FOLDER/data/corrected/*_1.corrected.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p metaspades \
  $PROJECT_FOLDER/data/assembled \
  $FR \
  $RR  \
  $PREFIX \
  -k 21,33,55,77
done
```

#### megahit
Several options are recommended for soil samples  
--k-min=27 (or higher)  
--kmin-1pass  
--k-min 27 --k-step 10 --k-max 87 (127)  
```shell
# using merged and unmerged pairs
for FR in $PROJECT_FOLDER/data/merged/*_1.unmerged.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  MR=$(sed 's/_1\.un/\./' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit \
  $PROJECT_FOLDER/data/assembled \
  $PREFIX \
  -r $MR,$FR,$RR\
  --k-min=27 --k-step 10 --k-max 127
done
```

```shell
# using unmerged reads
for FR in $PROJECT_FOLDER/data/corrected/*_1.corrected.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  MR=$(sed 's/_1\.un/\./' <<< $FR)
  PREFIX=$(grep -Po 'N[0-9]+.' <<<$FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c assemble -p megahit \
  $PROJECT_FOLDER/data/assembled \
  $PREFIX \
 -1 $FR -2 $RR -r $MR \
 --k-min=27 --k-step 10 --k-max 127
done
```   

## Taxonomy assignment
Metakallisto/kracken/centrifuge?

### Binning


#### metakallisto

