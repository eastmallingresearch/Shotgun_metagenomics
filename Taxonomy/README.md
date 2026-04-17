I have two pipeline kraken and kaiju for taxonomic profiling.  
Both are useful for soil metagenomics, and it is useful to run both

# Kraken pipeline
The CropDiversit HPC maintains several version of kraken databases  
Will check how to use these

## Database options

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

# prebuilt options
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20260226.tar.gz

```

## Run Kraken

```shell
for FR in $PROJECT_FOLDER/data/cleaned/*_1*.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  S=$(sed 's/\(.*\/\)\(.*_1\)\(\..*\)/\2/' <<< $FR)
  sbatch --mem=250000 -p himem -c 20 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_kraken.sh \
  $PROJECT_FOLDER/PATHTOKRAKENDB \
  $FR \
  $RR \
  ${S}.kraken.out \
  ${S}.kraken.report.out \
  $PROJECT_FOLDER/data/taxonomy/kraken/ \
  20
done
```
## Braken
Braken is used to estimate abundances from kraken reports.

```shell
for KR in $PROJECT_FOLDER/data/taxonomy/kraken/*.report.out; do
  S=$(sed 's/\(.*\/\)\(.*_1\)\(\..*\)/\2/' <<< $KR)
  sbatch --mem=20000 -p short $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_bracken.sh \
  $PROJECT_FOLDER/PATHTOKRAKENDB \
  $KR \
  200 \ # read length
  ${S} \
  $PROJECT_FOLDER/data/taxonomy/bracken/ 
done

```

# Kaiju pipeline
Kaiju needs a database - for taxonomy probably best to stick to the nr_euk database. Other, or custom options are available 



## Download and setup database  
```shell
kaiju-makedb -s nr_euk 
```

## Run kaiju
```shell
for FR in $PROJECT_FOLDER/data/cleaned/*_1*.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  S=$(sed 's/\(.*\/\)\(.*_1\)\(\..*\)/\2/' <<< $FR)
  sbatch --mem=120000 -p medium -c 20 $PROJECT_FOLDER/metagenomics_pipeline/scripts/   slurm/sub_kaiju.sh \
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
