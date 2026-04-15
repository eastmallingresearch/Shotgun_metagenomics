I have two pipeline kraken and kaiju for taxonomic profiling.  
Both are useful for soil metagenomics, and it is useful to run both

# Kraken pipeline



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
