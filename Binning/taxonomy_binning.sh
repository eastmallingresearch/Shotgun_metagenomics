# binning using metabat
# pipeline needs sorted bam files

# sort bam files
for f in *.bam; do
 PREFIX=$(echo $f|sed -e 's/\..*//')
 sbatch --mem-per-cpu 2000M -c 10 \
 ~/pipelines/metagenomics/scripts/slurm/sub_bam_sort.sh \
 10 /data/data2/scratch2/deakig/Oak/sorted $PREFIX /data/data2/scratch2/deakig/Oak/to_sort/$f
done

# get list of bam files for each assembly
A=$(for f in ../sorted/A*; do echo $f; done|tr  '\n' ' ')
G=$(for f in ../sorted/G*; do echo $f; done|tr  '\n' ' ')
L=$(for f in ../sorted/L*; do echo $f; done|tr  '\n' ' ')
W=$(for f in ../sorted/W*; do echo $f; done|tr  '\n' ' ')


# run metabat
# runMetaBat.sh -i assembly.fa.gz --unbinned -o assembly_name -m 1500 -x 0 --minCVSum 0.5 bam_files 

runMetaBat.sh  --unbinned -m 1500 -x 0 --minCVSum 0.5 \
~/projects/Oak_decline/metagenomics/data/assembled/megahit/ATTINGHAM/ATTINGHAM.contigs.fa $A &

runMetaBat.sh  --unbinned -m 1500 -x 0 --minCVSum 0.5 \
~/projects/Oak_decline/metagenomics/data/assembled/megahit/GTMONK/GTMONK.contigs.fa $G &

runMetaBat.sh  --unbinned -m 1500 -x 0 --minCVSum 0.5 \
~/projects/Oak_decline/metagenomics/data/assembled/megahit/LANGDALE/LANGDALE.contigs.fa $L &

runMetaBat.sh  --unbinned -m 1500 -x 0 --minCVSum 0.5 \
~/projects/Oak_decline/metagenomics/data/assembled/megahit/WINDING/WINDING.contigs.fa $W &

jgi_summarize_bam_contig_depths --outputDepth GTMONKdepth.txt ../sorted/G*
metabat2  --unbinned -m 1500 -x 0 --minCVSum 0.5 \
-i ~/projects/Oak_decline/metagenomics/data/assembled/megahit/GTMONK/GTMONK.contigs.fa -a GTMONKdepth.txt -o bins_dir/bin 


# Taxonomy assignment with kaiju
kaiju-makedb -s nr_euk 
# kaiju-makedb -s nr_euk --index-only -t 12 --no-download
kaiju -t nodes.dmp -f kaiju_db.fmi -i fortaxa.fa -o kaiju.out -z 10 -v

# bad method as has to read in index each time - slow
for f in bin.1.fa; do
  kaiju -t ../../../kaiju/nodes.dmp -f ../../../kaiju/nr_euk/kaiju_db_nr_euk.fmi -i $f -o $f.kaiju.out -z 10 -v
done

cat bin*.fa >ATTINGHAM.bins.fa
cat bin*.fa >GTMONK.bins.fa
cat bin*.fa >LANGDALE.bins.fa
cat bin*.fa >WINDING.bins.fa

kaiju -t ../../../kaiju/nodes.dmp -f ../../../kaiju/nr_euk/kaiju_db_nr_euk.fmi -i ATTINGHAM.bins.fa -o ATTINGHAM.kaiju.out -z 20 -v
kaiju -t ../../../kaiju/nodes.dmp -f ../../../kaiju/nr_euk/kaiju_db_nr_euk.fmi -i GTMONK.bins.fa -o GTMONK.kaiju.out -z 20 -v
kaiju -t ../../../kaiju/nodes.dmp -f ../../../kaiju/nr_euk/kaiju_db_nr_euk.fmi -i LANGDALE.bins.fa -o LANGDALE.kaiju.out -z 20 -v
kaiju -t ../../../kaiju/nodes.dmp -f ../../../kaiju/nr_euk/kaiju_db_nr_euk.fmi -i WINDING.bins.fa -o WINDING.kaiju.out -z 20 -v

# Total taxonomy counts
f=ATTINGHAM
kaiju2table -t ../../../kaiju/nodes.dmp -n ../../../kaiju/names.dmp -r phylum -o $f.phylum.tsv $f.kaiju.out &
kaiju2table -t ../../../kaiju/nodes.dmp -n ../../../kaiju/names.dmp -r class -o $f.class.tsv $f.kaiju.out &
kaiju2table -t ../../../kaiju/nodes.dmp -n ../../../kaiju/names.dmp -r order -o $f.order.tsv $f.kaiju.out &
kaiju2table -t ../../../kaiju/nodes.dmp -n ../../../kaiju/names.dmp -r family -o $f.family.tsv $f.kaiju.out &
kaiju2table -t ../../../kaiju/nodes.dmp -n ../../../kaiju/names.dmp -r genus -o $f.genus.tsv $f.kaiju.out &
kaiju2table -t ../../../kaiju/nodes.dmp -n ../../../kaiju/names.dmp -r species -l superkingdom,phylum,class,order,family,genus,species -o $f.species.tsv $f.kaiju.out &

# This is better - add taxon names to the output
kaiju-addTaxonNames -t ../../../kaiju/nodes.dmp -n ../../../kaiju/names.dmp -r superkingdom,phylum,class,order,family,genus,species -i ATTINGHAM.kaiju.out -o ATTINGHAM.names.out &
kaiju-addTaxonNames -t ../../../kaiju/nodes.dmp -n ../../../kaiju/names.dmp -r superkingdom,phylum,class,order,family,genus,species -i GTMONK.kaiju.out -o GTMONK.names.out &
kaiju-addTaxonNames -t ../../../kaiju/nodes.dmp -n ../../../kaiju/names.dmp -r superkingdom,phylum,class,order,family,genus,species -i LANGDALE.kaiju.out -o LANGDALE.names.out &
kaiju-addTaxonNames -t ../../../kaiju/nodes.dmp -n ../../../kaiju/names.dmp -r superkingdom,phylum,class,order,family,genus,species -i WINDING.kaiju.out -o WINDING.names.out &

# Get protein names from nr database (sqlite is fairly quick for this sort of query)
zgrep ">.*?\[" -oP nr.gz |sed 's/..$//'|sed 's/>//'|sed 's/MULTIGENE: //'|sed 's/ /|/' >nr.names
sqlite3 nr.db "CREATE TABLE nr(acc TEXT PRIMARY KEY, desc TEXT)"
sqlite3 -separator "|" nr.db ".import nr.names nr" 2>/dev/null

#echo "SELECT * FROM nr WHERE " >script.sql
#awk -F"\t" '{print $6}' OFS="," LANGDALE.kaiju.out|sed 's/.$//'|awk -F"," '{ for(i = 1; i <= NF; i++) { print "acc=\x27"$i"\x27 OR"; } }'|sed '$s/OR//' >>script.sql
#sqlite3 /data/data2/scratch2/deakig/kaiju/nr_euk/nr.db <script.sql > LANGDALE.prots.out
# nsqlite3 has a limit on the complexity of queries - 10000 OR statements
awk -F"\t" '{print $6}' OFS="," ATTINGHAM.kaiju.out|sed 's/.$//'|awk -F"," '{ for(i = 1; i <= NF; i++) { print "acc=\x27"$i"\x27 OR"; } }'|sed '$s/OR//'|split -l 9999
for f in x*; do 
 sed -i -e '$s/OR//' $f
 sed -i -e '1s/acc/SELECT * FROM nr WHERE acc/' $f
 sqlite3 /data/data2/scratch2/deakig/kaiju/nr_euk/nr.db <$f >> ATTINGHAM.prots.out
done

# then from R
#...R
library(tidyverse)
library(data.table)
dat <- fread("ATTINGHAM.names.out",fill=T,sep="\t")
dat[,acc:=sub(",.*","",V6)]
prot <- fread("ATTINGHAM.prots.out",header=F)
setnames(prot,c("acc","protein"))
prot <- unique(prot)
dat <- prot[dat,on=c("acc==acc")]




