# script to replace HirBin extractSequences module in clusterBinsToSubbins.py
# This script is not that speedy (several days/weeks/months (dunno, gave up waiting) faster than the hirbin version mind) and requires an arbitrary amount of memory dependent on chunk size and length of each protein
# Memory usage peaks at about 40G for chunk size of 100 and mc.cores 10 and 4
# Changed output names to match gff/cov/tab names to prevent overcounting

# this is no longer a low mem version as it is using parallel - replace mclapply with lapply for a far slower implementation (but less memory intensive)

# extract_subbin_domains.R should be run first to get a domain file in the correct format

# args[1] = domains file (.hmm.cut)
# args[2] = protein file (.pep)
# args[3] = output directory
# args[4] = chunk size

library(Biostrings)
library(data.table)
library(parallel)

# get command arguments
args <- commandArgs(TRUE)

# args <- c("domains_2","LANGDALE.pep", "LANGDALE_hirbin_output/forClustering/",100,"LANGDALE.gff")

# set number of records to read in each loop
chunks <- as.numeric(args[4])
records <- as.numeric(system2("grep",paste(" '>'" ,args[2],"|wc -l",sep=" "),stdout=T))
read_in <- round(records/chunks,0)

# load domains
domains  <- fread(args[1])
print("loaded domains")
gc()

# load proteins in chunks and process
proteins <- do.call(c,mclapply(seq(0,(chunks -1)),function(i) {
	p <- readAAStringSet(args[2],nrec=read_in,skip=(i*read_in))
	names(p) <- sub(" .*","",names(p))
	d <- domains[V2%in%names(p)]
	setorder(d,-V2)
	p <- p[d$V2]
	#p <- p[order(names(p),decreasing=T)]
	d_IR <- IRanges(start=d$V3,end=(d$V4-1),names=paste(d$V2,d$V1,d$start,d$end,sep=" "))
	p <- AAStringSet(p,start=start(d_IR),width=width(d_IR))
	names(p) <- names(d_IR)
	mcols(p) <- d$V1
	p
},mc.cores=10))


print("loaded proteins")
print("named proteins")
print("subsetted proteins") 
print("ordered proteins")
print("extracted proteins")
print("set metadata")

# get domain list
domain_list <- domains[,logical(1),keyby=V1]$V1

rm(domains)

# set the output directory (args[3]) - writeXStringSet won't run if it doesn't exist)
dir.create(file.path(args[3]), showWarnings = FALSE,recursive=T)

gc()

print("Writing domains - slow")
# loop through the domain list and get subset of proteins matching the domain, then write them to file (this is the slow part)
invisible(mclapply(domain_list,function(dom) {writeXStringSet(proteins[mcols(proteins)$X==dom],paste0(args[3],dom,".fasta"))},mc.cores=4))
