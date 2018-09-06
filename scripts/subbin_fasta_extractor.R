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
# args[5] = T/F (highmem)

library(Biostrings)
library(data.table)
library(parallel)


# get command arguments
args <- commandArgs(TRUE)
# args <- c("domains_2","LANGDALE.pep", "LANGDALE_hirbin_output/forClustering/",100,"T")


# set number of records to read in each loop
chunks <- as.numeric(args[4])
records <- as.numeric(system2("grep",paste(" '>'" ,args[2],"|wc -l",sep=" "),stdout=T))
read_in <- round(records/chunks,0)

# load domains
domains  <- fread(args[1],header=T)
print("loaded domains")
gc()

# load proteins in chunks and process (either using parallel, or not)
lfun <- function(i){
	p <- readAAStringSet(args[2],nrec=read_in,skip=(i*read_in))
	names(p) <- sub(" .*","",names(p))
	d <- domains[V2%in%names(p)]
	setorder(d,-V2)
	p <- p[d$V2]
	d_IR <- IRanges(start=d$V3,end=(d$V4-1),names=paste(d$V2,d$V1,d$start,d$end,sep=" "))
	p <- AAStringSet(p,start=start(d_IR),width=width(d_IR))
	names(p) <- names(d_IR)
	mcols(p) <- d$V1
	p
}

proteins <- if(args[5]){
	do.call(c,mclapply(seq(0,(chunks -1)),lfun,mc.cores=10))
} else {
	do.call(c,lapply(seq(0,(chunks -1)),lfun))
}

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
args[3] <- paste0(args[3],"/")
dir.create(file.path(args[3]), showWarnings = T,recursive=T)

gc()

# loop through the domain list and get subset of proteins matching the domain, then write them to file
# Using 4 cores here as proteins object is large -  
print("Writing domains")
wfun <- function(dom) {
	writeXStringSet(proteins[mcols(proteins)$X==dom],paste0(args[3],dom,".fasta"))
}
if(args[5]){
	invisible(mclapply(domain_list,wfun,mc.cores=4))
} else {
	invisible(lapply(domain_list,wfun))
}
