library(data.table)
library(parallel)

# args[1] = concatenated uc file (reduced.txt)
# args[2] = tab file regex
# args[3] = output file
# args <- c("reduced.txt","*.tab","countData.sub_bins")
args <- commandArgs(TRUE)

# get list of tab files
tab_files <- list.files(".",args[2],full.names=F,recursive=F) 

# read sub bins
sub_bins <- fread(args[1],header=F)
setnames(sub_bins,c("BIN_NAME","SUB_BIN_NAME"))

# function for joining tab files with the sub bins and getting sub bin counts 
qqfun <- function(X) {
  # read tab file
  DT <- fread(X,header=F)
  # set col names of tab file
  setnames(DT,c("BIN_NAME","DOM","START","END","count"))
  # left join sub bins on tab file
  DDT <- sub_bins[DT,on="BIN_NAME"]
  # group by sub bin and sum counts per sub bin	
  DDT <- DDT[,.(Count=sum(count)),.(SUB_BIN_NAME)]
  setnames(DDT,"Count",sub("(\\.tab)","",X))
  # garbage collection
  gc()
  DDT
}

# garbage collection before running parallel job
gc()

# run the qqfun function in parallel
qq <- mclapply(tab_files,qqfun,mc.cores=8)
# qq <- lapply(tab_files,qqfun,new.env())

# merge sub bin counts
count_table <- Reduce(function(...) {merge(..., all = TRUE)}, qq)

# write the final count table
fwrite(count_table,args[3],sep="\t",quote=F,na=0)
       
