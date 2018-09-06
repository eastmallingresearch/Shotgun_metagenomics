library(data.table)
# library(tidyverse)
library(parallel)


# args[1] = concatenated uc file (reduced.txt)
# args[2] = tab file regex
# args[3] = output file
#args <- commandArgs(TRUE)
 args <- c("reduced.txt","*.tab","countData.sub_bins")

# read in sub_bins
sub_bins   <- fread(args[1],header=F)
sub_bins[,"BIN_NAME":=paste(V1,V2,V4,V5,sep="_")]
sub_bins[,"SUB_BIN_NAME":=paste(V6,V7,V9,V10,sep="_")]
colsToDelete <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11")
sub_bins[, (colsToDelete) := NULL]


# read in tab files
tab_files <- list.files(".",args[2],full.names=F,recursive=F) 

qf <- function(X,env) {
	env$qq <- fread(X)
}
qq<-lapply(tab_files,qf,new.env())


# qq<-mclapply(tab_files,fread,mc.cores=8) - no speed advantage here

# get count file names and substitute to required format
names <- sub("(\\.tab)","",tab_files)

#### apply names to appropriate list columns and do some renaming  - splitting columns by string values (strsplit) is slow
gc()
colsToDelete <- c("TEMP","V1","DIR","BIN_ID")
mclapply(qq,function(DT) {
  DT[,c("BIN_ID","TEMP"):=tstrsplit(V1, "|",fixed=T)]
  DT[,c("DOM","START","END","DIR"):=tstrsplit(TEMP, ":",fixed=T)]
  setnames(DT,"V2","count")
  DT[,"BIN_NAME":=paste(BIN_ID,DOM,START,END,sep="_")]
  DT[, (colsToDelete) := NULL]
  setcolorder(DT,c("BIN_NAME","DOM","START","END","count"))
},mc.cores=4) # NOTE: this works as everything is being set by reference (DT references qq[[x]]), therefore no copies taken and original is modified


### sum_counts systimes: plyr;1563 DT;754 parallel_DT;150 

#sum_counts <- lapply(qq,function(bin_counts) left_join(bin_counts,sub_bins) %>% group_by(SUB_BIN_NAME) %>%  summarise(sum(count)))
#sum_counts <- lapply(1:length(sum_counts),function(i) {X<-sum_counts[[i]];colnames(X)[2]<- names[i];return(as.data.table(X))})
#count_table <- sum_counts %>% purrr::reduce(full_join,by="subbin") # dplyr method - much slower than using data table method here

# or data.table way - maybe (lots of copies, will it be faster?)
sum_counts <- mclapply(qq,function(DT) {
  DDT <- copy(sub_bins)
  DDT <- DDT[DT,on="BIN_NAME"] # this is not by reference
  DDT <- DDT[,.(Count=sum(count)),.(SUB_BIN_NAME)] # this is not by reference, but wayyyy faster than plyr
  DDT 
},mc.cores=12)
count_table <- Reduce(function(...) {merge(..., all = TRUE)}, sum_counts)
fwrite(count_table,"CHESTNUTS.countData.sub_bins",sep="\t")
       
countData[,"PFAM_NAME":=sub("(k[0-9]+_)([0-9]+_)(.*)(_[0-9]+_[0-9]+$)","\\3",countData$SUB_BIN_NAME)]