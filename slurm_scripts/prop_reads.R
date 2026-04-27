library(data.table)
library(tidyverse)

# get command arguments
args <- commandArgs(TRUE)
# args <- c("CN2_10_EKDN220043202-1A_HKJKGDSX5_L3_1","S.pcounts","S.ncounts","S.new_counts")

file_suffix <- args[1]#gsub("\\..*","",list.files(".",".*pcounts$",full.names=F,recursive=F))

# load files
pcounts    <- fread(args[2])#lapply(file_suffix,function(i) fread(paste0(i,".pcounts"))) 
ncounts    <- fread(args[3])#lapply(file_suffix,function(i) fread(paste0(i,".ncounts"))) 

# load new_counts
# new_counts <- mclapply(file_suffix,freadFun,mc.cores=10)
fwrite(as.data.table(t(as.data.table(paste0("V",1:23)))),"header.txt",sep=",")
#freadFun <- function(i) fread(cmd=paste0("cat header.txt ",i,".new_counts"),fill=T,header=T)
new_counts <- fread(cmd=paste0("cat header.txt ",args[4]),fill=T,header=T)#lapply(file_suffix,freadFun)


del.rows <- function(DT, del.idxs) {           # pls note 'del.idxs' vs. 'keep.idxs'
  keep.idxs <- setdiff(DT[, .I], del.idxs);  # select row indexes to keep
  cols = names(DT);
  DT.subset <- data.table(DT[[1]][keep.idxs]); # this is the subsetted table
  setnames(DT.subset, cols[1]);
  for (col in cols[2:length(cols)]) {
    DT.subset[, (col) := DT[[col]][keep.idxs]];
    DT[, (col) := NULL];  # delete
  }
  return(DT.subset);
}

new_counts <- del.rows(new_counts,c(1,which(new_counts[,1]<2)))
#lapply(new_counts,function(DT) del.rows(DT,c(1,which(DT[,1]<2))))
del_cols <- names(new_counts)[c(1,length(names(new_counts)))]#lapply(new_counts,function(DT)names(DT)[c(1,length(names(DT)))])
new_counts[,(del_cols):=NULL];
new_counts[,seq:=(1:nrow(new_counts))]

#invisible(lapply(seq_along(new_counts),function(i){new_counts[[i]][,(del_cols[[i]]):=NULL];new_counts[[i]][,seq:=(1:nrow(new_counts[[i]]))]}))

pcounts[,V1:=as.character(V1)] #lapply(pcounts,function(DT) DT[,V1:=as.character(V1)])
ncounts[,V1:=as.character(V1)] #lapply(ncounts,function(DT) DT[,V1:=as.character(V1)])

multi_hits <- melt(new_counts,id.vars="seq") #lapply(new_counts,melt,id.vars="seq")

#invisible(lapply(multi_hits,function(DT){DT[,variable:=NULL];DT[,value:=as.factor(value)];setnames(DT,"value","V1")}))
multi_hits[,variable:=NULL];
multi_hits[,value:=as.factor(value)]
setnames(multi_hits,"value","V1")

multi_hits <- del.rows(multi_hits,which(multi_hits[,2]==""))#lapply(multi_hits,function(DT) del.rows(DT,which(DT[,2]=="")))

final_counts <- pcounts[multi_hits,on="V1"]
final_counts[,prop:=V2/sum(V2),by=seq]
final_counts <- final_counts[,lapply(.SD,sum),by=V1,.SDcols="prop"]
final_counts <- merge(final_counts, ncounts, all=TRUE)
final_counts[,tot:=rowSums(.SD, na.rm = TRUE), .SDcols =c("V2","prop")]
 
fwrite(final_counts,"corrected.counts",sep="\t")
#lapply(seq_along(final_counts),function(i) fwrite(final_counts[[i]],paste0(file_suffix[[i]],".corrected.counts"),sep="\t"))