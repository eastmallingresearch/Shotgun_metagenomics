# Depricated use subbin_parser_v2.R in preference

library(data.table)

# get command arguments
args <- commandArgs(TRUE)

# location of files to load
tmpdir <- paste0(args[1],"/")

# load count files
qq <- lapply(list.files(tmpdir ,args[2],full.names=T),function(x) fread(x,sep="\t",select=c(9,10)))

# get the sample names  
names <- sub("_1\\.cov","",list.files(tmpdir ,args[2],full.names=F,recursive=F))

# aggregate by domain

qa <- lapply(qq,function(DT) DT[,sum(V10),by = V9])

# apply names to appropriate list columns (enables easy joining of all count tables)
qa <- lapply(seq(1:length(qa)),function(i) {X<-qa[[i]];colnames(X)[2] <- names[i];return(X)})

# merge count tables (full join)
countData <- Reduce(function(...) {merge(..., all = TRUE)}, qa)

# remove ID= from domains
countData[,V9:=sub("ID=","",countData[,V9])]

colnames(countData)[1] <- "Domain"

countData <- countData[,lapply(.SD, function(x) {x[is.na(x)] <- "0" ; x})]

write.table(countData,args[3],sep="\t",quote=F,row.names=F,col.names=T)
