library(data.table)

args <- commandArgs(TRUE)

domains  <- fread(args[1])

domains[,f := as.integer(substring(V2,nchar(V2)))]
domains[,t1:=V3*3-1]
domains[,t2:=V4*3-1]
domains[,t3:=(V5-V4)*3-3]
domains[,t4:=(V5-V3)*3-3]

qf <- function(X) {
	as.data.table(t(apply(X,1,function(l) {
		f <- l[5]
		r <- round(f/6,0)
		c(l[(r*2+1)] + f,l[(r*2+2)] + f)
	})))
}

domains[,c("start","end"):=qf(domains[,c(8:11,7)])]

colsToDelete <- c("t1","t2","t3","t4","V5","V6","f")
domains[,(colsToDelete):=NULL]

fwrite(domains,args[2],sep="\t",quote=F)


#domains[,"start":= apply(domains[,c(8:11,7)],1,function(l) {
#	f <- l[5]
#	r <- round(f/6,0)
#	l[(r*2+1)] + f
#})]
	
#domains[,"end":= apply(domains[,c(8:11,7)],1,function(l) {
#	f <- l[5]
#	r <- round(f/6,0)
#	l[(r*2+2)] + f
#})]

