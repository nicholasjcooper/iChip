


setwd("/chiswick/data/ncooper/iChipData")
library(reader)
doc.path <- "ImChip_T1D_20130913"
docs <- cat.path(doc.path,list.files(doc.path))
table1 <- reader(docs[1])
excl <- reader(docs[2])
prv(table1)
prv(excl)
table1a <- table1[-which(rownames(table1) %in% rownames(excl)),]
table1a <- table1a[order(table1a[,3]),]
table1a <- table1a[order(table1a[,2]),]
#table1a <- table1a[order(table1a[,11]),]
prv.large(table1a[,c(3,10,12,15,18)-1],rows=100,cols=7)
poz <- as.numeric(table1a[,3])
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")


# iChip regions
if(F) {
  print(load("/chiswick/data/ncooper/iChipData/dense.ic.regions.b36.RData"))
  ## gives  regions.gr
  cyto <- get.cyto(); cyto[["gene"]] <- rownames(cyto)
  cy.ranges <- toGenomeOrder(as(regions.gr,"RangedData"))
  cy.rangesX <- annot.cnv(cy.ranges,gs=cyto)  ### YES!!!! ###
  ichip.regions <- cy.rangesX
  save(ichip.regions,file="/chiswick/data/ncooper/iChipData/ichip.regions.RData")
} else {
  print(load("/chiswick/data/ncooper/iChipData/ichip.regions.RData"))
}


table.ranges <- RangedData(ranges=IRanges(start=poz,end=poz,names=table1a[,1]),
           space=table1a[,2],OR=table1a[,9],p.value=table1a[,11],fam.OR=table1a[,12],
           fam.p.value=table1a[,14],meta.OR=table1a[,15],meta.p.value=table1a[,17])


gs <- ichip.regions
table.ranges1 <- annot.cnv(table.ranges,gs=gs)
gs[["gene"]] <- paste("EXT",gs[["gene"]],sep="_")
start(gs) <- pmax(1,(start(gs)-50000))
end(gs) <- end(gs)+50000
table.ranges2 <- annot.cnv(table.ranges,gs=gs)
gs2 <- cyto
gs2[["gene"]] <- paste("OTHER",gs2[["gene"]],sep="_")
table.ranges3 <- annot.cnv(table.ranges,gs=gs2)
table.ranges1[["gene"]][which(table.ranges1$gene=="")] <- table.ranges2[["gene"]][which(table.ranges1$gene=="")]
table.ranges1[["gene"]][which(table.ranges1$gene=="")] <- table.ranges3[["gene"]][which(table.ranges1$gene=="")]

table.ranges <- table.ranges1
MHC <- RangedData(ranges=IRanges(start=25000000,end=35000000),space=6)
remv <- queryHits(findOverlaps(table.ranges[6],MHC))
#new6 <- table.ranges[6][-remv,]
st6 <- (chrIndices(table.ranges)[6,1]-1)
TR <- table.ranges[-(remv+st6),]
if(length(queryHits(findOverlaps(TR[6],MHC)))) { cat("MHC not removed\n") } else { cat("MHC successfully removed\n") }


# chr(cy.ranges) <- gsub("chr","",chr(cy.ranges))
# cy.chr <- gsub("chr","",chr(cy.ranges))
# annot.cnv(table.ranges,gs=cy.ranges)
# table.rangesX <- annot.cnv(table.ranges,gs=cy.ranges)


#table.ranges <- annot.cnv(table.ranges,gs=gs)
jj <- order(TR[["p.value"]])
tt <- as.data.frame(TR)[jj,]
kk <- which(tt[["p.value"]]<10^-5)
tt <- tt[kk,]
prv.large(tt,rows=100)


highlights <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { cat("more than one top hit: ",paste(wh,collapse=","),"\n") }
  next.row <- c(wh,top,length(X),length(which(X<3.23*(10^-7))))
  return(next.row) 
}


conditional <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { cat("more than one top hit: ",paste(wh,collapse=","),"\n") }
  next.row <- c(wh,top,length(X),length(which(X<3.23*(10^-7))))
  return(next.row) 
}


out.list <- tapply(tt$p.value,tt$gene,highlights)
out.snps <- tapply(tt$names,tt$gene,"[",1)
grp.snps <- tapply(tt$names,tt$gene,c)
out.frame <- cbind(sapply(out.list,"[",1),sapply(out.list,"[",2),sapply(out.list,"[",3),sapply(out.list,"[",4))
colnames(out.frame) <- c("whichSNP","best.p.value","region.hits","hits<bonferroni")
out.frame <- cbind(out.frame,out.frame[,3]-out.frame[,4])
colnames(out.frame)[5] <- c("hits>=bonferroni")
out.frame[,1] <- as.character(out.frame[,1]) 
out.frame[,1] <- out.snps
genes <- get.gene.annot()
bandz <- paste(chr(genes),genes[["band"]],sep="")
nmz <- rownames(out.frame)
nmz <- gsub("OTHER_","",nmz)
nmz <- gsub("EXT_","",nmz)
genz <- lapply(nmz,function(nmz) { genes[["gene"]][which(bandz %in% nmz)] })

ofn <- "/chiswick/data/ncooper/iChipData/temp.ichip-data.RData"
chr <- 12; st <- 1; en <- get.chr.lens()[12]
system(paste("~/github/iChip/load-ichip.R chr=",chr," file=",out.file=ofn,sep=""))
#system("~/github/iChip/load-ichip.R --chr",chr,"--start",st,"--end",en,"--file",out.file=ofn)

print(load(cat.path(fn=ofn,suf=12,ext="RData")))
#annotated.snp.support, t1d.data, t1d.support, control.data, control.support
Ph <- rep(c(1,0),times=c(nrow(t1d.data),nrow(control.data)))
myData <- rbind(t1d.data,control.data)
myDataFilt <- myData[,grp.snps[[4]]]
snp.rhs.tests(as.formula("Ph ~ rs4759229"),snp.data=myDataFilt,allow.missing=.1,tests=grp.snps[[4]])

tab <- all.summary[order(all.summary[,"CHR"]),]


colnames(table1a)
cn <- c("CHR","MAF_UKcontrol","OR_CC","P_CC","OR_Fam","P_Fam","OR_Meta","P_Meta")
snp.list <- outframe$whichSNP
more.data <- table1a[match(snp.list,table1a[,1]),cn]
colnames(more.data)[2] <- "MAF"
all.summary <- cbind(out.frame,more.data)
print(all.summary[order(all.summary[,"CHR"]),])

do.one.chr <- function(chr.tab) {
  grp <- tapply(chr.tab$gene,chr.tab$gene,function(X) X)
  maxs <- tapply(chr.tab$p.value,chr.tab$gene,min,na.rm=T)
  #names(maxs) <- grp
}

#chr.results.list <- lapply(tt,do.one.chr)
results.list <- do.one.chr(tt)
prv(results.list)
gen.list <- names(results.list)
gen.list.spl <- strsplit(gen.list,";",fixed=T)
prv(gs)
bands <- lapply(gen.list.spl,
                function(X) { 
                  ind <- match(X,gs$gene)
                  OO <- paste(unique(paste(chr(gs)[ind],gs$band[ind],sep="")),collapse=";")
                  return(OO)  })
res2 <- data.frame(genes=substr(gen.list,1,20),p.value=as.numeric(results.list),band=unlist(bands))
res2 <- res2[order(res2[,2]),]
prv.large(res2,rows=50,cols=3)
unique(res2$band)




#table1:
# 1.Chr 2.causal  3.best  4.r2  5.all.M>m  6.MAF  7.OR  8.P  9.prev  10.All  11.OR  12.P  13.Ref
# 1: band - function lookup from location (chr=3, pos =4)
# 2: gene names - function lookup from location (chr=3, pos =4)
# 3: rsid - snpid in col2, some of which will need rsid lookup
# 4: r^2: could look in dataset at col2, versus c#9
# 5: col 5 > col 6 ?
# 6: col 7
# 7: col 10
# 8: col 12
# 9:  prev snp?
# 10: prev snp alleles ?
# 11: prev snp OR ?
# 12: prev snp P ?
# 13: prev snp ref ?

