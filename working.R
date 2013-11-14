


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
print(load("/chiswick/data/ncooper/iChipData/dense.ic.regions.b36.RData"))
## gives  regions.gr

table.ranges <- RangedData(ranges=IRanges(start=poz,end=poz,names=table1a[,1]),
           space=table1a[,2],OR=table1a[,9],p.value=table1a[,11],fam.OR=table1a[,12],
           fam.p.value=table1a[,14],meta.OR=table1a[,15],meta.p.value=table1a[,17])
gs <- get.gene.annot()
table.ranges1 <- annot.cnv(table.ranges,gs=gs)
start(gs) <- pmax(1,(start(gs)-50000))
end(gs) <- end(gs)+50000
table.ranges2 <- annot.cnv(table.ranges,gs=gs)
start(gs) <- pmax(1,(start(gs)-50000))
end(gs) <- end(gs)+50000
table.ranges3 <- annot.cnv(table.ranges,gs=gs)
table.ranges1[["gene"]][which(table.ranges1$gene=="")] <- table.ranges2[["gene"]][which(table.ranges1$gene=="")]
table.ranges1[["gene"]][which(table.ranges1$gene=="")] <- table.ranges3[["gene"]][which(table.ranges1$gene=="")]

cyto <- get.cyto()
cy.ranges <- toGenomeOrder(as(regions.gr,"RangedData"))
cy.rangesX <- annot.cnv(cy.ranges,gs=cyto)

chr(cy.ranges) <- gsub("chr","",chr(cy.ranges))
cy.chr <- gsub("chr","",chr(cy.ranges))
annot.cnv(table.ranges,gs=cy.ranges)
table.rangesX <- annot.cnv(table.ranges,gs=cy.ranges)


#table.ranges <- annot.cnv(table.ranges,gs=gs)
table.ranges <- table.ranges1
jj <- order(table.ranges[["p.value"]])
tt <- as.data.frame(table.ranges)[jj,]
kk <- which(tt[["p.value"]]<10^-4)
tt <- tt[kk,]
prv.large(tt,rows=100)

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

