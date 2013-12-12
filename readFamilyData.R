dir1 <- "/chiswick/data/ncooper/iChipData/t1dgc/"
setwd(dir1)

source("~/github/plumbCNV/FunctionsCNVAnalysis.R")
library(NCmisc)
library(reader)
library(snpStats)


nchr <- 22

import.from.orig <- FALSE
if(import.from.orig) {
  chr <- 1:nchr
  cmd <- paste("plink --bfile t1dgc-asp-preqc --chr ",chr," --make-bed --out /chiswick/data/ncooper/iChipData/t1dgc/t1dgc-asp-preqc-chr",chr," --noweb",sep="")
  for (jj in 1:nchr) { system(cmd[jj]) }
}

if(import.from.orig) {
  fname.core <- "t1dgc-asp-preqc-chr"
  map.list <- fam.list <- vector("list",nchr)
  for (chr in 1:nchr) {
    bed <- cat.path(dir1,fname.core,suf=chr,ext="bed")
    bim <- cat.path(dir1,fname.core,suf=chr,ext="bim")
    fam <- cat.path(dir1,fname.core,suf=chr,ext="fam")
    snpDat <- read.plink(bed, bim, fam)
    snpMat <- snpDat$genotypes
    map.list[[chr]] <- snpDat$map
    fam.list[[chr]] <- snpDat$fam
    save(snpMat,file=cat.path(dir1,fn="SNPMAT",suf=chr,ext="RData"))
    loop.tracker(cc=chr,max=nchr)
  }
} 

for (chr in 1:nchr) {
  #snpMat <- load(cat.path(dir1,fn="SNPMAT",suf=chr,ext="RData"))
  snpMat <- get.SnpMatrix.in.file(cat.path(dir1,fn="SNPMAT",suf=chr,ext="RData"))
  loop.tracker(cc=chr,max=nchr)
}




## myMat <- read.plink(bed, bim, fam)  ## tooo big!

head(myMat$genotypes)
tail(myMat$genotypes)

sMat <- myMat$genotypes
map <- myMat$map
snp.info <- make.snp.info(map=myMat$map)
sample.info <- make.sample.info(dir=getwd(),id.list=rownames(myMat$fam),
                                phenotype=column.salvage(myMat$fam,"phenotype","affected"))
sel.22 <- which(colnames(sMat) %in% rownames(select.autosomes(snp.info)))
rs <- row.summary(sMat[,sel.22])
cs <- col.summary(sMat)

save(cs,rs,file="summaries.families.RData")
#save(cs,rs,cs2,file="summaries.finns.RData")
#save(cs,rs,cs2,snp.info,sample.info,file="summaries.finns.RData")

# assume chr ["X", "Y", "XY", "MT"] == [23, 24, 25, 26]
#load("summaries.finns.RData")

head(cs)
head(rs)

call.rate.summary(cs)
call.rate.summary(rs)

het.density.plots(rs,het.lo=.185,het.hi=.235)
excl.samp.sel <- hz.vs.callrate.plots(rs,callrate.samp.thr=.97,het.lo=.185,het.hi=.235,excl=T) 
excl.samp <- rownames(rs)[excl.samp.sel]
cs2 <- col.summary(sMat[-excl.samp.sel,])
writeLines(excl.samp,"callrateHzFailingSamples.txt")
system("cp /chiswick/data/Immunochip/project-added-by-Neil-2013-08-27/PLINK/finn-postqc-duplicate-samples.txt ~/Finnish/")
system("cp /chiswick/data/Immunochip/project-added-by-Neil-2013-08-27/PLINK/finn-postqc-monomorph-snps.txt  ~/Finnish/")
system("cp /chiswick/data/Immunochip/project-added-by-Neil-2013-08-27/PLINK/finn-postqc-sexcheck-samples.txt  ~/Finnish/")
	    
mendelz <- read.table(cat.path(dir1,"finn-d6-mendel.imendel"),header=T,comment="")
ME <- mendelz$N
hist(ME)
plot(density(ME[ME<200]))

length(which(ME>200))/length(ME)

hwe.density.plots(cs2)
hwe.vs.callrate.plots(cs2)


sample.info[["call.rate"]] <- rs$Call.rate
snp.info[["call.rate"]] <- cs$Call.rate
draw.density.plots(cs,snp.info=snp.info,sample.info=sample.info)

