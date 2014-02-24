#!/usr/bin/Rscript

## Chris W's load ichip script, slightly modified to use reader::parge.args() instead of chris' args function

library(methods)
library(reader)
#args <- getArgs(defaults=list(chr="19",snp=NULL, start=NULL, end=NULL, file="ichip-data.RData"))
args <- parse.args(com=c("chr","snp","start","end","file"),def=c("19","","","","ichip-data.RData"),list.out=TRUE)
#prv(args)
#args <- getArgs(defaults=list(chr="19",snp="rs34536443", start=NULL, end=NULL))
##print(args)
#args <- as.list(as.data.frame(t(argz)))
cat("Loading ichip data for chromosome",args$chr, "with optional limits:\n")
for (cc in 1:length(args)) { if(args[[cc]]=="") { args[[cc]] <- as.numeric(NULL) }  }
num.args <- c("start","end")
args[num.args] <- lapply(args[num.args], as.numeric)
#print(args[c("snp","start","end","file")])



## editted from Jason

library(snpStats, quietly=TRUE)

# define directory locations
Rscr1 <- "/ipswich/data/chrisw/IChip/data/SANGER"
Rscr2 <- "/ipswich/data/chrisw/IChip/data/UVA"
##Rscr3 <- "/ipswich/data/chrisw/IChip/data/CONTROLS"
Rscr4 <- "/ipswich/data/chrisw/IChip/data/T1D-patients"
Rscr5 <- "/ipswich/data/chrisw/IChip/data/BR4000-BR8000"
##Rscr6 <- "/ipswich/data/chrisw/IChip/data/RESULTS"

## snp support
cat("loading snp support\n")
(load(file = "/ipswich/data/chrisw/IChip/data/annotated-snp-support.RData")) # annotated.snp.support

annotated.snp.support <- annotated.snp.support[!is.na(annotated.snp.support$Chr),]
if(length(args$chr)>0)
  annotated.snp.support <- annotated.snp.support[annotated.snp.support$Chr == args$chr,]

if(length(args$start)>0)
  annotated.snp.support <- annotated.snp.support[annotated.snp.support$Pos > as.numeric(args$start),]

if(length(args$end)>0)
  annotated.snp.support <- annotated.snp.support[annotated.snp.support$Pos < as.numeric(args$end),]

if(length(args$snp)>0)
  annotated.snp.support <- annotated.snp.support[annotated.snp.support$dbSNP %in% args$snp |
                                                 rownames(annotated.snp.support) %in% args$snp ,]

#prv(annotated.snp.support)

annotated.snp.support <- annotated.snp.support[order(annotated.snp.support$Pos),]

#prv(annotated.snp.support)

target.snps <- rownames(annotated.snp.support)

#head(annotated.snp.support[,c("dbSNP","Pos","gene.annotation"),])
#tail(annotated.snp.support[,c("dbSNP","Pos","gene.annotation"),])


#
# sample supports
#

add.excl <- function(x) {
  excl.cols <- grep("^excl",colnames(x))
  x$excl <- apply(x[,excl.cols]==1,1,any, na.rm=TRUE)
  return(x)
}

load(file = paste(Rscr1,"/Robjects/SANGER-control-support.RData", sep = "")) # sample.support
sanger.support <- add.excl(sample.support)
#prv(sanger.support)

load(file = paste(Rscr2,"/Robjects/UVA-control-support.RData", sep = "")) # sample.support
uva.support <- add.excl(sample.support)

load(file = paste(Rscr4,"/Robjects/T1D-sample-support.RData", sep = "")) # sample.support
t1d.support <- add.excl(sample.support)

load(file = paste(Rscr5,"/Robjects/BR-sample-support.RData", sep = "")) # sample.support
br.support <- add.excl(sample.support)

cols <- intersect(intersect(colnames(sanger.support),colnames(br.support)),colnames(uva.support))
control.support <- rbind(sanger.support[,cols], uva.support[,cols], br.support[,cols])

rm(sample.support)

# load genotype data
fn <- paste(Rscr1,"/Robjects/SANGER-",args$chr,"-snp.RData", sep = "")
#prv(fn)
load(file = fn) # snp.data
#prv(target.snps)

if(!all(target.snps %in% colnames(snp.data))) {
  sd <- setdiff(target.snps, colnames(snp.data))
  cat(length(sd),"SNPs in target region according to snp.support file not found in SNP data and will be dropped:\n",sd,"\n")
  target.snps <- intersect(target.snps, colnames(snp.data))
  annotated.snp.support <- annotated.snp.support[target.snps,]
}
sanger.data <- snp.data[rownames(sanger.support),target.snps]

load(file = paste(Rscr2,"/Robjects/UVA-",args$chr,"-snp.RData", sep = "")) # snp.data
uva.data <- snp.data[rownames(uva.support),target.snps]

load(file = paste(Rscr5,"/Robjects/BR-",args$chr,"-snp.RData", sep = "")) # snp.data
br.data <- snp.data[rownames(br.support),target.snps]

load(file = paste(Rscr4,"/Robjects/T1D-",args$chr,"-snp.RData", sep = "")) # snp.data
t1d.data <- snp.data[rownames(t1d.support),target.snps]

rm(snp.data)

# apply allele switches and combine
load(file = paste(Rscr2,"/Robjects/UVA-",args$chr,"-snp-support.RData", sep = "")) # snp.support
sw.snps <- which(snp.support[target.snps, "switch"]==1)
uva.data <- switch.alleles(uva.data, sw.snps)

load(file = paste(Rscr5,"/Robjects/BR-",args$chr,"-snp-support.RData", sep = "")) # snp.support
sw.snps <- which(snp.support[target.snps, "switch"]==1)
br.data <- switch.alleles(br.data, sw.snps)

load(file = paste(Rscr4,"/Robjects/T1D-",args$chr,"-snp-support.RData", sep = "")) # snp.support
sw.snps <- which(snp.support[target.snps, "switch"]==1)
t1d.data <- switch.alleles(t1d.data, sw.snps)

# apply basic SNP QC
## mysumm <- function(x) {
##   cs <- col.summary(x)
##   cs$maf.NA <- is.na(cs$MAF)
##   cs$cr.less.99 <- (cs$Call.rate < 0.99)
##   cs$dev.hwe <- (abs(cs$z.HWE) > 3.09)
##   cs$excl <- cs$maf.NA | cs$cr.less.99 | cs$dev.hwe
##   return(cs)
## }

## cs.uva <- mysumm(uva.data)
## cs.sanger <- mysumm(sanger.data)
## cs.br <- mysumm(br.data)
## cs.t1d <- mysumm(t1d.data)

# Sanger ctrl v UVA ctrl exclusions - note that you really need to determine the threshold P-value
# and apply this to the ctrl analysis as well

# load ctrl v ctrl exclusions, none at the MAF > 1%
## excl <- union(rownames(cs.uva[cs.uva$excl,]), rownames(cs.sanger[cs.sanger$excl,]))
## excl <- union(excl, rownames(cs.br[cs.br$excl,]))
## m <- match(colnames(sanger.data), excl, nomatch = 0)
## length(keep <- colnames(sanger.data[,(m == 0)]))
## sanger.data <- sanger.data[,keep]

#prv(sanger.data); prv(uva.data); prv(br.data)

control.data <- rbind(sanger.data, uva.data, br.data)

## check
if(!identical(rownames(control.data),rownames(control.support)))
  stop("rbinding of control data failed to preserve sample order\n")

# check
par(mfrow=c(1,2))
cs.sanger <- col.summary(sanger.data)
cs.uva <- col.summary(uva.data)
cs.br <- col.summary(br.data)
snps <- intersect(colnames(sanger.data), colnames(uva.data))
plot(cs.sanger[snps, "P.AA"], cs.uva[snps, "P.AA"])

snps <- intersect(colnames(sanger.data), colnames(br.data))
plot(cs.sanger[snps, "P.AA"], cs.br[snps, "P.AA"])

# apply sample exclusions
cat("\nSample exclusions ....\n")

cat("controls")
print(table(control.support$excl))
control.data <- control.data[!control.support$excl,]
control.support <- control.support[!control.support$excl,]

cat("t1d")
print(table(t1d.support$excl))
t1d.data <- t1d.data[!t1d.support$excl,]
t1d.support <- t1d.support[!t1d.support$excl,]

cat("\nRemaining ....\n")
cat("SNPs:\t",ncol(t1d.data),"\n")
cat("T1D samples:\t",nrow(t1d.support),"\n")
cat("Control samples:\t",nrow(control.support),"\n")

## tidy a little
annotated.snp.support$SNP <- rownames(annotated.snp.support)
rownames(annotated.snp.support) <- make.names(rownames(annotated.snp.support))
colnames(t1d.data) <- make.names(colnames(t1d.data))
colnames(control.data) <- make.names(colnames(control.data))
rownames(t1d.support) <- make.names(rownames(t1d.support))
rownames(control.support) <- make.names(rownames(control.support))

#
# save
#
mysave <- function(..., file) {
  dir <- dirname(file)
  if(!file.exists(dir))
    dir.create(dir,recursive=TRUE)
  save(..., file=file)
}

#mysave(annotated.snp.support, t1d.data, t1d.support, control.data, control.support,
#       file = args$file)


save(annotated.snp.support, t1d.data, t1d.support, control.data, control.support,
       file = cat.path(dir="/chiswick/data/ncooper/iChipData/",fn=basename(args$file),suf=args$chr,ext="RData"))



