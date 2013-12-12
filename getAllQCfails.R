#R --no-save --slave < ~/github/iChip/conditionalAnalysis.R > condanalFri29.log

source("~/github/iChip/iFunctions.R")
work.dir <- "/chiswick/data/ncooper/iChipData/"
setwd(work.dir)
library(reader)
library(annotSnpStats)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")

if(!exists("covs")) { covs <- c(TRUE,FALSE)[2] }

load.uva.table <- TRUE  # results of UVA analyses for each SNP

### naughty SNPs ###
qc.cloud.fail <- c("imm_19_10324118","rs1016431","rs1160544","rs12743005","rs12927355","rs151234","rs1788097","rs202535","rs229533","rs2793108","rs3024493","rs34934650","rs3826110","rs4656327","rs4759229","rs4935063","rs56994090","rs5763751","rs5763779","rs61955086","rs666619","rs6691977","rs6710297","rs6907898","rs72729320","rs8074179","rs9585056","rs34536443","rs36038753") #,"rs61839660" last 2 new
ok <- c("rs3024493","rs61839660","rs56994090","rs12927355","rs151234","rs229533")
qc.cloud.fail <- qc.cloud.fail[!qc.cloud.fail %in% ok]


## load the results of the UVA meta, c/c and family analyses ##
if(load.uva.table) {
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  excl <- reader(docs[2])
}

######################
#### MAIN PROGRAM ####
######################

### GET LIST OF FILE NAMES FOR RDATA FILES FOR EACH CHROMOSOME ###
ofn <- "/chiswick/data/ncooper/iChipData/temp.ichip-data.RData"
chr.dat <- cat.path(fn="temp.ichip-data",suf=paste(1:22),ext="RData")
chr.dat <- as.list(chr.dat)

## do sample QC on each chrosome and recombine ##
if(!exists("ms")) { ms <- list.rowsummary(chr.dat) } # if re-running the script, save repeating the QC
chrz <- 1:22
## create filters for call rate and Heterozygosity ##
cr.filt <- ms$Call.rate>=0.953
hz.filt <- ms$Heterozygosity>=0.19 & ms$Heterozygosity<=0.235
sample.filt <- cr.filt & hz.filt  # logical filter
excl.ids <- rownames(get.SnpMatrix.in.file(chr.dat[[22]]))[!sample.filt]  # ID exclusion list
excl.snps <- rownames(excl) # from UVA, above
#not.top.snps = non significant ids from #clear.not.tops <- lapply(all.results,function(Y) { lapply(Y,clearly.suck,thresh=20) } )
excl.snps <- unique(c(excl.snps,qc.cloud.fail)) #  ,not.top.snps))
sample.excl1 <- reader("sample-exclusions-2011-08-09.tab")[,2]
sample.excl2 <- reader("unacknowledged-dups-2011-08-11.tab")[,5]
excl.ids <- unique(c(excl.ids,sample.excl1,sample.excl2))
###
SNPQC <- NULL
# run separately for each chromosome #
for(next.chr in chrz) {
  Header(paste("Chromosome",next.chr))
  chr <- next.chr
  print(load(cat.path(fn=ofn,suf=next.chr,ext="RData")))
  #loads: annotated.snp.support, t1d.data, t1d.support, control.data, control.support
  myData <- rbind(control.data,t1d.data) # combine the cases and controls
  colnames(myData) <- clean.snp.ids(colnames(myData))
  snp.qc <- col.summary(myData[sample.filt,])
  SNPQC <- rbind(SNPQC,snp.qc)
  maf <- snp.qc$MAF>0.005
  clr <- snp.qc$Call.rate>.99
  hwe <- abs(snp.qc$z.HWE)<3.75
  qc.excl.snps <- rownames(snp.qc)[which(!maf | !clr | !hwe)]
  qc.excl.snps <- qc.excl.snps[!qc.excl.snps %in% c("rs3842727","imm_11_2141424","rs12212193","imm_6_91053490")]
  excl.snps <- unique(c(excl.snps,qc.excl.snps))
  print(length(excl.snps))
}
writeLines(excl.snps,con="snpsExcluded2.txt")


