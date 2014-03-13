## This file is used to get the conditional results list 
# iterates through each table 1 hit and tests whether any additional SNPs are able to significantly contribute
# to the model


#R --no-save --slave < ~/github/iChip/conditionalAnalysis.R > condanalThu16.log

source("~/github/iChip/iFunctions.R")
work.dir <- "/chiswick/data/ncooper/iChipData/"
setwd(work.dir)
library(reader)
library(annotSnpStats)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")

### output: 'all.results'; a list containing the OR,OR_L,OR_H,P-value for each snp in 'topsnplist', by chromosome

#bonf <- .05/23236      #3.23*(10^-7) # bonferroni threshold #
covs <- TRUE
update.snp.qc <- FALSE #TRUE
#thresholds <- list(MAF=0.005,CR=.99,HWE=3.8905,SCR=0.953,HZlo=.19,HZhi=.235,
#                   bonf=3.68*(10^-7),bonfcond=.05/23236 )
thresholds <- list(MAF=0.005,CR=.90,HWE=6.46,SCR=0.953,HZlo=.19,HZhi=.235,
                   bonf=3.68*(10^-7),bonfcond=.05/23236 ) # relaxed as exclusions can be done later
current.qc.file <- "SNPQC_FEB14.RData"

if(!exists("covs")) { covs <- c(TRUE,FALSE)[2] }

load.uva.table <- TRUE  # results of UVA analyses for each SNP
use.imputation <- FALSE

# get hard coded lists of SNPs passing, failing QC, appearing in table 1, etc #tab1snps
source('~/github/iChip/hardCodedSnpLists.R', echo=FALSE)

## load the results of the UVA meta, c/c and family analyses ##
if(load.uva.table) {
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  excl <- reader(docs[2])
}

print(load(cat.path(work.dir,"all.support.RData")))   ## load object: all.support [snp support for whole chip]

######################
#### MAIN PROGRAM ####
######################

### GET LIST OF FILE NAMES FOR RDATA FILES FOR EACH CHROMOSOME ###
ofn <- "/chiswick/data/ncooper/iChipData/temp.ichip-data.RData"
chr.dat <- cat.path(fn="temp.ichip-data",suf=paste(1:22),ext="RData")
chr.dat <- c(chr.dat,cat.path(fn="temp.ichip-data",suf=c("X","MT"),ext="RData"))
chr.dat <- as.list(chr.dat)

## do sample QC on each chrosome and recombine ##
if(!exists("ms")) { ms <- list.rowsummary(chr.dat) } # if re-running the script, save repeating the QC
chrz <- c(1:22,"X")  #,"MT")
## create filter for call rate and Heterozygosity ##
sample.filt <- samp.summ(ms,CR=thresholds$SCR,HZlo=thresholds$HZlo,HZhi=thresholds$HZhi) # + print summary 
# prints SNPQC overall summary, but actually for analysis this is done chr by chr
ignore.for.now <- snp.summ(MAF=thresholds$MAF,CR=thresholds$CR,HWE=thresholds$HWE,qc.file=current.qc.file) 
excl.ids <- 
  suppressWarnings(rownames(get.SnpMatrix.in.file(chr.dat[[22]]))[!sample.filt])  # ID exclusion list
excl.snps <- clean.snp.ids(rownames(excl)) # from UVA, above
#not.top.snps = non significant ids from #clear.not.tops <- lapply(all.results,function(Y) { lapply(Y,clearly.suck,thresh=20) } )
excl.snps <- unique(c(excl.snps,qc.cloud.fail)) #  ,not.top.snps))
sample.excl1 <- reader("sample-exclusions-2011-08-09.tab")[,2]
sample.excl2 <- reader("unacknowledged-dups-2011-08-11.tab")[,5]
excl.ids <- unique(c(excl.ids,sample.excl1,sample.excl2,nonconsent))
###

all.results <- vector("list",length(chrz))
file.out <- cat.path(work.dir,"gwas.results",suf=gsub(":",".",gsub(" ","_",date())),ext="RData")
if(update.snp.qc) { SNPQC <- vector("list",22) }

# run separately for each chromosome #
for(next.chr in chrz) {
  Header(paste("Chromosome",next.chr))
  chr <- next.chr
  nchr<- if(chr=="X") { 23 } else { as.numeric(chr) }
  print(load(cat.path(fn=ofn,suf=next.chr,ext="RData")))
  #loads: annotated.snp.support, t1d.data, t1d.support, control.data, control.support
  
  # gather data and support for next chromosome
  rownames(control.support) <- remove.X(rownames(control.support))
  rownames(t1d.support) <- remove.X(rownames(t1d.support))
  
  myData <- rbind(control.data,t1d.data) # combine the cases and controls
  colnames(myData) <- clean.snp.ids(colnames(myData))
  if(paste(next.chr)=="11"){
    rsmtch <- which( colnames(control.data) %in% "imm_11_2138800")
    myData[,rsmtch] <- rep(as.raw("00"),times=nrow(myData)) # "imm_11_2138800" = "rs689"
    rs689.snpmatrix <- get(paste(load("taqman_rs689.RData"))) #get.taqman.snp.from.text.file()
    indzz <- match(rownames(rs689.snpmatrix),rownames(myData))
    myData[narm(indzz),rsmtch] <- rs689.snpmatrix[which(!is.na(indzz)),1]
  }
  
  snp.qc <- col.summary(myData[sample.filt,])
  if(update.snp.qc) { SNPQC[[nchr]] <- snp.qc }
  maf <- snp.qc$MAF>thresholds$MAF
  clr <- snp.qc$Call.rate>thresholds$CR
  hwe <- abs(snp.qc$z.HWE)<thresholds$HWE
  qc.excl.snps <- rownames(snp.qc)[which(!maf | !clr | !hwe)]
  excl.snps <- unique(c(excl.snps,qc.excl.snps))
  excl.snps <- excl.snps[!excl.snps %in% unique(rs.to.ic(c(unexclude,table1snpsfinaljan30)))]
  
  ## extract covariates as 'cov.dat' ##
  if(covs) {
    nms <- c(rownames(control.support),rownames(t1d.support))
    region13 <- c(control.support$b58cregion,t1d.support$b58cregion) 
    the.sex <- c(control.support$sex,t1d.support$sex)
    the.pheno <- make.pheno(myData,rownames(t1d.support),rownames(control.support))
    cov.dat <- data.frame(sex=the.sex,region=region13,pheno=the.pheno)
    rownames(cov.dat) <- nms
    ## cov.dat <- fix.rownames(cov.dat)
    cat(" using covariates for region and sex:\n")
    prv(cov.dat)
  }
  ## get the list of SNP/dbSNP ids that are in the current chromosome
  snpid.list <- clean.snp.ids(all.support$dbSNP[paste(all.support$Chr)==paste(next.chr)])
  snpic.list <- rs.to.ic(snpid.list) #snp.support$SNP[match(snpid.list,snp.support$dbSNP)]
  if(any(is.na(snpic.list))) { warning("missing dbSNP entries in list:",snpid.list[is.na(snpic.list)])}
  snpid.list <- snpid.list[!is.na(snpic.list)]
  snpic.list <- narm(snpic.list)
  if(length(snpic.list)<1) { next } # skip this chromosome if no SNPs were found to analyse
  ff <- !snpic.list %in% excl.snps
  snpic.list <- snpic.list[ff]
  snpid.list <- snpid.list[ff]
  
  snp.filt <- (match(snpic.list,colnames(myData)))
  if(any(is.na(snp.filt))) { 
    warning("These were missing: ",paste(snpic.list[is.na(snp.filt)],collapse=","),"\n") ; 
    X <- snpic.list[is.na(snp.filt)]
    Y <- clean.snp.ids(snpic.list[is.na(snp.filt)]); print(Y)
    snp.filt <- narm(snp.filt) 
  } 
  #1_113034078,1_113044737,1_113056157
  smp.filt <- which(!rownames(myData) %in% excl.ids)
  
  ## select desired SNPs and samples, imputed and make a dataframe and SnpMatrix version
  # DON'T IMPUTE FOR GWAS as we don't need to compare BICs #
  if(use.imputation) {
    imp.file <- cat.path(work.dir,pref="Imputed_chr",fn=next.chr,suf=paste("grp",grp,sep=""),ext="RData")
    print(load(ichip.imputation(myData, imp.file, smp.filt, snp.filt)))
  } else {
    myDatSnp <- myData[smp.filt,snp.filt]
    #myDat <- SnpMatrix.to.data.frame(myDatSnp) # only need snp matrix for plain gwas
  }
  # NB: old way was to source("~/github/iChip/imputationBit.R")
  
  ## summarise dimension of working dataset ##
  cat("Analysis dataset:\n")
  print(Dim(myDatSnp))
  
  ## Collect phenotype information  ##
  Pheno <- make.pheno(myDatSnp,rownames(t1d.support),rownames(control.support))
  
  if(any(is.na(Pheno))) { 
    stop("some samples were neither in the T1D nor Control dataset") 
  } else {
    cat("successfully matched",length(which(Pheno==1)),"cases and",length(which(Pheno==0)),"controls\n\n")
  }
  
  ## check covariates, and if ok/present, put into sex and region variables ##
  if(covs) { 
    which.covs <- match(rownames(myDatSnp),rownames(cov.dat))
    if(any(is.na(which.covs))) { 
      warning("missing covariates, ignoring covariates") ; covs <- F 
    } else {
      sex <- cov.dat$sex[which.covs]; region <- as.factor(cov.dat$region[which.covs])
    }
  } 
  
  ## MAIN ANALYSIS ##
  # run snprhs analysis for each SNP in the list on the current chromosome
  n.in.grp <- length(snpic.list)
  ## decide between model with and without covariates ##
  cat("running chromosome wide analysis on",n.in.grp,"SNPs in chromosome:",next.chr,"\n")
  lambda <- lambdas(myDatSnp,pheno=Pheno,cc1000=TRUE)
  rownames(lambda) <- colnames(myDatSnp)
  
  all.results[[nchr]] <- lambda
  save(all.results,file=file.out)
  cat("updated object in: ",file.out,"\n")
}  


# extract all as a table
print(load(file.out))

# uncomment rest of these line sif names screw up
#all.results <- (all.results[-1:-22])
c2 <- do.call("rbind",all.results)
#cn1 <- do.call("c",sapply(all.results,rownames))
#rownames(c1) <- cn1
print(load("/chiswick/data/ncooper/iChipData/compiledTableAllResultsPassingQC.RData"))
near.region <- grep("EXT",tt$gene)
out.region <- grep("OTHER",tt$gene)
in.region <- which(!1:nrow(tt) %in% c(near.region,out.region)) 

p.m.reg <- names(p.meta[names(p.meta) %in% rs.to.ic(tt$names[in.region])])
p.m.out <- names(p.meta[names(p.meta) %in% rs.to.ic(tt$names[out.region])])

save(c2,p.m.reg,p.m.out,file="lambdaTableVars.RData")

## A PROBLEM= OFTEN CC BEST != BEST - how did i deal with this before???

if(update.snp.qc) { save(SNPQC,file=cat.path(work.dir,fn=paste("SNPQC",simple.date(),sep="_"),ext="RData")) }
  
##################################################################################
if(F) {
  if(covs) { indx <- 13 } else { indx <- 1 }
  
  ORs <- lapply(all.results,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,1)) }) } )
  Ps <- lapply(all.results,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,4)) }) } )

  all.fn <- cat.path(work.dir,"totalGWAS.RData")
  save(ORs,Ps,all.results,file=all.fn)
  cat("wrote file",all.fn,"\n")
  
  indxr <- match(topsnplist,table1a$SNP_ID)
  
  
  ss <- cbind(unlist(ORs),unlist(Ps))
  rownames(ss) <- gsub(".OR","",rownames(ss))
  tt <- cbind(ss[topsnplist,],table1a[indxr,c("OR_CC","P_CC")])[,c(1,3,2,4)]
  tt[,1] <- or.conv(tt[,1])
  tt[,2] <- or.conv(tt[,2])
  colnames(tt) <- c("ncOR","uvaOR","ncPv","uvaPv")
  print(tt)
  
}

