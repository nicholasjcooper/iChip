# This script takes the tops SNPs from table 1 and for each tests every other SNP within 1 centimorgan  + 50kb
#  and then calculates bayes factors (log) relative to the best SNP (because the best snp is based on meta 
# analysis which includes family data, sometimes the 'best' is not the top SNP for a region, but BF's are
# always relative to the meta-analysis top SNP in the region). Returns result list.

# this is a newer version than indistinguishableAnalysis.R


#R --no-save --slave < ~/github/iChip/conditionalAnalysis.R > condanalFri29.log

source("~/github/iChip/iFunctions.R")
work.dir <- "/chiswick/data/ncooper/iChipData/"
setwd(work.dir)
library(reader)
library(annotSnpStats)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")

### output: 'all.results'; a list containing the OR,OR_L,OR_H,P-value for each snp in 'topsnplist', by chromosome

bonf <- .05/135836  #3.64  #3.23*(10^-7) # bonferroni threshold #
thresholds <- list(MAF=0.005,CR=.99,HWE=3.8905,SCR=0.953,HZlo=.19,HZhi=.235,
                                     bonf=3.68*(10^-7),bonfcond=.05/23236 )
current.qc.file <- "SNPQC_FEB14.RData"
cm.window <- .1; bp.ext <- 0
covs <- TRUE

if(!exists("covs")) { covs <- c(TRUE,FALSE)[2] }

load.uva.table <- TRUE  # results of UVA analyses for each SNP
load.ichip.regions <- TRUE # ichip regions as RangedData
use.imputation <- TRUE

# get hard coded lists of SNPs passing, failing QC, appearing in table 1, etc #tab1snps now superceded
source('~/github/iChip/hardCodedSnpLists.R', echo=FALSE)

# get resulting Snps from 'getMetaTable1.R'
print(load("/chiswick/data/ncooper/iChipData/finalMetaTopHitsFEB17.RData"))

### get list of snps to test:
#### should i be reading bonf.snps from /chiswick/data/ncooper/iChipData/finalMetaTopHits15JAN.RData ??
#topsnplist <- tab1snps # table1snpsfinal ##NEW!  i think the former is right, although latter was last to be run
topsnplist <- bonf.snps
topsnplist[topsnplist %in% "rs3842727"] <- "rs689"
topsnplist <- topsnplist[!topsnplist %in% qc.cloud.fail]
topsnplist <- ids.by.pos(topsnplist)

## not sure whether this is necessary???
if(F) {
gg <- grep("rs689",topsnplist)
if(length(gg)>0) { topsnplist <- topsnplist[-gg] } # need to remove this SNP as not actually on iChip, just Taqman
}

print(load(cat.path(work.dir,"all.support.RData")))   ## load object: all.support [snp support for whole chip]
topsnplistic <- all.support$SNP[match(topsnplist,all.support$dbSNP)]

## load the results of the UVA meta, c/c and family analyses ##
if(load.uva.table) {
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
  #prv.large(table1a[,c(3,10,12,15,18)-1],rows=100,cols=7)
  poz <- as.numeric(table1a[,3])
  #table1a$SNP_ID[match("imm_19_10324118",table1a$SNP_ID)] <- "rs34536443"
}

#stop()

## load the ichip dense mapping regions  
if(load.ichip.regions) {
  # iChip regions
  print(load("/chiswick/data/ncooper/iChipData/ichip.regions.RData"))
  if(!exists("cyto"))  { cyto <- get.cyto(); cyto[["gene"]] <- rownames(cyto) }
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
  ## create filter for call rate and Heterozygosity ##
  sample.filt <- samp.summ(ms,CR=thresholds$SCR,HZlo=thresholds$HZlo,HZhi=thresholds$HZhi) # + print summary 
  # prints SNPQC overall summary, but actually for analysis this is done chr by chr
  ignore.for.now <- snp.summ(MAF=thresholds$MAF,CR=thresholds$CR,HWE=thresholds$HWE,qc.file=current.qc.file)
  excl.ids <- suppressWarnings(rownames(get.SnpMatrix.in.file(chr.dat[[22]]))[!sample.filt])  # ID exclusion list
  excl.snps <- clean.snp.ids(rownames(excl)) # from UVA, above
  #not.top.snps = non significant ids from #clear.not.tops <- lapply(all.results,function(Y) { lapply(Y,clearly.suck,thresh=20) } )
  excl.snps <- unique(c(excl.snps,qc.cloud.fail)) #,not.top.snps))
  sample.excl1 <- reader("sample-exclusions-2011-08-09.tab")[,2]
  sample.excl2 <- reader("unacknowledged-dups-2011-08-11.tab")[,5]
  excl.ids <- unique(c(excl.ids,sample.excl1,sample.excl2,nonconsent))
  if(T) {
    cat(out.of(length(excl.ids),nrow(ms))," samples fail on Hz, Callrate)\n",sep="")
    cat(out.of(length(sample.excl1),nrow(ms))," samples fail on DIL exclusions (MI?)\n",sep="")
    cat(out.of(length(sample.excl2),nrow(ms))," samples fail on duplications (DIL))\n",sep="")
    cat(out.of(length(excl.ids),nrow(ms))," samples fail for all reasons)\n",sep="")
  }
  ###

  all.results <- vector("list",22)
  file.out <- cat.path(work.dir,"all.results",suf=gsub(":",".",gsub(" ","_",date())),ext="RData")

  # run separately for each chromosome #
  for(next.chr in chrz[chrz>0 & chrz<23]) {
    Header(paste("Chromosome",next.chr))
    chr <- next.chr
    print(load(cat.path(fn=ofn,suf=next.chr,ext="RData")))
    #loads: annotated.snp.support, t1d.data, t1d.support, control.data, control.support

    # gather data and support for next chromosome
    rownames(control.support) <- remove.X(rownames(control.support))
    rownames(t1d.support) <- remove.X(rownames(t1d.support))
    ####moved to all.support
    ##snp.support <- clean.snp.support(annotated.snp.support) # fix little issues with SNP vs dbSNP names
    ##snp.support <- clean.snp.support(snp.support)
    myData <- rbind(control.data,t1d.data) # combine the cases and controls
    colnames(myData) <- clean.snp.ids(colnames(myData))
    if(next.chr==11){
      rsmtch <- which( colnames(control.data) %in% "imm_11_2138800")
      myData[,rsmtch] <- rep(as.raw("00"),times=nrow(myData)) # "imm_11_2138800" = "rs689"
      rs689.snpmatrix <- get(paste(load("taqman_rs689.RData"))) #get.taqman.snp.from.text.file()
      indzz <- match(rownames(rs689.snpmatrix),rownames(myData))
      myData[narm(indzz),rsmtch] <- rs689.snpmatrix[which(!is.na(indzz)),1]
    }
    snp.qc <- col.summary(myData[sample.filt,])
    maf <- snp.qc$MAF>thresholds$MAF
    clr <- snp.qc$Call.rate>thresholds$CR
    hwe <- abs(snp.qc$z.HWE)<thresholds$HWE
    qc.excl.snps <- rownames(snp.qc)[which(!maf | !clr | !hwe)]
      
    #qc.excl.snps <- qc.excl.snps[!qc.excl.snps %in% unexclude]
    excl.snps <- unique(c(excl.snps,qc.excl.snps))
    excl.snps <- excl.snps[!excl.snps %in% unique(c(unexclude,topsnplist))]
    
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
    snpid.list <- topsnplist[topsnplist %in% all.support$dbSNP[all.support$Chr==next.chr]]
    snpic.list <- rs.to.ic(snpid.list) #snp.support$SNP[match(snpid.list,snp.support$dbSNP)]
    if(any(is.na(snpic.list))) { warning("missing dbSNP entries in list:",snpid.list[is.na(snpic.list)])}
    snpid.list <- snpid.list[!is.na(snpic.list)]
    snpic.list <- narm(snpic.list)
    if(length(snpic.list)<1) { next } # skip this chromosome if no SNPs were found to analyse  

    ### get snps surrounding top snps by 0.1cM ##
    grp.labs <- get.nearby.snp.lists(snpid.list,cM=cm.window,bp.ext=bp.ext,build=37,excl.snps=excl.snps)
    
    n.grps <- length(grp.labs) ; chr.results <- vector("list",n.grps)
    if(length(snpic.list)!=n.grps) { 
      stop("why is length(snpic.list)=",length(snpic.list),"but length(grp.labs)=",n.grps) }
    #### MAIN LOOP OF TOP-SNP GROUPINGS ####
    for (grp in 1:n.grps) {
      snp.filt <- (match(grp.labs[[grp]],colnames(myData)))
      if(any(is.na(snp.filt))) { 
        warning("These were missing:",paste(grp.labs[[grp]][is.na(snp.filt)],collapse=","),"\n") ; 
        X <- grp.labs[[grp]][is.na(snp.filt)]
        Y <- clean.snp.ids(grp.labs[[grp]][is.na(snp.filt)]); print(Y)
        snp.filt <- narm(snp.filt) 
      } 
      #1_113034078,1_113044737,1_113056157
      smp.filt <- which(!rownames(myData) %in% excl.ids)
      ## select desired SNPs and samples, imputed and make a dataframe and SnpMatrix version
      if(use.imputation) {
        imp.file <- cat.path(work.dir,pref="Imputed_chr",fn=next.chr,suf=paste("grp",grp,sep=""),ext="RData")
        print(load(ichip.imputation(myData, imp.file, smp.filt, snp.filt)))
      } else {
        myDatSnp <- myData[smp.filt,snp.filt]
        myDat <- SnpMatrix.to.data.frame(myDatSnp)
      }
      # NB: old way was to source("~/github/iChip/imputationBit.R")
      
      ## summarise dimension of working dataset ##
      cat("Analysis dataset:\n")
      print(Dim(myDat))
      
      ## Collect phenotype information  ##
      Pheno <- make.pheno(myDat,rownames(t1d.support),rownames(control.support))

      ## do test of difference in genotype distribution between sanger and uva controls for each SNP ##
      sanger.not.uva <- substr(rownames(myDat[Pheno==0,]),1,1)=="7"
      ppp <- numeric(ncol(myDat)); names(ppp) <- colnames(myDat)
      for (cc in 1:ncol(myDat)) {
        ppp[cc] <- suppressWarnings(chisq.test(table(sanger.not.uva[Pheno==0],round(myDat[Pheno==0,cc])))$p.value)
      }
      baduns <- which(ppp<2.15e-06)  #(0.05/length(ppp)))
      #"imm_6_91053490"
      toppers <- (names(ppp)[baduns] %in% rs.to.ic(topsnplist))
      if(any(toppers)) { warning("top snp ",names(ppp)[baduns][toppers]," failed sanger/uva qc") ; baduns <- baduns[-which(toppers)] }
      grp.labs[[grp]] <- grp.labs[[grp]][!grp.labs[[grp]] %in% rs.to.ic(names(ppp)[baduns])]
      cat(length(baduns),"/",length(ppp)," SNPs failed on UVA vs Sanger genotype frequencies\n",sep="")
      #if(length(baduns)>0) { cat(paste(names(ppp)[baduns],collapse=","),"\n") }      
      if(length(baduns)==length(ppp)) { next; print("all failed") }
      if(any(is.na(Pheno))) { 
        stop("some samples were neither in the T1D nor Control dataset") 
      } else {
        cat("successfully matched",length(which(Pheno==1)),"cases and",length(which(Pheno==0)),"controls\n\n")
      }
      
      ## check covariates, and if ok/present, put into sex and region variables ##
      if(covs) { 
        which.covs <- match(rownames(myDat),rownames(cov.dat))
        if(any(is.na(which.covs))) { 
          warning("missing covariates, ignoring covariates") ; covs <- F 
        } else {
          sex <- cov.dat$sex[which.covs]; region <- as.factor(cov.dat$region[which.covs])
        }
      } 
      
      ## MAIN ANALYSIS ##
      # run GLM analysis for each SNP in the list on the current chromosome
      n.in.grp <- length(grp.labs[[grp]])
      #bic <- numeric(n.in.grp)
      ## decide between model with and without covariates ##
      cat("running GWAS on",n.in.grp,"SNPs in surrounding",cm.window,"centimorgan region of marker:",snpid.list[grp],"\n")
      if(any(!grp.labs[[grp]] %in% colnames(myDat))) {
        warning("missing ",length(which((!grp.labs[[grp]] %in% colnames(myDat))))," SNPs in myDat")
        print(grp.labs[[grp]][which((!grp.labs[[grp]] %in% colnames(myDat)))])
        grp.labs[[grp]] <- grp.labs[[grp]][grp.labs[[grp]] %in% colnames(myDat)]
        n.in.grp <- length(grp.labs[[grp]])
      }
      result <- vector("list",n.in.grp)
      bic <- numeric(n.in.grp)
      for(ss in 1:n.in.grp) {
        if(!covs) {
          fm <- paste("Pheno ~ ",grp.labs[[grp]][ss])
        } else {
          fm <- paste("Pheno ~ sex + region +",grp.labs[[grp]][ss])
        }
        nxt <- glm(as.formula(fm), family = binomial(logit),data=myDat)
        result[[ss]] <- mysumfun(nxt,p.digits=250,ci=T)[[1]]
        bic[ss] <- BIC(nxt)
        loop.tracker(ss,n.in.grp)
      }
      #names(result) <- paste(grp.labs[[grp]],grp.snps[[grp]],sep="/")
      #print(result) # print results with dbSNP and SNP id
      #names(result) <- names(bic) <- grp.snps[[grp]]
      names(result) <- names(bic) <- ic.to.rs(grp.labs[[grp]])
      chr.results[[grp]] <- list(GLM=result,BIC=bic,QC=ppp)  # add result, just with dbSNP name
    }
    all.results[[next.chr]] <- chr.results
    save(all.results,file=file.out)
    cat("updated object in: ",file.out,"\n")
  }


  ##
if(F) {
  summ <- bonfs.filt[,c(1:3,1:4)]
  summ[,4] <- Chr(summ[,1])
  summ[,5] <- Pos(summ[,1])
  summ[,6:7] <- AB(summ[,1])
  colnames(summ)[4:7] <- c("Chr","Pos","Allele.A","Allele.B")
  summ <- summ[order(as.numeric(summ[,4])),]
  rs689.swap <- which(summ[,1]=="rs3842727")
  summ[rs689.swap,] <- c("rs689",2.43E-201,161,11,2182224,"A","T")
  write.csv(summ,"table1summary2.csv")
  
  summ <- cbind(summ,table1a[match(summ[,1],table1a$SNP_ID),]) # everything but rs689
  write.csv(summ,"table1summaryanno.csv")
  
  #exonic  INS	11	2182224	T?	A?	0.3	T?	2.389	0.0020408	2.434966E-201	2.117	.04?	5.3E-47	2.29	.0325?	4.86E-130

}
##
