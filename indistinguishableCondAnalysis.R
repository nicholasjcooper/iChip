## This file is used to get the indistinguishable list 
# (bayes factors of the difference to the best) for each conditional hit
# same as the 'indistinguishable' analysis, but limited to the conditionals found to be significant


#R --no-save --slave < ~/github/iChip/conditionalAnalysis.R > condanalFri29.log

source("~/github/iChip/iFunctions.R")
work.dir <- "/chiswick/data/ncooper/iChipData/"
setwd(work.dir)
library(reader)
library(annotSnpStats)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")

### output: 'all.results'; a list containing the OR,OR_L,OR_H,P-value for each snp in 'topsnplist', by chromosome

bonf <- 3.23*(10^-7) # bonferroni threshold #
covs <- TRUE

if(!exists("covs")) { covs <- c(TRUE,FALSE)[2] }

load.uva.table <- TRUE  # results of UVA analyses for each SNP
load.ichip.regions <- TRUE # ichip regions as RangedData
restart.mid <- FALSE  #T
n.done <- 0  #16
use.imputation <- TRUE

# get hard coded lists of SNPs passing, failing QC, appearing in table 1, etc #tab1snps
source('~/github/iChip/hardCodedSnpLists.R', echo=FALSE)

print(load("conditionalSummaryObjects.RData"))
# get grpz,cond,condlist,glm.results 


ct.fn <- "conditionalTests.csv"
if(file.exists(ct.fn)) {
  ct <- reader(ct.fn,stringsAsFactors=F)
  ct[,1] <- ic.to.rs(ct[,1])
  ct[,2] <- ic.to.rs(ct[,2])
  if(any(is.na(ct[,1]))) { stop("invalid entry in table",ct.fn) }
  condlist <- apply(ct,1,function(X) { paste(c(narm(X[1:2]))) })
  groupz <- ct$GRP
  chrz <- ct$CHR
  cat("\n\nRunning conditional analyses from:",ct.fn,"\n")
  prv.large(ct,rows=20,cols=6)
} else {
 # e.g, define manually
 #condlist is loaded
 chrz <- Chr(sapply(condlist,"[",1))
 condlist <- c("rs2111485","rs10795791","rs66718203","rs757411","rs12720356")  #"rs2111485","rs10795791", #16 "rs66718203","rs12927355",
 chrz <- c(2,10,16,17,19)
 groupz <- c(1,1,1,1,1)
}

if(!restart.mid) { n.done <- 0 ; cond.count <- 0 }
# get resulting Snps from 'getMetaTable1.R'
print(load("finalMetaTopHits6FEB.RData"))

### get list of snps to test
##topsnplist <- bonf.snps
topsnplist <- unique(unlist(condlist))
#topsnplist <- condlist #c("rs2111485","rs10795791","rs757411","rs12720356")  #16 "rs66718203","rs12927355",

#gg <- grep("rs689",topsnplist)
#if(length(gg)>0) { topsnplist <- topsnplist[-gg] } # need to remove this SNP as not actually on iChip, just Taqman
topsnplist[topsnplist %in% "rs3842727"] <- "rs689" # do now with rs689 instead
topsnplist <- topsnplist[!topsnplist %in% not.top.snps]

print(load(cat.path(work.dir,"all.support.RData")))   ## load object: all.support [snp support for whole chip]
topsnplistic <- all.support$SNP[match(topsnplist,all.support$dbSNP)]

## load the results of the UVA meta, c/c and family analyses ##
if(load.uva.table & !restart.mid) {
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
  table1a$SNP_ID[match("imm_19_10324118",table1a$SNP_ID)] <- "rs34536443"
}



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
  if(!exists("ms")) { ms <- list.rowsummary(chr.dat) } else { cat("found existing whole genome QC\n") } 
  # ^ ie., if re-running the script, save repeating the QC
  #chrz <- 1:22
  ## create filter for call rate and Heterozygosity ##
  sample.filt <- samp.summ(ms,CR=0.953,HZlo=0.19,HZhi=0.235) # + print summary 
  # prints SNPQC overall summary, but actually for analysis this is done chr by chr
  ignore.for.now <- snp.summ(MAF=0.005,CR=.99,HWE=3.8905,qc.file="snpqc.RData") 
  excl.ids <- suppressWarnings(rownames(get.SnpMatrix.in.file(chr.dat[[22]]))[!sample.filt])  # ID exclusion list
  excl.snps <- clean.snp.ids(rownames(excl)) # from UVA, above
  #not.top.snps = non significant ids from #clear.not.tops <- lapply(all.results,function(Y) { lapply(Y,clearly.suck,thresh=20) } )
  excl.snps <- unique(c(excl.snps,qc.cloud.fail)) #,not.top.snps))
  sample.excl1 <- reader("sample-exclusions-2011-08-09.tab")[,2]
  sample.excl2 <- reader("unacknowledged-dups-2011-08-11.tab")[,5]
  excl.ids <- unique(c(excl.ids,sample.excl1,sample.excl2,nonconsent))
  ###
  if(!restart.mid) {
    all.results <- vector("list",length(groupz))
    file.out <- cat.path(work.dir,"indcond.results",suf=gsub(":",".",gsub(" ","_",date())),ext="RData")
  }
  # run separately for each chromosome #
  for(next.chr in unique(chrz[chrz>n.done])) {
    Header(paste("Chromosome",next.chr))
    chr <- next.chr
    print(load(cat.path(fn=ofn,suf=next.chr,ext="RData")))
    #loads: annotated.snp.support, t1d.data, t1d.support, control.data, control.support

    # gather data and support for next chromosome
    rownames(control.support) <- remove.X(rownames(control.support))
    rownames(t1d.support) <- remove.X(rownames(t1d.support))
    #snp.support <- clean.snp.support(annotated.snp.support) # fix little issues with SNP vs dbSNP names
    #snp.support <- clean.snp.support(snp.support)
    myData <- rbind(control.data,t1d.data) # combine the cases and controls
    colnames(myData) <- clean.snp.ids(colnames(myData))
    if(next.chr==11){
      rsmtch <- which( colnames(control.data) %in% "imm_11_2138800")
      if(length(rsmtch)==0)
      myData[,rsmtch] <- rep(as.raw("00"),times=nrow(myData)) # "imm_11_2138800" = "rs689"
      rs689.snpmatrix <- get(paste(load("taqman_rs689.RData"))) #get.taqman.snp.from.text.file()
      indzz <- match(rownames(rs689.snpmatrix),rownames(myData))
      myData[narm(indzz),rsmtch] <- rs689.snpmatrix[which(!is.na(indzz)),1]
    }
    cat("Applying QC thresholds to current chromosome dataset..")
    snp.qc <- col.summary(myData[sample.filt,])
    maf <- snp.qc$MAF>0.005
    clr <- snp.qc$Call.rate>.99
    hwe <- abs(snp.qc$z.HWE)<3.8905
    qc.excl.snps <- rownames(snp.qc)[which(!maf | !clr | !hwe)]
    excl.snps <- unique(c(excl.snps,qc.excl.snps))
    excl.snps <- excl.snps[!excl.snps %in% unique(c(unexclude,topsnplist))]
    cat("complete\n")
    ## extract covariates as 'cov.dat' ##
    if(covs) {
      cat("extracting covariates for current chromosome dataset...")
      nms <- c(rownames(control.support),rownames(t1d.support))
      region13 <- c(control.support$b58cregion,t1d.support$b58cregion) 
      the.sex <- c(control.support$sex,t1d.support$sex)
      the.pheno <- make.pheno(myData,rownames(t1d.support),rownames(control.support))
      cov.dat <- data.frame(sex=the.sex,region=region13,pheno=the.pheno)
      rownames(cov.dat) <- nms
     ## cov.dat <- fix.rownames(cov.dat)
      #cat(" using covariates for region and sex:\n")
      cat("complete\n")
      prv(cov.dat)
    }

    ## get the list of SNP/dbSNP ids that are in the current chromosome
    snpid.list <- topsnplist[topsnplist %in% all.support$dbSNP[all.support$Chr==next.chr]]
    snpic.list <- rs.to.ic(snpid.list) #snp.support$SNP[match(snpid.list,snp.support$dbSNP)]
    if(any(is.na(snpic.list))) { warning("missing dbSNP entries in list:",snpid.list[is.na(snpic.list)])}
    snpid.list <- snpid.list[!is.na(snpic.list)]
    snpic.list <- narm(snpic.list)
    if(length(snpic.list)<1) { cat("found no SNPs in this chromsome!\n"); next } # skip this chromosome if no SNPs were found to analyse  
   
    ### READY TO ADD RECWINDOW FROM 'WORKING.R' ##
    #recwindow(chr,st,en=st,window=0.1,bp=0)
    if(T) {
      which.snps <- match(snpid.list,all.support$dbSNP)
      if(any(is.na(which.snps))) { stop(paste("NAs in dbSNP match:",paste(snpid.list[is.na(which.snps)],collapse=","))) }
      snps.locs <- Pos(snpid.list) #snp.support$Pos[which.snps]
      if(!all(all.support$SNP[which.snps]==snpic.list)) { warning("bad rs-ic conversion!\n") }
      #warning("dup SNPs:",all.support$SNP[which.snps][duplicated(all.support$SNP[which.snps])],"\n")
      snp.rd <- RangedData(ranges=IRanges(start=snps.locs,
                           end=snps.locs,names=all.support$SNP[which.snps]),
                           space=rep(next.chr,length(snps.locs)))
      snp.rd <- annot.cnv(snp.rd,gs=cyto); colnames(snp.rd) <- "band"
      bands <- snp.rd$band
      nxt.window <- lapply(snps.locs, function(X,...) { recwindow(st=X,...) },chr=next.chr,window=1)
      if(next.chr==16) {        
        cat("changing window[[2]] from",paste(nxt.window[[2]],collapse=","),
         "to",paste(c(10663100,11601037),collapse=","),"\n")
        nxt.window[[2]] <- c(10663100,11601037) 
      }
      print(nxt.window) 
      st.window <- lapply(nxt.window, "[",1)
      en.window <- lapply(nxt.window, "[",2)
      n.snps <- vector("list",length(st.window))
      for(cc in 1:length(st.window)) {
        n.snps[[cc]] <- which(all.support$Chr==next.chr &
                              all.support$Pos>=st.window[cc] & 
                              all.support$Pos<=en.window[cc] &
                            (!all.support$SNP %in% excl.snps) &
                            (!all.support$dbSNP %in% excl.snps) 
        )
      }
      grp.labs <- lapply(n.snps,function(X) { all.support$SNP[X] })
      grp.snps <- lapply(n.snps,function(X) { all.support$dbSNP[X] })
      for (cc in 1:length(grp.labs)) { 
        if(any(is.na(grp.snps))) { stop(paste("NAs in grp.labs match:",paste(grp.snps[is.na(grp.snps)],collapse=","))) }
        grp.snps[[cc]][is.na(grp.snps[[cc]])] <- grp.labs[[cc]][is.na(grp.snps[[cc]])]
        grp.snps[[cc]][duplicated(grp.snps[[cc]])] <- grp.labs[[cc]][duplicated(grp.snps[[cc]])]
      }
      if(length(unique(bands))!=length(bands)) { warning("these bands are not unique ==> ",paste(bands[duplicated(bands)],collapse=",")) }
      grpz2 <- 1:length(bands)
      names(grp.snps) <- names(grp.labs) <- paste(grpz2,bands,sep=":")
    } 
    # ? ? ? grp.labs[[1]] <- NULL  
    n.grps <- length(grp.labs) ; chr.results <- vector("list",n.grps)
    #  if(length(snpic.list)!=n.grps) { 
    #    stop("why is length(snpic.list)=",length(snpic.list),"but length(grp.labs)=",n.grps) }
    #### MAIN LOOP OF TOP-SNP GROUPINGS ####
    grpz <- groupz[chrz==next.chr]
    if(any(!grpz %in% 1:length(grp.labs))) { stop("invalid group number(s) produced ",paste(grpz,collapse=",")) }
    for (grp in unique(grpz)) {
    #for (snp.cc in which(snpic.list %in% snpic.list[ct$GRP[ct$CHR==next.chr]==grp])) {
    #  cond.count <- cond.count+1
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
      toppers <- (names(ppp)[baduns] %in% rs.to.ic(topsnplist))
      if(any(toppers)) { warning("top snp ",names(ppp)[baduns][toppers]," failed sanger/uva qc") ; baduns <- baduns[-which(toppers)] }
      grp.labs[[grp]] <- grp.labs[[grp]][!grp.labs[[grp]] %in% rs.to.ic(names(ppp)[baduns])]
      grp.snps[[grp]] <- grp.snps[[grp]][!grp.snps[[grp]] %in% ic.to.rs(names(ppp)[baduns])]
      cat(length(baduns),"/",length(ppp)," SNPs failed on UVA vs Sanger genotype frequencies\n",sep="")
      if(length(baduns)>0 & length(baduns)<50) { cat(paste(names(ppp)[baduns],collapse=","),"\n") }      
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
    ## multi in grp loop
    rssubs <- condlist[which(ct$GRP==grp & ct$CHR==next.chr)]
    subgrps <- lapply(rssubs,rs.to.ic)
    subtxt <- unlist(lapply(rssubs,paste,collapse=","))
    cat("Looking at subgroups",paste(subtxt,collapse="; "),"\n")
    for (snp.cc in 1:length(subgrps)) {
      cond.count <- cond.count+1
      ## MAIN ANALYSIS ##
      # run GLM analysis for each SNP in the list on the current chromosome
      n.in.grp <- length(grp.labs[[grp]])
      #n.in.grp <- min(200,n.in.grp)
      ## decide between model with and without covariates ##
      cat("running GWAS on",n.in.grp,"SNPs in surrounding 1 centimorgan + 50Kb region of marker:",snpid.list[grp],"\n")
      cat("conditioning on: ",paste(condlist[[cond.count]],collapse=","),"\n")
      condn <- all.support$SNP[match(condlist[[cond.count]],all.support$dbSNP)] 
      if(any(!grp.labs[[grp]] %in% colnames(myDat))) { 
        warning("missing",length(which((!grp.labs[[grp]] %in% colnames(myDat)))),"SNPs in myDat")
        print(grp.labs[[grp]][which((!grp.labs[[grp]] %in% colnames(myDat)))])
        grp.snps[[grp]] <- grp.snps[[grp]][grp.labs[[grp]] %in% colnames(myDat)]
        grp.labs[[grp]] <- grp.labs[[grp]][grp.labs[[grp]] %in% colnames(myDat)]
        n.in.grp <- length(grp.labs[[grp]])
      }
      result <- vector("list",n.in.grp)
      bic <- numeric(n.in.grp)
      for(ss in 1:n.in.grp) {
        if(!covs) {
          fm <- paste("Pheno ~ ",paste(condn,collapse=" + "),"+",grp.labs[[grp]][ss])
        } else {
          fm <- paste("Pheno ~ sex + region +",paste(condn,collapse=" + "),"+",grp.labs[[grp]][ss])
        }
        nxt <- glm(as.formula(fm), family = binomial(logit),data=myDat)
        result[[ss]] <- mysumfun(nxt,p.digits=250,ci=T)[[1]]
        bic[ss] <- BIC(nxt)
        loop.tracker(ss,n.in.grp)
      }
      #names(result) <- paste(grp.labs[[grp]],grp.snps[[grp]],sep="/")
      #print(result) # print results with dbSNP and SNP id
      names(result) <- names(bic) <- ic.to.rs(grp.labs[[grp]])
      tort <- which(names(result) %in% cond)
      if(length(tort)>0) { result <- result[-tort] }
      #names(result) <- names(bic) <- grp.snps[[grp]]
      all.results[[cond.count]] <- list(GLM=result,BIC=bic)  # add result, just with dbSNP name
    } }
    #all.results[[next.chr]] <- chr.results
    save(all.results,file=file.out)
    cat("updated object in: ",file.out,"\n")
  }


  ##




if(F) {
  ### get LD for rs689 and rs3842753 ###
  
  #n DIL974 DIL969
  O0 <- as.raw("01") # i.e , letter O instead of zero!
  O1 <- as.raw("02")
  O2 <- as.raw("03")
  
  raw.mat <- rbind(
    cbind(rep(O0,180) , rep(O0,180) ),      #A  A   A   A
    cbind(rep(O0,120) , rep(O1,120) ),      #120 A  A   A   T
    cbind(rep(O0,341) , rep(O2,341) ),      #341 A  A   T   T
    cbind(rep(O1,14) , rep(O0,14) ),        #14 A  C   A   A
    cbind(rep(O1,1712) , rep(O1,1712) ),    #1712 A  C   A   T
    cbind(rep(O1,1) , rep(O2,1) ),          #1 A  C   T   T
    cbind(rep(O2,2095) , rep(O0,2095) ),    #2095 C  C   A   A
    cbind(rep(O2,2) , rep(O1,2)  ))         #2 C  C   A   T
  snpMat <- as(raw.mat,"SnpMatrix")
  colnames(snpMat) <- c("rs3842753","rs689")
  ld(snpMat,stats=c("D.prime","R.squared"),depth=1)
}
