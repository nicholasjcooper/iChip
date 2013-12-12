source("~/github/iChip/iFunctions.R")
setwd("/chiswick/data/ncooper/iChipData")
library(reader)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")

### output: 'all.results'; a list containing the OR,OR_L,OR_H,P-value for each snp in 'topsnplist', by chromosome

bonf <- 3.23*(10^-7) # bonferroni threshold #
if(!exists("first.deg")) { first.deg <- c("yes","no","none")[2] }
if(!exists("only.sibs")) { only.sibs <- c("only","no","any")[2] }
if(!exists("covs")) { covs <- c(TRUE,FALSE)[2] }

load.fam.info <- TRUE  # variable for affected family members
load.uva.table <- FALSE  # results of UVA analyses for each SNP
load.ichip.regions <- FALSE # ichip regions as RangedData

### get list of snps to test
topsnplist <- reader(fn="/chiswick/data/ncooper/iChipData/topsnplist.txt")
gg <- grep("rs689",topsnplist)
if(length(gg)>0) { topsnplist <- topsnplist[-gg] } # need to remove this SNP as not actually on iChip, just Taqman


## load affected status##
if(load.fam.info) {
  ### load a variable for each ID that says "yes" or "no" for whether there is an affected immediate fam / sib * 2 types
  id.switcher <- reader("~/Documents/necessaryfilesICHIP/alt.id.lookup.txt")
  first.degree <- reader("/ipswich/data/Immunochip/support/casecontrol/subject-lookup-2013-11-08.tab")
  
  sib.only <- (c("no","yes")[(1+as.numeric(first.degree$t1d_sib=="yes" & first.degree$t1d_mother=="no" & first.degree$t1d_father=="no"))])
  sib.incl <- first.degree$t1d_sib
  
  fsi <- match(first.degree$X.subjectid,id.switcher$subjectid)
  new.id <- id.switcher$sampleid[fsi]
  if(any(is.na(new.id))) { cat(length(which(is.na(new.id))),"/",length(new.id)," IDs were not found for first degree status\n")}
  FD <- first.degree$t1d_first_degree[!is.na(new.id)]
  FDSO <- sib.only[!is.na(new.id)]
  FDSI <- sib.incl[!is.na(new.id)]
  names(FD) <- names(FDSO) <- names(FDSI) <- narm(new.id)
  cat("summary of first degree status:\n")
  print(table(FD))
  cat("\n")
}


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
  prv.large(table1a[,c(3,10,12,15,18)-1],rows=100,cols=7)
  poz <- as.numeric(table1a[,3])
}
  

## load the ichip dense mapping regions  
if(load.ichip.regions) {
  # iChip regions
  print(load("/chiswick/data/ncooper/iChipData/ichip.regions.RData"))
  #  cyto <- get.cyto(); cyto[["gene"]] <- rownames(cyto)
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
  ###

  all.results <- vector("list",22)

  # run separately for each chromosome #
  for(next.chr in chrz) {
    Header(paste("Chromosome",next.chr))
    chr <- next.chr
    print(load(cat.path(fn=ofn,suf=next.chr,ext="RData")))
    #loads: annotated.snp.support, t1d.data, t1d.support, control.data, control.support
  
    # gather data and support for next chromosome
    snp.support <- clean.snp.support(annotated.snp.support) # fix little issues with SNP vs dbSNP names
    myData <- rbind(control.data,t1d.data) # combine the cases and controls
    
    ## extract covariates as 'cov.dat' ##
    if(covs) {
      nms <- c(rownames(control.support),rownames(t1d.support) )
      bdz <- substr(nms,1,1)
      nms[bdz=="X"] <- substr(nms,2,100000)[bdz=="X"]
      region13 <- c(control.support$b58cregion,t1d.support$b58cregion) 
      the.sex <- c(control.support$sex,t1d.support$sex)
      cov.dat <- data.frame(sex=the.sex,region=region13)
      rownames(cov.dat) <- nms
     ## cov.dat <- fix.rownames(cov.dat)
      cat(" using covariates for region and sex:\n")
      prv(cov.dat)
    }
    
    ## get the list of SNP/dbSNP ids that are in the current chromosome
    snpid.list <- topsnplist[topsnplist %in% snp.support$dbSNP]
    snpic.list <- snp.support$SNP[match(snpid.list,snp.support$dbSNP)]
    if(any(is.na(snpic.list))) { warning("missing dbSNP entries in list:",snpid.list[is.na(snpic.list)])}
    snpid.list <- snpid.list[!is.na(snpic.list)]
    snpic.list <- narm(snpic.list)
    if(length(snpic.list)<1) { next } # skip this chromosome if no SNPs were found to analyse
    
    ## set up affected / case filters #
    if(only.sibs!="no") { if(only.sibs=="any") { FD <- FDSI } else { FD <- FDSO } } # set which level of relatives to test
    # get list of YES and NO case IDs
    t1y <- rownames(t1d.data)[rownames(t1d.data) %in% names(FD[FD=="yes"])] # cases with first degree relative / sib
    t1n <- rownames(t1d.data)[rownames(t1d.data) %in% names(FD[FD=="no"])] # cases without first degree relative / sib
    ## filter to select only the iChip IDs from the dataset
    snp.filt <- narm(match(snpic.list,colnames(myData)))
    
    ### decide whether to select only cases with affective relatives, adjust sample filter accordingly
    if(first.deg!="none") {
      if(first.deg=="yes") {
        smp.filt <- which((!rownames(myData) %in% excl.ids) & (!rownames(myData) %in% t1n))
      } else {
        smp.filt <- which((!rownames(myData) %in% excl.ids) & (!rownames(myData) %in% t1y))
      }
    } else {
      smp.filt <- which(!rownames(myData) %in% excl.ids)
    }
    
    ## select desired SNPs and samples
    myDataFilt <- myData[smp.filt,snp.filt]
    myDat <- SnpMatrix.to.data.frame(myDataFilt)
    
    ## summarise dimension of working dataset ##
    cat("Analysis dataset:\n")
    print(Dim(myDat))
    
    
    ## Collect phenotype information  ##
    Pheno <- rep(NA,nrow(myDataFilt)) # missing (default)
    Pheno[rownames(myDataFilt) %in% rownames(t1d.data)] <- 1 # CASE if in the T1d dataset row (id) names
    Pheno[rownames(myDataFilt) %in% rownames(control.data)] <- 0 # CONTROL if in the Controls dataset row (id) names
    if(any(is.na(Pheno))) { 
      stop("some samples were neither in the T1D nor Control dataset") 
    } else {
      cat("successfully matched",length(which(Pheno==1)),"cases and",length(which(Pheno==0)),"controls\n\n")
      if(first.deg!="none") { cat("Cases with first degree relative:",length(t1y),"without:",length(t1n),"\n") }
    }
    
    ## check covariates, and if ok/present, put into sex and region variables ##
    if(covs) { 
      which.covs <- match(rownames(myDataFilt),rownames(cov.dat))
      if(any(is.na(which.covs))) { 
        warning("missing covariates, ignoring covariates") ; covs <- F 
      } else {
        sex <- cov.dat$sex[which.covs]; region <- as.factor(cov.dat$region[which.covs])
      }
    } 
    
    ## MAIN ANALYSIS ##
    # run GLM analysis for each SNP in the list on the current chromosome
    result <- vector("list",length(snpid.list))
    
    ## decide between model with and without covariates ##
    for(ss in 1:length(snpic.list)) {
      if(!covs) {
        fm <- paste("Pheno ~ ",snpic.list[ss])
      } else {
        fm <- paste("Pheno ~ sex + factor(region) +",snpic.list[ss])
      }
      nxt <- glm(as.formula(fm), family = binomial(logit),data=myDat)
      result[[ss]] <- mysumfun(nxt,p.digits=250,ci=T)[[1]]
    }
    names(result) <- paste(snpic.list,snpid.list,sep="/")
    print(result) # print results with dbSNP and SNP id
    names(result) <- snpid.list
    all.results[[next.chr]] <- result  # add result, just with dbSNP name
  }

  ##



### UNCOMMENT HERE TO MANUALLY SAVE THE RESULTS DESIRED ###
 #cov.results <- all.results
 #plain.results =  case/control analysis, univariate
 #cov.results =  case/control analysis with region and sex
 #yes.results =  case/control analysis, univariate, only affect sibs
 #cov.yes.results =  case/control analysis with region and sex, only affect sibs
 #no.results =  case/control analysis, univariate, no affect sibs
 #cov.no.results =  case/control analysis with region and sex, no affect sibs
 #save(plain.results,cov.results,yes.results,cov.yes.results,no.results,cov.no.results,file="/chiswick/data/ncooper/iChipData/results.RData")
 #save(yes.results,no.results,file="/chiswick/data/ncooper/iChipData/results2.RData")


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
