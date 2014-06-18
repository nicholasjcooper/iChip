
get.indel.from.text.file <- function(fn="cc-rs333.tab",label="rs333",dir="/chiswick/data/ncooper/iChipData/",as.code=TRUE) {
  cur.dir <- getwd()
  lookup.tab <- read.delim("/ipswich/data/Immunochip/support/casecontrol/sample-lookup-2011-08-09.tab")
  setwd(dir)
  geno.dat <- read.delim(fn)
  new.ids <- paste(lookup.tab[,2][match(geno.dat[,1],lookup.tab[,4])])
  bad.uns <- (new.ids=="NA" | is.na(new.ids))
  geno.dat[!bad.uns,1] <- new.ids[!bad.uns]
  geno.dat[["SNP"]] <- label
  geno.dat[["A1"]] <- gsub("INS","+1",gsub("DEL","-1",geno.dat[,4])) 
  geno.dat[["A2"]] <- gsub("INS","+1",gsub("DEL","-1",geno.dat[,5]))
  setwd(cur.dir)
  if(as.code) { 
   a1 <- as.numeric(geno.dat$A1)
   a2 <- as.numeric(geno.dat$A2)
   code <- rep(NA,times=length(a1))
   code[(a1==-1 & a2==-1)] <- 0
   code[(a1==-1 & a2==1)] <- 1
   code[(a1==1 & a2==-1)] <- 1
   code[(a1==1 & a2==1)] <- 2
   names(code) <- geno.dat[,1]
   return(code)
  } else {
   return(geno.dat)
  }
}

get.taqman.snp.from.text.file <- function(fn="rs689.tab",label="rs689",dir="/chiswick/data/ncooper/INS/")
{
  # get a taqman snp from a text file into snpMatrix format
  ## The SNP data file 'rs689.tab' was obtained as follows:
  #cd /ipswich/data/shared/neil/ncooper
  #head * | cut -f 1-7 
  #grep DIL969 *support*
  #cat cc-genotype.tab  | cut -f 1-3,68-69 > /chiswick/data/ncooper/INS/rs689.tab
  cur.dir <- getwd()
  ## alternative ID lookup data for all samples
  lookup.tab <- read.delim("/ipswich/data/Immunochip/support/casecontrol/sample-lookup-2011-08-09.tab")
  #X.sampleid X.dil_subjectid
  #sanger-controls 5412557223_R06C01 UKBS-547071
  #t1d-cases 53591710  101321602A
  #uva-controls 7852856  B58C48_7G
  
  setwd(dir)
  geno.dat <- read.delim(fn)
  new.ids <- paste(lookup.tab[,2][match(geno.dat[,1],lookup.tab[,4])])
  bad.uns <- (new.ids=="NA")
  geno.dat[!bad.uns,1] <- new.ids[!bad.uns]
  geno.dat[["SNP"]] <- label
  # write data to a text file so can use the 'read.long' function to convert to snpMatrix
  write.table(geno.dat[!bad.uns,],file="geno.dat",sep="\t",col.names=F,row.names=F,quote=F)
  mysnp <- read.long("geno.dat", samples=16453, snps=1,
                     fields = c(snp = 6, sample = 1, genotype = NA, confidence = NA,allele.A = 4, allele.B = 5),
                     no.call="N")
  paste(mysnp$genotypes[1:5,1])
  
  Y <- mysnp$genotypes
  setwd(cur.dir)
  return(Y)
}


get.snp.from.T1DGC <- function(rs.id="rs2334499")
{
  # returns the data (in a vector of snp.matrix format)
  # for a single snp present in T1d/controls datasets in T1DGC files
  # function is very specific to this directory, etc
  ## eg. GET DATA FOR rs2334499 ##
  require(snpStats)
  dir <- "/ipswich/data/T1DGC/R-objects/genotypes"
  fnms <- paste(dir,c("case.snp.chr11.RData","control.snp.chr11.RData",
                      "sanger.snp.chr11.RData"),sep="/")
  
  ## alternative ID lookup data for all samples
  lookup.tab <- read.delim("/ipswich/data/Immunochip/support/casecontrol/sample-lookup-2011-08-09.tab")
  #X.sampleid X.dil_subjectid
  #sanger-controls 5412557223_R06C01 UKBS-547071
  #t1d-cases 53591710  101321602A
  #uva-controls 7852856  B58C48_7G
  mini.mat <- vector("list",length(fnms))
  for (cc in 1:length(fnms)) {
    next.file <- get(paste(load(fnms[cc])))
    print(fnms[cc])
    col.loc <- (which(colnames(next.file) %in% rs.id))
    print(col.loc)
    mini.mat[[cc]] <- next.file[,col.loc]
  }
  samp.sup <- "/ipswich/data/T1DGC/R-objects/support/WTCCC-sample-support.RData"
  samp1 <- get(paste(load(samp.sup))) #sanger.id WTCCC123.., dil.id# UKBS-5...
  samp.sup <- "/ipswich/data/T1DGC/R-objects/support/T1DGC-sample-support.RData"
  samp2 <- get(paste(load(samp.sup))) # subject 101010...C #sample  WG0012-DNA_..
  
  samp1.ids <- lookup.tab[match(samp1$dil.id,lookup.tab$X.dil_subjectid),"X.sampleid"]
  samp2.ids <- lookup.tab[match(samp2$subject,lookup.tab$X.dil_subjectid),"X.sampleid"]
  
  minimat <- mini.mat
  s1nr_of_mini3 <- match(names(mini.mat[[3]]),samp1$sanger.id)
  s2nr_of_mini1 <- match(names(mini.mat[[1]]),samp2$sample)
  names(minimat[[3]]) <- samp1.ids[s1nr_of_mini3]
  names(minimat[[1]]) <- samp2.ids[s2nr_of_mini1]
  
  case.new.snp <- minimat[[1]][!is.na(names(minimat[[1]]))]
  cont.new.snp <- minimat[[3]][!is.na(names(minimat[[3]]))]
  X <- c(case.new.snp,cont.new.snp)
  print(is(X)); print(length(X))
  rN <- names(X)  
  Y <- new("SnpMatrix",(X))
  print(dim(Y))
  rownames(Y) <- rN; colnames(Y) <- rs.id
  return(Y)
}



mysumfun <- function(glmr,o.digits=3,p.digits=6,lab=T)
{
  co <- summary(glmr)$coefficients
  predz <- rownames(co)[-1]
  label <- paste(summary(glmr)$call)[2]
  p <- co[2:nrow(co),4]; #print(p)
  o.r <- exp(co[2:nrow(co),1]); #print(o.r)
  p <- round(p,p.digits)
  out.fr <- cbind(round(o.r,o.digits),p)
 # outlist <- list(round(o.r,o.digits),p)
 # names(outlist) <- c("OR","p-value")
  colnames(out.fr) <- c("OR","p-value")
  if(lab) { out.fr <- list(out.fr); names(out.fr) <- label }
  return(out.fr)
}


impute.rs7111341 <- function(myData) 
{
  # function does nothing but attempt to impute this snp using the others
  newData <- myData[(!is.na(rowMeans(myData))),]
  for (cc in 1:ncol(myData)) {
    newData[,cc] <- as.factor(newData[,cc])
  }
  model.snp <- polr(rs7111341 ~ rs689 + rs11043056 + rs6357 + rs76275158 + rs3842727,  data=newData)
  tt <- table(predict(model.snp),newData$rs7111341)
  print(tt)
  impute.pc <- sum(diag(tt))/sum(tt)
  return(impute.pc)
}


# DIL specific, extract the ichip data for a specified gene
extract.ichip.data.for.gene <- function(chr.n=11,gn.nm="INS",ichipLabels,rs.ids,dil.grps=c("T1D","BR","UVA","SANGER"))
{
  # define directories
  home.dir <-  "/ipswich/data/jason/ImmunoChip/"
  grp.nms <- c("T1D-patients","BR4000-BR8000","UVA","SANGER")
  data.dir <- paste(grp.nms,"/Robjects/",sep="")
  ## data and support file defn ##
  prefixes <- dil.grps # RData file prefixes for each grp
  main.sup.fn <- paste(home.dir,"annotated-snp-support.RData",sep="")
  dat.fn <- paste(home.dir,data.dir,prefixes,"-",chr.n,"-snp.RData",sep="")         
  sup.fn <- paste(home.dir,data.dir,prefixes,"-",chr.n,"-snp-support.RData",sep="")
  ## load in RData files
  # main support info 'annotated.snp.support'
  print(paste("Loaded object:",load(main.sup.fn)))
  ## Limit data and main support to target chromosome and sort by position
  chr.support <- (annotated.snp.support[annotated.snp.support$Chr==chr.n,])
  chr.support <- chr.support[order(as.numeric(chr.support$Pos)),]
  ## subset tagged with target gene name
  annotated.subset <- grep(gn.nm, paste(chr.support$gene.annotation))
  # select all Snps with position between max/min of range
  annotated.subset <- c(min(annotated.subset):max(annotated.subset))
  new.support <- chr.support[annotated.subset,]
  # read in the four datasets (data and support) to lists
  # select those only in relevant gene range
  Snp.data <- Snp.support <- list()
  grp.lens <- numeric()
  for (cc in 1:length(dat.fn)) {
    # chromosome data for target gene 'snp.data'
    print(paste("Loaded object:",load(dat.fn[cc])))
    Snp.data[[cc]] <- (snp.data)[,match(rownames(new.support),colnames(snp.data))]
    # support for chromosome 'snp.support'
    print(paste("Loaded object:",load(sup.fn[cc]))) 
    nr <- nrow(new.support)
    Snp.support[[cc]] <- data.frame(allele.A=character(nr),allele.B=character(nr),
                                    chr=integer(nr),switch=logical(nr))
    col.match <- match(colnames(Snp.support[[cc]]),colnames(snp.support))
    Snp.support[[cc]][,!is.na(col.match)] <- (snp.support)[match(rownames(new.support),
                                                                 rownames(snp.support)),col.match[!is.na(col.match)]]
    rownames(Snp.support[[cc]]) <- rownames(snp.support)[match(rownames(new.support),
                                                               rownames(snp.support))]
  }
  if(!all(sapply(Snp.data,is)[1,]=="SnpMatrix")) { stop("SnpMatrix objects not correctly imported") }
  if(!all(sapply(Snp.support,is)[1,]=="data.frame")) { stop("Snpsupport data.frame objects not correctly imported") }
  # merge from list to single object
  grp.lens <- sapply(Snp.data,nrow)
  grp.lab <- as.factor(paste(unlist(mapply(rep,grp.nms,grp.lens))))
  SNP.data <- do.call("snp.rbind",Snp.data)
  SNP.support <- do.call("cbind",Snp.support)
  out <- list(SNP.data,SNP.support,grp.lab)
  names(out) <- c("SNP.data","SNP.support","grp.lab")
  return(out)
}


# easier LD function for snpStats::SnpMatrix objects
# targ.cols allows specification of which snps to calculate LD for either by name, column # or logical vector,
#  but actually it's probably easier to just force evaluation of the whole frame
# r2 is an easy way of specifying whether to use r2 (default) or d.prime
# targ.rownames allows specification of custom rownames for the resulting table
get.LD.mat <- function(SNP.data,r2=TRUE,targ.cols=c(1:ncol(SNP.data)),targ.rownames=NULL)
{
  if(all(targ.cols %in% colnames(SNP.data))) { 
    targ.cols <- match(targ.cols, colnames(SNP.data)) 
  } else {
    if(any(targ.cols %in% colnames(SNP.data))) { 
      warning("some values of targ.cols were not found in the colnames of SNP.data, will use all columns")
    }
  }
  if(is.logical(targ.cols)) { targ.cols <- which(targ.cols) }
  if(any(!targ.cols %in% 1:ncol(SNP.data))) { 
    warning("invalid targ.cols entered, reverting to all"); targ.cols <- 1:ncol(SNP.data) }
  dm <- length(targ.cols)
  if(length(targ.rownames)!=dm) { targ.rownames <- NULL }
  if(is.null(targ.rownames)) { targ.rownames <- colnames(SNP.data)[targ.cols] }
  result <- matrix(nrow=dm,ncol=dm)
  colnames(result) <- rownames(result) <- targ.rownames
  if(r2) { statz <- "R.squared" } else { statz <- "D.prime" }
  for (cc in 1:dm){
    for (dd in 1:dm) {
      result[cc,dd] <- ld(SNP.data[,targ.cols[cc]],SNP.data[,targ.cols[dd]],stats=statz)
    }
  }
  diag(result) <- NA
  return(result)
}



add.new.snp <- function(snp,snp.mat,case.code=NULL)
{
  # add a snp (in snpMatrix format) to an existing snpMatrix
  # at the same time, if case.code!=NULL, update the phenotypes
  # to reflect only samples common to both snpMatrix objects
  which.in.snp <- which(rownames(snp.mat) %in% rownames(snp))
  xx <- snp.mat[which.in.snp,]
  indx <- match(rownames(xx),rownames(snp))
  YY <- snp.cbind(xx,snp[indx,])
  if(is.null(case.code)) {
    return(YY)
  } else {
    case.code2 <- case.code[which.in.snp]  
    return(list(YY,case.code2))
  }
}


code.snps.to.binary.vars <- function(snp.mat,cod.sc=1,global=F)
{
  # Creating a separate variable for each SNP in the global environment,  
  # with coding:  0/1/2 for genotypes, NA for missing
  all.rs <- colnames(snp.mat)
  if(!global) { out.dat <- as.data.frame(snp.mat)}
  coding.schemes <- list(c(NA,0,0,1),c(NA,0,1,1),c(NA,0,1,0),c(NA,0,1,2))
  codez <- c(-1,0,1,2)
  for(cc in 1:length(all.rs)) {
    next.snp <- all.rs[cc]
    vec <- as.numeric(paste(snp.mat[,next.snp]))-1
    orig <- vec
    ncodes <- length(coding.schemes[[cod.sc]])
    for (dd in 1:ncodes) { vec[orig==codez[dd]] <- coding.schemes[[cod.sc]][dd]  }
    if(global) {
      assign(next.snp,vec,pos=1)
    } else {
      out.dat[,next.snp] <- vec
    }
  }
  if(!global) { return(out.dat) }
}

