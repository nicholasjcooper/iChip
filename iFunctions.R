if(Sys.info()[["user"]]=="ncooper")
{
  source("~/github/iChip/firstscriptFunctions.R")
  source("~/github/plumbCNV/geneticsFunctions.R")
}

require(snpStats)
require(reader)
require(genoset)


# stats function convenience wrappers
p.to.Z <- function(p) { O <- qnorm(1-(p/2),F); O[!is.finite(O)] <- NA; return(O) }
Z.to.p <- function(Z) { O <- 2*pnorm(-abs(Z)); O[!is.finite(O)] <- NA; return(O) }
l10 <- function(x) { O <- log10(x); O[!is.finite(O)] <- NA; return(O) }
Cor <- function(...) { cor(...,use="pairwise.complete") }
pt2 <- function(q, df, log.p=FALSE) {  2*pt(-abs(q), df, log.p=log.p) }



# function specific to getMetatable.R script
#  gets a list of equivalent SNPs to the current snp

get.equivs <- function(id,iden.list,na.fail=T) {
  id <- ic.to.rs(id) 
  do.one <- function(id,iden.list,na.fail) {
    if(id %in% ic.to.rs(unlist(iden.list))) {
      lnum <- which(sapply(iden.list,function(X) { id %in% ic.to.rs(X) }))
      idens <- ic.to.rs(iden.list[[lnum]])
      return(idens[-which(idens %in% id)])
    } else {
      if(na.fail) { return(NA) } else { return(id) }
    }
  }
  return(sapply(id,do.one,iden.list=iden.list,na.fail=na.fail))
}


# retrieve a simple representation of date_time or just date, for time-stamping file names, etc
simple.date <- function(sep="_",long=FALSE,time=TRUE) {
  myt <- format(Sys.time(), "%a %b %d %X %Y")
  if(long) {return(gsub(":",".",gsub(" ",sep,myt))) }
  dt <- strsplit(myt,":",fixed=TRUE)[[1]][1]
  splt <- strsplit(dt," ")[[1]]
  tm <- as.numeric(tail(splt,1))
  pr.tm <- head(splt,length(splt)-1)
  pr.tm[2] <- toupper(pr.tm[2])
  ampm <- {if(as.numeric(tm)>11) {"PM"} else {"AM"}}
  tm <- {if(tm>12) { tm-12 } else { if(tm<1) { tm+12 } else { tm } }}
  if(nchar(paste(pr.tm[3]))==1) { pr.tm[3] <- paste("0",pr.tm[3],sep="" ) }
  if(!time) { out <- paste(pr.tm[-1],collapse="") } else {
    out <- paste(paste(pr.tm[-1],collapse=""),sep,
               tm, ampm ,sep="") }
  return(out)
}

# for a subset 'n' and total 'N' nicely prints text n/N and/or percentage%
out.of <- function(n,N=100,digits=2,pc=TRUE,oo=TRUE) {
  pct <- 100*(n/N)
  outof <- paste(n,"/",N,sep="")
  percent <- paste(round(pct,digits),"%",sep="")
  if(pc & oo) {
    outof <- paste(outof," (",percent,")",sep="")
  } else {
    if(pc) { outof <- percent }
  }
  return(outof)
}

### GLOBAL QC STATS ###
samp.summ <- function(ms,CR=0.953,HZlo=0.19,HZhi=0.235) {
  nsamp <- nrow(ms)
  cr.filt <- ms$Call.rate>=CR
  hz.filt <- ms$Heterozygosity>=HZlo & ms$Heterozygosity<=HZhi
  sample.filt <- cr.filt & hz.filt
  cat(out.of(length(which(!cr.filt)),nsamp)," samples fail on call rate < ",CR,"\n",sep="")
  cat(out.of(length(which(!hz.filt)),nsamp)," samples fail on ",HZlo,"<Hz<",HZhi,"\n",sep="")
  cat(out.of(length(which(!sample.filt)),nsamp)," samples fail on call rate or HWE\n",sep="")
  return(sample.filt)
}

snp.summ <- function(MAF=0.005,CR=.99,HWE=3.8905,qc.file="snpqc.RData") {
  if(file.exists(qc.file)) { SNPQC <- reader(qc.file) } else { stop("couldn't find",qc.file) }
  mini.snp.qc <- function(SNPQC,MAF=0.005,CR=.99,HWE=3.8905) {
    maf <- SNPQC$MAF>MAF
    clr <- SNPQC$Call.rate>CR
    hwe <- abs(SNPQC$z.HWE)<HWE
    mono <- SNPQC$MAF==0 | SNPQC$MAF<0.0005
    maf.not.mono <- SNPQC$MAF>=0.0005 & SNPQC$MAF<MAF
    nsnp <- nrow(SNPQC)
    excl.rules <- rownames(SNPQC)[which(!maf | !clr | !hwe)]
    cat(out.of(length(which(!maf)),nsnp)," snps fail on MAF < ",MAF," ... \n",sep="")
    cat("... of which ",length(which(mono))," were monomorphic and ",length(which(maf.not.mono))," were just rare\n",sep="")
    cat(out.of(length(which(!clr)),nsnp)," snps fail on call rate < ",CR,"\n",sep="")
    cat(out.of(length(which(!hwe)),nsnp)," snps fail on HWE < ",HWE," [p=",round(2*(1-pnorm(HWE)),8),"]\n",sep="")
    cat(out.of(length(excl.rules),nsnp)," snps fail on MAF, Callrate or HWE \n",sep="")
    return(excl.rules)
  }
  nsnp <- nrow(SNPQC)
  Header("DIL")
  mono <- SNPQC$MAF==0 | SNPQC$MAF<0.0005
  #excl.rules <- mini.snp.qc(SNPQC[!mono,],MAF=MAF,CR=CR,HWE=HWE)
  excl.rules <- mini.snp.qc(SNPQC,MAF=MAF,CR=CR,HWE=HWE)
  num.excl.rules <- length(excl.rules)
  if(!file.exists("snpsExcluded.txt")) { warning("couldn't find exclusions file") } else {
    which.cut <- readLines(cat.path(getwd(),"snpsExcluded.txt"))
    total.excl <- length((unique(which.cut)))
    cat(out.of(total.excl-num.excl.rules,nsnp),"snps fail for other reasons\n")
    cat(out.of(total.excl,nsnp),"snps failed altogether\n")
    doc.path <- "ImChip_T1D_20130913"
    docs <- cat.path(doc.path,list.files(doc.path))
    if(file.exists(docs[2])) {  
      Header("UVA")
      excl.uva <- rownames(reader(docs[2])) 
      cat(out.of(length(excl.uva),nsnp),"snps failed according to UVA list\n")
      cat(out.of(length(excl.uva %in% which.cut),length(excl.uva),oo=F),"of these also failed according to DIL list\n")
      uva.lookup <- match(rs.to.ic(excl.uva),rs.to.ic(rownames(SNPQC)))
      print(length(which(is.na(uva.lookup))))
      excl.rules.uva <- mini.snp.qc(SNPQC[uva.lookup,],MAF=.005,CR=.95,HWE=4.89)
      num.excl.rules.uva <- length(excl.rules.uva)
      cat(out.of(length(excl.uva)-num.excl.rules.uva,length(excl.uva)),"UVA failing snps fail for other reasons\n")
    }
  }
  return(excl.rules)
}
####################

# with input as a glm model, returns a nice table of coefficients, p values, confidence intervals and SEs
mysumfun <- function(glmr,o.digits=3,p.digits=6,lab=TRUE,ci=FALSE)
{
  co <- summary(glmr)$coefficients
  predz <- rownames(co)[-1]
  label <- paste(summary(glmr)$call)[2]
  p <- co[2:nrow(co),4]; #print(p)
  o.r <- exp(co[2:nrow(co),1]); #print(o.r)
  p <- round(p,p.digits)
  
  # outlist <- list(round(o.r,o.digits),p)
  # names(outlist) <- c("OR","p-value")
  if(ci) {
    co1 <- co[2:nrow(co),1]-(1.96*co[2:nrow(co),2])
    co2 <- co[2:nrow(co),1]+(1.96*co[2:nrow(co),2])
    if(sum(co1)<=sum(co2)) { co.l <- co1; co.h <- co2 } else {  co.h <- co1; co.l <- co2 }
    co.l <- exp(co.l); co.h <- exp(co.h)
    out.fr <- cbind(round(o.r,o.digits),round(co.l,o.digits),round(co.h,o.digits),p)
    colnames(out.fr) <- c("OR","OR-low","OR-hi","p-value")
  } else {
    out.fr <- cbind(round(o.r,o.digits),p)
    colnames(out.fr) <- c("OR","p-value")
  }
  if(lab) { out.fr <- list(out.fr); names(out.fr) <- label }
  return(out.fr)
}

# flip odds ratios to always be positive (e.g, to change whether with respect to minor/major alleles)
or.pos <- function(X) { 
  x <- X; sel <- !is.na(X)
  X <- X[sel]
  X[X<1] <- 1/(X[X<1])
  x[sel] <- X
  return(x)
}

# flip odds ratios - pretty pointless really as it's just 1/x
or.flip <- function(X){
  x <- X; sel <- !is.na(X)
  X <- X[sel]
  X <- 1/(X)
  x[sel] <- X
  return(x)
}



# retrieve t1d previous SNP hits from t1dbase in either build hg18/hg19 coords (b36/37)
get.t1d.snps <- function(ucsc="hg19") {
  filenm <- "t1dhits.tab"
  cat("attempting to download t1d hits from t1dbase\n")
  url36 <- "http://www.t1dbase.org/webservice/RegionDownloads/=/model/variantsTAB/species=Human&disease_id=1&type=assoc&build=GRCh36"
  url37 <- "http://www.t1dbase.org/webservice/RegionDownloads/=/model/variantsTAB/species=Human&disease_id=1&type=assoc&build=GRCh37"
  urL <- switch(ucsc, hg18=url36,  hg19=url37)
  success <- T
  success <- tryCatch(download.file(urL ,filenm ,quiet=T),error=function(e) { F } )
  if(!is.logical(success)) { success <- T }
  if(success) {
    print("download successful")
    t1dh <- readLines(filenm,)
    firsts <- substr(t1dh,1,2)
    t1dh <- t1dh[firsts!="##"]
    len.lst <- strsplit(t1dh,"\t")
    rsids <- sapply(len.lst,"[",3)
    if(substr(rsids[1],1,2)!="rs") { rsids <- rsids[-1] }
    return(unique(rsids))
  } else {
    stop("couldn't reach t1dbase website at: ",urL)
  }
}

# create a phenotype vector for a dataframe with rownames that are subject ids, where
# cases and controls are text vectors of which IDs are that category; ctrls coded 0, cases 1
make.pheno <- function(X,cases,controls) {
  Pheno <- rep(NA,nrow(X)) # missing (default)
  Pheno[rownames(X) %in% cases] <- 1 # CASE if in the T1d dataset row (id) names
  Pheno[rownames(X) %in% controls] <- 0 # CONTROL if in the Controls dataset row (id) names
  return(Pheno)
}


# remove leading X from variable names (e.g, if original name started with a number and changed by make.names)
remove.X <- function(str) {
  bdz <- substr(str,1,1)
  str[bdz=="X"] <- substr(str,2,100000)[bdz=="X"]
  return(str)
}


# add X to the start of any string which has a digit as the first character
add.x <- function(str) {
  bdz <- substr(str,1,1)
  numy <-(paste(bdz) %in% paste(c(0:9)))
  str[numy] <- paste("X",str[numy],sep="") 
  return(str)
}



# convert from immunochip ids to rs-ids
ic.to.rs <- function(ic.ids) {
  ic.ids <- clean.snp.ids(ic.ids)
  if(!exists("all.support")) { print(load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
  outlist <- all.support$dbSNP[match(ic.ids,all.support$SNP)]
  outlist2 <- all.support$dbSNP[match(ic.ids,all.support$dbSNP)]
  outlist[is.na(outlist)] <- outlist2[is.na(outlist)]
  return(outlist)
}


# convert from rs-ids to immunochip ids
rs.to.ic <- function(rs.ids) {
  rs.ids <- clean.snp.ids(rs.ids)
  if(!exists("all.support")) { print(load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
  outlist <- all.support$SNP[match(rs.ids,all.support$dbSNP)]
  outlist2 <- all.support$SNP[match(rs.ids,all.support$SNP)]
  outlist[is.na(outlist)] <- outlist2[is.na(outlist)]
  return(outlist)
}



# for an immunochip or rs-id returns the chromosome it is a member of
Chr <- function(id) {
  ic.chr <- function(ic.ids) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { print(load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
    outlist <- all.support$Chr[match(ic.ids,all.support$SNP)]
    return(outlist)
  }
  ic <- ic.chr(rs.to.ic(id))
  return(ic)
}


# for an immunochip or rs-id returns the genome position (build = 36 or 37)
Pos <- function(id,build=36,warn.build=FALSE) {
  ic.pos <- function(ic.ids,build=36,warn.build=FALSE) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { print(load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
    if(build==37 & ("Pos37" %in% colnames(all.support))) {
      if(warn.build) { cat("Build 37/hg19 coordinates:\n") }
      outlist <- all.support$Pos37[match(ic.ids,all.support$SNP)]
    } else {
      if(warn.build) { cat("Build 36/hg18 coordinates:\n") }
      outlist <- all.support$Pos[match(ic.ids,all.support$SNP)]
    }
    return(outlist)
  }
  ic <- ic.pos(rs.to.ic(id),build=as.numeric(build[1]),warn.build=warn.build)
  return(ic)
}


# returns a table of the annotated allele 1 and allele 2 for a set of rs or ichip ids
AB <- function(id) {
  ic.ab <- function(ic.ids) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { print(load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
    outlist <- cbind(all.support$A1[match(ic.ids,all.support$SNP)],all.support$A2[match(ic.ids,all.support$SNP)])
    return(outlist)
  }
  out <- matrix(nrow=length(id),ncol=2)
  ic.ab.id <- ic.ab(rs.to.ic(id))
  return(ic.ab.id)
}


# this finds the full set of rows in a dataframe that share a duplicated column value (coln)
# i.e, differs to duplicated() which only returns 2..n duplicates and ignores the first members of the sets
# appends row numbers to allow to easily keep track of where these were in the original dataframe
dup.tracker <- function(X,col="dbSNP") {  
  keeps <- (X[(which(duplicated(X[[col]],fromLast=T))),])
  kills <- (X[(which(duplicated(X[[col]],fromLast=F))),])
  dupz <- rbind(keeps,kills)
  wkp <- which(duplicated(X[[col]],fromLast=F))
  wkl <- which(duplicated(X[[col]],fromLast=T))
  rn <- c(wkp,wkl)
  dupz[["rn"]] <- rn
  return(dupz[order(dupz[[col]]),])
}


## read in downloaded annotation from immunobase and return implied support file as dataframe
immunobase.support <- function(src="/chiswick/data/ncooper/iChipData/iChipAnnot.txt") {
  ici <- reader(src)
  ici[["rsid"]] <- rownames(ici)  
  ici[["A1"]] <- substr(ici[,"Alleles"],1,1)
  ici[["A2"]] <- substr(ici[,"Alleles"],3,3)
  info <- ici[,c("rsid","Immunochip.ID","Chromosome","Start","A1","A2")]
  colnames(info) <- c("rsID","iChipID","Chr","Pos","A1","A2")
  if(!anyDuplicated(info[,2])) { rownames(info) <- info[,2] } else { warning("ichip ids were not unique") }
  return(info)
}


# make sure support file naming and symbol convention match internal conventions, see clean.snp.ids
clean.snp.support <- function(X) {
  X[["dbSNP"]][is.na(X$dbSNP)] <- X$SNP[is.na(X$dbSNP)]
  X[["SNP"]] <- clean.snp.ids(X$SNP)
  X[["dbSNP"]] <- clean.snp.ids(X$dbSNP)
  rownames(X) <- clean.snp.ids(rownames(X))
  return(X)
}



# standardize snp ids so they would always appear the same, all _,.;, etc replaced with _
# all names leading with a number preceeded with X. mainly 'make.names' standard R-conventions
clean.snp.ids <- function(snpid.list) {
  snpid.list <- make.names(snpid.list)
  snpid.list <- gsub(".","_",snpid.list,fixed=T)
  return(snpid.list)
}


## convert from build 36 to build 37 coordinates. should work for GRanges or RangedData or Chr,Pos
# and should return modified GRanges, or otherwise a dataframe of rownames,Chr,Pos.
# works for regions or SNPs. if using GRanges, can get an output of a different length to input
# internals from Ollie Burren
conv.36.37 <- function(ranged=NULL,chr=NULL,pos=NULL,chain.file="/home/oliver/R/stuff/hg18ToHg19.over.chain") {
  require(GenomicRanges); require(rtracklayer); require(genoset)
  if(!file.exists(chain.file)) { stop("couldn't find ollie's script in: ",chain.file) }
  chn <- import.chain(chain.file)
  toranged <- F
  if(!is.null(chr) & !is.null(pos)) { ranged <- data.frame.to.ranged(cbind(chr,pos),start="pos",end="pos") }
  if(is(ranged)[1]=="RangedData") {
    wd <- width(ranged)
    if(all(wd==1)) { SNPs <- TRUE } else { SNPs <- FALSE }
    ranged[["XMYINDEXX"]] <- orn <- rownames(ranged)
    ranged[["XMYCHRXX"]] <- ocr <- chr(ranged)
    ranged <- set.chr.to.char(ranged)
    #print(head(ranged))
    ranged.gr <- as(ranged,"GRanges"); toranged <- T
  } else { 
    if(!is(ranged)[1]=="GRanges") { stop("need GRanges or RangedData 'ranged' object") } 
    wd <- width(ranged)
    if(all(wd==1)) { SNPs <- TRUE } else { SNPs <- FALSE }
  }
  ranged.gr.37<-liftOver(ranged.gr,chn)
  myfun <- function(x) { 
    data.frame(start=min(start(x),na.rm=T),end=max(end(x),na.rm=T)) 
  }
  if(!SNPs) {
    new.coords.df <- do.call("rbind",lapply(ranged.gr.37,myfun))
    ranged.gr.37<-ranged.gr
    ranges(ranged.gr.37)<-with(new.coords.df,IRanges(start=start,end=end))
    seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
  } else {
    seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
    #new.coords.df <- as.data.frame(ranged.gr.37)
  }
  if(!toranged) { return(ranged.gr.37) }
  ranged.gr.37 <- toGenomeOrder(as(as(ranged.gr.37,"IRangesList"),"RangedData"))
  if(all(c("XMYINDEXX","XMYCHRXX") %in% colnames(ranged.gr.37))) {
    RN <- ranged.gr.37[["XMYINDEXX"]]
    nr <- nrow(ranged.gr.37)
    if(length(orn)>length(RN)) { 
      cat("conversion failed for",length(orn[!orn %in% RN]),"rows, NAs produced\n") ;  
      print(head(orn[!orn %in% RN],20) )
      ln <- orn[!orn %in% RN]
      #return(ranged)
      newchr <- gsub("chr","",chr(ranged[match(ln,ranged$XMYINDEXX),]))
      extra <- data.frame(Chr=newchr,Pos=rep(NA,times=length(ln)))
      rownames(extra) <- ln
    }
    out <- data.frame(Chr=ranged.gr.37[["XMYCHRXX"]],Pos=start(ranged.gr.37))
    rownames(out) <- RN
    if(length(orn)>length(RN)) {
      out <- rbind(out,extra)
      out <- out[orn,]
    }
    return(out)
  } else { warning("missing key columns for chr, snp-name") }
  return(ranged.gr.37)
}



# convert snpStats:SnpMatrix object nicely to a dataframe where coding becomes 0,1,2,NA
SnpMatrix.to.data.frame <- function(SnpMat) {
  cov.data <- as.data.frame(SnpMat)
  for(jj in 1:ncol(cov.data)) { 
    nuxt <- as.numeric(cov.data[,jj])-1
    nuxt[nuxt<0] <- NA
    cov.data[,jj] <- nuxt
    # assign(colnames(cov.data)[jj], nuxt)
  }
  return(cov.data)
}


# convert from a dataframe to a SnpMatrix
data.frame.to.SnpMatrix <- function(X){
  mxx <- max(X,na.rm=TRUE)
  if(mxx>3) { warning("Dataframe does not appear to contain allele codes") }
  X <- round(X)
  if(mxx==3) { X <- X-1 ; X[X<0] <- NA }
  NN <- as.matrix(X)
  #NN <- round(NN)
  SS <- as(NN,"SnpMatrix")
  return(SS)
}


## for a snp dataset, determine whether reference allele is the major or minor allele 
# Where 2 copies coded as highest value = reference, e.g, if AA=0, AB=1, BB=2, then B is reference.
# Combines this with frequencies of the alleles to evaluate whether 'BB' is major or minor
# takes SnpMatrix or data.frame coded 0,1,2,NA as input
majmin <- function(X,checks=TRUE) {
  ## workhorse internal function ##
  do.mm <- function(x,snpmat=NULL) { 
    if(!is.null(snpmat)) { 
      if(!snpmat) { 
        tt <- table(round(as.numeric(x))) 
      }
    }
    if(is.null(snpmat)) {  
      tt <- table(round(as.numeric(x)))
      if("3" %in% names(tt)) { snpmat <- T } else { snpmat <- F }
    }
    if(snpmat) {
      xx <- round(as.numeric(x))-1; xx[xx==-1] <- NA
      tt <- table(xx)
    }
    if(checks) {
      # checks will slow down running so can optionally turn them off
      if(length(tt)<1) { warning("all allele counts were empty!"); return(NA) }
      if(!all(names(tt) %in% paste(c(0,1,2)))) { warning("invalid genotype values in X: ",paste(head(unique(names(tt))),collapse=",")) ; return(NA)  }
    }
    type <- "unknown"
    if("0" %in% names(tt)) { A0 <- tt[["0"]] } else { A0 <- 0 }
    if("1" %in% names(tt)) { A1 <- tt[["1"]] } else { A1 <- 0 }
    if("2" %in% names(tt)) { A2 <- tt[["2"]] } else { A2 <- 0 } 
    a0 <- (A0*2)+A1; a2 <- (A2*2)+A1
    if(a0 > a2) { type <- "minor" }
    if(a2 > a0) { type <- "major" }
    if(a0==a2) { type <- "neutral"  }
    if( length(which(c(A0,A1,A2)==0))==2 ) { type <- "monomorph" }
    return(type)
  }
  ## main code ##
  if(length(Dim(X))!=2) { 
    if(length(Dim(X))==1) { return(do.mm(as.numeric(X))) } else {
      warning("invalid object for major/minor allele testing"); return(NA)
    }
  }
  snpmat <- F
  if(is(X)[1] %in% "SnpMatrix") { snpmat <- T } else {
    tt1 <- table(round(as.numeric(X[,1])));  if("3" %in% names(tt1)) { snpmat <- T }
  }
  all.typz <- apply(X,2,do.mm,snpmat=snpmat)
  return(factor(all.typz))
}


## function that calculates which way around the data is for case vs control (pheno) for a GWAS test
# will indicate with respect to cases whether they have more reference allele, or less, 
# or also if the heterozygous is the affected genotype
# works on a SnpMatrix or dataframe coded 0,1,2,NA (autodetects which)
caseway <- function(X, pheno, checks=TRUE, long=FALSE) {
  ## workhorse internal function ##
  do.cw <- function(x,ph,snpmat=NULL) { 
    if(!is.null(snpmat)) { 
      if(!snpmat) { 
        tt <- table(round(as.numeric(x)),ph) 
      }
    }
    if(is.null(snpmat)) {  
      tt <- table(round(as.numeric(x)),ph)
      if("3" %in% rownames(tt)) { snpmat <- T } else { snpmat <- F }
    }
    if(snpmat) {
      xx <- round(as.numeric(x))-1; xx[xx==-1] <- NA
      tt <- table(xx,ph)
    }
    #tt <- narm(tt)
    if(checks) {
      # checks will slow down running so can optionally turn them off
      if(length(tt)<1) { warning("all allele counts were empty!"); return(NA) }
      if(!all(rownames(tt) %in% paste(c(0,1,2)))) { warning("invalid genotype values in X: ",paste(head(unique(names(tt))),collapse=",")) ; return(NA)  }
      if(!all(colnames(tt) %in% paste(c(0,1,2)))) { warning("invalid genotype values in X: ",paste(head(unique(names(tt))),collapse=",")) ; return(NA)  }
    }
    if("0" %in% rownames(tt)) { Ctrl0 <- tt["0","0"]; Case0 <- tt["0","1"] } else { Ctrl0 <- Case0 <- 0 }
    if("1" %in% rownames(tt)) { Ctrl1 <- tt["1","0"]; Case1 <- tt["1","1"] } else { Ctrl1 <- Case1 <- 0 }
    if("2" %in% rownames(tt)) { Ctrl2 <- tt["2","0"]; Case2 <- tt["2","1"] } else { Ctrl2 <- Case2 <- 0 } 
    Ctrl0[is.na(Ctrl0)] <- 0; Ctrl1[is.na(Ctrl1)] <- 0; Ctrl2[is.na(Ctrl2)] <- 0
    Case0[is.na(Case0)] <- 0; Case1[is.na(Case1)] <- 0; Case2[is.na(Case2)] <- 0
    ctrl.pc0 <- Ctrl0/sum(Ctrl0,Ctrl1,Ctrl2); case.pc0 <- Case0/sum(Case0,Case1,Case2)
    ctrl.pc1 <- Ctrl1/sum(Ctrl0,Ctrl1,Ctrl2); case.pc1 <- Case1/sum(Case0,Case1,Case2)
    ctrl.pc2 <- Ctrl2/sum(Ctrl0,Ctrl1,Ctrl2); case.pc2 <- Case2/sum(Case0,Case1,Case2)
    #prv(ctrl.pc0,ctrl.pc1,ctrl.pc2,case.pc0,case.pc1,case.pc2)
    if(long) { res <- "unclear results" } else { res <- "???" }
    if(long) { r1 <- "cases have more 1, less 0,2" } else { r1 <- "CasesHet+" }
    if(long) { r2 <- "cases have less 1, more 0,2" } else { r2 <- "CasesHet-" }
    if(long) { r3 <- "cases have more 0, less 2" } else { r3 <- "CasesRef-" }
    if(long) { r4 <- "cases have more 2, less 0" } else { r4 <- "CasesRef+" }
    if((case.pc0 < ctrl.pc0) & (case.pc2 < ctrl.pc2)) { res <-r1  }
    if((case.pc0 > ctrl.pc0) & (case.pc2 > ctrl.pc2)) { res <-r2  }
    if((case.pc0 > ctrl.pc0) & (case.pc2 < ctrl.pc2)) { res <-r3  }
    if((case.pc0 < ctrl.pc0) & (case.pc2 > ctrl.pc2)) { res <-r4  }
    if((case.pc0 == ctrl.pc0) & (case.pc2 < ctrl.pc2)) { res <-r3  }
    if((case.pc0 == ctrl.pc0) & (case.pc2 > ctrl.pc2)) { res <-r4  }
    if((case.pc0 > ctrl.pc0) & (case.pc2 == ctrl.pc2)) { res <-r3  }
    if((case.pc0 < ctrl.pc0) & (case.pc2 == ctrl.pc2)) { res <-r4  }
    return(res)
  }
  ## main code ##
  if(!max(Dim(pheno)) %in% Dim(X)) { warning("Phenotype data different size to dataset X"); return(NA)}
  if(all(pheno %in% c(1,2))) { pheno <- pheno-1 }
  if(!all(pheno %in% c(0,1))) { warning("Phenotype must be coded as controls,cases=0,1; or =1,2"); return(NA) }
  if(length(Dim(X))!=2) { 
    if(length(Dim(X))==1) { return(do.cw(as.numeric(X),ph=pheno)) } else {
      warning("invalid object for case/control orientation testing"); return(NA)
    }
  }
  snpmat <- F
  if(is(X)[1] %in% "SnpMatrix") { snpmat <- T } else {
    tt.temp <- table(round(as.numeric(X[,1])));  if("3" %in% names(tt.temp)) { snpmat <- T }
  }
  all.res <- apply(X,2,do.cw,ph=pheno,snpmat=snpmat)
  return(factor(all.res))
}



##' log sum function @author Claudia Giambartolomei
logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}



# 
# if(any(colnames(regions.gr.37) %in% "XMYINDEXX")) {
#   rownames(regions.gr.37) <- RN <- regions.gr.37[["XMYINDEXX"]]
#   regions.gr.37 <- regions.gr.37[,-which(colnames(regions.gr.37) %in% "XMYINDEXX")]
# } else { warning("not sure what happened missing index column") ; return(regions.gr.37) }
# if(any(colnames(regions.gr.37) %in% "XMYCHRXX")) {
#   out <- as.data.frame(regions.gr.37)
#   rownames(out) <- RN
#   out <- out[,c("XMYCHRXX","start")]
#   colnames(out) <- c("Chr","Pos")
#   return(out)
# }



## functions from Chris W - adapted


# create a factor that can split a group of 'size' entries into categories size 'by'
# merge last two groups if final group size is less than min.pc of 'by'
# if fac.out is false, return start/end ranges rather than a grouping factor
# (Nick)
make.split <- function(size,by,fac.out=T,min.pc=0.5) {
  stepz <- round(seq(from=1,to=(size+1),by=by))
  if((tail(stepz,1)) != (size+1)) { stepz <- c(stepz,(size+1)) }
  split.to <- length(stepz)-1
  out1 <- cbind(stepz[1:split.to],(stepz[2:(split.to+1)]-1))
  repz <- 1+apply(out1,1,diff)
  out <- rep(1:split.to,repz)
  lrg <- max(out)
  if(!length(which(out==lrg))>(by*min.pc)) {
    out[out==lrg] <- lrg-1
  }
  if(!fac.out) {
    out <- cbind(start(Rle(out)),end(Rle(out)))
  }
  return(out)
}


# snpStats imputation only works if there are correlated SNPs with non-missing values
# that can be used to interpolate missing SNPs. If any correlated SNPs are missing
# 'impute.missing' will leave these blank. This function mops up the remainder
# by randomly inserting values consistent with the minor allele frequency of each SNP
# (Nick)
randomize.missing <- function(X) {
  miss1 <- function(x) { 
    TX <- table(x)
    naz <- which(is.na(x))
    if(length(naz)>0 & length(TX)>0) {
      x[naz] <- sample(as.numeric(names(TX)),size=length(naz),
                       replace=T,prob=as.numeric(TX))
    }
    return(x)
  }
  # randomly generate replacements for missing values using current distribution for each column of X
  if(is.null(dim(X))) { warning("not a matrix/data.frame") ; return(X) }
  count.miss <- function(x) { length(which(is.na(x))) }
  nmiss <- apply(X,2,count.miss)
  FF <- nrow(X)
  select <- nmiss>0
  if(length(which(select))<1) { return(X) }
  if(length(which(select))==1) { X[,select] <- miss1(X[,select]); return(X) }
  X[,select] <- apply(X[,select],2,miss1)
  return(X)
}


## wrapper for the snpStats::snp.imputation function
# allows for stratified imputation, and the parameter 'by' allows
# speeding up of imputation for large sets by running the imputation
# in smaller batches (as large ones are really slow)
# (chris)
impute.missing <- function (X, bp = 1:ncol(X), strata = NULL, numeric = FALSE, verbose=FALSE, by=NULL, ...) {
  N <- as(X, "numeric")
  if(any(Dim(X)!=Dim(N))) { stop("'as numeric' lost desired dimensionality") }
  if (!is.null(strata)) {
    strata <- as.factor(strata)
    if (length(levels(strata)) > 10) 
      stop("too many levels in strata\n")
    for (i in levels(strata)) {
      cat("\nstrata", i, "\n")
      wh <- which(strata == i)
      N[wh, ] <- impute.missing(X[wh, , drop = FALSE], 
                                bp, numeric = TRUE, ...)
    }
  }
  if (is.null(strata) & is.numeric(by)) {
    strata <- as.factor(make.split(ncol(X),by=by))
    # added by nick
    for (i in levels(strata)) {
      cat("\nsplit", i, "\n")
      wh <- which(strata == i)
      N[, wh] <- impute.missing(X[,wh , drop = FALSE], 
                                bp, numeric = TRUE, ...)
    }
  }
  else {
    csumm <- col.summary(X)
    use <- csumm[, "Certain.calls"] == 1
    X2 <- X[, use]
    bp <- bp[use]
    imp <- (csumm[, "Call.rate"] < 1 & !is.na(csumm[, "Call.rate"]))[use]
    cat(sum(imp,na.rm=T), "to impute\n")
    mx <- max(which(imp),na.rm=T)
    for (i in which(imp)) {
      loop.tracker(i,mx)
      supres <- capture.output(rule <- snp.imputation(X2[, -i, drop = FALSE], X2[, 
           i, drop = FALSE], bp[-i], bp[i]))
      if (is.null(rule@.Data[[1]])) 
        next
      imp <- impute.snps(rules = rule, snps = X2[, rule@.Data[[1]]$snps, drop = FALSE], ...) 
      wh.na <- which(is.na(N[, i]))
      N[wh.na, colnames(X2)[i]] <- imp[wh.na]
    }
    cat("\n")
  }
  if (numeric) {
    return(as.data.frame(N))
  }
  else {
    print(Dim(N))
    print(Dim(X@snps))
    return(new("aSnpMatrix", .Data = new("SnpMatrix", data = (round(N) + 1),
                                         , nrow = nrow(N), ncol = ncol(N), dimnames = dimnames(N)), 
               snps = X@snps, samples = X@samples, phenotype = X@phenotype, 
               alleles = X@alleles))
  }
}

## imputation subscript ## 
# Nick
# used by conditional analysis and indistinguishable analyses scripts, does the imputation,
# saving imputed results as we go so don't need to keep recalculating the same SNPs
# will automatically proceed in the most efficient way possible
ichip.imputation <- function(myData, imp.file, smp.filt=NULL, snp.filt=NULL) {
  if(is.null(smp.filt)) { smp.filt <- 1:nrow(myData) }
  if(is.null(snp.filt)) { snp.filt <- 1:ncol(myData) }
  if(!file.exists(imp.file)) {
    myDataFilt <- myData[smp.filt,snp.filt]
    myDat <- impute.missing(myDataFilt,numeric=T)
    myDat <- randomize.missing(myDat)
    myDatSnp <- data.frame.to.SnpMatrix(myDat)
    save(myDat,myDatSnp,cov.dat,file=imp.file)
    cat("wrote imputed to:",imp.file,"\n")
  } else { 
    cat("loaded imputed data from:",imp.file,"\n")
    print(load(imp.file))
    ### note this row/col check is not 100% foolproof!
    if(nrow(myDat)!=length(smp.filt)) { stop("mismatching number of samples in loaded file") }
    if(ncol(myDat)!=length(snp.filt)) { 
      ## if this has changed a bit, try to ressurrect without recalculating the whole thing
      cat("mismatching number of snps in loaded file\n") 
      if(!exists("bigDat")) { print(load("allImputed.RData")) }
      targs <- colnames(myData)[snp.filt]
      gotem <- narm(match(targs,colnames(bigDat)))
      aintgotem <- targs[!targs %in% colnames(bigDat)]
      if(length(aintgotem)>0) {
        if(length(gotem)>0) {
          cat("combining",length(gotem),"previously imputed with",length(aintgotem),"from scratch\n")
          myDat.part1 <- impute.missing(myData[smp.filt,match(aintgotem,colnames(myData))],numeric=T)
          if(!exists("bigDat")) { print(load("allImputed.RData")) }
          myDat.part2 <- bigDat[,gotem]
          if(any(rownames(myDat.part1)!=rownames(myDat.part2))) { stop("samples did not match up") }
          myDat.cbind <- cbind(myDat.part1,myDat.part2)
          indz <- match(targs,colnames(myDat.cbind))
          if(any(is.na(indz))) { stop("combined data still missing target SNPs") }
          myDat <- myDat.cbind[,indz]
        } else {
          ## do all from scratch
          cat("none were previously imputed, calculating from scratch\n")
          myDataFilt <- myData[smp.filt,snp.filt]
          myDat <- impute.missing(myDataFilt,numeric=T)
        }
      } else { 
        ## have all the snps, but need to prune some
        cat("trimming loaded data to subset needed\n")
        indz <- match(targs,colnames(bigDat))
        if(any(is.na(indz))) { stop("not sure why loaded data is missing target SNPs") }
        myDat <- bigDat[,indz]
      }
      myDatSnp <- data.frame.to.SnpMatrix(myDat)
      myDat <- impute.missing(myDatSnp,numeric=T)
      myDat <- randomize.missing(myDat)
      myDatSnp <- data.frame.to.SnpMatrix(myDat)
      save(myDat,myDatSnp,cov.dat,file=imp.file)
      cat("wrote imputed to:",imp.file,"\n")
    } else { 
      myDat <- randomize.missing(myDat)
      myDatSnp <- data.frame.to.SnpMatrix(myDat)
      save(myDat,myDatSnp,cov.dat,file=imp.file)
      cat("loaded file is a good match to requested SNP-set\n")
    }
  }
  return(imp.file)
}




##'Calculate approximate Bayes factors from p values and MAF
##'
##' this is a function to calculate approximate Bayes factors from p
##' values and MAF - for reference see Wakefield, J (2009) Bayes
##' factors for genome-wide association studies: comparison with
##' p-values. Genetic Epidemiology 33: 79–86.
##' @title abf
##' @param p p value
##' @param maf minor allele frequency
##' @param n0 number of controls
##' @param n1 number of cases
##' @param scale0 by default, =n0
##' @param scale1 by default, =n1
##' @return ABF
##' @export
##' @author Chris Wallace
abf <- function(p,maf, n0=9500, n1=6670, scale0=n0, scale1=n1) { 
  # evidence for null - ie low ABF => support for alternative
  z <- qnorm(p/2, lower.tail=FALSE)
  x0 <- 0; x1 <- 1; x2 <- 2 # multiplicative model
  d2 <- (1-maf)^2 * x0^2 + 2*maf*(1-maf)*x1 + maf^2 * x2^2
  d1 <- (1-maf)^2 * x0 + 2*maf*(1-maf)*x1 + maf^2 * x2
  V <- (n0 + n1) / ( n0 * n1 * (d2-d1) )
  ## scale
  scale <- ((n0 + n1)/(scale0 + scale1)) * (scale0/n0) * (scale1/n1)
  V <- V * scale
  W <- ( log(1.5)/qnorm( 0.99, lower.tail=FALSE) )^2
  VW <- V+W
  2 * log(sqrt(VW/V) * exp( - z^2 * W / (2 * VW) ))
}


# fix.rownames <- function(data) {
#   nms <- rownames(data)
#   nms[substr(nms,1,2)=="NA"] <- nms[substr(nms,2,100000)=="NA"] # there are NAxxxx names that should be Axxxx
#   nms[nms=="R58201702_C02"] <- "58201702_C02" # this name has a preceeding R in one place, not in another
#   rownames(data) <- nms.
#   return(data)
# }

## Given a region, add window cM either side and return.
## Everything is HapMap v2, release 22, build 36.

## Options to plot the rates and window exist:
##   do.plot=TRUE -> produce plot
##   add.plot=TRUE -> add lines to existing plot
##   do.lines=TRUE -> use lines to indicate positions of window
## NB, add.plot and do.lines are ignored, unless do.plot=TRUE.
## function from David.  Used because when running this on the queue,
## sometimes need multiple connections to open the damn files.
multitry <- function(expr, times=5, silent=FALSE, message=""){
  warn <- options()$warn
  options(warn=-1)
  for (i in 1:times) {
    res <- try(expr, silent=TRUE)
    if (inherits(res, "try-error")){
      if (i==times) {
        options(show.error.messages = TRUE)
        stop(geterrmessage(), " (", times, " failed attempts)", message)
      }
      next
    }
    break
  }
  options(warn=warn)
  if (!silent && (i>1)) {
    warning(geterrmessage(), " (", times-1, " failed attempts)")
  }
  res
}


# chris' function to get a centimorgan window from intervals
# vector input
recwindow <- function(chr,st,en=st,window=0.1, # cM either side
                      do.plot=FALSE, # if wanted to plot
                      add.plot=FALSE,do.lines=TRUE,...) {
  rates <- read.table(gzfile(sprintf("/dunwich/scratch/chrisw/HapMap/rates_rel22/genetic_map_chr%s_b36.txt.gz",chr)),
                      header=TRUE)
  cm.st <- rates[which.min(abs(rates$position-st)),3]
  cm.en <- rates[which.min(abs(rates$position-en)),3]
  
  mx <- max(window,1)
  kk <- rates[which.min(abs(rates[,3]-(cm.st-window))) : which.min(abs(rates[,3]-(cm.en+window))),]
  cat("n hapmap snps in window =",nrow(kk),"\n")
  from <- min(kk[,1])
  to <- max(kk[,1])
  
  if(do.plot) {
    kk <- rates[abs(rates[,3]-cm.st)<mx | abs(rates[,3]-cm.en)<mx,]
    if(add.plot) {
      lines(kk[,1:2])
    } else {
      plot(kk[,1:2],type="l",main=paste("Recombination rates on chr",chr),
           xlab="chromosome position (bp)",ylab="rec rate (cM/Mb)",...)
    }
    if(window>0 & do.lines) {
      abline(v=c(from,to),col="red")
      abline(v=c(st,en),col="blue")
      legend("topleft",lty=c(1,1),col=c("red","blue"),legend=c("window","target"))
    }
  }
  
  cat("window size is\nleft: ",(st-from)/1000,"kb\tright: ",(to-en)/1000,"kb\ttotal: ",(to-from)/1000,"kb\n",sep="")
  return(c(from,to))
}



### this is used by getMetaTable1.R to get the top snp in each region-group ###
highlights <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { cat("more than one top hit: ",paste(wh,collapse=","),"\n")  }
  # if(length(wh)==2) {  wh <- wh[-1] } ## wh[-1] is hack to adjust for 1p13.2, although don't know why i bothered as it doesn't carry through
  next.row <- c(wh[1],top,length(X),length(which(X<bonf)))
  return(next.row) 
}

multihit.check <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { return(T) } else { return(F) }
}

multihit.return <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  return(wh)
}


## this function gets marginal likelihoods from the BIC, similar to how do.bic.max works on all.results
do.bic.marg <- function(X,dif=3) {
  if(!"BIC" %in% names(X)) { return(NA) }
  bic <- X$BIC;
  if(is.list(bic)) { bic <- unlist(bic) }
  #which.min <- (which(bic==min(bic)))
  denom <- logsum(-bic/2)
  #logsum((-bic/2)[-which.min])
  next.marg <- next.exp <- numeric(length(bic))
  for (cc in 1:length(bic)) {
    which.min <- which(bic==sort(bic)[cc])[1]
    next.min <- (-bic/2)[which.min]
    next.marg[cc] <- next.min - denom
    next.exp[cc] <- (exp(next.marg[cc])); #print(next.exp[cc])
  }
  names(next.marg) <- names(next.exp) <- names(bic)[order(bic)]
  #print(paste("SNP:",names(bic)[top.grp]))
  #lapply(X$GLM[top.grp],function(Z) { print(tail(Z,1)[4]) })
  return(list(MAR=next.marg,EXP=next.exp))
}


## this is the function that extracts 'BIC' from the all.results list and converts to relative Bayes factors
# against the lowest BIC in the list
do.bic.max <- function(X,dif=3) {
  if(!"BIC" %in% names(X)) { return(NA) }
  bic <- X$BIC;
  if(is.list(bic)) { bic <- unlist(bic) }
  max.bic <- -min(bic/2)
  dif.bic <- (-bic/2)-max.bic
  return(rev(sort(dif.bic)))
}


# extract the BIC information from 'all.results' list which is returned by the indistinguishableAnalysis.R file
suck.bic <- function(X,dif=3) {
  bic <- X$BIC;
}

# not used in the end?
clearly.suck <- function(X,thresh=100) {
  if(!"BIC" %in% names(X)) { return(NA) }
  bic <- X$BIC;
  if(is.list(bic)) { bic <- unlist(bic) }
  max.bic <- -min(bic/2)
  dif.bic <- (-bic/2)-max.bic
  failers <- (rev(sort(dif.bic)))
  failers <- failers[abs(failers)>=thresh]
  return(names(failers))
}





## this is very specific ##
plot.one.chr <- function(result,chr=1,new.plot=TRUE,text=TRUE,x.off=0,line.col="blue",pt.col="red",text.pos=1,label.offset=3, text.cex=1) {
  wone <- result
  tab <- sapply(wone[[chr]],"[",13)  # or no number
  #sorter <- order(tab[1,])
  sorter <- order(tab)
  wone[[chr]] <- wone[[chr]][sorter]
  n.in.chr <- length(wone[[chr]])
  xo <- x.off
  l1 <- numeric(n.in.chr)
  if(new.plot) {
    plot(c(0,(n.in.chr+1)),c(0,2.5),xlab="marker(s)",ylab="Odds Ratio",col="white",main=paste("Chr",chr),bty="l",xaxt="n")
    abline(h=1,col="grey",lty="dashed",lwd=1.25)
  }
  for (jj in (1:n.in.chr)) {
    l1[jj] <- y <- tail(wone[[chr]][[jj]][,1],1)
    yl <- tail(wone[[chr]][[jj]][,2],1)
    yh <- tail(wone[[chr]][[jj]][,3],1)
    symbols(jj+xo, y, circles= 0.05,add=T,fg="black",inches=F,bg=pt.col)
    arrows(jj+xo,y,jj+xo,yl,angle=90,length=0.05)
    arrows(jj+xo,y,jj+xo,yh,angle=90,length=0.05)
    if(text) {
      text(jj,yl,labels=names(wone[[chr]])[jj],pos=text.pos,cex=text.cex,offset=label.offset)
    }
  }
  lines((1:n.in.chr)+xo,l1,col=line.col)
  return()
}




conditional <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { cat("more than one top hit: ",paste(wh,collapse=","),"\n") }
  next.row <- c(wh,top,length(X),length(which(X<bonf)))
  return(next.row) 
}




data.frame.to.ranged <- function(dat,ids=NULL,start="start",end="end",width=NULL,
                                 chr="chr",exclude=NULL,ucsc="hg18") 
{
  ## convert any data frame with chr,start,end, or pos data into a RangedData object
  # not case sensitive
  ## abandon longer names as they clash with function names
  st <- paste(start); en <- paste(end); ch <- paste(chr); wd <- paste(width)
  must.use.package(c("genoset","IRanges"),T)
  if(is.matrix(dat)) { dat <- as.data.frame(dat) }
  if(!is.data.frame(dat)) { stop("Error: not a dataframe")}
  key.nms <- c(ids,st,en,ch,wd)
  tries <- 0
  while(!all(key.nms %in% colnames(dat))) { 
    colnames(dat) <- tolower(colnames(dat)); key.nms <- tolower(key.nms)
    st <- tolower(st); en <- tolower(en); ch <- tolower(ch); wd <- tolower(wd)
    if(tries>1) {
      if((tolower(st)=="pos" | tolower(en)=="pos") & !(tolower(st)=="pos" & tolower(en)=="pos")) {
        st <- en <- "pos"
      } else {
        if(tolower(st)=="start" & tolower(en)=="end") { st <- en <- "pos" }
      }
    }
    key.nms <- c(ids,st,en,ch,wd)
    tries <- tries+1
    if(tries > 3) { warning("columns not found") ; break }
  }
  if(!is.null(ids)) { 
    if(anyDuplicated(dat[[ids]])==0) { 
      id <- dat[[ids]] 
    } else { 
      key.nms <- key.nms[-match(ids,key.nms)] # allow non-unique ids as regular
      ids <- NULL
      warning("id must be unique to form rownames, will insert as a separate column") 
    }
  }
  if(is.null(ids)) { 
    if(!is.null(rownames(dat)) & all(rownames(dat)!=paste(1:nrow(dat)))) { 
      id <- rownames(dat)
    } else { 
      id <- paste(1:nrow(dat)) 
    }
  }
  if(length(ch)>0) { ch <- gsub("chr","",dat[[ch]],ignore.case=T) } else { ch <- NULL }
  if(length(st)>0) { st <- dat[[st]] } else { st <- NULL }
  if(length(en)>0) { en <- dat[[en]] } else { en <- NULL }
  if(length(wd)>0) { wd <- dat[[wd]] } else { wd <- NULL }
  outData <- RangedData(ranges=IRanges(start=st,end=en,names=id),space=ch,universe=ucsc[1])
  outData <- toGenomeOrder(outData,strict=T)
  # note when adding data subsequently that 'RangedData' sorts by genome order, so need
  # to resort any new data before adding.
  if(is.null(rownames(outData))) { rownames(outData) <- paste(1:nrow(outData)) }
  reorder <- match(rownames(outData),id)
  more.cols <- colnames(dat)[!colnames(dat) %in% key.nms]
  more.cols <- more.cols[!more.cols %in% exclude]
  if(length(more.cols)>0) {
    for (cc in 1:length(more.cols)) {
      outData[[more.cols[cc]]] <- dat[[more.cols[cc]]][reorder]
    }
  }
  return(outData)
}




