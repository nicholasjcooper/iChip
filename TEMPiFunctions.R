if(Sys.info()[["user"]]=="ncooper")
{
  source("~/github/iChip/firstscriptFunctions.R") # only needed for internal analyses
  source("~/github/plumbCNV/generalCNVFunctions.R")
}

require(snpStats)
require(reader)
require(genoset)


##file includes the generally useful functions: simple.date, out.of, randomize.missing



# stats function convenience wrappers
p.to.Z <- function(p) { O <- qnorm(1-(p/2),F); O[!is.finite(O)] <- NA; return(O) }
Z.to.p <- function(Z) { O <- 2*pnorm(-abs(Z)); O[!is.finite(O)] <- NA; return(O) }
l10 <- function(x) { O <- log10(x); O[!is.finite(O)] <- NA; return(O) }
Cor <- function(...) { cor(...,use="pairwise.complete") }
pt2 <- function(q, df, log.p=FALSE) {  2*pt(-abs(q), df, log.p=log.p) }


if(F) {
  LL <- length(which(!is.na(qchisq(p.cc,1)))) ; plot(qchisq(1-((1:LL)/LL),1),qchisq(1-sort(p.cc),1),xlim=c(0,6),type="l",main="Case-control Analysis (ChiSq)",xlab="expected",ylab="observed"); lines(x=qchisq(1-((1:LL)/LL),1),y=qchisq(1-((1:LL)/LL),1),col="red",lty="dotted")
  ww <- (which(!is.na(-log10(p.cc)))); LL <- length(ww)
  plot(-log10(((1:LL)/LL)),-log10(sort(p.cc)),xlim=c(0,6),type="l",main="Case-control Analysis -Log10(p)",xlab="expected",ylab="observed"); lines(x=-log10(((1:LL)/LL)),y=-log10(((1:LL)/LL)),col="red",lty="dotted")
  
  do.the.meta.ones <- function(p.meta) {
    LL <- length(which(!is.na(qchisq(p.meta,1)))) ; plot(qchisq(1-((1:LL)/LL),1),qchisq(1-sort(p.meta),1),xlim=c(0,6),type="l",main="Meta Analysis (ChiSq)",xlab="expected",ylab="observed"); lines(x=qchisq(1-((1:LL)/LL),1),y=qchisq(1-((1:LL)/LL),1),col="red",lty="dotted")
    ww <- (which(!is.na(-log10(p.meta)))); LL <- length(ww)
    plot(-log10(((1:LL)/LL)),-log10(sort(p.meta)),xlim=c(0,6),type="l",main="Meta Analysis -Log10(p)",xlab="expected",ylab="observed"); lines(x=-log10(((1:LL)/LL)),y=-log10(((1:LL)/LL)),col="red",lty="dotted")
  }
  
  print(load("/chiswick/data/ncooper/iChipData/compiledTableAllResultsPassingQC.RData"))
  near.region <- grep("EXT",tt$gene)
  out.region <- grep("OTHER",tt$gene)
  in.region <- which(!1:nrow(tt) %in% c(near.region,out.region))
  
  print(load("pvaluesforqq.RData"))
  p.m.reg <- p.meta[names(p.meta) %in% rs.to.ic(tt$names[in.region])]
  p.m.out <- p.meta[names(p.meta) %in% rs.to.ic(tt$names[out.region])]
  do.the.meta.ones(p.m.reg)
  do.the.meta.ones(p.m.out)
}



pduplicated <- function(X) {
  if(length(Dim(X))>1) {  stop("can only enter a vector into this function") }
  return((duplicated(X,fromLast=T) | duplicated(X,fromLast=F)))
}


comma <- function(...) {
  paste(...,collapse=",")
}


## the following are internal functions for the 'lambdas' function below

get.allele.counts <- function(myData,cc1000=FALSE) {
  ii <-  col.summary(myData)
  ii[["majmin"]] <- c("minor","major")[as.numeric(round(ii$RAF,3)!=round(ii$MAF,3))+1]
  if(cc1000) { ii$Calls <- rep(1000,nrow(ii)) }
  aa <- aA <- AA <- rep(0,nrow(ii))
  aa[which(ii$majmin=="minor")] <- (ii$P.BB*ii$Calls)[which(ii$majmin=="minor")]
  aa[which(ii$majmin=="major")] <- (ii$P.AA*ii$Calls)[which(ii$majmin=="major")]
  AA[which(ii$majmin=="minor")] <- (ii$P.AA*ii$Calls)[which(ii$majmin=="minor")] 
  AA[which(ii$majmin=="major")] <- (ii$P.BB*ii$Calls)[which(ii$majmin=="major")]
  aA <- ii$P.AB*ii$Calls
  ii[["aa"]] <- aa
  ii[["aA"]] <- aA
  ii[["AA"]] <- AA
  colnames(ii)[1] <- "TOTAL"
  return(ii[c("aa","aA","AA","TOTAL")])
}


# convert inflation factors to lamba1000
lambda_nm <- function(Lnm,n=1000,m=1000,nr,mr) { 1 + ((Lnm-1)*(((1/nr)+(1/mr))/((1/n)+(1/m)))) }
# likelihood for one marker
Likelihood_Lj <- function(c,Lnm) { rchisq(c/Lnm,df=1)/Lnm }
# total likelihood across all K markers  : http://www.nature.com/ng/journal/v36/n4/full/ng1333.html

LLikelihood_L <- function(Cj,LNMj) {
  tot <- 0
  for (cc in 1:length(Cj)) { 
    tot <- tot + log(Likelihood_Lj(Cj[cc],LNMj[cc])) 
  }
  return(tot) 
}

Y_2 <- function(r1,r2,n1,n2,N,R) { (N*((N*(r1+(2*r2)))-(R*(n1+(2*n2))))^2) / ((R*(N-R))*((N*(n1+(4*n2)))-(n1+(2*n2))^2)) }

X_2 <- function(r1,r2,n1,n2,N,R) { (2*N*((2*N*(r1+(2*r2)))-(R*(n1+(2*n2))))^2) / ((4*R*(N-R))*((2*N*(n1+(2*n2)))-(n1+(2*n2))^2)) }

#Y2 ~ L*X2

#Case     r0  r1  r2  R
#Control  s0  s1  s2  S
#Total    n0  n1  n2  N

#http://en.wikipedia.org/wiki/Population_stratification


## function that calculates SNP-WISE Lambda and Lambda1000 statistics for inflation due to population structure
# works on a SnpMatrix or dataframe coded 0,1,2,NA (autodetects which)
lambdas <- function(X, pheno, checks=TRUE, cc1000=FALSE) {
  ## workhorse internal function ##
  do.lambda <- function(x,ph,snpmat=NULL) {
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
      if(!all(rownames(tt) %in% paste(c(0,1,2)))) { 
         warning("invalid genotype values in X: ",paste(head(unique(names(tt))),collapse=",")) ; return(NA) 
      }
      if(!all(colnames(tt) %in% paste(c(0,1,2)))) {
         warning("invalid genotype values in X: ",paste(head(unique(names(tt))),collapse=",")) ; return(NA)  
      }
    }

    if("0" %in% rownames(tt)) { Ctrl0 <- tt["0","0"]; Case0 <- tt["0","1"] } else { Ctrl0 <- Case0 <- 0 }
    if("1" %in% rownames(tt)) { Ctrl1 <- tt["1","0"]; Case1 <- tt["1","1"] } else { Ctrl1 <- Case1 <- 0 }
    if("2" %in% rownames(tt)) { Ctrl2 <- tt["2","0"]; Case2 <- tt["2","1"] } else { Ctrl2 <- Case2 <- 0 }
    Ctrl0[is.na(Ctrl0)] <- 0; Ctrl1[is.na(Ctrl1)] <- 0; Ctrl2[is.na(Ctrl2)] <- 0
    Case0[is.na(Case0)] <- 0; Case1[is.na(Case1)] <- 0; Case2[is.na(Case2)] <- 0

    A0 <- Case0+Ctrl0; A1 <- Case1+Ctrl1; A2 <- Case2+Ctrl2
    a0 <- (A0*2)+A1; a2 <- (A2*2)+A1
    if(a0 > a2) { type <- "minor" }
    if(a2 > a0) { type <- "major" }
    if(a0==a2) { type <- "neutral"  }
    if( length(which(c(A0,A1,A2)==0))==2 ) { type <- "monomorph" }
    #cat(type,"\n")
    if(type=="minor") {
      #Case
      r0 <- Case2  ; r1 <- Case1  ; r2 <- Case0  ; R <- Case0+Case1+Case2
      #Control  s0  s1  s2  S
      s0 <- Ctrl2  ; s1 <- Ctrl1  ; s2 <- Ctrl0  ; S <- Ctrl0+Ctrl1+Ctrl2
    } else {
      #Case
      r0 <- Case2  ; r1 <- Case1  ; r2 <- Case0  ; R <- Case0+Case1+Case2
      #Control  s0  s1  s2  S
      s0 <- Ctrl2  ; s1 <- Ctrl1  ; s2 <- Ctrl0  ; S <- Ctrl0+Ctrl1+Ctrl2
    }
    return(c(r0,r1,r2,R,s0,s1,s2,S))  #,n0,n1,n2,N,xx2,yy2,LL,L1000))
    #return(c(LL,L1000))
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
  if(snpmat) {
    cnts0 <- get.allele.counts(X[pheno==0,],cc1000=cc1000)
    cnts1 <- get.allele.counts(X[pheno==1,],cc1000=cc1000)
    cnts <- cbind(cnts1,cnts0)
  } else {
    cnts <- apply(X,2,do.lambda,ph=pheno,snpmat=snpmat)
    cnts <- t(cnts)
  }
  #Total
  colnames(cnts) <- c("r0","r1","r2","R","s0","s1","s2","S")
  if(cc1000) { cnts <- round(cnts) }
  tryCatch(detach("cnts"),error=function(e) {NULL}) # in case erroneously attached from earlier
  attach(cnts)
  n0 <- r0 + s0 ; n1 <- r1 + s1 ; n2 <- r2 + s2 ; N <- R + S
  xx2 <- X_2(r1,r2,n1,n2,N,R)
  yy2 <- Y_2(r1,r2,n1,n2,N,R)
  LL <- yy2/xx2
  L1000 <- lambda_nm(Lnm=LL,n=1000,m=1000,nr=R,mr=S)
  detach(cnts)
  all.res <- cbind(cnts,xx2,yy2,LL,L1000)
  colnames(all.res)[9:12] <- c("X2","Y2","Lambda","L1000")
  return(all.res)
}


# calculate and print the overall lambdas reflecting inflation
calculate.overall.lambdas <- function(surround=FALSE) {
  print(load("lambdaTableVars.RData")) #"c2"      "p.m.out" "p.m.reg"
  # c2 is a table of results including armitage trend test 'Y2' column, calc derived using lambda.R
  work.dir <- "/chiswick/data/ncooper/iChipData/"
  print(load("/chiswick/data/ncooper/iChipData/finalMetaTopHitsFEB17.RData"))
  # ^file contains objects from getMetaTable1.R
  # "bonfs.filt"     "non.bonfs.filt" "bonf.snps"      "non.bonf.snps" 
  print(load("/chiswick/data/ncooper/iChipData/compiledTableAllResultsPassingQC.RData"))
  # "tt"
  print(load("pvaluesforqq.RData"))
  # "p.meta" "p.cc"
  near.region <- grep("EXT",tt$gene)
  out.region <- grep("OTHER",tt$gene)
  in.region <- which(!1:nrow(tt) %in% c(near.region,out.region))
  potential.novel.region <- which(1:nrow(tt) %in% c(near.region,out.region))
  p.m.reg <- names(p.meta[names(p.meta) %in% rs.to.ic(tt$names[in.region])])
  p.m.out <- names(p.meta[names(p.meta) %in% rs.to.ic(tt$names[out.region])])
  p.m.nov <- names(p.meta[names(p.meta) %in% rs.to.ic(tt$names[potential.novel.region])])
  
  if(surround) {
    t1dsnps <- calibrate.cond.bonf(topsnplist,cm.window=0.2,bp.ext=0,build=37,qclist="snpsExcluded.txt",ret.snps=T)
  }
  t1region <- tt$gene %in% rownames(bonfs.filt)
  t1region.novel <- (tt$gene %in% rownames(bonfs.filt)) & (rs.to.ic(tt$names) %in% rs.to.ic(p.m.nov))
  t1regsnps <- rs.to.ic(tt$names[t1region])
  t1regsnps.novel <- rs.to.ic(tt$names[t1region.novel])
  p.m.out.not.t1 <- p.m.out[!p.m.out %in% t1regsnps.novel]
  cat("there were",length(t1regsnps.novel),"SNPs in novel regions\n")
  if(surround) { cat("raw lambda for the cM regions surrounding t1dhits:",median(c2[t1dsnps,"Y2"],na.rm=T)/.456,"\n") }
  cat("raw lambda for the regions containing t1dhits:",median(c2[t1regsnps,"Y2"],na.rm=T)/.456,"\n")
  lam0 <- median(c2[,"Y2"],na.rm=T)/.456
  cat("raw lambda for all SNPs passing QC:",lam0,"\n")  
  lam1000.0 <- lambda_nm(lam0,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T))
  cat("lambda_1000 for all SNPs passing QC:",lam1000.0,"\n")
  if(surround) { lam1 <- median(c2[-which(rownames(c2) %in% t1dsnps),"Y2"],na.rm=T)/.456 }
  lam2 <- median(c2[-which(rownames(c2) %in% t1regsnps),"Y2"],na.rm=T)/.456
  if(surround) { cat("raw lambda excluding the 0.2 cM regions surrounding t1dhits:",lam1,"\n") }
  cat("raw lambda excluding the dense regions containing t1dhits or other table-1 novel regions:",lam1,"\n")
  if(surround) { lam1000.1 <- lambda_nm(lam1,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T)) }
  lam1000.2 <- lambda_nm(lam2,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T))
  if(surround) { cat("lambda_1000 excluding the cM regions surrounding t1dhits:",lam1000.1,"\n") }
  cat("lambda_1000 excluding the dense regions containing t1dhits or other table-1 novel regions:",lam1000.2,"\n")
  lam3 <- median(c2[p.m.out,"Y2"],na.rm=T)/.456
  lam4 <- median(c2[p.m.reg,"Y2"],na.rm=T)/.456
  lam5 <- median(c2[p.m.out.not.t1,"Y2"],na.rm=T)/.456
  cat("raw lambda for outside ichip dense regions:",lam3,"\n")
  cat("raw lambda for inside ichip dense regions:",lam4,"\n")
  cat("raw lambda for outside ichip dense regions and excluding table1 novel regions:",lam5,"\n")
  lam1000.3 <- lambda_nm(lam3,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T))
  lam1000.4 <- lambda_nm(lam4,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T))
  lam1000.5 <- lambda_nm(lam5,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T))
  cat("lambda_1000 for outside ichip dense regions:",lam1000.3,"\n")
  cat("lambda_1000 for inside ichip dense regions:",lam1000.4,"\n")
  cat("lambda_1000 for outside ichip dense regions and excluding table1 novel regions:",lam1000.5,"\n")
  
  #L = 1.641179
  #L1000 = 1.074309
  do.the.meta.ones <- function(p.meta,suf="") {
    pdf(cat.path("","chiqq",ext="pdf",suf=suf))
    #prv(p.meta)
    #p.meta <- narm(p.meta)
    txt <- paste("Meta Analysis (ChiSq) [",toheader(paste(gsub("."," ",suf,fixed=T),"regions")),"]",sep="")
    LL <- length(which(!is.na(qchisq(p.meta,1))))
    xx <- qchisq(1-((1:LL)/LL),1); yy <- qchisq(1-sort(p.meta),1)
    cond <- (is.na(xx) | is.na(yy) | !is.finite(yy) | !is.finite(xx))
    xx <- xx[!cond]; yy <- yy[!cond]
    three4 <- function(x) { mean(c(min(x,na.rm=T), rep(max(x,na.rm=T),3))) }
    plot(xx,yy,type="l",main=txt,xlab="expected",ylab="observed",xlim=c(0,20))
    text(three4(xx),three4(yy),paste("slope =",round(coefficients(lm(yy~xx))[2],3)))
    lines(x=qchisq(1-((1:LL)/LL),1),y=qchisq(1-((1:LL)/LL),1),col="red",lty="dotted")
    dev.off()
  }
  
  # QQ for dense versus non dense regions #
  p.meta.reg <- p.meta[names(p.meta) %in% rs.to.ic(tt$names[in.region])]
  #p.meta.out <- p.meta[names(p.meta) %in% rs.to.ic(tt$names[out.region])]
  p.meta.out <- p.meta[names(p.meta) %in% rs.to.ic(p.m.out.not.t1)]
  do.the.meta.ones(p.meta.reg,"inside.dense")
  do.the.meta.ones(p.meta.out,"outside.dense")
  # QQ for t1d versus non dense regions #
  p.meta.reg <- p.meta[names(p.meta) %in% rs.to.ic(t1regsnps)]
  p.meta.out <- p.meta[!names(p.meta) %in% rs.to.ic(t1regsnps)]
  do.the.meta.ones(p.meta.reg,"inside.t1d")
  do.the.meta.ones(p.meta.out,"outside.t1d")
  
  return(c(lam1000.1,lam1000.2,lam1000.3,lam1000.4,lam1000.5))
}


#send cw: bonferonnis, numbers of snps


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



### GLOBAL QC STATS ###
samp.summ <- function(ms,CR=0.953,HZlo=0.19,HZhi=0.235,by.pheno=FALSE) {
  if(by.pheno) {
    if(exists("control.data") & exists("t1d.data")) {
      ph <- make.pheno(ms,rownames(t1d.data),rownames(control.data))
      cat("Sample QC for CASES (phenotype=1)\n")
      sample.filt1 <- samp.summ(ms=ms[ph==1,],CR=CR,HZlo=HZlo,HZhi=HZhi,by.pheno=FALSE)
      cat("Sample QC for CONTROLS (phenotype=0)\n")
      sample.filt0 <- samp.summ(ms=ms[ph==0,],CR=CR,HZlo=HZlo,HZhi=HZhi,by.pheno=FALSE)
      return(c(sample.filt0,sample.filt1))
    }
  }
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
    cat(out.of(length(which(!clr & maf)),nsnp)," snps fail on call rate < ",CR,"\n",sep="")
    cat(out.of(length(which(!hwe & clr & maf)),nsnp)," snps fail on HWE < ",HWE," [p=",round(2*(1-pnorm(HWE)),8),"]\n",sep="")
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


# flip odds ratios 
# true/false vector of whether the odds ratio should be greater than 1
# GWAS convention is that ORs>1 when cases have more of the minor allele than controls
or.flip <- function(X,greater.than.1=rep(T,length(X))) {
  gt1 <- greater.than.1 # save space
  if(!is.numeric(X)) { stop("X should be numeric (ie, odds ratios)") }
  if(!is.logical(gt1) | (length(gt1)!=length(X) & length(gt1)!=1) ) { 
    stop("'greater.than.1' should be logical (ie, TRUE/FALSE for each member of X)") }
  xx <- X; sel <- !is.na(gt1) & !is.na(X)
  X <- X[sel]; gtone <- gt1[sel]
  gt1s <- X>1
  lt1s <- X<1
  X[gt1s & !gtone] <- 1/X[gt1s & !gtone]
  X[lt1s & gtone] <- 1/X[lt1s & gtone]
  xx[sel] <- X
  return(xx)
}



# to extract info from 22 separate chromosomes
# pheno <- make.pheno(rawdata,rownames(t1d.data),rownames(control.data)) # assume all are same length
# for (cc in 1:22) {
#   
#   load(paste("temp.ichip-data",cc,".RData",sep=""))
#   rawdata <- rbind(control.data,t1d.data)
#   tt <- fix.OR.directions(tt, OR.col=c(4,6,8), snp.data=rawdata, pheno=pheno, partial=T,verbose=F)
#   loop.tracker(cc,22)
# }
### specific to ichip paper analysis
fix.OR.directions <- function(results, OR.col=1, snp.data=NULL, pheno=NULL, partial=FALSE,
                              case.list=NULL, control.list=NULL, effect.labels=TRUE, verbose=TRUE, alt.return=FALSE) {
  nameslist <- colnames(snp.data)
  if(length(OR.col)>0) {
    if(is.character(OR.col)) {
      if(any(!OR.col %in% colnames(results))) { stop("OR.col must be column names or numbers of 'results'") } else { do.or.col <- T }
    }
    if(is.numeric(OR.col)) {
      if(any(OR.col > ncol(results))) { stop("OR.col was greater than the number of columns in 'results'") } else { do.or.col <- T }
    } 
  } else { effect.labels <- TRUE }
  rr_in_ss <- rownames(results) %in% nameslist
  if(!all(rr_in_ss)) { 
    if(!partial) { warning(length(which(!rr_in_ss))," rownames from results not found in colnames of snp.data") };
    incomplete <- T  
  } else { incomplete <- F }
  if(!any(rr_in_ss)) { stop("no rownames from results were found in colnames of snp.data") }
  if(is.null(pheno) & !is.character(case.list) & is.character(control.list)) {  
    pheno <- make.pheno(snp.data,case.list,control.list)
  }
  if(length(pheno)!=nrow(snp.data)) { 
    stop("pheno (or case.list+control.list) must be ",
         "the same length as the number of rows in snp.data") }
  indx <- match(rownames(results)[rr_in_ss],nameslist)
  majminlist <- majmin(snp.data[,indx])
  reflist <- caseway(snp.data[,indx],pheno)
  if(alt.return) { return(list(majminlist=majminlist,reflist=reflist,
                               nameslist=nameslist,indx=indx,rr_in_ss=rr_in_ss)) }
  # get logical vectors
  l1 <- (majminlist=="minor" & reflist=="CasesRef-")
  l2 <- (majminlist=="minor" & reflist=="CasesRef+")
  l3 <- (majminlist=="major" & reflist=="CasesRef-")
  l4 <- (majminlist=="major" & reflist=="CasesRef+")
  OR.dir <- rep(as.logical(NA),length(indx)); 
  caseshave <- rep("nodiff",length(indx))
  # assign non-default values based on logical vectors above
  OR.dir[l1] <- FALSE; caseshave[l1] <- "MoreMajor"
  OR.dir[l2] <- TRUE;  caseshave[l2] <- "MoreMinor"
  OR.dir[l3] <- TRUE;  caseshave[l3] <- "MoreMinor"
  OR.dir[l4] <- FALSE; caseshave[l4] <- "MoreMajor"
  if(do.or.col) {
    for(cc in 1:length(OR.col)) { 
      results[rr_in_ss,OR.col[cc]] <- or.flip(results[rr_in_ss,OR.col[cc]],OR.dir)
      if(verbose) { cat("inverted any inconsistent odds-ratios in column",OR.col[cc],"\n") }
    }
  }
  if(effect.labels) {
   # if(incomplete) { results[["CasesHave"]][!rr_in_ss] <- NA }
    results[["CasesHave"]][rr_in_ss] <- caseshave
  }
  return(results)
}

  
# to extract info from 22 separate chromosomes
#  pheno <- make.pheno(rawdata,rownames(t1d.data),rownames(control.data)) # assume all are same length
#  for (cc in 1:22) {
#    load(paste("temp.ichip-data",cc,".RData",sep=""))
#    rawdata <- rbind(control.data,t1d.data)
#    tt2 <- add.allele.to.result(tt2, OR.col=8, snp.data=rawdata, pheno=pheno, partial=T)
#    loop.tracker(cc,22)
#  }

add.allele.to.result <- function(results, OR.col=1, snp.data=NULL, pheno=NULL, partial=FALSE,
                              case.list=NULL, control.list=NULL) {
  lll <- fix.OR.directions(results, OR.col=OR.col[1], snp.data=snp.data, pheno=pheno, partial=partial,verbose=FALSE,
                          case.list=case.list, control.list=control.list, effect.labels=FALSE,alt.return=TRUE)
  list.to.env(lll) # gets: majminlist,reflist,nameslist,indx,rr_in_ss
  al.names <- c("allele.A","allele.B")
  if(all(al.names %in% colnames(results))) {
    tt <- results
  } else {
    tt <- cbind(results,AB(rownames(results)))
    colnames(tt)[(ncol(tt)-c(1,0))] <- al.names
  }
  ## T>A   #  cases h ave more T, less A
  ## T is ref because ref is allele.B
  ## if cases ref+ then  allele.B >  allele.A
  ## if cases ref- then  allele.A >  allele.B
  ## else allele.A ~ allele.B
  
  # get logical vectors
  l1 <- (majminlist=="minor")
  l2 <- (majminlist=="major")
  relship <- paste((tt$allele.A[rr_in_ss]),"~",(tt$allele.B[rr_in_ss]))
  relship[l1] <- paste((tt$allele.A[rr_in_ss]),">",(tt$allele.B[rr_in_ss]))[l1]
  relship[l2] <- paste((tt$allele.B[rr_in_ss]),">",(tt$allele.A[rr_in_ss]))[l2]
  if(!"effect" %in% colnames(tt)) { tt[["effect"]] <- NA }
  #prv(reflist,l1,l2,relship,indx,rr_in_ss)
  tt[rr_in_ss,"effect"] <- relship
  return(tt)
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
# or can also enter X as a vector of all ids, e.g, X = rownames(someDataFrame)
make.pheno <- function(X,cases,controls) {
  if(length(Dim(X))==1) {
    Pheno <- rep(NA,length(X)) # missing (default)
    Pheno[X %in% cases] <- 1 # CASE if in the T1d dataset row (id) names
    Pheno[X %in% controls] <- 0 # CONTROL if in the Controls dataset row (id) names
  } else {
    Pheno <- rep(NA,nrow(X)) # missing (default)
    Pheno[rownames(X) %in% cases] <- 1 # CASE if in the T1d dataset row (id) names
    Pheno[rownames(X) %in% controls] <- 0 # CONTROL if in the Controls dataset row (id) names
  }  
  return(Pheno)
}


# remove leading X from variable names (e.g, if original name started with a number and changed by make.names)
remove.X <- function(str,char="X") {
  bdz <- substr(str,1,1)
  str[bdz==char] <- substr(str,2,100000)[bdz==char]
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
ic.to.rs <- function(ic.ids,dir=NULL) {
  if(!exists("work.dir")) { if(is.null(dir)) { work.dir <- getwd() } else { work.dir <- dir } }
  ic.ids <- clean.snp.ids(ic.ids)
  if(!exists("all.support")) { (load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
  outlist <- all.support$dbSNP[match(ic.ids,all.support$SNP)]
  outlist2 <- all.support$dbSNP[match(ic.ids,all.support$dbSNP)]
  outlist[is.na(outlist)] <- outlist2[is.na(outlist)]
  return(outlist)
}


# convert from rs-ids to immunochip ids
rs.to.ic <- function(rs.ids,dir=NULL) {
  if(!exists("work.dir")) { if(is.null(dir)) { work.dir <- getwd() } else { work.dir <- dir } }
  rs.ids <- clean.snp.ids(rs.ids)
  if(!exists("all.support")) { (load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
  outlist <- all.support$SNP[match(rs.ids,all.support$dbSNP)]
  outlist2 <- all.support$SNP[match(rs.ids,all.support$SNP)]
  outlist[is.na(outlist)] <- outlist2[is.na(outlist)]
  return(outlist)
}



# for an immunochip or rs-id returns the chromosome it is a member of
Chr <- function(id,dir) {
  if(!exists("work.dir")) { if(is.null(dir)) { work.dir <- getwd() } else { work.dir <- dir } }
  ic.chr <- function(ic.ids) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { (load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
    outlist <- all.support$Chr[match(ic.ids,all.support$SNP)]
    return(outlist)
  }
  ic <- ic.chr(rs.to.ic(id))
  return(ic)
}


# for an immunochip or rs-id returns the genome position (build = 36 or 37)
Pos <- function(id,build=36,warn.build=FALSE,dir=NULL,snps.only=FALSE) {
  if(!exists("work.dir")) { if(is.null(dir)) { work.dir <- getwd() } else { work.dir <- dir } }
  ic.pos <- function(ic.ids,build=36,warn.build=FALSE) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { (load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
    if(build==37 & ("Pos37" %in% colnames(all.support))) {
      if(warn.build) { cat("Build 37/hg19 coordinates:\n") }
      outlist <- all.support$Pos37[match(ic.ids,all.support$SNP)]
    } else {
      if(warn.build) { cat("Build 36/hg18 coordinates:\n") }
      outlist <- all.support$Pos[match(ic.ids,all.support$SNP)]
    }
    return(outlist)
  }
  query <- rs.to.ic(id)
  if(!snps.only & all(is.na(query))) { 
    ## unless the 'snps.only' function is set, then if it looks like we have not been handed
    ## snp ids, then check for band ids or gene ids instead
    numpqs <- (length(grep("q",id))+length(grep("p",id)))
    if(numpqs==length(id)) { try.band <- T } else { try.band <- F }
    if(try.band) {
      suppressWarnings(test <- Pos.band(id,dir=dir,build=build,warn.build=warn.build))
      if(!is.null(test)) { return(test) }
    }
    suppressWarnings(test <- Pos.gene(id,dir=dir,build=build,warn.build=warn.build))
    if(!is.null(test)) { return(test) } 
  } 
  ic <- ic.pos(query,build=as.numeric(build[1]),warn.build=warn.build)
  return(ic)
}

# for an immunochip or rs-id returns the genome position (build = 36 or 37)
# bioC - whether to return ranges or dataframe
# band - whether to include band/stripe in returned object
# dir - putting a dir will speedup future lookups as annotation will be save here
# unique - some genes have split ranges, this merges to give only 1 range per gene
Pos.gene <- function(genes,build=36,warn.build=FALSE,dir=NULL,bioC=FALSE,band=FALSE,unique=TRUE,map.cox.to.6=TRUE) {
  ucsc <- ucsc.sanitizer(build)
  char.lim <- 100
  ga <- get.gene.annot(dir=dir,ucsc=ucsc,range.out=bioC,unique=unique,map.cox.to.6=map.cox.to.6)
  if(warn.build) {
    if(build!="hg18") {    cat("Build 37/hg19 coordinates:\n") 
    } else {    cat("Build 36/hg18 coordinates:\n") }
  }
  mt <- match(genes,ga$gene)
  failz <- paste(genes[is.na(mt)],collapse=", "); if(nchar(failz)>char.lim) { failz <- paste(substr(failz,1,char.lim),",...",sep="") }
  if(length(mt)<1 | all(is.na(mt))) { 
    warning("did not find any 'genes' features: ",failz) ; return(NULL) }
  if(any(is.na(mt))) { 
    cnt <- length(which(is.na(mt)))
    warning("did not find the following ",cnt," 'genes' features: ",failz) 
  }
  outlist <- ga[sort(mt[!is.na(mt)]),]
  if(!band) { outlist <- outlist[,-which(colnames(outlist) %in% "band")] }
  if(unique & ("gene" %in% colnames(outlist))) {
    rownames(outlist) <- outlist[["gene"]]
    outlist <- outlist[,-which(colnames(outlist) %in% "gene")]
    if(all(genes %in% rownames(outlist))) {
      outlist <- outlist[genes,]
    }
  }
  return(outlist)
}

Pos.band <- function(bands,build=36,warn.build=FALSE,dir=NULL,bioC=FALSE) {
  ucsc <- ucsc.sanitizer(build)
  char.lim <- 100
  ga <- get.cyto(ucsc=ucsc,bioC=bioC,dir=dir)
  if(warn.build) {
    if(build!="hg18") {    cat("Build 37/hg19 coordinates:\n") 
    } else {    cat("Build 36/hg18 coordinates:\n") }
  }
  mt <- match(bands,rownames(ga))
  failz <- paste(bands[is.na(mt)],collapse=", "); if(nchar(failz)>char.lim) { failz <- paste(substr(failz,1,char.lim),",...",sep="") }
  msg <- ("format for bands is: chromosome[p/q]xx.xx ; e.g, 13q21.31, Yq11.221, 6p23, etc")
  if(length(mt)<1 | all(is.na(mt))) { 
    warning("did not find any 'bands' features: ",failz) ; warning(msg); return(NULL) }
  if(any(is.na(mt))) { 
    cat("format for bands is: chromosome[p/q]xx.xx ; e.g, 13q21.31, Xq27.1, 6p23, etc")
    cnt <- length(which(is.na(mt)))
    warning("did not find the following ",cnt," 'bands' features: ",failz,"...") ; warning(msg)
  }
  outlist <- ga[sort(mt[!is.na(mt)]),]
  if(any(colnames(outlist) %in% "negpos")) { outlist <- outlist[,-which(colnames(outlist) %in% "negpos")] }
  if(all(bands %in% rownames(outlist)) & !bioC) {
    outlist <- outlist[bands,]
  }
  return(outlist)
}

Band <- function(genes=NULL,chr=NULL,ranges=NULL,build=36,dir=getwd(),...) {
  if(all(is.character(genes))) { 
    Band.gene(genes=genes,build=build,dir=dir,...)
  } else {
    Band.pos(chr=chr,ranges=ranges,build=build,dir=dir,...)
  }
}

Band.gene <- function(genes,build=36,dir=getwd(),append.chr=TRUE,data.frame=FALSE) {
  pg <- Pos.gene(genes,build=build,warn.build=FALSE,dir=dir,bioC=F,band=TRUE,unique=TRUE)
  if(data.frame) {
    if(all(c("start","end") %in% colnames(pg))) { pg <- pg[,-which(colnames(pg) %in% c("start","end"))] }
    if(all(c("gene") %in% colnames(pg))) { rownames(pg) <- pg[["gene"]] ; pg <- pg[,-which(colnames(pg) %in% c("gene"))] }
    out <- pg
  } else {
    if(append.chr) {
      out <- paste(pg[["chr"]],pg[["band"]],sep="")
    } else {
      out <- pg[["band"]]
    }
  }
  return(out)    
}



# will return a vector if concat=T, else a list (allowing for multiple hits per range)
# if bioC=T, if concat=F, will return a RangedData object which will have a column 'rangeindex' showing
# for each row, which query location index the band/pos information originated (to allow
# for multiple bands for ranges). Otherwise if concat=T, will just concatenate multi hits into single rows
# to match dimension of the query set
# can input chr, pos, or chr,start,end, or ranges=RangedData
# concat will return bands separated by semi colons when bioC=FALSE if there are multiple in 1 range
Gene.pos <- function(chr,pos=NA,start=NA,end=NA,ranges=NULL,build=36,dir=NULL,bioC=FALSE,concat=TRUE) {
  ucsc <- ucsc.sanitizer(build)
  if(is(ranges)[1]!="RangedData") {
    if(any(!is.na(pos))) { start <- pos; end <- pos }
    if(length(chr)==1 & length(start)>1) { chr <- rep(chr,times=length(start)) }
    if(length(chr)!=length(start)) { stop("chr vector must have same length as pos or start/end") }
    if(any(is.na(chr))) { stop("cannot have chr=NA") }
    Pos <- matrix(ncol=2,nrow=length(start))
    for (cc in 1:length(start)) {
      Pos[cc,] <- force.chr.pos(Pos=c(start[cc],end[cc]),Chr=chr[cc],dir=dir,ucsc=ucsc)
    }
    #if(any(tolower(substr(chr,1,3))!="chr")) { chr <- gsub("chr",chr,sep="") }
    #chr <- gsub("chrchr","chr",chr)
    testData <- RangedData(ranges=IRanges(start=Pos[,1],end=Pos[,2]),space=chr,index=1:length(chr),universe=ucsc[1])
    testData <- toGenomeOrder2(testData,strict=T)
  } else {
    testData <- ranges # set.chr.to.char(ranges)
    if("index" %in% colnames(testData)) { warning("'index' is a reserved column name for ranges objects passed to this function so will be replaced. Consider renaming this column if this is undesired") }
    testData[["index"]] <- 1:nrow(testData)
  }
  testData <- set.chr.to.numeric(testData,keep=T)
  #return(testData)
  ga <- get.gene.annot(ucsc=ucsc,dir=dir)
  #ga <- set.chr.to.numeric(ga,keep=F)
  newDataList <- vector("list",nrow(testData))
  overlaps <- findOverlaps(testData,ga)
  genez <- ga$gene[subjectHits(overlaps)]
  indexz <- queryHits(overlaps)
  #prv(indexz,genez)
  if(length(indexz)<1) { return(NA) }
  if(!concat) {
    newData <- testData[queryHits(overlaps),]
    newData[["gnm.index"]] <- indexz
    newData[["gene"]] <- genez
  } else {
    out <- tapply(genez,factor(indexz),c,simplify=FALSE)
    out <- sapply(out,function(X) { paste(X,collapse=";") })
    newData <- testData
    newData[["gene"]] <- rep("intergenic",nrow(newData))
    if(!is.null(names(out))) {
      newData[["gene"]][as.numeric(names(out))] <- out
    } else {
      newData[["gene"]] <- out
    }
  }
  if(bioC) {
    return(newData)
  } else {
    return(newData[["gene"]][order(newData[["index"]])])
  }
  #if(!all(chr %in% chr2(ga))) { stop("invalid chromosome(s) entered") } # redundant i think
}



# if bioC=F
# will return a vector if concat=T, else a list (allowing for multiple hits per range)
# if bioC=T, if concat=F, will return a RangedData object which will have a column 'rangeindex' showing
# for each row, which query location index the band/pos information originated (to allow
# for multiple bands for ranges). Otherwise if concat=T, will just concatenate multi hits into single rows
# to match dimension of the query set
# can input chr, pos, or chr,start,end, or ranges=RangedData
# concat will return bands separated by semi colons when bioC=FALSE if there are multiple in 1 range
Band.pos <- function(chr,pos=NA,start=NA,end=NA,ranges=NULL,build=36,dir=NULL,bioC=FALSE,concat=TRUE) {
  ucsc <- ucsc.sanitizer(build)
  if(is(ranges)[1]!="RangedData") {
    if(any(!is.na(pos))) { start <- pos; end <- pos }
    if(length(chr)==1 & length(start)>1) { chr <- rep(chr,times=length(start)) }
    if(length(chr)!=length(start)) { stop("chr vector must have same length as pos or start/end") }
    if(any(is.na(chr))) { stop("cannot have chr=NA") }
    Pos <- matrix(ncol=2,nrow=length(start))
    for (cc in 1:length(start)) {
      Pos[cc,] <- force.chr.pos(Pos=c(start[cc],end[cc]),Chr=chr[cc],dir=dir,ucsc=ucsc)
    }
    #if(any(tolower(substr(chr,1,3))!="chr")) { chr <- gsub("chr",chr,sep="") }
    #chr <- gsub("chrchr","chr",chr)
    testData <- RangedData(ranges=IRanges(start=Pos[,1],end=Pos[,2]),space=chr,index=1:length(chr),universe=ucsc[1])
    testData <- toGenomeOrder2(testData,strict=T)
  } else {
    testData <- ranges # set.chr.to.char(ranges)
    if("index" %in% colnames(testData)) { warning("'index' is a reserved column name for ranges objects passed to this function so will be replaced. Consider renaming this column if this is undesired") }
    testData[["index"]] <- 1:nrow(testData)
  }
  testData <- set.chr.to.numeric(testData)
  #return(testData)
  cyto <- get.cyto(ucsc=ucsc,bioC=TRUE,dir=dir)
  cyto <- set.chr.to.numeric(cyto,keep=F)
  #if(any(colnames(cyto) %in% "negpos")) { cyto <- cyto[,-which(colnames(cyto) %in% "negpos")] }
  newDataList <- vector("list",nrow(testData))
  overlaps <- findOverlaps(testData,cyto)
  bandz <- rownames(cyto)[subjectHits(overlaps)]
  indexz <- queryHits(overlaps)
  if(length(indexz)<1) { return(NA) }
  if(!concat) {
    newData <- testData[queryHits(overlaps),]
    newData[["gnm.index"]] <- indexz
    newData[["band"]] <- bandz
  } else {
    out <- tapply(bandz,factor(indexz),c,simplify=FALSE)
    if(concat) { 
      out <- sapply(out,function(X) { paste(X,collapse=";") }) 
      newData <- testData
      #newData[["band"]] <- out 
      newData[["band"]] <- rep("",nrow(newData))
      if(!is.null(names(out))) {
        newData[["band"]][as.numeric(names(out))] <- out
      } else {
        newData[["band"]] <- out
      }
    }
  }
  if(bioC) {
    return(newData)
  } else {
    return(newData[["band"]][order(newData[["index"]])])
  }
  #if(!all(chr %in% chr2(ga))) { stop("invalid chromosome(s) entered") } # redundant i think
}


#' Retrieve SNP ids or positions in specified range
#' 
#' @param chr integer, chromosome, should be a number from 1 to 25, where 23,24,25 are X,Y,MT
#' Alternatively chr can be a RangedData or GRanges object in which case SNP lists will be returned
#' in a list for each row of the ranges object.
#' @param start integer, genomic position to define the start of the range to look for SNPs,
#'  should be between 1 and the length of the chromosome 'chr'
#' @param end integer, genomic position to define the end of the range to look for SNPs,
#'  should be between 1 and the length of the chromosome 'chr', and >= start
#' @param ids logical, if TRUE will return snp ids (chip ids, for rs-ids, use ic.to.rs on the output), 
#' or if FALSE will return the chromosome positions of the SNPs.
#' @param build integer whether to use build 36/37 parameters, 36/37 is preferred, but can enter
#' using any form recognised by ucsc.sanitizer()
#' @param dir string, specify the directory to store or look for annotation (defaults to current)
#' @return Set of SNP ids (when ids=TRUE), or otherwise genomic positions within chromosome 'chr', that
#' fall within the genomic range described by the chr, start, and end parameters. Alternatively, if
#' chr is a RangedData or GRanges object then multiple SNP lists will be returned
#' in a list for each row of the ranges object.
#' @examples
#' snps.in.range(1,9000000,10000000)
#' snps.in.range(10,19000000,20000000,ids=T,36)
#' snps.in.range(10,19000000,20000000,ids=F,36)
snps.in.range <- function(chr, start=NA, end=start, ids=TRUE,build=37,dir=NULL) { 
  # ids - whether to return ichip SNP ids or positions
  if(is(chr)[1]=="RangedData" | is(chr)[1]=="GRanges") {
    chrz <- chr2(chr); stz <- start(chr); enz <- end(chr)
    output <- vector("list",nrow(chr))
    for(cc in 1:nrow(chr)) {
      output[[cc]] <- snps.in.range(chrz[cc],stz[cc],enz[cc],ids=ids,build=build,dir=dir)
    }
    names(output) <- rownames(chr)
    return(output)
  }
  ucsc <- ucsc.sanitizer(build)
  if(is.null(dir)) { dir <- getwd() }
  if(length(chr)>1) { warning("chr should be length 1, using only first entry"); chr <- chr[1] }
  if(length(start)>1) { warning("start should be length 1, using only first entry"); start <- start[1] }
  if(length(end)>1) { warning("end should be length 1, using only first entry"); end <- end[1] }
  if(start>end) { warning("start was higher than end, so switching") }
  the.range <- sort(c(start,end))
  if(!exists("work.dir")) { if(is.null(dir)) { work.dir <- getwd() } else { work.dir <- dir } }
  if(!exists("all.support")) { (load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
  all.chr <- all.support$Chr
  if(ucsc=="hg19" & ("Pos37" %in% colnames(all.support))) {
    all.pos <- all.support$Pos37[all.chr %in% chr]
  } else {
    all.pos <- all.support$Pos[all.chr %in% chr]
  }
  if(length(all.pos)<1) { warning("no positions found for 'chr' specified"); return(NULL) }
  validz <- which(all.pos>=the.range[1] & all.pos<=the.range[2])
  if(ids) {
    out <- rownames(all.support)[all.chr %in% chr][validz]
  } else {
    if(ucsc=="hg19" & ("Pos37" %in% colnames(all.support))) {
      out <- all.support$Pos37[(all.chr %in% chr)][validz]
    } else {
      out <- all.support$Pos[(all.chr %in% chr)][validz]
    }
  }
  return(out)
}


#' Retrieve 'n' closest SNP ids or positions near specified locus
#' 
#' @param chr integer, chromosome, should be a number from 1 to 25, where 23,24,25 are X,Y,MT
#' @param pos integer, genomic position, should be between 1 and the length of the chromosome 'chr'
#' @param n integer, the number of nearest SNPs to seek, if there aren't enough in the annotation
#' then NAs will fill the gaps to force the return value length to equal 'n'
#' @param side character, can be 'either', 'left' or 'right' and specifies which side of the 'pos'
#' to look for nearest snps (where left is decreasing genomic position and right is increasing)
#' @param ids logical, if TRUE will return snp ids (chip ids, for rs-ids, use ic.to.rs on the output), 
#' or if FALSE will return the chromosome positions of the SNPs.
#' @param limit integer, a limit on the maximum distance from the position 'pos' can be specified
#' @param build integer whether to use build 36/37 parameters, 36/37 is preferred, but can enter
#' using any form recognised by ucsc.sanitizer()
#' @param dir string, specify the directory to store or look for annotation (defaults to current)
#' @return Set of SNP ids (when ids=TRUE), or otherwise genomic positions within chromosome 'chr'.
#' If the number of SNPs on the chromosome or the bounds of the 'side' and 'limit' parameters
#' restrict the number returned to less than 'n' then the return value will be padded with NAs.
#' @examples
#' nearest.snp(1,159000000,n=10)
#' nearest.snp(1,159000000,n=10,build=37)
#' nearest.snp(1,159000000,n=10,build=36,ids=F)
#' nearest.snp(1,159000000,n=10,build=37,ids=F)
#' nearest.snp(6,25000000,n=10,build=37,ids=F,side="left")
#' nearest.snp(6,25000000,n=10,build=37,ids=F,side="right")
nearest.snp <- function(chr, pos, n=1, side=c("either","left","right"),ids=TRUE,
                        limit=NULL,build=37,dir=NULL) { 
  # ids - whether to return ichip SNP ids or positions
  ucsc <- ucsc.sanitizer(build)
  if(is.null(dir)) { dir <- getwd() }
  if(length(chr)>1) { warning("chr should be length 1, using only first entry"); chr <- chr[1] }
  if(length(pos)>1) { warning("pos should be length 1, using only first entry"); pos <- pos[1] }
  if(!exists("work.dir")) { if(is.null(dir)) { work.dir <- getwd() } else { work.dir <- dir } }
  if(!exists("all.support")) { (load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
  side <- tolower(side[1]); 
  if(!side %in% c("either","left","right")) {
    side <- "either"; warning("invalid side argument, defaulting to 'either'") }
  if(!is.null(limit)) { if(!is.numeric(limit)) { limit <- NULL; warning("invalid limit argument, defaulting to NULL") } }
  all.chr <- all.support$Chr
  if(ucsc=="hg19" & ("Pos37" %in% colnames(all.support))) {
    all.pos <- all.support$Pos37[all.chr %in% chr]
  } else {
    all.pos <- all.support$Pos[all.chr %in% chr]
  }
  if(length(all.pos)<1) { warning("no positions found for 'chr' specified"); return(NULL) }
  difz <- pos-all.pos
  all.true <- difz==difz
  if(is.null(limit)) { lfilt <- all.true } else { lfilt <- abs(difz)<=limit }
  if(side=="left") { filt <- difz>0 & lfilt }
  if(side=="right") { filt <- difz<0 & lfilt }
  if(side=="either") { filt <- all.true & lfilt }
  Difz <- abs(difz[filt])
  if(length(Difz)<n)  { warning("fewer than ",n," positions found for 'chr' specified (within 'limit'), NAs returned") }
  indx <- order(Difz)[1:n]
  if(ids) {
    out <- rownames(all.support)[all.chr %in% chr][filt][indx]
  } else {
    if(ucsc=="hg19" & ("Pos37" %in% colnames(all.support))) {
      out <- all.support$Pos37[(all.chr %in% chr)][filt][indx]
    } else {
      out <- all.support$Pos[(all.chr %in% chr)][filt][indx]
    }
  }
  return(out)
}

# returns a table of the annotated allele 1 and allele 2 for a set of rs or ichip ids
AB <- function(id) {
  if(!exists("work.dir")) { work.dir <- getwd() }
  ic.ab <- function(ic.ids) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { (load(cat.path(work.dir,"all.support.RData"))) }  ## load object: all.support [snp support for whole chip]
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


## FIX UVA XCHR name screw up ##
rmv.uva.dup.rows <- function(table1) {
  print(Dim(table1))
  bad.tab <- table1[which(is.na(ic.to.rs(rownames(table1)))),]
  nr <- substr(rownames(bad.tab),1,nchar(rownames(bad.tab))-1)
  table1 <- table1[which(!rownames(table1) %in% paste(nr,"1",sep="")),] # remove duplicate rows
  print(Dim(table1))
  return(table1)
}


# standardize snp ids so they would always appear the same, all _,.;, etc replaced with _
# all names leading with a number preceeded with X. mainly 'make.names' standard R-conventions
clean.snp.ids <- function(snpid.list) {
  snpid.list <- make.names(snpid.list)
  snpid.list <- gsub(".","_",snpid.list,fixed=T)
  return(snpid.list)
}


ensemblify <- function(X) {
  X <- paste(X)
  if(!is.character(X)) { stop("X (ID list) must be character") }
  if(!any(nchar(X)==15)) { 
    X <- pad.left(c("00000000000",X),"0")
    X <- paste("ENSG",X,sep="")
    X <- X[-1]
  }
  if(!all(nchar(X)==15)) {
    warning("not all ENSEMBL ids had 15 characters, input may be incorrect")
  }
  if(any(substr(X,1,4)!="ENSG")) {
    warning("it looks like at least 1 X element had an invalid prefix, should be: ENSG00xxxxxx") 
  }
  prv(X)
  return(X)
}

## convert ensembl ids to gene ids - ... are args passed to get.gene.annot()
ENS.to.GENE <- function(id.list,ucsc="hg18",name.dups=TRUE,name.missing=TRUE,to.gene=TRUE) {
  must.use.package(c("biomaRt","genoset","gage"),T)
  ucsc <- ucsc.sanitizer(ucsc)
  if(ucsc=="hg18") {
    ens <- useMart("ENSEMBL_MART_ENSEMBL",
                   dataset="hsapiens_gene_ensembl",
                   host="may2009.archive.ensembl.org",
                   path="/biomart/martservice",
                   archive=FALSE)
  } else {
    ens <- useMart("ensembl")
  }
  ens <- useDataset("hsapiens_gene_ensembl",mart=ens)
  data(egSymb)
  egSymb[,"eg"] <- ensemblify(egSymb[,"eg"])
  ### now have the gene data with the ensembl ids ##
  if(to.gene) {
    id.list <- ensemblify(id.list)
    indx <- match(id.list,egSymb[,"eg"])
    txt <- "ENSEMBL"; col <- "sym"
  } else {
    indx <- match(id.list,egSymb[,"sym"])
    txt <- "Gene"; col <- "eg"
  }
  missin <- length(which(is.na(indx))); valid <- length(indx)-missin
  if(valid<1) {  warning("did not find any ",txt," ids from id.list in the bioMart human gene reference"); return(NULL) }
  if(missin>0) {  warning(out.of(missin,(valid+missin))," of id.list did not match any ",txt," ids in the bioMart human gene reference") }
  #prv(id.list,egSymb);
  outData <- egSymb[indx,col]
  if(name.missing & any(is.na(outData))) {
    outData[is.na(outData)] <- paste("MISSING",pad.left(1:length(which(is.na(outData))),"0"),sep="_")
  }
  if(any(duplicated(outData))) { 
    if(name.dups) { 
      cnt <- 2
      while(any(duplicated(outData))) { 
        if(cnt==2) {
          outData[duplicated(outData)] <- paste(outData[duplicated(outData)],cnt,sep="_")
        } else {
          outData[duplicated(outData)] <- gsub(paste("_",cnt-1,sep=""),paste("_",cnt,sep=""),outData[duplicated(outData)])
        }
        cnt <- cnt + 1
      }
    } else { 
      warning("duplicated gene names produced, select 'name.dups=TRUE' to append numbers to make these unique")
    }
  }
  return(outData)
}

## convert gene ids to ensembl ids - ... are args passed to get.gene.annot()
GENE.to.ENS <- function(id.list,...,to.gene=FALSE) {
  if(to.gene) { to.gene <- FALSE ; warning("to.gene is always false for GENE.to.ENS") }
  return(ENS.to.GENE(id.list,...,to.gene=FALSE))
}


# convenience function to reverse the conversion
conv.37.36 <- function(ranged=NULL,chr=NULL,pos=NULL,chain.file="/home/oliver/R/stuff/hg19ToHg18.over.chain") {
  return(conv.36.37(ranged=ranged,chr=chr,pos=pos,chain.file=chain.file))
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
  if(!is.null(chr) & !is.null(pos)) { ranged <- data.frame.to.ranges(cbind(chr,pos),start="pos",end="pos") }
  #return(ranged)
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
    ranged.gr <- ranged
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
  if(!toranged | T) { return(ranged.gr.37) }
  ranged.gr.37 <- toGenomeOrder2(as(ranged.gr.37,"RangedData"))
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
  if(is.data.frame(X)) {
    if(any(sapply(lapply(X,is),"[",1) %in% c("character","factor"))) {
      for(cc in 1:ncol(X)) {
        X[[cc]] <- as.numeric(X[[cc]])
      }
    }
  }
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
  if(snpmat) {  
    ii <-  col.summary(X)
    all.typz <- c("minor","major")[as.numeric(round(ii$RAF,3)!=round(ii$MAF,3))+1]
  } else {
    all.typz <- apply(X,2,do.mm,snpmat=snpmat)
  }
  return(factor(all.typz))
}


## function that calculates which way around the data is for case vs control (pheno) for a GWAS test
# will indicate with respect to cases whether they have more reference allele, or less, 
# or also if the heterozygous is the affected genotype
# works on a SnpMatrix or dataframe coded 0,1,2,NA (autodetects which)
# using SnpMatrix with het.effects=FALSE can be much faster (5-100x) than other options
caseway <- function(X, pheno, checks=TRUE, long=FALSE, het.effects=FALSE) {
  # coding of output based on long=T/F
  if(long) { r1 <- "cases have more 1, less 0,2" } else { r1 <- "CasesHet+" }
  if(long) { r2 <- "cases have less 1, more 0,2" } else { r2 <- "CasesHet-" }
  if(long) { r3 <- "cases have more 0, less 2" } else { r3 <- "CasesRef-" }
  if(long) { r4 <- "cases have more 2, less 0" } else { r4 <- "CasesRef+" }
  ## workhorse internal function ##
  do.cw <- function(x,ph,snpmat=NULL,r1,r2,r3,r4) { 
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
    if((case.pc0 < ctrl.pc0) & (case.pc2 < ctrl.pc2)) { res <- r1  }
    if((case.pc0 > ctrl.pc0) & (case.pc2 > ctrl.pc2)) { res <- r2  }
    if((case.pc0 > ctrl.pc0) & (case.pc2 < ctrl.pc2)) { res <- r3  }
    if((case.pc0 < ctrl.pc0) & (case.pc2 > ctrl.pc2)) { res <- r4  }
    if((case.pc0 == ctrl.pc0) & (case.pc2 < ctrl.pc2)) { res <- r3  }
    if((case.pc0 == ctrl.pc0) & (case.pc2 > ctrl.pc2)) { res <- r4  }
    if((case.pc0 > ctrl.pc0) & (case.pc2 == ctrl.pc2)) { res <- r3  }
    if((case.pc0 < ctrl.pc0) & (case.pc2 == ctrl.pc2)) { res <- r4  }
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
  if(is(X)[1] %in% "SnpMatrix") { 
    snpmat <- T
    if(!het.effects) {
      SSTS <- single.snp.tests(pheno, snp.data=X, score=T)
      direc <- effect.sign(SSTS)
      if(long) { r3 <- "cases have more of the allele coded '0'" } else { r3 <- "CasesRef-" }
      if(long) { r4 <- "cases have more of the allele coded '2'" } else { r4 <- "CasesRef+" }
      all.res <- rep("???",length(direc))
      all.res[direc==1] <- r4
      all.res[direc==-1] <- r3
      return(factor(all.res))
    }
  } else {
    tt.temp <- table(round(as.numeric(X[,1])));  if("3" %in% names(tt.temp)) { snpmat <- T }
  }
  all.res <- apply(X,2,do.cw,ph=pheno,snpmat=snpmat,r1=r1,r2=r2,r3=r3,r4=r4)
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
  sample.names <- rownames(myData)[smp.filt]
  #print("1");print(length(smp.filt)); print(Dim(myData))
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
    #print("2");print(length(smp.filt)); print(Dim(myDat))
    if(nrow(myDat)!=length(smp.filt)) { 
      if(!all(sample.names %in% rownames(myDat))) {
        use.exist <- FALSE
        stop("missing samples in loaded file, suggest deleting existing imputation and re-run") 
      } else {
        cat("re-arranging samples in file to match smp.filt\n")
        indxz <- match(sample.names,rownames(myDat))
        if(length(which(is.na(indxz)))<1) {  myDat <- myDat[indxz,] } else { stop("rearranging failed") }
        if(nrow(myDat)!=length(smp.filt)) { stop("something went wrong with sample matching") }
        if(any(rownames(myDat)!=sample.names)) { stop("selection went wrong with sample matching") }
      }
    }
    if(ncol(myDat)!=length(snp.filt)) { 
      ## if this has changed a bit, try to ressurrect without recalculating the whole thing
      cat("mismatching number of snps in loaded file\n") 
      if(!exists("bigDat")) { print(load("allImputed.RData")) }
      targs <- colnames(myData)[snp.filt]
      gotem <- narm(match(targs,colnames(bigDat)))
      aintgotem <- targs[!targs %in% colnames(bigDat)]
      if(length(aintgotem)>0) {
        use.exist <- length(gotem)>0
        if(use.exist) {
          cat("combining",length(gotem),"previously imputed with",length(aintgotem),"from scratch\n")
          myDat.part1 <- impute.missing(myData[smp.filt,match(aintgotem,colnames(myData))],numeric=T)
          if(!exists("bigDat")) { print(load("allImputed.RData")) }
          myDat.part2 <- bigDat[,gotem]
          if(any(rownames(myDat.part1)!=rownames(myDat.part2))) { 
            if(!all(rownames(myDat.part1) %in% rownames(myDat.part2))) {
              warning("samples missing from existing dataset") 
              cat("samples were missing from existing dataset, ")
              use.exist <- FALSE
            } else {
              indzx <- match(rownames(myDat.part1),rownames(myDat.part2))
              if(length(which(is.na(indzx)))<1) {  myDat.part2 <- myDat.part2[indzx,] } else { stop("rearranging samps failed") }
            }
          }
        }
        if(use.exist) {
          # ie, if still true [becomes false if samples were missing]
          myDat.cbind <- cbind(myDat.part1,myDat.part2)
          indz <- match(targs,colnames(myDat.cbind))
          if(any(is.na(indz))) { stop("combined data still missing target SNPs") }
          myDat <- myDat.cbind[,indz]
        } else { cat("no existing valid imputation, ") }
        if(!use.exist) {
          ## do all from scratch
          cat("imputing from scratch\n")
          myDataFilt <- myData[smp.filt,snp.filt]
          myDat <- impute.missing(myDataFilt,numeric=T)
        }
      } else { 
        ## have all the snps, but need to prune some
        cat("trimming loaded data to subset needed\n")
        indz <- match(targs,colnames(bigDat))
        if(any(is.na(indz))) { stop("not sure why loaded data is missing target SNPs") }
        if(length(indz)==0) { stop("no target snps in list") }
        myDat <- bigDat[,indz]
      }
      #prv(myDat)
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



find.overlapping.regions  <- function(ranged) {
  if(is(ranged)[1]==("GRanges")) { ranged <- as(ranged,"RangedData") }
  if(!is(ranged)[1]==("RangedData")) { stop("Need RangedData or GRanges object to proceed")} else {
    ranged <- toGenomeOrder(ranged,strict=T) }  
  rr <- NULL
  chrz <- narm(rownames(chrInfo(ranged)))
  LenC <- length(chrz)
  if(LenC>1) { 
    ret <- vector("list",LenC)
    for (ccc in 1:LenC) { 
      cat("Chr",chr(ranged[ccc])[1],":\n")
      ret[[ccc]] <- find.overlapping.regions(ranged[ccc]) 
      if(!is.null(ret[[ccc]])) { 
        print(paste("chr",ccc))
        uu <- ranged[ccc][sort(unique(as.numeric(ret[[ccc]]))),]
        print(uu); if(ccc==1) { rr <- uu } else { rr <- rbind(rr,uu) }
      }
    }
    names(ret) <- paste(chrz)[1:LenC]
    isnull <- sapply(ret,function(X) { all(is.na(X)) })
    return(rr) #ret[!isnull])
  }
  ov <- findOverlaps(ranged)
  qh <- queryHits((ov))
  sh <- subjectHits((ov))
  ovs <- which(qh!=sh)
  if(length(ovs)>0) {
    overlaps <- cbind(qh,sh)[ovs,]
    overlaps <- t(apply(overlaps,1,sort))
   #s return(overlaps)
    otxt <- apply(overlaps,1,function(x) { paste(x,collapse=",") })
    overlaps <- overlaps[!duplicated(otxt),]
    otxt <- otxt[!duplicated(otxt)]
    otxt <- paste("overlap found for regions: ",otxt,"\n")
    cat(otxt,"\n",sep="")
    if(length(Dim(overlaps))==1) { dim(overlaps) <- c(1,2) }
    if(!is.null(dim(overlaps))) {
      colnames(overlaps) <- c("region.a","region.b")
    }
    return(overlaps)
  } else {
    cat("no overlaps\n")
    return(NA)
  }
}

# chris' function to get a centimorgan window from intervals
# vector input
recomWindow <- function(ranged=NULL,chr=NA,st=NA,en=st,window=0.1, # cM either side
                      do.plot=FALSE, # if wanted to plot
                      add.plot=FALSE,do.lines=TRUE,bp.ext=0,...) {
  if(is(ranged)[1] %in% c("RangedData","GRanges")) { 
    if(is(ranged)[1]=="GRanges") { ranged <- as(ranged,"RangedData") }
    ranged <- toGenomeOrder(ranged,strict=T)
    ss <- start(ranged); ee <- end(ranged); cc <- chr(ranged)
    out <- recomWindow(chr=cc,st=ss,en=ee,window=window,bp.ext=bp.ext,...)
    outData <- RangedData(ranges=IRanges(start=out[,1],end=out[,2],names=rownames(ranged)),space=cc)
    outData <- toGenomeOrder(outData,strict=TRUE)
    for (zz in 1:ncol(ranged)) { outData[[colnames(ranged)[zz]]] <- ranged[[colnames(ranged)[zz]]]  }
    if(is(ranged)[1]=="GRanges") { outData <- as(toGenomeOrder(outData,strict=T),"GRanges") }
    return(outData)
  } else {
    if(all(!is.na(chr)) & all(!is.na(st)) & all(!is.na(en))) {
      if(length(chr)==length(st) & length(st)==length(en)) {
        if(length(chr)>1) {
          # run for a vector
          out <- matrix(ncol=2,nrow=length(chr)); colnames(out) <- c("start","end")
          for (dd in 1:length(chr)) {
            out[dd,] <- recomWindow(chr=chr[dd],st=st[dd],en=en[dd],window=window,bp.ext=bp.ext,...)
          }
          return(out)
        } else {
          ## continue as normal, just a single coordinate/range to process
        }
      } else {
        stop("invalid input, st, en and chr need to be the same length")
      }
    } else {
      stop("invalid input, either use a RangedData object, or else chr, st and en")
    }
  }
  rate.fn <- sprintf("/dunwich/scratch/chrisw/HapMap/rates_rel22/genetic_map_chr%s_b36.txt.gz",chr)
  #print(rate.fn)
  rates <- read.table(gzfile(rate.fn),header=TRUE)
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
  cat("new window size is\nleft: ",(st-from+bp.ext)/1000,"kb\tright: ",(to-en+bp.ext)/1000,"kb\ttotal: ",(to-from+(2*bp.ext))/1000,"kb\n",sep="")
  if(bp.ext>0) { cat("in addition to cM distance, window was extended by",bp.ext,"base pairs on either side\n")} 
  from <- max(c(0,(from-bp.ext)))
  to <- min(c((to+bp.ext),get.chr.lens()[chr][1]),na.rm=T)
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




data.frame.to.ranges <- function(dat,ids=NULL,start="start",end="end",width=NULL,
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
  if(length(ch)>0) { ch1 <- gsub("X","chrX",gsub("chr","",dat[[ch]],ignore.case=T)) } else { ch1 <- NULL }
  if(length(st)>0) { st1 <- as.numeric(dat[[st]]) } else { st1 <- NULL }
  if(length(en)>0) { en1 <- as.numeric(dat[[en]]) } else { en1 <- NULL }
  if(length(wd)>0) { en1 <- st1+as.numeric(dat[[wd]]) } # { en1 <- st1+dat[[wd]] }
  #print(length(st1)); print(head(st1))
  outData <- RangedData(ranges=IRanges(start=st1,end=en1,names=id),space=ch1,universe=ucsc[1])
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


## usage 
# export.all.or.pv(TR,qc.excluded.snps,qc.cloud.fail,fn="forolly.RData")
# function takes the object 'TR' tidies and makes a clean dataframe of results ready for export
export.all.or.pv <- function(TR,qc.excluded.snps,qc.cloud.fail,fn="forolly.RData") {
  jj <- order(TR[["meta.p.value"]])
  tt <- as.data.frame(TR)[jj,]
  kk2 <- which(!(rs.to.ic(tt[,5]) %in% unique(rs.to.ic(c(qc.excluded.snps,qc.cloud.fail)))))
  tt <- tt[kk2,]
  tt <- tt[,-3:-4]
  colnames(tt)[1:2] <- c("Chr","Pos")
  aabb <- AB(tt$names)
  forolly <- cbind(tt,aabb)
  colnames(forolly)[11:12] <- c("allele.A","allele.B")
  colnames(forolly)[10] <- "band"
  colnames(forolly)[3] <- "rsid"
  rownames(forolly) <- rs.to.ic(forolly[,3])
  
  save(forolly,file=fn)
  return(forolly)
}


# Create a table comparing the meta and case-control p-values for the old version of the table (uva)
# versus the new version (DIL, includes extra ~2,500 CBR samples) 
# e.g, 
#  > compare.to.previous.results(bonfs.filt,first=T)
#  > compare.to.previous.results(bonfs.filt,first=F,true.reps=c(6,20,23,39))
compare.to.previous.results <- function(bonfs.filt,first=T,true.reps=NULL,
                                        prv.fn="finalMetaTopHits7FEB.RData") {
  new.tab <- bonfs.filt
  print(load(prv.fn))
  old.tab <- bonfs.filt
  all.reg <- unique(c(rownames(new.tab),rownames(old.tab)))
  nr <- length(all.reg); nc <- 5; com.tab <- as.data.frame(matrix(nrow=nr,ncol=nc))
  rownames(com.tab) <- all.reg
  colnames(com.tab) <- c("oldSnp","old_p","newSnp","new_p","equivalent.to.NEW")
  sel1 <- which(rownames(com.tab) %in% rownames(old.tab))
  sel2 <- which(rownames(com.tab) %in% rownames(new.tab))
  com.tab[sel1,1:2] <- old.tab[match(rownames(com.tab)[sel1],rownames(old.tab)),1:2]
  com.tab[sel2,3:4] <- new.tab[match(rownames(com.tab)[sel2],rownames(new.tab)),1:2]
  com.tab[sel2,5] <- sapply(get.equivs(com.tab[sel2,3],iden.list),paste,collapse=",")
  ## by inspection first time, checked whether any old snps in the equivalent list (shouldn't be):
  # if they are, put the row numbers into the vector 'true.reps'
  if(first) { 
    iii <- (which("NA"!=(com.tab[,5]) & (com.tab[,1]!=com.tab[,3])))
    print(com.tab[iii,]) ; stop() 
  }
  # then once vector is checked, first ==FALSE
  if(length(true.reps)>0) {  com.tab[true.reps,3] <- com.tab[true.reps,1] }
  com.tab[sel2,5] <- sapply(get.equivs(com.tab[sel2,3],iden.list),paste,collapse=",")
  
  com.tab[,4] <- paste(substr(com.tab[,4],1,6),substr(com.tab[,4],nchar(com.tab[,4])-3,nchar(com.tab[,4])),sep="")
  chrzz <- apply(cbind(Chr(com.tab[,1]),Chr(com.tab[,3])),1,mean,na.rm=T)
  com.tab <- com.tab[order(chrzz),]
  
  com.tab[["equiv.SNPs"]] <- com.tab[,5]
  com.tab[,5] <- Pos(com.tab[,3],build=37)

  colnames(com.tab)[5] <- "New.snp.pos"

  return(com.tab)
}



## my GWAS with SNPstats will not have ORs/betas with the right directions to be consistent 
# with UVA's TDT family analysis. In order to ensure the OR/betas are the right way around,
# flip each if necessary to have log(sign) consistent with the UVA case-control analysis
# e.g. 
#    new.table.fn <- "nick.meta.table2.RData"
#    table.nick <- convert.OR.directions(new.table.fn)
#    save(table.nick,file="nick.meta.table3.RData")
convert.OR.directions <- function(new.table.fn) {
  nick.table <- reader(new.table.fn)
  print(Dim(nick.table))
  rownames(nick.table) <- clean.snp.ids(rownames(nick.table))
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  excl <- reader(docs[2])
  uva.table <- reader(docs[1])
  rownames(uva.table) <- clean.snp.ids(rownames(uva.table))
  uva.table <- uva.table[rownames(nick.table),]
  wwu <- sign(log(uva.table[,"OR_CC"])) 
  names(wwu) <- rownames(uva.table)
  wwn <- sign(log(nick.table[,"OR_CC"])) 
  names(wwn) <- rownames(nick.table)
  wwun <- wwu[match(names(wwn),names(wwu))]
  length(which(is.na(wwun)))
  wwun[is.na(wwun)] <- wwn[is.na(wwun)] 
  print(head(cbind(wwn,wwun)))
  print(length(which(wwn==wwun)))
  print(length(which(wwn!=wwun)))
  print(Dim(nick.table))
  nick.table[["OR_CC_Raw"]] <- nick.table[,"OR_CC"]
  nick.table[,"OR_CC"][wwn!=wwun] <- 1/(nick.table[,"OR_CC"][wwn!=wwun])
  cat("removing columns",paste(colnames(nick.table)[15:17],collapse=","),"\n")
  cat("replacing with:",paste(colnames(meta.me(nick.table)),collapse=","),"\n")
  table.nick <- cbind(nick.table[,-15:-17],(meta.me(nick.table)))
  colnames(table.nick)[17:ncol(table.nick)] <- c("OR_Meta","b_Meta","SE_Meta","Z_Meta","P_Meta")
  return(table.nick)
}


# returns list of snps from old table, not in new
# prints summary of new vs old pvalues and odds ratios for those snps
# e.g, examine.no.longer.t1.snps(table1snpsfinaljan30,bonf.snps,table1a,identicals,T,T)
examine.no.longer.t1.snps <- function(table1snpsfinaljan30,bonf.snps,table1a,identicals,do.OR=FALSE,do.p=TRUE) {
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  excl <- reader(docs[2])
  uva.table <- reader(docs[1])
  why <- table1snpsfinaljan30[!ic.to.rs(table1snpsfinaljan30) %in% ic.to.rs(bonf.snps)]
  ct.fn <- "conditionalTests.csv"
  ct <- reader(ct.fn,stringsAsFactors=F)
  conditionals <- ct$TABLE1[ct$COND!=0]
  if(length(why)>0 & exists("conditionals")) {  why <- why[!ic.to.rs(why) %in% ic.to.rs(conditionals)] }
  if(length(why)>0) {   why <- why[!ic.to.rs(why) %in% ic.to.rs(identicals)] }
  if(length(why)==0) { return("none found a mystery")}  
  
  if(do.p) {
    # compare pvalues new vs old for SNPs in the old table 1, no longer in the new table 1
    ww <- cbind(table1a[rs.to.ic(why),"P_Meta"],
                uva.table[rs.to.ic(why),"P_Meta"],Chr(why),Pos(why))[order(Chr(why)),]
    xx <- cbind(table1a[rs.to.ic(why),"P_CC"],
                uva.table[rs.to.ic(why),"P_CC"])[order(Chr(why)),]
    zz<-cbind(xx,ww)
    colnames(zz) <- c("New_CC","Old_CC","New_Meta","Old_Meta","Chr","Pos")
    zz <- zz[,c(5,6,1,2,3,4)]
    rownames(zz) <- why[order(Chr(why))]
    print(zz,digits=4)
  }
  if(do.OR) {
    # compare odds ratios new vs old for SNPs in the old table 1, no longer in the new table 1
    ww2 <- cbind(table1a[rs.to.ic(why),"OR_Fam"],
                 uva.table[rs.to.ic(why),"OR_Fam"],Chr(why),Pos(why))[order(Chr(why)),]
    xx2 <- cbind(table1a[rs.to.ic(why),"OR_CC"],
                 uva.table[rs.to.ic(why),"OR_CC"])[order(Chr(why)),]
    zz2 <-cbind(xx2,ww2)
    colnames(zz2) <- c("New_CC","Old_CC","New_Fam","Old_Fam","Chr","Pos")
    zz2 <- zz2[,c(5,6,1,2,3,4)]
    rownames(zz2) <- why[order(Chr(why))]
    print(zz2,digits=4)  
  }
  return(why) # returns list of snps from old table, not in new
}


# T/F see whether snps in a list have been checked for signal clouds previous
unchecked <- function(snps) {
  source("~/github/iChip/hardCodedSnpLists.R",local=TRUE,echo=FALSE)
  list.to <- ((!snps %in% qc.cloud.fail) & (!snps %in% ok) & (!snps %in% good.snps.checked))
  return(list.to)
}

#' Posterior probability of association function
#'
#' @param p p-value you want to test [p<0.367]
#' @param prior prior odds for the hypothesis (Ha) being tested
#' @return prints calculations, then returns the posterior 
#' probability of association given the observed p-value 
#' under the specified prior
#' @references
#' Equations 1, 2 from
#' http://www.readcube.com/articles/10.1038/nrg2615
#' Equations 2, 3 from
#' http://www.tandfonline.com/doi/pdf/10.1198/000313001300339950
#' @examples
#' ps <- rep(c(.05,.01),3)
#' prs <- rep(c(.05,.50,.90),each=2)
#' mapply(ps,prs,FUN=ppa)  # replicate Nuzzo 2014 table
#' # try with bayes factors
#' ppa(BF=3,prior=.9)
#' ppa(BF=10,prior=.5)
ppa <- function(p=.05, prior=.5, BF=NULL, quiet=TRUE) {
  if(any(p<=0 | p>=(1/exp(1)))) { stop("invalid p value") }
  if(any(prior<=0 | prior>=(1))) { stop("invalid prior") }
  if(is.null(BF)) { 
    # calculate bayes factors from p, if BF not entered
    if(!quiet) { cat("\np value:",p,"with prior:",prior,"\n") }
    BF <- (-exp(1)*(p)*log(p) )^(-1)
    # NB: ^invert BF so in terms of % support for Ha 
  } else { 
    if(!quiet) { cat("\nprior:",prior,"with ") }
    if(any(BF<0)) { stop("invalid bayes factor (BF)") }
  }
  if(!quiet) { cat("bayes factor:",BF,"\n") }
  P0 <- (prior/(1-prior)) * (BF) 
  if(!quiet) { cat("posterior odds = bayes factor * H1/H0 prior:",P0,"\n") }
  ppa <- (P0/(1+P0)) 
  if(!quiet) { cat("posterior probability of association:",ppa,"\n") }
  return(ppa)
}

# to get meta analysis parameters from table containing case-control and family data beta, se values
meta.me <- function(X) {
  OR.CC <- X[,"OR_CC"]
  beta.CC  <- log(X[,"OR_CC"])
  se.CC <- X[,"SE_CC"]
  OR.family <- X[,"OR_Fam"]
  beta.family  <- log(X[,"OR_Fam"])
  se.family <- X[,"SE_Fam"]
  z.CC <- beta.CC/se.CC
  z.family <- beta.family/se.family
  
  inv.CC <- 1 / (se.CC^2)
  inv.family <- 1 / (se.family^2)
  var.meta <- 1 / (inv.CC+inv.family)
  weight.CC <- inv.CC * var.meta
  weight.family <- inv.family * var.meta
  
  famN <- 3819*2  #3509*2   #  3819*2   #  10796
  ccN <-  6683+12173 # including CBR, or for UVA analyses use instead: 9416+6670
  WeightFam = sqrt(famN)/(sqrt(famN)+sqrt(ccN))
  #WeightFam = wf
  WeightCC <- 1-WeightFam
  
  beta.meta <- round((weight.CC * beta.CC) + (weight.family * beta.family),digit=3)
  z.metaW1 <- round((weight.CC * z.CC) + (weight.family * z.family),digit=6)
  z.metaW2 <- round((WeightCC * z.CC) + (WeightFam * z.family),digit=6)
  se.meta <- round(sqrt(var.meta), digit=3)
  z.meta <- beta.meta/se.meta
  OR.meta <- exp(beta.meta)
  p.meta <- 2*pnorm(-abs(z.meta))
  p.metaW1 <- 2*pnorm(-abs(z.metaW1))
  p.metaW2 <- 2*pnorm(-abs(z.metaW2))
  out <- (cbind(OR.meta,beta.meta,se.meta,z.meta,p.meta)) #,z.metaW1,p.metaW1,z.metaW2,p.metaW2))
  rownames(out) <- rownames(X)
  return(out)
}



# turn the 'condit.res' object returned by conditionalAnalysis.R into a 
# clean table ready to paste into table 1
condit.to.res <- function(condit.res) {
  newt1rown <- unlist(lapply(condit.res,names))
  newt1p <- unlist(lapply(condit.res,function(X) { lapply(X,"[",1) }))
  newt1simp <- unlist(unlist(lapply(condit.res,function(X) { lapply(X,"[","glm") }),recursive=F),recursive=F)
  OR_SE <- lapply(newt1simp,function(X) { Y <- tail(X,1); return(Y[c(1,5)]) })
  out <- cbind(newt1rown,newt1p,sapply(OR_SE,"[",1),sapply(OR_SE,"[",2))
  colnames(out) <- c("rsid","aov.P","CC_OR","CC_SE.beta")
  out <- cbind(out,Chr(out[,"rsid"]),Pos(out[,"rsid"],37))
  colnames(out)[5:6] <- c("Chr","Pos")
  out <- as.data.frame(shift.rownames(out))
  rn <- rownames(out)
  out <- as.data.frame(lapply(out,as,"numeric"))
  rownames(out) <- rn
  out[["allele.A"]] <- out[["allele.B"]] <- rep("A",nrow(out))
  out[,6:7] <- AB(rownames(out))
  return(out)
}


## for a list of snp-ids from iChip, obtain the nearby SNP-lists within 0.1cm, etc
# do.bands labels each sublist by the band name, but faster not to do this
get.nearby.snp.lists <- function(snpid.list,cM=0.1,bp.ext=0,build=37,excl.snps=NULL,do.bands=TRUE) {
  #if(!exists("all.support")) { print(load("all.support.RData")) }
  snpic.list <- rs.to.ic(snpid.list)
  cyto <- get.cyto(dir=getwd()); cyto[["gene"]] <- rownames(cyto)
  #which.snps <- match(snpid.list,all.support$dbSNP)
  #if(any(is.na(which.snps))) { stop(paste("NAs in dbSNP match:",paste(snpid.list[is.na(which.snps)],collapse=","))) }
  snps.locs36 <- Pos(snpid.list,36) #snp.support$Pos[which.snps]
  snps.locs37 <- Pos(snpid.list,37)
  next.chr <- unique(Chr(snpid.list)); if(length(next.chr)>1) { stop("enter snpids from only 1 chromosome at a time!") }
  if(build==36) { snps.locs <- snps.locs36 } else { snps.locs <- snps.locs37 }
  if(any(snps.locs!=sort(snps.locs))) { 
    warning("snp-ids not in position order, rearrangement is preferred but will attempt to continue")
    sort.back <- match(snps.locs,sort(snps.locs))
  } else { sort.back <- 1:length(snps.locs) }
  ddz <- snpic.list[duplicated(snpic.list)]
  if(length(ddz)>0) { warning("dup SNPs:",ddz,"\n") }
  snp.rd <- RangedData(ranges=IRanges(startSnps.locs,endSnps.locs,names=snpic.list),
                       space=rep(next.chr,length(snps.locs)))
  snp.rd <- toGenomeOrder(snp.rd,strict=T) # think it autosorts anyway, but just in case
  if(do.bands) {
    snp.rd <- annot.cnv(snp.rd,gs=cyto); colnames(snp.rd) <- "band"
    bands <- snp.rd$band
  }
  ## recomWindow uses build36 only, so convert back afterwards
  nxt.window <- lapply(snps.locs36, function(X,...) { recomWindow(st=X,...) },chr=next.chr,window=cM,bp.ext=bp.ext)
  if(build==36) {
    st.window <- sapply(nxt.window, "[",1)
    en.window <- sapply(nxt.window, "[",2)
    pozz <- all.support$Pos
  } else {
    st.window <- conv.36.37(chr=next.chr,pos=sapply(nxt.window, "[",1))$Pos
    en.window <- conv.36.37(chr=next.chr,pos=sapply(nxt.window, "[",2))$Pos
    pozz <- all.support$Pos37
  }
  n.snps <- vector("list",length(st.window))
  for(cc in 1:length(st.window)) {
    n.snps[[cc]] <- which(all.support$Chr==next.chr &
                            pozz>=st.window[cc] & 
                            pozz<=en.window[cc] &
                            (!all.support$SNP %in% excl.snps) &
                            (!all.support$dbSNP %in% excl.snps) 
    )
  }
  grp.labs <- lapply(n.snps,function(X) { all.support$SNP[X] })
  if(do.bands) {
    if(length(unique(bands))!=length(bands)) { warning("these bands are not unique ==> ",paste(bands[duplicated(bands)],collapse=",")) }
    grpz <- 1:length(bands)
    names(grp.labs) <- paste(grpz,bands,sep=":")
  }
  grp.labs <- grp.labs[sort.back]
  return(grp.labs)
}


# order rs-ids or ichip ids by genome chr, position
ids.by.pos <- function(ids) {
  pp <- Pos(ids)
  if(any(is.na(pp))) { stop("invalid id list, could not find position for all") }
  ids <- ids[order(pp)]
  cc <- Chr(ids)
  ids <- ids[order(cc)]
  return(ids)
}


# calibrate the bonferroni threshold for conditional analyses
calibrate.cond.bonf <- function(snplist,cm.window=0.1,bp.ext=0,build=37,qclist="snpsExcluded.txt",ret.snps=FALSE) {
  all.snps.tested <- NULL
  qc.excluded.snps <- reader(qclist)
  qc.excluded.snps <- qc.excluded.snps[!rs.to.ic(qc.excluded.snps) %in% rs.to.ic(snplist)]
  
  for (cc in 1:22) {
    cat("chr",cc,"\n")
    snpid.list <- snplist[Chr(snplist) %in% cc]
    if(length(snpid.list)<1) { next }
    grp.labs <- get.nearby.snp.lists(snpid.list,cM=cm.window,bp.ext=bp.ext,build=37,excl.snps=qc.excluded.snps,do.bands=FALSE)
    all.snps.tested <- c(all.snps.tested,unlist(grp.labs))  
  }
  
  ast <- all.snps.tested[!rs.to.ic(all.snps.tested) %in% rs.to.ic(c(qc.excluded.snps,qc.cloud.fail))]
  bcf <- length(ast)
  bcfu <- length(unique(ast))
  cat("implied bonferroni (count all tests) threshold is:",.05/bcf,"\n")
  cat("implied bonferroni (count unique snps only) threshold is:",.05/bcfu,"\n")
  if(ret.snps) {
    return(unique(ast))
  } else {
    return(list(all.tests=bcf,unique.snps=bcfu))
  }
}


  
