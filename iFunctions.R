###NAMESPACE ADDITIONS###
# Depends: R (>= 2.10), grDevices, graphics, stats, utils, genoset, reader, NCmisc, snpStats
# Imports: 
# Suggests:
# importFrom(proftools, readProfileData, flatProfile)
# importFrom(tools, toHTML)
# importFrom(BiocInstaller, biocVersion)
# import(grDevices, graphics, stats, utils)
###END NAMESPACE###



## options ##
options(chip.info="/chiswick/data/ncooper/iChipData/ImmunoChip_ChipInfo_New.RData") # if you can access the file, you won't need to change this path
options(ucsc="hg19") # depends on which analysis, need to set something though
options(save.annot.in.current=1)  # 1 = TRUE, stores annotation in current folder to speed up subsequent lookups


if(Sys.info()[["user"]]=="ncooper")
{
  #source("~/github/plumbCNV/generalCNVFunctions.R", echo=FALSE)
  if(getwd()!= "/home/ncooper"){
    source("/home/ncooper/github/iChip/ChipInfoClass.R", echo=FALSE)
  } else {
    warning("You need to source the file 'ChipInfoClass.R to get some of these functions to work properly")
  }
  internal.analyses <- FALSE
  if(internal.analyses) {
    source("/home/ncooper/github/iChip/firstscriptFunctions.R", echo=FALSE) # only needed for internal analyses
    source('/home/ncooper/github/iChip/specificIFunctions.R', echo=FALSE) # only needed for internal analyses
  }
}

if(getwd()!= "/home/ncooper"){
  must.use.package("snpStats",T)
  must.use.package("reader")
  must.use.package("genoset",T)
}
#


##file includes the generally useful functions: simple.date, out.of, randomize.missing
 # simple.date, out.of should now be in NCmisc

'internals'

finitize <- function(X) {
  if(is.data.frame(X)) { X <- as.matrix(X) }
  return(X[is.finite(X)])
}

minna <- function(...) {
  if(length(list(...))==1) { 
    min(finitize(...),na.rm=TRUE)
  } else {
    min(...,na.rm=TRUE)
  }
}
maxna <- function(...) {
  if(length(list(...))==1) { 
    max(finitize(...),na.rm=TRUE)
  } else {
    max(...,na.rm=TRUE)
  }
}
meanna <- function(...) {
  if(length(list(...))==1) { 
    mean(finitize(...),na.rm=TRUE)
  } else {
    mean(...,na.rm=TRUE)
  }
}
medianna <- function(...) {
  if(length(list(...))==1) { 
    median(finitize(...),na.rm=TRUE)
  } else {
    median(...,na.rm=TRUE)
  }
}
sdna <- function(...) {
  if(length(list(...))==1) { 
    sd(finitize(...),na.rm=TRUE)
  } else {
    sd(...,na.rm=TRUE)
  }
}
sumna <- function(...) {
  if(length(list(...))==1) { 
    sum(finitize(...),na.rm=TRUE)
  } else {
    sum(...,na.rm=TRUE)
  }
}
sortna <- function(...) {
  sort(..., na.last=TRUE)
}

# stats function convenience wrappers

#' Convert p-values to Z-scores
#' 
#' Simple conversion of two-tailed p-values to Z-scores. Written
#' in a way that allows maximum precision for small p-values.
#' @param p p-values (between 0 and 1), numeric, scalar, vector or matrix, 
#' or other types coercible using as.numeric()
#' @return Z scores with the same dimension as the input
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Z.to.p
#' @examples
#' p.to.Z(0.0001)
#' p.to.Z("5E-8")
#' p.to.Z(c(".05",".01",".005"))
#' p.to.Z(matrix(runif(16),nrow=4))
p.to.Z <- function(p) { 
  if(!is.numeric(p)) { p <- as.numeric(p) }
  if(!is.numeric(p)) { stop("p was not coercible to numeric type") }
  ll <- length(which(p<0 | p>1))
  if(ll>0) { warning(ll, " invalid p-values set to NA"); p[p<0 | p>1] <- NA }
  O <- qnorm((p/2),F)
  O[!is.finite(O)] <- NA
  return(-O) 
}

#' Convert Z-scores to p-values
#' 
#' Simple conversion of Z-scores to two-tailed p-values. Written
#' in a way that allows maximum precision for small p-values.
#' @param Z Z score, numeric, scalar, vector or matrix, or other types coercible
#'  using as.numeric()
#' @return p-valuues with the same dimension as the input
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso p.to.Z
#' @examples
#' Z.to.p("1.96")
#' Z.to.p(p.to.Z(0.0001))
#' Z.to.p(37, T)
#' Z.to.p(39, T) # maximum precision exceeded, warnings on
#' Z.to.p(39) # maximum precision exceeded, warnings off
Z.to.p <- function(Z, warn=FALSE) {
  if(!is.numeric(Z)) { Z <- as.numeric(Z) }
  if(!is.numeric(Z)) { stop("Z was not coercible to numeric type") }
  if(any(abs(Z)>=38) & warn) { warning("maximum precision exceeded, p < 10^-300") }
  O <- 2*pnorm(-abs(Z))
  O[!is.finite(O)] <- NA
  return(O) 
}

# internal
l10 <- function(x) { O <- log10(x); O[!is.finite(O)] <- NA; return(O) }
# internal
Cor <- function(...) { cor(...,use="pairwise.complete") }
# internal
pt2 <- function(q, df, log.p=FALSE) {  2*pt(-abs(q), df, log.p=log.p) }





#' Posterior probability of association function
#'
#' @param p p-value you want to test [p<0.367]
#' @param prior prior odds for the hypothesis (Ha) being tested
#' @return prints calculations, then returns the posterior 
#' probability of association given the observed p-value 
#' under the specified prior
#' @export
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

#' Function to add commas for large numbers
#' 
#' Often for nice presentation of genomic locations it is helpful
#' to insert commas every 3 digits when numbers are large. This function
#' makes it simple and allows specification of digits if a decimal number
#' is in use.
#' @param x a vector of numbers, either as character, integer or numeric form
#' @param digits integer, if decimal numbers are in use, how many digits to display, 
#' same as input to base::round()
#' @return returns a character vector with commas inserted every 3 digits
#' @export
#' @examples
#' comify("23432")
#' comify(x=c(1,25,306,999,1000,43434,732454,65372345326))
#' comify(23432.123456)
#' comify(23432.123456,digits=0)
comify <- function(x,digits=2) {
  if(length(Dim(x))>1) { stop("x must be a vector") }
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           if(length(x)>1) { return(sapply(x, comify, digits=digits)) }
  x <- round(as.numeric(x),digits=digits)
  x <- (paste(x)); dec <- ""
  if(any(grep(".",x))) {
    x.plus.dp <- strsplit(x,".",fixed=TRUE)[[1]]
    if(length(x.plus.dp)>2) { stop("x contained invalid decimal point(s)") }
    xx <- x.plus.dp[1]
    if(length(x.plus.dp)==2) { dec <- paste(".",x.plus.dp[2],sep="") }
  } else { xx <- x }
  splt <- strsplit(xx,"")[[1]]
  nm <- rev(splt)
  cnt <- 0; new <- NULL
  LL <- length(nm)
  for (cc in 1:LL) {
    new <- c(nm[cc],new)
    cnt <- cnt+1
    if(cnt>2 & cc!=LL) { new <- c(",",new); cnt <- 0 }
  }
  return(paste(paste(new,collapse=""),dec,sep=""))
}



# internal function to allow flexible input for the build parameter
ucsc.sanitizer <- function(build,allow.multiple=FALSE,show.valid=FALSE) {
  build.alt <- c("hg15","hg20","hg17","hg18","hg19","hg38",17,18,19,20,35,36,37,38,
                "build35","build36","build37","build38","b35","b36","b37","b38")
  build.new <- c("hg15","hg20",rep(c("hg17","hg18","hg19","hg38"),times=5))
  if(show.valid) { return(cbind(valid=build.alt,mapsTo=build.new)) }
  build <- build.new[match(tolower(build),build.alt)]
  if(any(is.na(build))) { 
    warning("Illegal build parameter '",build[1],"', defaulting to hg18") 
    build[is.na(build)] <- "hg18" 
  }
  if(allow.multiple) {
    return(build)
  } else {
    return(build[1])
  }
}


# ok as long as at least one non-missing snp in the summary
#' See snpStats::col.summary. Same in every way, except for the undesirable
#' behaviour of snpStats when a SNP has 100% missing values it is ignored
#' in the call-rate summary (rather than given a zero). This can unintentionally
#' mean that call-rate filters do not filter SNPs with 100% missing values.
#' This function is simply a wrapper that cleans up this problem.
# col.summary2 <- function(object,...) {
#   if(!is(object)[1]=="SnpMatrix")   { stop("'object' must be a SnpMatrix (snpStats package)") } 
#   if(any(!names(list(...)) %in% c("rules","uncertain"))) { 
#     stop("... contained invalid arguments to snpStats::col.summary") }
#   
# }

#internal
pduplicated <- function(X) {
  if(length(Dim(X))>1) {  stop("can only enter a vector into this function") }
  return((duplicated(X,fromLast=T) | duplicated(X,fromLast=F)))
}


#internal
comma <- function(...) {
  paste(...,collapse=",")
}


## internal function for the 'lambdas' function below
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



#' Normalize Lambda inflation factors to specific case-control count
#' 
#' Lambda inflation statistics are influenced by the size of the generating datasets. To facilitate
#' comparison to other studies, this function converts a given lambda from nr cases and mr controls,
#' to n cases and m controls, where m=n=1000 is the most common normalization. All values other than
#' 'Lnm' are forced within this function to be positive integers.
#' @param Lnm numeric, a raw Lambda inflation statistic, generated from nr cases and mr controls
#' @param n integer, desired number of 'virtual' cases to normalise to
#' @param m integer, desired number of 'virtual' controls to normalise to
#' @param nr integer, original number of cases that Lnm was derived from
#' @param mr integer, original number of controls that Lnm was derived from
#' @return A normalized Lambda coefficient
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @references Freedman M.L., et al. Assessing the impact of population stratification
#'  on genetic association studies. Nat. Genet. 2004;36:388-393.
#' @seealso lambdas
#' @examples
#' require(snpStats) ; data(testdata)
#' pheno <- rep(0,nrow(Autosomes))
#' pheno[subject.data$cc=="case"] <- 1
#' raw.L <- lambdas(Autosomes,pheno,output="lambda")
#' L1000 <- lambda_nm(raw.L,1000,1000,200,200)
#' raw.L; L1000
#' lambda_nm(1.56,1000,1000,6500,9300) # lambda1000<lambda when n's>1000,1000
lambda_nm <- function(Lnm,n=1000,m=1000,nr,mr) { 
  if(!is.numeric(Lnm)) { stop("Lnm must be numeric") }
  if(!is.numeric(n)) { stop("n must be numeric") } else { n <- abs(round(n)) }
  if(!is.numeric(m)) { stop("m must be numeric") } else { m <- abs(round(m)) }
  if(!is.numeric(nr)) { stop("nr must be numeric") } else { nr <- abs(round(nr)) }
  if(!is.numeric(mr)) { stop("mr must be numeric") } else { mr <- abs(round(mr)) }
  return(1 + ((Lnm-1)*(((1/nr)+(1/mr))/((1/n)+(1/m)))) )
}


# internal functions for lambdas #
Y_2 <- function(r1,r2,n1,n2,N,R) { (N*((N*(r1+(2*r2)))-(R*(n1+(2*n2))))^2) / ((R*(N-R))*((N*(n1+(4*n2)))-(n1+(2*n2))^2)) }
X_2 <- function(r1,r2,n1,n2,N,R) { (2*N*((2*N*(r1+(2*r2)))-(R*(n1+(2*n2))))^2) / ((4*R*(N-R))*((2*N*(n1+(2*n2)))-(n1+(2*n2))^2)) }
Likelihood_Lj <- function(c,Lnm) { rchisq(c/Lnm,df=1)/Lnm } # likelihood for one marker
LLikelihood_L <- function(Cj,LNMj) {
  # total likelihood across all K markers  : http://www.nature.com/ng/journal/v36/n4/full/ng1333.html
  tot <- 0
  for (cc in 1:length(Cj)) { 
    tot <- tot + log(Likelihood_Lj(Cj[cc],LNMj[cc])) 
  }
  return(tot) 
}



#' Calculate Lambda inflation factors for SNP dataset
#' 
#' This function calculates SNP-wise or overall Lambda and Lambda_1000 statistics for inflation due
#' to population structure. It works on a SnpMatrix object or dataframe coded 0,1,2,NA (autodetects which).
#' @param x numeric, it's the thing that goes in
#' @return if snp.wise is false, the scalar value(s) specified by 'output', or otherwise a matrix
#' of parameters for each SNP, optionally limited to just the lambda column(s) by the value of 'output'.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @references Freedman M.L., et al. Assessing the impact of population stratification
#'  on genetic association studies. Nat. Genet. 2004;36:388-393.
#' @seealso lambda_nm
#' @examples
#' # http://en.wikipedia.org/wiki/Population_stratification
#' # Note that Y2 ~ L*X2, where allele counts are symbolized:
#' # Case     r0  r1  r2  R
#' # Control  s0  s1  s2  S
#' # Total    n0  n1  n2  N
#' require(snpStats) ; data(testdata)
#' pheno <- rep(0,nrow(Autosomes))
#' pheno[subject.data$cc=="case"] <- 1
#' lambdas(Autosomes,pheno)
#' lambdas(Autosomes[,1:5],pheno,snp.wise=T) # list everything snp-wise for the first 5 SNPs
#' lambdas(Autosomes[,1:5],pheno) # just for the first 5 SNPs
lambdas <- function(X, pheno, checks=TRUE, 
                    output=c("all","lamba","l1000","both"),snp.wise=FALSE) {
  ## i don't think this should ever be used: ?? cc1000=TRUE, ??
  cc1000 <- FALSE
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
  if(!snp.wise) { cc1000 <- FALSE }
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
  output <- tolower(paste(output[1]))
  output <- gsub("lamda","lambda",output); output <- gsub("lamba","lambda",output)
  if(!output %in% c("lambda","l1000","both")) { output <- "all" }
  if(!snp.wise) {
    # return overall scalar result(s) across all SNPs (median based)
    lam <- medianna(all.res[,"Y2"])/.456
    if(output=="all") { output <- "both" }
    if(output %in% c("both","l1000")) { 
      lam1000 <- lambda_nm(lam,1000,1000,medianna(all.res$R),medianna(all.res$S))
    }
    out <- switch(output,lambda=lam,l1000=lam1000, both=c(Lambda=lam,L1000=lam1000))
  } else {
    # return separate result(s) for each SNP
    out <- switch(output,all=all.res,lambda=all.res[["Lambda"]],
                l1000=all.res[["L1000"]], both=all.res[,c("Lambda","L1000")])
  }
  return(out)
}




#' Function to produce a clean table from logistic regression done via GLM
#' 
#' With input as a glm model result object, returns a clean dataframe with 
#' coefficients, p values, confidence intervals and standard errors. Multiple
#' options for which columns to return and digits to display
#'
#' @param glm.result should be type 'glm', the result of a call to glm() where
#' the family parameter specifies logistic regression, i.e, family=binomial("logit")
#' @param intercept logical, whether to display the intercept as the first row of results
#' @param ci logical, whether to include upper and lower confidence limits in the table of results
#' @param se logical, whether to include standard error (SE) in the table of results
#' @param alpha percentage, the type 1 error rate to apply to calculation of confidence
#' limits, e.g, alpha=.01 results in a 99% confidence interval.
#' @param o.digits, integer, number of digits to use displaying odds-ratios and standard errors
#' @param p.digits, integer, number of digits to round p-values to, not that a high number ensures
#' that very small p-value are truncated to be displayed as zero.
#' @param label, logical, whether to label output matrices as a list, with name corresponding to 
#' the formula for the model
#' @return the output returned
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' ph <- c(0,0,0,1,1,1) # phenotype
#' X <- rnorm(6) # independent variable
#' test.glm <- glm(ph ~ X,family=binomial('logit'))
#' logistic.summary(test.glm) # very simple example
#' XX <- sim.cor(10,3)
#' X1 <- XX[,1]; X2 <- XX[,2]; X3 <- XX[,3]
#' ph <-c(rep(0,5),rep(1,5))
#' test.glm <- glm(ph ~ X1 + X2 + X3,family=binomial('logit'))
#' logistic.summary(test.glm,intercept=TRUE) # include intercept
#' logistic.summary(test.glm,ci=TRUE,se=FALSE) # show 95% confidence intervals, hide standard error
#' logistic.summary(test.glm,ci=TRUE,alpha=.01) # get 99% confidence intervals
logistic.summary <- function(glm.result,intercept=FALSE,ci=FALSE,se=TRUE,alpha=.05,o.digits=3,p.digits=299,label=FALSE)
{
  if(!"glm" %in% is(glm.result)) { stop("'glm.result' was not a glm result (of type 'glm'") }
  if(family(glm.result)[[1]]!="binomial" | family(glm.result)[[2]]!="logit") { 
    warning("this function is designed to work for glm logistic regression, i.e; ",
            "glm(family=binomial('logit'))", " invalid output is likely")
  }
  co <- summary(glm.result)$coefficients
  alpha <- force.percentage(alpha,default=.05)
  if(rownames(co)[1]!="(Intercept)") { intercept <- TRUE }
  predz <- rownames(co)
  if(!intercept) { predz <- predz[-1] }
  lab <- paste(summary(glm.result)$call)[2]
  #prv(co)
  if(length(Dim(co))<2) {
    warning("glm result was empty of parameters")
    return(NULL)
  }
  if(!intercept & nrow(co)<2) { 
    warning("intercept was not in use, and there were no other parameters in the models") 
    return(NULL) 
  }
  if(!intercept) { fst <- 2 } else { fst <- 1 }
  p <- co[fst:nrow(co),4]; #print(p)
  o.r <- exp(co[fst:nrow(co),1]); #print(o.r)
  p <- round(p,p.digits)
  
  # outlist <- list(round(o.r,o.digits),p)
  # names(outlist) <- c("OR","p-value")
  SE <- co[fst:nrow(co),2]
  if(ci) {
    #prv(co)
    co1 <- co[fst:nrow(co),1]-(p.to.Z(alpha)*SE)
    co2 <- co[fst:nrow(co),1]+(p.to.Z(alpha)*SE)
    if(sum(co1)<=sum(co2)) { co.l <- co1; co.h <- co2 } else {  co.h <- co1; co.l <- co2 }
    co.l <- exp(co.l); co.h <- exp(co.h)
    out.fr <- cbind(round(o.r,o.digits),round(SE,o.digits),round(co.l,o.digits),round(co.h,o.digits),p)
    colnames(out.fr) <- c("OR","SE","OR-low","OR-hi","p-value")
  } else {
    out.fr <- cbind(round(o.r,o.digits),round(SE,o.digits),p)
    colnames(out.fr) <- c("OR","SE","p-value")
  }
  if(is.null(rownames(out.fr)) & nrow(out.fr)==1) {
    rownames(out.fr) <- predz
  }
  if(!se) {
    if(length(Dim(out.fr))==2) {
      out.fr <- out.fr[,-which(colnames(out.fr) %in% "SE")]
    } else {
      if(length(Dim(out.fr))==1) {
        out.fr <- out.fr[-which(colnames(out.fr) %in% "SE")]
      }
    }
  }
  if(label) { out.fr <- list(out.fr); names(out.fr) <- lab }
  return(out.fr)
}



#' Download GWAS hits from t1dbase.org
#' 
#' Retrieve human disease top GWAS hits from t1dbase in either build hg18 or hg19 coords (b36/37).
#' 28 Diseases currently available
#' @param disease integer (1-28), or character (abbreviation), or full name of one of the listed
#' diseases. A full list of options can be obtained by setting show.codes=TRUE.
#' @param build character, either "hg18" or "hg19". Will also accept build number, 36 or 37.
#' @param show.codes logical, if set to TRUE, instead of looking up t1dbase, will simply return
#' a table of available diseases with their index numbers and abbreviations.
#' @return A character vector of SNP rs-ids
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @references PMID: 20937630
#' @examples
#' get.t1dbase.snps(disease="CEL",build=36) # get SNP ids for celiac disease in build-36/hg18
#' get.t1dbase.snps(disease="AS") # get SNP ids for Ankylosing Spondylitis in build-37/hg19
#' get.t1dbase.snps(show.codes=TRUE) # show codes/diseases available to download
#' get.t1dbase.snps(disease=27) # get SNP ids for Alopecia Areata
#' get.t1dbase.snps("Vitiligo")
get.t1dbase.snps <- function(disease="T1D",build=NULL,show.codes=FALSE) {
  disease.codes <- c("Type 1 Diabetes", "Crohns Disease","Rheumatoid Arthritis",
  "Systemic Scleroderma",  "Ulcerative Colitis","Inflammatory Bowel Disease",  "Multiple Sclerosis",
  "Bipolar Disorder",  "Diabetes Mellitus",  "Coronary Artery Disease",  "Hypertension",  "Celiac Disease",
  "Systemic Lupus Erythematosus",  "Ankylosing Spondylitis",  "Type 2 Diabetes",  "Sjogren Syndrome",
  "Graves' Disease",  "Juvenile Rheumatoid Arthritis",  "Vitiligo",  "Primary Biliary Cirrhosis",
  "Psoriasis",  "Idiopathic Membranous Nephropathy",  "Immunoglobulin A Deficiency",
  "Autoimmune Thyroid Disease",  "Juvenile Idiopathic Arthritis",  "Narcolepsy",  "Alopecia Areata",
  "Alzheimer's Disease")
  abbr <- c("T1D","CD","RA","SCL","UC","IBD","MS","BD","DM","CAD","HYP",
            "CEL","SLE","AS","T2D","SS","GD","JRA","VIT",
            "PBC","PSO","IMN","IGA","ATD","JIA","NAR","AA","AD")
  code.table <- cbind(Abbreviation=abbr,FullNames=disease.codes)
  if(show.codes) { cat("values for the 'disease' parameter can be specified by the following index numbers or abbreviations:\n")
                   print(code.table,quote=F) ; return() }
  disease <- disease[1]
  if(toupper(disease) %in% abbr) { 
    disN <- match(toupper(disease),toupper(abbr)) 
  } else {
    if(toupper(disease) %in% toupper(disease.codes)) {
      disN <- match(toupper(disease),toupper(disease.codes))
    } else {
      if(as.numeric(disease) %in% 1:length(disease.codes)) { 
        disN <- as.numeric(disease)
      } else {
        stop("Invalid input for 'disease', use show.codes=TRUE to see list of codes/abbreviations")
      }
    }
  }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  filenm <- paste(dir=getwd(),pref=tolower(abbr[disN]),"hits",suf=build,ext="tab")
  cat("attempting to download",abbr[disN],"hits from t1dbase\n")
  url36 <- paste("http://www.t1dbase.org/webservice/RegionDownloads/=/model/variantsTAB/species=Human&disease_id=",disN,"&type=assoc&build=GRCh36",sep="")
  url37 <- paste("http://www.t1dbase.org/webservice/RegionDownloads/=/model/variantsTAB/species=Human&disease_id=",disN,"&type=assoc&build=GRCh37",sep="")
  urL <- switch(build, hg18=url36,  hg19=url37)
  success <- T
  success <- tryCatch(download.file(urL ,filenm ,quiet=T),error=function(e) { F } )
  if(!is.logical(success)) { success <- T }
  if(success) {
    t1dh <- readLines(filenm,)
    firsts <- substr(t1dh,1,2)
    t1dh <- t1dh[firsts!="##"]
    len.lst <- strsplit(t1dh,"\t")
    rsids <- sapply(len.lst,"[",3)
    if(substr(rsids[1],1,2)!="rs") { rsids <- rsids[-1] }
    if(length(rsids)<1) { 
      cat("download successful but the list of hits for",disease.codes[disN],"was empty\n") 
    } else {
      cat("download successful for",disease.codes[disN],"\n")
    }
    return(unique(rsids))
  } else {
    stop("couldn't reach t1dbase website at: ",urL)
  }
}



#' Retrieve current ChipInfo annotation object
#' 
#' This function returns the current 'ChipInfo' annotation object, containing chromosome,
#' id, position, strand, 'rs' id, allele 1, allele 2 for each SNP of a microarray chip,
#' in either hg18 or hg19 (build 36/37) coordinates.
#' This package makes extension use of this class of annotation object for the working
#' microarray chip, e.g, default is ImmunoChip, but Metabochip is also built-in,
#' and you can also load your own annotation if using a different chip. The class
#' of the object used is 'ChipInfo' which is a GRanges object, modified to always
#' have columns for A1, A2 (alleles), rs.id, and a quality control flag. The
#' default display is tidier than GRanges, it has nice coersion to and frame data.frame
#' and indexing by chromosome using [[n]] has been added, in addition to normal [i,j]
#' indexing native to GRanges.
#' @param build character, either "hg18" or "hg19". Will also accept build number, 36 or 37.
#' @param refresh logical, TRUE to just load whatever object is already in memory (except
#' when first using a function in this package, there should be a ChipInfo object loaded),
#' or FALSE to reload from the original source. For instance you may wish to do this when 
#' you want to use a different chip, different build, or if the annotation has been modifed
#' via a manual correction).
#' @param alternate.file character, name of an alternative RData file containing a ChipInfo
#' object to use instead of the object found in getOption("chip.info"). Such files may
#' contain two objects, one for build 36 and one for build 37; when doing this, the object
#' names should contain the characters '36' and '37' respectively.
#' @return returns the current ChipInfo object [S4]
#' @seealso ChipInfo, build, rs.id, QCfail, convTo36, convTo37, A1, A2
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' get.support() # shows the current ChipInfo object (default is 'ImmunoChip' build 36)
get.support <- function(build=NULL,refresh=FALSE,alternate.file=NULL) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  ## NEED TO ADD SUPPORT HERE TO USE CORRECT OUT OF 36/37
  if(!exists("all.support",envir=globalenv())) { refresh <- T }
  if(refresh) {
    use.options <- TRUE
    if(is.character(alternate.file)) {
      if(file.exists(alternate.file)) {
        use.options <- FALSE
      }
    }
    if(use.options) {
      file <- getOption("chip.info")
    } else {
      file <- alternate.file
    }
    all.support <- reader(file)
    if(length(all.support)>1) { typz <- sapply(lapply(all.support,is),"[",1) } else { typz <- is(all.support)[1] }
    if(all(typz=="ChipInfo")) { 
      assign("all.support",value=all.support,envir=globalenv())
    } else {
      stop("object (all.support) in the file",file,
           "should have type ChipInfo, or else object all.support in the global environment has been modified")
    }
  }
  if(!exists("all.support",envir=globalenv())) { stop("ChipInfo data object 'all.support' not found") }  
  all.support <- get("all.support",envir=globalenv())
  nnn <- (names(all.support))
  if(length(nnn)>1) { 
    # multiple objects in file, probably different builds
    any37 <- c(grep("37",nnn),grep("hg19",nnn))
    any36 <- c(grep("36",nnn),grep("hg18",nnn))
    if(build=="hg19" & length(any37)>0) { all.support <- all.support[[nnn[any37[1]]]]}
    if(build=="hg18" & length(any36)>0) { all.support <- all.support[[nnn[any36[1]]]]}
  }
  return(all.support)
}

#' Convert from chip ID labels to dbSNP rs-ids
#' 
#' Most SNPs will have an 'rs-id' from dbSNP/HapMap, and these are often the standard for reporting or
#' annotation lookup. These can differ from the IDs used on the chip. This functions looks at the current
#' snp support (ChipInfo object) and returns rs-ids in place of chip IDs.
#' @param ids character, meant to be a list of chip ids, but if rs-ids are present they will not be altered.
#' @return A character vector of SNP rs-ids, where the input was chip ids, rs-ids or a mixture, any text
#' other than this will result in NA values being returned in the character vector output.
#' @export
#' @seealso rs.to.id, GENE.to.ENS, ENS.to.GENE
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' id.to.rs(c("imm_11_2138800","rs9467354","vh_1_1108138")) # middle one is already a rs.id
id.to.rs <- function(ids) {
  ids <- clean.snp.ids(ids)
  all.support <- get.support()
  if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found") }  ## load object: all.support [snp support for whole chip]
  rsvec <- mcols(all.support)$rs.id[match(ids,rownames(all.support))]
  rsvec2 <- mcols(all.support)$rs.id[match(ids,mcols(all.support)$rs.id)]
  rsvec[is.na(rsvec)] <- rsvec2[is.na(rsvec)]
  return(rsvec)
}


#' Convert from dbSNP rs-ids to chip ID labels
#' 
#' Most SNPs will have an 'rs-id' from dbSNP/HapMap, and these are often the standard for reporting or
#' annotation lookup. These can differ from the IDs used on the chip. This functions looks at the current
#' snp support (ChipInfo object) and looks up chip IDs based on rs-ids.
#' @param ids character, meant to be a list of rs-ids, but if chip-ids are present they will not be altered.
#' @param multi.list logical, some rs-ids could map to multiple chip ids. It is recommended that if that is
#' the case then a letter should be appended to duplicate rs-ids to make them unique in the ChipInfo object,
#' e.g, rs1234, rs1234b, rs1234c, etc. If multi.list is TRUE, then the id list will be returned as a list,
#' and any time an rs-id is entered without a letter suffix, all possible corresponding chip ids will be
#' listed.
#' @return A character vector of SNP chip-ids, where the input was rs-ids, chip-ids or a mixture, any text
#' other than this will result in NA values being returned in the character vector output. Or, if multi-list
#' is true, then returns a list instead, which takes more than 1 value where there are multiple chip-ids 
#' with the same rs-id; if there are no such rs-id duplicates the result will still be a list.
#' @export
#' @seealso id.to.rs, GENE.to.ENS, ENS.to.GENE
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' rs.to.id(c("rs689","rs9467354","rs61733845"))  # middle one has no chip id
#' test.ids <- c("rs61733845","rs2227313","rs11577783","rs3748816","rs12131065","rs3790567","rs2270614")
#' rs.to.id(test.ids, multi.list=TRUE) # list with duplicates
rs.to.id <- function(rs.ids,multi.list=FALSE) {
  rs.ids <- clean.snp.ids(rs.ids)
  all.support <- get.support()
  if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found") }  ## load object: all.support [snp support for whole chip]
  if(multi.list) {
    X0 <- rmv.trail(rs.ids)
    X1 <- paste0(X0,"b"); X2 <- paste0(X0,"c"); X3 <- paste0(X0,"d"); X4 <- paste0(X0,"a")
    X <- cbind(rs.to.id(X0),rs.to.id(X1),rs.to.id(X2),rs.to.id(X3),rs.to.id(X4))
    idvec <- apply(X,1,function(x) { unique(narm(x)) })
    if(!is.list(idvec)) { idvec <- as.list(idvec) }
      #warning("'multi.list' option was used, but no duplicate rs-ids found, so returning a vector, not a list") }
  } else {
    idvec <- rownames(all.support)[match(rs.ids,mcols(all.support)$rs.id)]
    idvec2 <- rownames(all.support)[match(rs.ids,rownames(all.support))]
    idvec[is.na(idvec)] <- idvec2[is.na(idvec)]
  }
  return(idvec)
}



#' Find chromosome for SNP ids, gene name or band
#' 
#' Allows retrieval of the chromosome associated with a SNP-id, HGNC gene label, karyotype band,
#' or vector of such ids. For SNPs the ids can be either chip ids, or rs-ids, but must be contained
#'  in the current annotation. Default behaviour is to assume 'id' are SNP ids, but if none are
#'  found in the SNP annotation, the id's will be passed to functions Pos.gene() and Pos.band() to
#'  see whether a result is found. This latter step will only happen if no SNP ids are retreived in
#'  the first instance, and if snps.only=TRUE, then genes and bands will not be searched and NA's 
#'  returned. If you are repeatedly searching for chromosomes for genes/bands, using the dedicated 
#'  Pos.gene and Pos.band functions would be slightly faster than relying on the fallback behaviour
#'  of the Chr() function.  See documentation for these functions for more information. The build
#'  used will be that in the current ChipInfo object.
#' @param ids character, a vector of rs-ids or chip-ids representing SNPs in the current ChipInfo
#'  annotation, or gene ids, or karyotype bands
#' @param dir character, only relevant when gene or band ids are entered, in this case 'dir' is the location
#' to download gene and cytoband information; if left as NULL, depending on the value of 
#' getOption("save.annot.in.current"), the annotation will either be saved in the working directory to 
#' speed-up subsequent lookups, or deleted after use.
#' @param snps.only logical, if TRUE, only search SNP ids, ignore the possibility of genes/cytobands.
#' @return A character vector of Chromosomes for each ids, with NA values where no result was found.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Pos
#' @examples
#' Chr(c("rs689","rs9467354","rs61733845"))
#' Chr("CTLA4")
#' Chr("13q21.31")
#' Chr(c("CTLA4","PTPN22"),snps.only=TRUE) # fails as these are genes
#' Chr(c("rs689","PTPN22","13q21.31")) # mixed input, will default to SNPs, as at least 1 was found
Chr <- function(ids,dir=NULL,snps.only=FALSE) {
  ic.chr <- function(ids) {
    ic.ids <- clean.snp.ids(ids)
    all.support <- get.support()
    if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found")  }  ## load object: all.support [snp support for whole chip]
    outlist <- chr(all.support)[match(ic.ids,rownames(all.support))]
    return(outlist)
  }
  query <- rs.to.id(ids)
  if(!snps.only & all(is.na(query))) { 
    ## unless the 'snps.only' function is set, then if it looks like we have not been handed
    ## snp ids, then check for band ids or gene ids instead
    chr.ify <- function(X) { tolower(names(X)); return(gsub("chr","",X["chr"],ignore.case=TRUE)) }
    numpqs <- (length(grep("q",ids))+length(grep("p",ids)))
    if(numpqs==length(ids)) { try.band <- T } else { try.band <- F }
    if(try.band) {
      suppressWarnings(test <- Pos.band(ids,dir=dir))
      if(!is.null(test)) { return(chr.ify(test)) }
    }
    suppressWarnings(test <- Pos.gene(ids,dir=dir))
    if(!is.null(test)) { return(chr.ify(test)) } 
  } 
  ic <- ic.chr(query)
  return(ic)
}


#' Find the chromosome position for SNP ids, gene name or band
#' 
#' Allows retrieval of the the chromosome position associated with a SNP-id, HGNC gene label, 
#'  karyotype band, or vector of such ids. For SNPs the ids can be either chip ids, or rs-ids,
#'  but must be contained in the current annotation. Default behaviour is to assume 'id' are 
#'  SNP ids, but if none are found in the SNP annotation, the id's will be passed to functions
#'  Pos.gene() and Pos.band() to see whether a result is found. This latter step will only happen
#'  if no SNP ids are retreived in the first instance, and if snps.only=TRUE, then genes and bands
#'  will not be searched and NA's returned. If you are repeatedly searching for positions for 
#'  genes/bands, using the dedicated Pos.gene() and Pos.band() functions would be slightly faster
#'  than relying on the fallback behaviour of the Pos() function. Note that the position for
#'  genes and bands are not a single point, so the result will be a range with start and end, 
#'  see 'values' below. See documentation for these functions for more information.
#' @param ids character, a vector of rs-ids or chip-ids representing SNPs in the current ChipInfo annotation,
#'  or gene ids, or karyotype bands
#' @param dir character, only relevant when gene or band ids are entered, in this case 'dir' is the location
#' to download gene and cytoband information; if left as NULL, depending on the value of 
#' getOption("save.annot.in.current"), the annotation will either be saved in the working directory to 
#' speed-up subsequent lookups, or deleted after use.
#' @param snps.only logical, if TRUE, only search SNP ids, ignore the possibility of genes/cytobands.
#' @return When ids are SNP ids, returns a numeric vector of positions for each id, with NA values
#'  where no result was found. When ids are genes or karyotype bands, will return a data.frame with
#'  columns 'chr' [chromosome], 'start' [starting position of feature], 'end' [end position of feature], 
#'  and the band without the chromosome prefix, if ids are bands. Note that this function cannot
#'  retrieve multiple ranges for a single gene (e.g, OR2A1), which means you'd need to use Pos.gene().
#'  The coordinates used will be of version getOption(ucsc="hg18"), or build(get.support()), which
#'  should be equivalent.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr
#' @examples
#' Pos(c("rs689","rs9467354","rs61733845"))
#' Pos("CTLA4") # returns a range
#' Pos("13q21.31") # returns a range
#' Pos(c("CTLA4","PTPN22"),snps.only=TRUE) # fails as these are genes
#' Pos(c("rs689","PTPN22","13q21.31")) # mixed input, will default to SNPs, as at least 1 was found
Pos <- function(ids,dir=NULL,snps.only=FALSE) {
  all.support <- get.support()
  ic.pos <- function(ic.ids) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found") }  ## load object: all.support [snp support for whole chip]
    outlist <- start(all.support)[match(ic.ids,rownames(all.support))]
    return(outlist)
  }
  query <- rs.to.id(ids)
  if(!snps.only & all(is.na(query))) { 
    ## unless the 'snps.only' function is set, then if it looks like we have not been handed
    ## snp ids, then check for band ids or gene ids instead
    numpqs <- (length(grep("q",ids))+length(grep("p",ids)))
    if(numpqs==length(ids)) { try.band <- T } else { try.band <- F }
    if(try.band) {
      suppressWarnings(test <- Pos.band(ids,dir=dir))
      if(!is.null(test)) { return(test) }
    }
    suppressWarnings(test <- Pos.gene(ids,dir=dir))
    if(!is.null(test)) { return(test) } 
  } 
  ic <- ic.pos(query)
  return(ic)
}



#' Order rs-ids or ichip ids by chrosome and position
#' 
#' Simple function to sort a character list of SNP ids into genome order.
#' @param ids character, vector of SNP rs-ids or chip-ids, see rs.to.id()
#' @return the same vector 'ids', sorted by genome position
#' @seealso rs.to.id, id.to.rs, Chr, Pos
#' @export
#' @examples
#' snp.ids <- c("rs3842724","imm_11_2147527","rs689","rs9467354","rs61733845")
#' Chr(snp.ids) # shows each is on a different chromosome
#' Pos(snp.ids)
#' ids.by.pos(snp.ids)
#' Chr(ids.by.pos(snp.ids))
#' Pos(ids.by.pos(snp.ids))
ids.by.pos <- function(ids) {
  if(!is.character(ids)) { stop("ids must be a character vector") }
  pp <- Pos(ids)
  if(any(is.na(pp))) { stop("invalid id list, 'ids' must all be valid, with position information in the ChipInfo object") }
  ids <- ids[order(pp)]
  cc <- Chr(ids)
  ids <- ids[order.chr(cc)]
  return(ids)
}


#internal function to properly sort chromosome labels as text
order.chr <- function(chrs) {
  # sort chr nms
  if(is.numeric(chrs)) { chrs <- paste(chrs) }
  if(!is.character(chrs)) { stop("chrs should be a character or integer vector") }
  asn <- function(X) { suppressWarnings(as.numeric(X)) }
  textz <- is.na(asn(chrs))
  nums <- chrs[!textz]
  txts <- chrs[textz]
  ns <- which(!textz)[order(asn(nums))]
  #print(max(ns,na.rm=TRUE))
  #print(order(txts)); print(txts)
  #print(which(textz))
  ts <- which(textz)[order(txts)]
  out <- c(ns,ts)
  return(out)
}

#internal
sort.chr <- function(chr) { chr[order.chr(chr)] }



#' Find the chromosome, start and end position for gene names
#' 
#' Allows retrieval of the the chromosome position associated with a HGNC gene label, 
#'  or vector of such labels. Note that the position returned for genes is not a 
#'  single point as for SNPs, so the result will be a chromosome, then a position range with
#'  start and end.
#' @param genes character, a vector of gene ids
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download gene annotation information to; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, else 
#' a data.frame
#' @param band logical, whether to include band/stripe in returned object
#' @param one.to.one logical, some genes have split ranges, TRUE merges these to give only 1 range 
#' per gene, NB: this is the default behaviour when using the more general Pos() function
#' @param remap.extra logical, if TRUE genes with chromosome annotation 'c6_cox' and 'c6_QBL' will
#'  be mapped to chromosome 6, and 'NT_xxxx' chromosome labels will all be mapped to 'Z_NT', etc
#' @param discard.extra logical, if TRUE then any gene hit with chromosome not in 1:22, X, Y, XY, MT, 
#' will be discarded.
#' @param warnings logical, whether to show warnings when some/all ids are not matched to the 
#' reference
#' @return Returns a data.frame with columns 'chr' [chromosome], 'start' [starting position of the
#'  gene],'end' [end position of the gene], or if bioC=TRUE, then returns a GRanges object with
#'  equivalent information, and if band=TRUE, then an extra column is added with band information
#'  If returning a data.frame, then it will be in the same order as 'genes'. If bioC=TRUE, then
#'  the result will be in genome order, regardless of the order of 'genes'.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.band, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' Pos.gene(c("CTLA4","PTPN22"))
#' Pos.gene("MICA",build=36)
#' Pos.gene("MICA",build=37)
#' Pos.gene(c("CTLA4","PTPN22"),bioC=TRUE,band=TRUE)
#' Pos.gene(c("CTLA4","OR2A1"),one.to.one=TRUE) # OR2A1 is split over two ranges
#' Pos.gene(c("CTLA4","OR2A1"),one.to.one=FALSE)
#' Pos.gene("RNU2-1",one.to.one=FALSE,bioC=T,remap.extra=FALSE) # RNU2-1 is split over several ranges
Pos.gene <- function(genes,build=NULL,dir=NULL,bioC=FALSE,band=FALSE,one.to.one=TRUE,
                     remap.extra=FALSE,discard.extra=TRUE,warnings=TRUE) {
  if(!is.character(genes)) { stop("'genes' must be a character vector of gene names") }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  char.lim <- 100
  ga <- get.gene.annot(dir=dir,build=build,bioC=bioC,one.to.one=one.to.one,
                       remap.extra=remap.extra,discard.extra=discard.extra,GRanges=TRUE)
  typ <- is(ga)[1]
  if(typ=="GRanges") {  mt <- match(genes,mcols(ga)$gene) } else { mt <- match(genes,ga$gene) }
  failz <- paste(genes[is.na(mt)],collapse=", "); if(nchar(failz)>char.lim) { failz <- paste(substr(failz,1,char.lim),",...",sep="") }
  if(length(mt)<1 | all(is.na(mt))) { 
    if(warnings) { warning("did not find any 'genes' features: ",failz) }; return(NULL) }
  if(any(is.na(mt))) { 
    cnt <- length(which(is.na(mt)))
    if(warnings) { warning("did not find the following ",cnt," 'genes' features: ",failz) }
  }
  if(!one.to.one) { 
    order.important <- FALSE  # not currently implemented
    if(order.important) {
     if(length(genes)>0) {
       mt <- NULL
       for(gg in 1:length(genes)) {
         if(typ=="GRanges") { 
           mt <- c(mt,which(mcols(ga)$gene %in% genes[gg]))
         } else {
           mt <- c(mt,which(ga$gene %in% genes[gg]))
         }
       }
     } else { warning("length of genes entered was zero?") }     
    } else {
      if(typ=="GRanges") { 
        mt <- which(mcols(ga)$gene %in% genes)
      } else {
        mt <- which(ga$gene %in% genes)
      }
    }
  }
  outlist <- ga[(narm(mt)),]
  if(typ=="GRanges") {
    if(!band) { mcols(outlist) <- mcols(outlist)[,-which(colnames(mcols(outlist)) %in% "band")]  }
    cnn <- colnames(mcols(outlist)); rnn <- mcols(outlist)$gene
  } else {
    if(!band) { outlist <- outlist[,-which(colnames(outlist) %in% "band")] }
    cnn <- colnames(outlist); rnn <- outlist$gene
  }
  if(one.to.one & ("gene" %in% cnn) & !anyDuplicated(genes)) {
    rownames(outlist) <- rnn
    if(all(genes %in% rownames(outlist))) {
      outlist <- outlist[match(genes,rnn),]
    }
    # outlist <- outlist[,-which(colnames(outlist) %in% "gene")]
  }
  if(bioC) { outlist <- as(outlist,"GRanges") }
  return(outlist)
}


#' Find the chromosome, start and end position for cytoband names
#' 
#' Allows retrieval of the the chromosome position of a karyotype/cytoband label, 
#'  or vector of such labels. Note that the position returned for bands is not a 
#'  single point as for SNPs, so the result will be a chromosome, then a position range with
#'  start and end, and lastly the band without the chromosome prefix
#' @param bands character, a vector of cytoband labels, chromosome[p/q]xx.xx ; 
#' e.g, 13q21.31, Yq11.221, 6p23, etc
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download cyto annotation information; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, else 
#' a data.frame
#' @return Returns a data.frame with columns 'chr' [chromosome], 'start' [starting position of the
#'  gene],'end' [end position of the gene] and 'band' [band without the chromosome prefix],
#'  or if bioC=TRUE, then returns a GRanges object with equivalent information.
#'  If returning a data.frame, then it will be in the same order as 'bands'. If bioC=TRUE, then
#'  the result will be in genome order, regardless of the order of 'bands'.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.gene, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' Pos.band("1p13.2")
#' Pos.band("Yq11.221",build=36)
#' Pos.band("Yq11.221",build=37)
#' Pos.band(c("13q21.31","1p13.2","2q33.2","6p23"),bioC=TRUE)
Pos.band <- function(bands,build=NULL,dir=NULL,bioC=FALSE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  char.lim <- 100
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.cyto(build=build,bioC=bioC,dir=dir,GRanges=FALSE)
  mt <- match(bands,rownames(ga))
  failz <- paste(bands[is.na(mt)],collapse=", "); if(nchar(failz)>char.lim) { failz <- paste(substr(failz,1,char.lim),",...",sep="") }
  msg <- ("format for bands is: chromosome[p/q]xx.xx ; e.g, 13q21.31, Yq11.221, 6p23, etc")
  if(length(mt)<1 | all(is.na(mt))) { 
    warning("did not find any 'bands' features: ",failz) ; warning(msg); return(NULL) }
  if(any(is.na(mt))) { 
    cat("format for bands is: chromosome[p/q]xx.xx ; e.g, 13q21.31, Xq27.1, 6p23, etc")
    cnt <- length(which(is.na(mt)))
    warning("did not find the following ",cnt," 'bands' features: ",failz) ; warning(msg)
  }
  outlist <- ga[sort(mt[!is.na(mt)]),]
  if(any(colnames(outlist) %in% "negpos")) { outlist <- outlist[,-which(colnames(outlist) %in% "negpos")] }
  if(all(bands %in% rownames(outlist)) & !bioC) {
    outlist <- outlist[bands,]
  }
  return(outlist)
}


#' Retrieve the cytoband(s) for snp ids, genes or locations
#' 
#' Allows retrieval of the the cytoband/karyotype label, based on multiple
#'  possible input featues, including SNP chip or rs-ids, HGNC gene labels, GRanges or
#'  RangedData object, chromosome and position vectors. The most robust way to use the
#'  function is to use the parameter names to imply the type of input, e.g, use the 'genes'
#'  parameter to input gene labels, the 'snps' parameter to enter SNP ids, etc. However,
#'  if you enter the first argument as a GRanges or RangedData object instead of using the
#'  'ranges' argument, this will be detected and automatically moved to the 'ranges' parameter.
#' @param genes character, an optional vector of gene ids, or RangedData/GRanges object
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions to retrieve the band from
#' @param ranges optional GRanges or RangedData object describing positions for which we want bands
#' @param snps optional SNP ids, e.g, chip ids or rs-ids, to retrieve the band they fall within
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download cyto annotation information; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param ... further arguments to Band.gene if entering gene names, or further arguments to 
#' Band.pos if entering ranges, or chr, pos/start/end
#' @return Returns a vector of bands, if any entries span more than one band, the bands will be
#' concatenated as character type, delimited by semicolons (;)
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.gene, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' Band(chr=1,pos=1234567) # using chr,pos vectors
#' rd <- RangedData(ranges=IRanges(start=87654321,end=87654321),space=1)
#' gr <- as(rd,"GRanges")
#' Band(rd)    # using RangedData, autodetects this parameter should be 'ranges' not 'genes'
#' Band(ranges=gr) # using GRanges
#' Band("SLC6A4")  # serotonin gene [5-HTT]
#' a.few.snps <- c("rs3842724","imm_11_2147527","rs9467354")
#' Band(a.few.snps) # using SNP ids in the 'genes' parameter (still works!)
#' Band(snps=a.few.snps) # using SNP ids with the dedicated 'snps' parameter is quicker
#' Band(chr="X",pos=8000000)
#' # Band() with longer ranges  #
#' Band(chr=12,start=40000000,end=50000000,build="hg19") # concatenates if range spans multiple bands
#' Band(chr=12,start=40000000,end=50000000,build="hg18") # one extra band in the older annotation
Band <- function(genes=NULL,chr=NULL,ranges=NULL,snps=NULL,build=NULL,dir=NULL,...) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(all(is.character(genes))) { 
    Band.gene(genes=genes,build=build,dir=dir,...,warnings=FALSE)
  } else {
    if(!is.null(snps)) {
      pp <- Pos(snps); ch <- Chr(snps)
      Band.pos(chr=ch,pos=pp,build=build,dir=dir)
    } else {
      if((is(genes)[1] %in% c("RangedData","GRanges")) & is.null(ranges) ) { ranges <- genes }
      Band.pos(chr=chr,ranges=ranges,build=build,dir=dir,...)
    }
  }
}


#' Retrieve the cytoband(s) for genes labels
#' 
#' Allows retrieval of the the cytoband/karyotype label for HGNC gene labels.
#' @param genes character, an optional vector of gene ids, or RangedData/GRanges object
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download cyto annotation information; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param append.chr logical, it is typical that the chromosome character preceeds cytoband labels,
#'  but if this parameter is set to FALSE, it will be left off.
#' @param data.frame logical, if data.frame is true, instead of returning a vector of full cytoband
#' labels, a data.frame will be returned.
#' @param warnings logical, if warnings=FALSE and SNP ids are entered instead of Gene labels,
#' then the function will automatically detect this and return the result of Band(snps='genes')
#' @return Returns a vector of bands, if any entries span more than one band, the bands will be
#' concatenated as character type, delimited by semicolons (;). If data.frame is true, instead of 
#' returning a vector of full cytoband labels, a data.frame will be returned with a 'chr'
#' [chromosome] column, 'band' cytoband label  without the chromosome prefix, and rownames 
#' equal to 'genes'
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.gene, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' a.few.snps <- c("rs3842724","imm_11_2147527","rs9467354")
#' Band.gene("HLA-C") # using chr,pos vectors
#' Band.gene(a.few.snps)  # fails with warning as these are SNPs, not genes
#' Band.gene(a.few.snps,warnings=FALSE) # with warnings=FALSE this recognises snps were entered and carries on
#' Band.gene("SLC6A4")  # serotonin gene [5-HTT]
#' Band.gene("SLC6A4",append.chr=F)
#' Band.gene("SLC6A4",data.frame=T)
Band.gene <- function(genes,build=NULL,dir=getwd(),append.chr=TRUE,data.frame=FALSE,warnings=TRUE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  pg <- Pos.gene(genes=genes,build=build,dir=dir,bioC=F,band=TRUE,one.to.one=TRUE,warnings=warnings)
  ## a bit contradictory what i've done here with 'warnings' :
  if(length(pg)<1 & !warnings) { 
    #if(warnings) { warning("perhaps parameter 'genes' should have been 'snps'?") } 
    return(Band(snps=genes)) 
  }
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



#' Find the gene(s) overlapping a chromosome location
#' 
#' Allows retrieval of genes intersected by a chromosome and position, which can be entered
#' using chr, pos/start/end vectors, or a RangedData or GRanges object
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions to retrieve the possible overlapping gene(s)
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' start or end if this is entered, and vice-versa
#' @param start integer, an optional vector of start points for chromosome ranges
#' @param end integer, an optional vector of end points for chromosome ranges
#' @param ranges optional GRanges or RangedData object describing positions for which we want genes,
#' removing the need to enter chr, pos, start or end
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download gene annotation information to; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, or
#' RangedData if 'ranges' is RangedData, else a data.frame
#' @param one.to.one logical, whether to concatenate multiple hits for the same range into one result,
#' or spread the result over multiple lines, one for each gene overlapped
#' @return Returns a set of genes separated by semicolons (if more than one) for each range entered.
#' If bioC=TRUE, returns the equivalent as a GRanges object, unless a RangedData object was used
#' for the ranges parameter, in which case a RangedData object would be returned. If one.to.one is
#' FALSE, then instead of concatenating multiple genes into one line per range, each is listed 
#' separately as a new row, with an index added to correspond to the original input order of ranges,
#' if bioC=TRUE; or just adds additional elements to the resulting vector if bioC=FALSE.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.band, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' Gene.pos(chr=6, start=31459636, end=31462760)
#' Gene.pos(chr=22, pos=3452345) # no gene here
#' Gene.pos(Chr("rs689"),Pos("rs689")) # combine with Chr() and Pos() to find gene(s) for a the rs689 SNP
#' Gene.pos(chr=1,start=114000000,end=115000000,build="hg19") # multiple genes in range
#' Gene.pos(chr=1,start=114000000,end=115000000,one.to.one=FALSE) # list separately
#' Gene.pos(Pos.gene(c("CTLA4","PTPN22"),bioC=TRUE)) # use the ranges object returned by Pos.gene()
Gene.pos <- function(chr=NA,pos=NA,start=NA,end=NA,ranges=NULL,build=NULL,dir=NULL,bioC=FALSE,one.to.one=TRUE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is(chr)[1] %in% c("RangedData","GRanges"))  { ranges <- chr } # in case first parameter is used
  typ <- is(ranges)[1]
  if(!typ %in% c("RangedData","GRanges")) {
    if(any(!is.na(pos))) { start <- pos; end <- pos }
    if(length(chr)==1 & length(start)>1) { chr <- rep(chr,times=length(start)) }
    if(length(chr)!=length(start)) { stop("chr vector must have same length as pos or start/end") }
    if(any(is.na(chr))) { stop("cannot have chr=NA") }
    Pos <- matrix(ncol=2,nrow=length(start))
    for (cc in 1:length(start)) {
      Pos[cc,] <- force.chr.pos(Pos=c(start[cc],end[cc]),Chr=chr[cc],dir=dir,build=build)
    }
    #if(any(tolower(substr(chr,1,3))!="chr")) { chr <- gsub("chr",chr,sep="") }
    #chr <- gsub("chrchr","chr",chr)
    testData <- RangedData(ranges=IRanges(start=Pos[,1],end=Pos[,2]),space=chr,index=1:length(chr),universe=build[1])
    testData <- toGenomeOrder2(testData,strict=T)
  } else {
    testData <- ranges # set.chr.to.char(ranges)
    if(typ=="GRanges") { testData <- as(ranges,"RangedData") }
    #chr <- chr(testData)
    if(!bioC) { warning("bioC was set false, but ranges argument in use so overriding") ; bioC <- T }
    if("index" %in% colnames(testData)) { 
      warning("'index' is a reserved column name for ranges objects passed to this function so will be replaced. Consider renaming this column if this is undesired") }
    testData[["index"]] <- 1:nrow(testData)
  }
  testData <- set.chr.to.numeric(testData,keep=T)
  #return(testData)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.gene.annot(build=build,dir=dir,GRanges=FALSE)
  #ga <- set.chr.to.numeric(ga,keep=F)
  newDataList <- vector("list",nrow(testData))
  overlaps <- findOverlaps(testData,ga)
  genez <- ga$gene[subjectHits(overlaps)]
  indexz <- queryHits(overlaps)
  #prv(indexz,genez)
  if(length(indexz)<1) { return(NA) }
  if(!one.to.one) {
    newData <- testData[queryHits(overlaps),]
    newData[["gnm.index"]] <- indexz
    newData[["gene"]] <- genez
  } else {
    out <- tapply(genez,factor(indexz),c,simplify=FALSE)
    out <- sapply(out,function(X) { X <- narm(unique(X)); X <- X[X!=""] ; paste(X,collapse=";") })
    newData <- testData
    newData[["gene"]] <- rep("intergenic",nrow(newData))
    if(!is.null(names(out))) {
      newData[["gene"]][as.numeric(names(out))] <- out
    } else {
      newData[["gene"]] <- out
    }
  }
  if(bioC) {
    if(typ!="RangedData") { newData <- as(newData,"GRanges") }
    return(newData)
  } else {
    OO <- newData[["gene"]][order(newData[["index"]])]
    OO <- narm(unique(OO)); OO <- OO[OO!=""]
    return(OO)
  }
  #if(!all(chr %in% chr2(ga))) { stop("invalid chromosome(s) entered") } # redundant i think
}





#' Find the cytoband(s) overlapping a chromosome location
#' 
#' Allows retrieval of cytobands/karyotypes intersected by a chromosome and position, which can be 
#' entered using chr, pos/start/end vectors, or a RangedData or GRanges object
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions to retrieve the possible overlapping cytoband(s)
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' start or end if this is entered, and vice-versa
#' @param start integer, an optional vector of start points for chromosome ranges
#' @param end integer, an optional vector of end points for chromosome ranges
#' @param ranges optional GRanges or RangedData object describing positions for which we want bands,
#' removing the need to enter chr, pos, start or end
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download gene annotation information to; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, or
#' RangedData if 'ranges' is RangedData, else a data.frame
#' @param one.to.one logical, whether to concatenate multiple hits for the same range into one result,
#' or spread the result over multiple lines, one for each cytoband overlapped
#' @return Returns a set of cytobands separated by semicolons (if more than one) for each range entered.
#' If bioC=TRUE, returns the equivalent as a GRanges object, unless a RangedData object was used
#' for the ranges parameter, in which case a RangedData object would be returned. If one.to.one is
#' FALSE, then instead of concatenating multiple cytobands into one line per range, each is listed 
#' separately as a new row, with an index added to correspond to the original input order of ranges,
#' if bioC=TRUE; or just adds additional elements to the resulting vector if bioC=FALSE.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.band, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' Band.pos(chr=6, start=31459636, end=31462760)
#' Band.pos(chr=22, pos=3452345) 
#' Band.pos(Chr("rs689"),Pos("rs689")) # combine with Chr() and Pos() to find the cytoband for a the rs689 SNP
#' Band.pos(chr=1,start=110000000,end=120000000,build="hg19") # multiple cytobands in range
#' Band.pos(chr=1,start=110000000,end=120000000,one.to.one=FALSE) # list separately
#' Band.pos(Pos.band(c("13q21.31","1p13.2"),bioC=TRUE)) # use the ranges object returned by Pos.band()
#' # note that three ranges are returned for each entry as the start/end overlap the adjacent ranges
Band.pos <- function(chr=NA,pos=NA,start=NA,end=NA,ranges=NULL,build=NULL,dir=NULL,bioC=FALSE,one.to.one=TRUE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is(chr)[1] %in% c("RangedData","GRanges"))  { ranges <- chr } # in case first parameter is used
  typ <- is(ranges)[1]
  if(!typ %in% c("RangedData","GRanges")) {
    if(any(!is.na(pos))) { start <- pos; end <- pos }
    if(length(chr)==1 & length(start)>1) { chr <- rep(chr,times=length(start)) }
    if(length(chr)!=length(start)) { stop("chr vector must have same length as pos or start/end") }
    if(any(is.na(chr))) { stop("cannot have chr=NA") }
    Pos <- matrix(ncol=2,nrow=length(start))
    for (cc in 1:length(start)) {
      Pos[cc,] <- force.chr.pos(Pos=c(start[cc],end[cc]),Chr=chr[cc],dir=dir,build=build)
    }
    #if(any(tolower(substr(chr,1,3))!="chr")) { chr <- gsub("chr",chr,sep="") }
    #chr <- gsub("chrchr","chr",chr)
    testData <- RangedData(ranges=IRanges(start=Pos[,1],end=Pos[,2]),space=chr,index=1:length(chr),universe=build[1])
    testData <- toGenomeOrder2(testData,strict=T)
  } else {
    if(typ=="GRanges") { ranges <- as(ranges,"RangedData") }
    testData <- ranges # set.chr.to.char(ranges)
    if("index" %in% colnames(testData)) { warning("'index' is a reserved column name for ranges objects passed to this function so will be replaced. Consider renaming this column if this is undesired") }
    testData[["index"]] <- 1:nrow(testData)
  }
  testData <- set.chr.to.numeric(testData)
  #return(testData)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  cyto <- get.cyto(build=build,bioC=TRUE,dir=dir,GRanges=FALSE)
  cyto <- set.chr.to.numeric(cyto,keep=F)
  #if(any(colnames(cyto) %in% "negpos")) { cyto <- cyto[,-which(colnames(cyto) %in% "negpos")] }
  newDataList <- vector("list",nrow(testData))
  overlaps <- findOverlaps(testData,cyto)
  bandz <- rownames(cyto)[subjectHits(overlaps)]
  indexz <- queryHits(overlaps)
  if(length(indexz)<1) { return(NA) }
  if(!one.to.one) {
    newData <- testData[queryHits(overlaps),]
    newData[["gnm.index"]] <- indexz
    newData[["band"]] <- bandz
  } else {
    out <- tapply(bandz,factor(indexz),c,simplify=FALSE)
    if(one.to.one) { 
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
    if(typ!="RangedData") { newData <- as(newData,"GRanges") }
    return(newData)
  } else {
    return(newData[["band"]][order(newData[["index"]])])
  }
  #if(!all(chr %in% chr2(ga))) { stop("invalid chromosome(s) entered") } # redundant i think
}



#' Returns the A and B allele for SNP ids
#' 
#' For a set of chip ids or rs ids, returns a two column matrix containing the A and B allele. 
#' For snpStats objects the default is that A,B are coded in alphabetical order, so A,C; A,T; 
#' C,T; C,G are possible A,B pairs
#' @param ids character, a list of chip ids or rs-ids as contained in the current ChipInfo object
#' @return Returns a two column matrix containing the A and B allele.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.band, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' snp.ids <- c("rs3842724","rs9729550","rs1815606","rs114582555","rs1240708","rs6603785")
#' AB(snp.ids) 
AB <- function(ids) {
  all.support <- get.support()
  ic.ab <- function(ic.ids) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { all.support <- get.support() }  ## load object: all.support [snp support for whole chip]
    outlist <- cbind(A1(all.support)[match(ic.ids,rownames(all.support))],A2(all.support)[match(ic.ids,rownames(all.support))])
    return(outlist)
  }
  out <- matrix(nrow=length(ids),ncol=2)
  ic.ab.id <- ic.ab(rs.to.id(ids))
  return(ic.ab.id)
}


#' Convert ensembl ids to HGNC gene ids 
#' 
#' Retrieve the gene IDs (HGNC) corresponding to a list of ensembl gene ids.
#' Note that this will not find all IDs found on ensembl.org, as it uses bioMart which
#' seems to be incomplete, but this only pertains to a small minority of genes, so this
#' function should have general utility for most applications. This is of course the case
#' at the time of writing - bioMart is likely to be updated at some point.
#' @param ens character, a list of ensembl gene ids, of the form ENSG00xxxxxxxxx
#' @param ... further arguments to get.gene.annot()
#' @param dir character, 'dir' is the location to download gene and cytoband information; if
#' left as NULL, depending on the value of getOption("save.annot.in.current"), the annotation
#' will either be saved in the working directory to speed-up subsequent lookups, or deleted 
#' after use.
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param name.dups logical, if TRUE then duplicates will have a suffix appended to force the
#' list to be unique (e.g, so it would be usable as rownames, or in a lookup table). Otherwise
#' duplicate entries will just appear in the list multiple times
#' @param name.missing logical, if TRUE then missing values will be named as MISSING_n (n=1
#'  to # of missing), ensuring a valid unique name if the results are to be used as rownames,
#' etc. If FALSE then these will be left as NA. 
#' @return Returns a vector of HGNC gene ids corresponding to the 'ens' ensembl ids entered,
#' any ids not found will be returned as MISSING_n (n=1 to # of missing), if name.missing=TRUE.
#' If name.missing is FALSE then missing will be set to NA. Similarly with 'name.dups', if
#' duplicates are found and name.dups is true, each will be appended with suffix _n; else
#' their names will be left as is.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso GENE.to.ENS, rs.to.id, id.to.rs
#' @examples
#' ENS.ids <- c("ENSG00000183214", "ENSG00000163599", "ENSG00000175354", "ENSG00000134460")
#' ENS.to.GENE(ENS.ids)
#' gene.ids <- c("HLA-B","IFIH1","fake_gene!","FUT2")
#' ENS.to.GENE(GENE.to.ENS(gene.ids)) # lookup fails for the fake id, gives warning
ENS.to.GENE <- function(ens,...,dir=NULL,build=NULL,name.dups=FALSE,name.missing=TRUE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  must.use.package(c("biomaRt","genoset","gage"),T)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.gene.annot(...,dir=dir,bioC=FALSE,ens.id=TRUE,GRanges=FALSE)
  ### now have the gene data with the ensembl ids ##
  #print(head(ga))
  indx <- match(ens,ga$ens.id)
  missin <- length(which(is.na(indx))); valid <- length(indx)-missin
  if(valid<1) { warning("did not find any ENSEMBL ids from 'ens' in the bioMart human gene reference"); return(NULL) }
  if(missin>0) { warning(out.of(missin,(valid+missin))," of 'ens' did not match any ENSEMBL ids in the bioMart human gene reference") }
  outData <- ga$gene[indx]
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


#' Convert gene ids to ensembl ids
#' 
#' Retrieve the ensembl IDs corresponding to a list of common gene names (HGNC format).
#' @param genes character, gene labels, e.g, "APOE"
#' @param ... further arguments to get.gene.annot()
#' @param dir character, 'dir' is the location to download gene and cytoband information; if
#' left as NULL, depending on the value of getOption("save.annot.in.current"), the annotation
#' will either be saved in the working directory to speed-up subsequent lookups, or deleted 
#' after use.
#' @return Returns a vector of HGNC gene ids corresponding to the 'ens' ensembl ids entered,
#' any ids not found will be returned as NA.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso GENE.to.ENS, rs.to.id, id.to.rs
#' @examples
#' gene.ids <- c("MICA","PTPN2","IL2RA","APOE")
#' GENE.to.ENS(gene.ids)
GENE.to.ENS <- function(genes,...,dir=NULL) {
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.gene.annot(...,dir=dir,bioC=FALSE,ens.id=TRUE,GRanges=FALSE)
  ### now have the gene data with the ensembl ids ##
  indx <- match(genes,ga$gene)
  missin <- length(which(is.na(indx))); valid <- length(indx)-missin
  if(valid<1) { warning("did not find any gene ids from 'genes' in the bioMart human gene reference"); return(NULL) }
  if(missin>0) { warning("at least one of 'genes' did not match any gene ids in the bioMart human gene reference") }
  outData <- ga$ens.id[indx]
  return(outData)
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
#' @param ids logical, if TRUE will return snp ids (chip ids, for rs-ids, use id.to.rs on the output), 
#' or if FALSE will return the chromosome positions of the SNPs.
#' @param build integer whether to use build 36/37 parameters, 36/37 is preferred, but can enter
#' using any form recognised by ucsc.sanitizer()
#' @param dir string, specify the directory to store or look for annotation (defaults to current)
#' @export
#' @return Set of SNP ids (when ids=TRUE), or otherwise genomic positions within chromosome 'chr', that
#' fall within the genomic range described by the chr, start, and end parameters. Alternatively, if
#' chr is a RangedData or GRanges object then multiple SNP lists will be returned
#' in a list for each row of the ranges object.
#' @examples
#' snps.in.range(1,9000000,10000000)
#' snps.in.range(10,19000000,20000000,ids=T)
#' snps.in.range(10,19000000,20000000,ids=F) # return positions instead of rs-ids
snps.in.range <- function(chr, start=NA, end=start, ids=TRUE) { 
  # ids - whether to return ichip SNP ids or positions
  if(is(chr)[1]=="RangedData" | is(chr)[1]=="GRanges") {
    chrz <- chr2(chr); stz <- start(chr); enz <- end(chr)
    output <- vector("list",nrow(chr))
    for(cc in 1:nrow(chr)) {
      output[[cc]] <- snps.in.range(chrz[cc],stz[cc],enz[cc],ids=ids)
    }
    names(output) <- rownames(chr)
    return(output)
  }
  if(length(chr)>1) { warning("chr should be length 1, using only first entry"); chr <- chr[1] }
  if(length(start)>1) { warning("start should be length 1, using only first entry"); start <- start[1] }
  if(length(end)>1) { warning("end should be length 1, using only first entry"); end <- end[1] }
  if(start>end) { warning("start was higher than end, so switching") }
  the.range <- sort(c(start,end))
  all.support <- get.support()
  #if(!exists("work.dir")) { if(is.null(dir)) { work.dir <- getwd() } else { work.dir <- dir } }
  if(!exists("all.support")) { all.support <- get.support() }  ## load object: all.support [snp support for whole chip]
  all.chr <- chr(all.support)
  all.pos <- start(all.support)[all.chr %in% chr]
  if(length(all.pos)<1) { warning("no positions found for 'chr' specified"); return(NULL) }
  validz <- which(all.pos>=the.range[1] & all.pos<=the.range[2])
  if(ids) {
    out <- rownames(all.support)[all.chr %in% chr][validz]
  } else {
    out <- start(all.support)[(all.chr %in% chr)][validz]
  }
  return(out)
}


#' Retrieve the 'n' closest SNP ids or positions near specified locus
#' 
#' @param chr integer, chromosome, should be a number from 1 to 25, where 23,24,25 are X,Y,MT
#' @param pos integer, genomic position, should be between 1 and the length of the chromosome 'chr'
#' @param n integer, the number of nearest SNPs to seek, if there aren't enough in the annotation
#' then NAs will fill the gaps to force the return value length to equal 'n'
#' @param side character, can be 'either', 'left' or 'right' and specifies which side of the 'pos'
#' to look for nearest snps (where left is decreasing genomic position and right is increasing)
#' @param ids logical, if TRUE will return snp ids (chip ids, for rs-ids, use id.to.rs on the output), 
#' or if FALSE will return the chromosome positions of the SNPs.
#' @param limit integer, a limit on the maximum distance from the position 'pos' can be specified
#' @param build integer whether to use build 36/37 parameters, 36/37 is preferred, but can enter
#' using any form recognised by ucsc.sanitizer()
#' @param dir string, specify the directory to store or look for annotation (defaults to current)
#' @export
#' @seealso exp.window.nsnp, nearest.gene
#' @return Set of SNP ids (when ids=TRUE), or otherwise genomic positions within chromosome 'chr'.
#' If the number of SNPs on the chromosome or the bounds of the 'side' and 'limit' parameters
#' restrict the number returned to less than 'n' then the return value will be padded with NAs.
#' @examples
#' nearest.snp(1,159000000,n=10) # return ids
#' nearest.snp(1,159000000,n=10,build=37)
#' nearest.snp(1,159000000,n=10,build=36,ids=F) # return positions
#' nearest.snp(1,159000000,n=10,build=37,ids=F)
#' nearest.snp(6,25000000,n=10,build=37,ids=F,side="left")  # only SNPs to the left of the locus
#' nearest.snp(6,25000000,n=10,build=37,ids=F,side="right") # only SNPs to the right of the locus
nearest.snp <- function(chr, pos, n=1, side=c("either","left","right"),ids=TRUE,limit=NULL,build=NULL) { 
  # ids - whether to return ichip SNP ids or positions
  if(length(chr)>1) { warning("chr should be length 1, using only first entry"); chr <- chr[1] }
  if(length(pos)>1) { warning("pos should be length 1, using only first entry"); pos <- pos[1] }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  all.support <- get.support(build=build)
  if(!exists("all.support")) { all.support <- get.support() }  ## load object: all.support [snp support for whole chip]
  side <- tolower(side[1]); 
  if(!side %in% c("either","left","right")) {
    side <- "either"; warning("invalid side argument, defaulting to 'either'") }
  if(!is.null(limit)) { if(!is.numeric(limit)) { limit <- NULL; warning("invalid limit argument, defaulting to NULL") } }
  all.chr <- chr(all.support)
  all.pos <- start(all.support)[all.chr %in% chr]
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
    out <- start(all.support)[(all.chr %in% chr)][filt][indx]
  }
  return(out)
}


#' Retrieve the 'n' closest GENE labels or positions near specified locus
#' 
#' @param chr integer, chromosome, should be a number from 1 to 25, where 23,24,25 are X,Y,MT
#' @param pos integer, genomic position, should be between 1 and the length of the chromosome 'chr'
#' @param n integer, the number of nearest GENEs to seek, if there aren't enough in the annotation
#' then NAs will fill the gaps to force the return value length to equal 'n'
#' @param side character, can be 'either', 'left' or 'right' and specifies which side of the 'pos'
#' to look for nearest genes (where left is decreasing genomic position and right is increasing)
#' @param ids logical, if TRUE will return GENE labels, 
#' or if FALSE will return the chromosome positions of the genes
#' @param limit integer, a limit on the maximum distance from the position 'pos' can be specified
#' @param build integer whether to use build 36/37 parameters, 36/37 is preferred, but can enter
#' using any form recognised by ucsc.sanitizer()
#' @param dir string, specify the directory to store or look for annotation (defaults to current)
#' @param ga RangedData object, e.g, result of get.gene.annot(); gene annotation to save download
#' time if repeatedly calling this function
#' @export
#' @seealso exp.window.nsnp, nearest.snp, get.gene.annot
#' @return Set of GENE ids (when ids=TRUE), or otherwise genomic positions within chromosome 'chr'.
#' If the number of gemes on the chromosome or the bounds of the 'side' and 'limit' parameters
#' restrict the number returned to less than 'n' then the return value will be padded with NAs.
#' @examples
#' nearest.gene(1,159000000,n=10) # return ids
#' nearest.gene(1,159000000,n=10,build=37)
#' nearest.gene(1,159000000,n=10,build=36,ids=F) # return positions
#' nearest.gene(1,159000000,n=10,build=37,ids=F)
#' nearest.gene(6,25000000,n=10,build=37,ids=F,side="left")  # only genes to the left of the locus
#' nearest.gene(6,25000000,n=10,build=37,ids=F,side="right") # only genes to the right of the locus
nearest.gene <- function(chr, pos, n=1, side=c("either","left","right"),ids=TRUE,limit=NULL,build=NULL, ga=NULL) { 
  # ids - whether to return ichip SNP ids or positions
  if(length(chr)>1) { warning("chr should be length 1, using only first entry"); chr <- chr[1] }
  if(length(pos)>1) { warning("pos should be length 1, using only first entry"); pos <- pos[1] }
  if(is.null(build)) { build <- getOption("ucsc") }
  chrom <- paste(chr)
  build <- ucsc.sanitizer(build)
  if(is(get.gene.annot)[1]=="RangedData") { 
    if(!"gene" %in% colnames(ga)) { ga <- NULL }
  } else { ga <- NULL }
  if(is.null(ga)) {  
    ga <- get.gene.annot(build=build,GRanges=F) 
    if(!exists("ga")) { stop("couldn't find gene annotation") }  ## load object: ga [gene database]
    ga <- ga[ga$gene!="",]
  }
  side <- tolower(side[1]); 
  if(!side %in% c("either","left","right")) {
    side <- "either"; warning("invalid side argument, defaulting to 'either'") }
  if(!is.null(limit)) { if(!is.numeric(limit)) { limit <- NULL; warning("invalid limit argument, defaulting to NULL") } }
  all.chr <- paste(chr2(ga))
  all.st <- start(ga)[all.chr %in% chrom]
  all.en <- end(ga)[all.chr %in% chrom]
  #prv(all.st,all.en)
  if(length(all.st)<1) { warning("no positions found for 'chr' specified"); return(NULL) }
  difzS <- pos-all.st
  difzE <- pos-all.en
  all.true <- difzS==difzS
  if(is.null(limit)) { lfilt <- all.true } else { lfilt <- abs(difzS)<=limit | abs(difzE)<=limit }
  within <- difzS>0 & difzE<0
  if(side=="left") { filt <- ( difzE>0 & lfilt ) | within }
  if(side=="right") { filt <- ( difzS<0 & lfilt ) | within }
  if(side=="either") { filt <- ( all.true & lfilt ) | within }
  #print(length(which(filt)))
  tab <- rbind(abs(difzS),abs(difzE))
  minz <- apply(tab,2,min,na.rm=T)
  Difz <- abs(minz[filt])
  if(length(Difz)<n)  { warning("fewer than ",n," genes found for 'chr' specified (within 'limit'), NAs returned") }
  indx <- order(Difz)[1:n]
 # prv(minz,Difz,filt,indx)
  subi <- ga[["gene"]][all.chr %in% chrom][filt]
  #prv(ga,subi,all.chr,chrom)
  if(ids) {
    out <- ga[["gene"]][all.chr %in% chrom][filt][indx]
  } else {
    out <- start(ga)[(all.chr %in% chrom)][filt][indx]
  }
  return(out)
}


#internal
# standardize snp ids so they would always appear the same, all _,.;, etc replaced with _
# all names leading with a number preceeded with X. mainly 'make.names' standard R-conventions
clean.snp.ids <- function(snpid.list) {
  snpid.list <- make.names(snpid.list)
  snpid.list <- gsub(".","_",snpid.list,fixed=T)
  return(snpid.list)
}






#' Wrapper to construct GRanges object from chr,pos or chr,start,end
#' 
#' Slightly simplifies the creation of a GRanges object, allowing flexible input of
#' chr, pos, or chr,start,end, and specification of rownames and the 'genome' parameter
#' for specifying the build/coordinate type, e.g, hg18, build 37, etc. Designed for
#' a simplified GRanges object without metadata, and where the 'strand' data is of
#' no interest, so if strand/metadata is to be used, use the original GRanges() constructor.
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions for the GRanges object
#' @param pos integer/numeric, for SNPs, can enter positions just once in 'pos' instead of entering the same value
#' for start and end
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param start integer/numeric, specify the start position of ranges to encode in the new
#' GRanges object, alongside 'end' (do not use 'pos' if using start+end)
#' @param end integer/numeric, specify the end position of ranges to encode in the new
#' GRanges object, alongside 'start' (do not use 'pos' if using start+end)
#' @return Returns a GRanges object with the ranges, build and rownames specified. Rownames
#' will be 1:nrow if the 'row.names' parameter is empty. The strand information will default
#' to '+' for all entries, and the metadata will be empty (this function is only for creation
#' of a very basic GRanges object).
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.gene, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' g1 <- make.granges(chr=c(1,4,"X"),pos=c(132432,434342,232222))
#' g2 <- make.granges(chr=c(22,21,21),start=c(1,1,1),end=c(1000,10000,100000),row.names=c("1K","10K","100K"))
#' g1 ; g2
make.granges <- function(chr,pos=NULL,start=NULL,end=NULL,row.names=NULL,build=NULL) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is.null(start) & is.null(end) & !is.null(pos)) {
    dF <- cbind(paste(chr),round(as.numeric(pos)))
  } else {
    if(!is.null(start) & !is.null(end) & is.null(pos)) {
      pos <- cbind(round(as.numeric(start)),round(as.numeric(end)))
      dF <- cbind(paste(chr),pos) 
    } else {
      stop("must use either 'pos' or 'start' and 'end'")
    }
  }
  #prv(paste(row.names)); prv(dF)
  if(is.character(row.names)) { if(length(row.names)==nrow(dF)) { rownames(dF) <- row.names } else { warning("row.names had an incorrect length")} }
  if(!any(Dim(pos)==length(chr))) { stop("chr and pos must be of the same length") }
  if(length(Dim(pos))>1) { 
    if(!ncol(pos) %in% c(1,2)) { 
      stop("pos must be a vector of SNP locations, or a 2-column object with start and end coordinates") 
    } else {
      if(ncol(pos)==2) {
        if(any(pos[,2]<pos[,1])) { warning("end coordinates should be equal or greater than start coordinates") }
      }
    }
  }
  if(ncol(dF)==3) { 
    colnames(dF) <- c("chr","start","end") } else { colnames(dF) <- c("chr","pos") }
  #return(dF)
  ranged <- data.frame.to.granges(dF,start=colnames(dF)[2],end=tail(colnames(dF),1),build=build) 
  if(!any(rownames(ranged) %in% row.names)) {
    if(is.character(row.names)){ if(length(row.names)==nrow(ranged)) { rownames(ranged) <- row.names }}
  }
  return(ranged)
}


## internal function with extra mapping hits for immunochip that aren't in the chain file for 36-37
hard.coded.conv <- function() {
  chrzM <- c("7","7","9","5","7","14","17","4","8","8","15","7","6","6","2","2","4","17","19")
  pos36M <- c("142154515","142160115","132183222","17767156","141943232","27591752","41560151",
              "103951975","17510484","17501697","81350958","141911612","74644736","74644390",
              "1203295","21043693","4020119","78644427","52569727")
  pos37M <- c("142474939","142480539","135153668","17731427","142224511","28521898","44204373",
              "103732866","17466212","17457420","83559954","142108941","74588007","74587661",
              "1213294","21190209","3969218","81051007","47877928")
  rsidM <- c("rs10952532","rs10952534","rs11243704","rs11953245","rs17274","rs1952843",
             "rs2016730","rs223413","rs2427715","rs2517168","rs2621228",
             "rs2855938","rs2917890","rs2917891","rs4971417","rs6547409","rs6842556",
             "rs7502442","rs755327")
  chrzI <- c("3","3","3","3","3","6","7","7","8","8","8","8","17","17","X")
  pos36I <- c("50875374","50882163","50885514","50908888","195567372","119257505",
              "50323690","67383261","10961083","10961130","10975096","10975127",
              "21628754","59781521","75211826")
  pos37I <- c("50900354","50907147","50910499","50908888","194086083","119150813",
              "50353144","67745402","10923673","10923720","10937686","10937717",
              "21704627","59781521","75295444")
  rsidI <- c("rs12639243","rs62717061","rs4346541","imm_3_50908888",
             "rs4974514","rs284919","rs7804185","rs3113138",
             "rs2898255","rs2409687","rs7827367","rs6601557",
             "rs17052332","rs1131012","rs929032")
  chrz <- c(chrzI,chrzM)
  pos36 <- c(pos36I,pos36M)
  pos37 <- c(pos37I,pos37M)
  rsid <- c(rsidI,rsidM)
  return(list(chr=chrz,pos36=pos36,pos37=pos37,rs.id=rsid))
}



#' Convert from build 37 to build 36 SNP coordinates
#' 
#' Convert range or SNP coordinates between builds using a chain file. Depending on the chain file
#' this can do any conversion, but the default will use the hg19 to hg18 (37-->36) chain file
#' built into this package. The positions to convert can be entered using using chr, pos vectors,
#'  or a RangedData or GRanges object. This function is a wrapper for liftOver() from rtracklayer,
#' providing more control of input and output and 'defensive' preservation of order and length
#' of the output versus the input ranges/SNPs. 
#' @param chr character, an optional vector of chromosomes to combine with 'pos' to describe
#'  positions to convert from build hg19 to hg18
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' a ranges object if this is provided along with 'chr'
#' @param ranges optional GRanges or RangedData object describing positions for which conversion
#' should be performed. No need to enter chr, pos if using ranges
#' @param ids if the ranges have ids (e.g, SNP ids, CNV ids), then by including this parameter
#' when using chr, pos input, the output object will have these ids as rownames. For ranges input
#' these ids would already be in the rownames of the GRanges or RangedData object, so use of
#' this parameter should be unnecessary
#' @param ... additional arguments to make.granges(), so in other words, can use 'start' and
#' 'end' to specify ranges instead of 'pos'.
#' @return Returns positions converted from build 37 to 36. If using the 'ranges' parameter 
#' for position input, the object returned will be of the same format. If using chr and pos 
#' to input, then the object returned will be a data.frame with columns, chr and pos with 
#' rownames 'ids'. Output will be the same length as the input, which is not necessarily the
#'  case for liftOver() which does the core part of this conversion. Using vector or GRanges 
#'  input will give a resulting data.frame or GRanges object respectively that has the same
#'  order of rownames as the original input. Using RangedData will result in an output that
#'   is sorted by genome order, regardless of the original order.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso conv.36.37, convTo37, convTo36
#' @examples
#' gene.labs <- c("CTLA4","IL2RA","HLA-C")
#' pp <- Pos.gene(gene.labs,build=37)
#' gg <- GRanges(ranges=IRanges(start=pp$start,end=pp$end),seqnames=pp$chr)
#' conv.37.36(gg) # order of output is preserved   ### HERE!!! ###
#' rr <- as(gg,"RangedData")
#' conv.37.36(rr) # note the result is same as GRanges, but in genome order
conv.37.36 <- function(ranges=NULL,chr=NULL,pos=NULL,...,ids=NULL) {
  chain.file <- "/home/oliver/R/stuff/hg19ToHg18.over.chain"
  return(conv.36.37(ranges=ranges,chr=chr,pos=pos,...,ids=ids,chain.file=chain.file))
}



#' Convert from build 36 to build 37 SNP coordinates
#' 
#' Convert range or SNP coordinates between builds using a chain file. Depending on the chain file
#' this can do any conversion, but the default will use the hg18 to hg19 (36-->37) chain file
#' built into this package. The positions to convert can be entered using using chr, pos vectors,
#'  or a RangedData or GRanges object. This function is a wrapper for liftOver() from rtracklayer,
#' providing more control of input and output and 'defensive' preservation of order and length
#' of the output versus the input ranges/SNPs. 
#' @param chr character, an optional vector of chromosomes to combine with 'pos' to describe
#'  positions to convert to an alternative build
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' a ranges object if this is provided along with 'chr'
#' @param ranges optional GRanges or RangedData object describing positions for which conversion
#' should be performed. No need to enter chr, pos if using ranges
#' @param ids if the ranges have ids (e.g, SNP ids, CNV ids), then by including this parameter
#' when using chr, pos input, the output object will have these ids as rownames. For ranges input
#' these ids would already be in the rownames of the GRanges or RangedData object, so use of
#' this parameter should be unnecessary
#' @param chain.file character, a file location for the liftOver chain file to use for the
#' conversion. If this argument is left blank the default UCSC file that converts from hg18
#' to hg19 will be used. However chain files for other conversions are available from
#' http://crossmap.sourceforge.net/, and you could also customize these or create your own.
#' So this function can be used for conversion between any in-out build combination, using
#' this argument, not just 36--<37.
#' @param ... additional arguments to make.granges(), so in other words, can use 'start' and
#' 'end' to specify ranges instead of 'pos'.
#' @return Returns positions converted from build 36 to 37 (or equivalent for alternative chain 
#' files). If using the 'ranges' parameter for position input, the object returned will be of
#' the same format. If using chr and pos to input, then the object returned will be a data.frame
#' with columns, chr and pos with rownames 'ids'. Output will be the same length as the input,
#' which is not necessarily the case for liftOver() which does the core part of this conversion.
#' Using vector or GRanges input will give a resulting data.frame or GRanges object respectively
#' that has the same order of rownames as the original input. Using RangedData will result in an
#' output that is sorted by genome order, regardless of the original order. If ranges has no
#' rownames, or if 'ids' is blank when using chr, pos, ids of the form rngXXXX will be generated
#' in order to preserve the original ordering of locations.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso conv.37.36, convTo37, convTo36
#' @references http://crossmap.sourceforge.net/
#' @examples
#' # various chain files downloadable from http://crossmap.sourceforge.net/ #
#' gene.labs <- c("CTLA4","IL2RA","HLA-C")
#' snp.ids <- c("rs3842724","rs9729550","rs1815606","rs114582555","rs1240708","rs6603785")
#' pp <- Pos(snp.ids); cc <- Chr(snp.ids)
#' conv.36.37(chr=cc,pos=pp,ids=snp.ids)
#' pp <- Pos(gene.labs)
#' gg <- GRanges(ranges=IRanges(start=pp$start,end=pp$end),seqnames=pp$chr)
#' conv.36.37(gg) # order of output is preserved
#' rr <- as(gg,"RangedData")
#' conv.36.37(rr) # note the result is same as GRanges, but in genome order
conv.36.37 <- function(ranges=NULL,chr=NULL,pos=NULL,...,ids=NULL,chain.file="/home/oliver/R/stuff/hg18ToHg19.over.chain") {
  require(GenomicRanges); require(rtracklayer); require(genoset)
  if(!file.exists(chain.file)) { stop("couldn't find chain file: ",chain.file) }
  chn <- import.chain(chain.file)
  #toranged <- F
  outType <- is(ranges)[1]
     
  if(!is.null(chr) & (!is.null(pos) | all(c("start","end") %in% names(list(...))))) {
    if(is.null(pos) & length(chr)==1) {
      if(length(list(...)$start)>1) {
        warning("when using start/end, 'chr' must have the same length as 'start'") 
      }
    }
    if(is.null(ids)) { ids <- paste0("rng",1:(max(length(chr),length(pos)))) } 
    ranges <- make.granges(chr=chr,pos=pos,row.names=ids,...)
    orn <- ids
  } else {
    if(is.null(rownames(ranges))) { rownames(ranges) <- paste0("rng",1:nrow(ranges)) }    
    orn <- rownames(ranges)
  }
  # return(ranges)
  if(is(ranges)[1]=="RangedData") { ranges <- as(ranges, "GRanges") }
  if(is(ranges)[1] %in% c("RangedData","GRanges")) {
    wd <- width(ranges)
    if(all(wd==1)) { SNPs <- TRUE } else { SNPs <- FALSE }
    mcols(ranges)[["XMYINDEXX"]] <- rownames(ranges)
    mcols(ranges)[["XMYCHRXX"]] <- ocr <- chr2(ranges)
    #prv(orn,ocr)
    opos <- start(ranges)
    ranges <- set.chr.to.char(ranges)
    #print(head(ranged))
    ranged.gr <- ranges # as(ranges,"GRanges"); #toranged <- T
  } else {
    stop("input specified resulted in an invalid GRanges/RangedData 'ranged' object, type ",is(ranges)[1]) 
  } 
  # change CHR-XY to CHR-X prior to liftOver, then change back #
  xy.ind <- grep("XY",seqnames(ranged.gr))
  if(length(xy.ind)>0) {
    found.xy <- TRUE
    xy.id <- rownames(ranged.gr)[xy.ind]
    if(!"chrX" %in% seqlevels(ranged.gr)) { seqlevels(ranged.gr) <- c(seqlevels(ranged.gr),"chrX") }
    seqnames(ranged.gr)[xy.ind] <- "chrX"
  } else { found.xy <- FALSE }
  ranged.gr.37 <- liftOver(ranged.gr,chn)
  myfun <- function(x) { 
    data.frame(start=minna(start(x)),end=maxna(end(x))) 
  }
  if(!SNPs ) {
    new.coords.df <- do.call("rbind",lapply(ranged.gr.37,myfun))
    ranged.gr.37<-ranged.gr
    ranges(ranged.gr.37)<-with(new.coords.df,IRanges(start=start,end=end))
    #seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
    seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
    out <- ranged.gr.37
  } else {
    seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
    #seqnames(ranged.gr.37)<-gsub("chr","",seqnames(ranged.gr.37))
    out <- as(ranged.gr.37,"IRangesList")
    #seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
    #new.coords.df <- as.data.frame(ranged.gr.37)
  }
  # seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
  # out <- as(ranged.gr.37,"IRangesList")
  out <- as(out,"RangedData")
  #return(ranged.gr.37)
  #ranged.gr.37 <- set.chr.to.numeric(ranged.gr.37)
  #if(!toranged | T) { return(ranged.gr.37) }
  ranged.gr.37 <- out #toGenomeOrder2(out)
  #return(ranged.gr.37)
  if(all(c("XMYINDEXX","XMYCHRXX") %in% colnames(ranged.gr.37))) {
    RN <- ranged.gr.37[["XMYINDEXX"]]
    nr <- nrow(ranged.gr.37)
    MAXDISPLAY <- 50
    if(length(orn)>length(RN)) { 
      cat("conversion failed for",length(orn[!orn %in% RN]),"rows, original positions kept:\n") ;  
      failz <- orn[!orn %in% RN]
      cat(comma(head(failz,MAXDISPLAY))) 
      if(length(failz)>MAXDISPLAY) { cat(", ... and",length(failz)-MAXDISPLAY,"more\n")  } else { cat("\n") }
      ln <- orn[!orn %in% RN]
      #return(ranges)
      newchr <- gsub("chr","",chr2(ranges[match(ln,ranges$XMYINDEXX),]))
      noopos <- start(ranges[match(ln,ranges$XMYINDEXX),])
      hcc <- hard.coded.conv()
      ifNAthen0 <- function(X) { X[is.na(X)] <- 0; return(X) }
      h36 <- which(noopos %in% hcc$pos36 & newchr==ifNAthen0(hcc$chr[match(noopos,hcc$pos36)])); l36 <- length(h36)
      h37 <- which(noopos %in% hcc$pos37 & newchr==ifNAthen0(hcc$chr[match(noopos,hcc$pos37)])); l37 <- length(h37)
      if(l36>=l37 & l36>0) {
        noopos[h36] <- hcc$pos37[match(noopos[h36],hcc$pos36)]
        cat("found",l36,"of the missing SNP hg18-hg19 lookups in an internal table\n")
      } else {
        if(l36<l37 & l37>0) {
          noopos[h37] <- hcc$pos36[match(noopos[h37],hcc$pos37)]
          cat("found",l37,"of the missing SNP hg19-hg18 lookups in an internal table\n")
        } else {
          ## no matches to extras table
        }
      }
      extra <- data.frame(Chr=newchr,Start=noopos,End=noopos)
      rownames(extra) <- ln
    } #else { cat("length is already the same\n") }
    Ind <- match(ranged.gr.37[["XMYINDEXX"]],orn)
    out <- data.frame(Chr=ranged.gr.37[["XMYCHRXX"]],Start=start(ranged.gr.37),End=end(ranged.gr.37),ind=Ind)
    rownames(out) <- RN
    if(length(orn)>length(RN)) {
      #prv(out,extra)
      out <- out[,-4] # 4 is the 'ind' column
      out <- rbind(out,extra)
      out <- out[orn,]
    } #else { cat("length is now the same\n") }
    #return(out) 
  } else { warning("missing key columns for chr, snp-name")  }
  #print(outType)
  #return(out)
  #prv(out)
  ranged.rd <- toGenomeOrder2(data.frame.to.ranged(out))
  #print(colnames(ranged.rd))
  ranged.gr.37 <- as(ranged.rd,"GRanges")
  #print(colnames(ranged.gr.37))
  if(found.xy) {
    xy.ind <- match(xy.id,rownames(ranged.gr.37))
    lmis <- length(which(is.na(xy.ind)))
    if(lmis>0) { warning('liftOver function removed ",lmis," chrX/chrY ranges'); xy.ind <- narm(xy.ind) }
    if(!"XY" %in% seqlevels(ranged.gr.37)) { seqlevels(ranged.gr.37) <- c(seqlevels(ranged.gr.37),"XY") }
    seqnames(ranged.gr.37)[xy.ind] <- "XY"
  }
  if(outType=="GRanges") { 
    #return(ranged.gr.37)
    cn37 <- colnames(mcols(ranged.gr.37))
    if("ind" %in% cn37) { 
      mind <- as.numeric(mcols(ranged.gr.37)[["ind"]])
      #prv(ranged.gr.37,mind)
      ranged.gr.37 <- ranged.gr.37[order(mind),]
      mcols(ranged.gr.37) <- mcols(ranged.gr.37)[,-which(cn37 %in% "ind")] 
    } # else { warning("couldn't find index column, GRanges object not sorted in original order") }
    if(all(rownames(ranged.gr.37) %in% orn)) { ranged.gr.37 <- ranged.gr.37[orn,] }
    return(ranged.gr.37)
  } else {
    if(outType=="RangedData") {
      #if("ind" %in% colnames(ranged.gr.37)) { ranged.gr.37 <- ranged.gr.37[,-which(colnames(ranged.gr.37) %in% "ind")] }
      return(toGenomeOrder2(as(ranged.gr.37,"RangedData")))
    } else {
      #prv(ranged.gr.37)
      out <- ranged.to.data.frame(ranged.gr.37,include.cols=FALSE,use.names=TRUE)
      cn37 <- colnames(mcols(ranged.gr.37))
      if("ind" %in% cn37) { 
        mind <- as.numeric(mcols(ranged.gr.37)[["ind"]])
        out <- out[order(mind),,drop=FALSE]
        if("ind" %in% colnames(out)) { out <- out[,-which(colnames(out) %in% "ind")] }
      }
      if(is.null(dim(out))) { dim(out) <- c(length(out)/3,3) }
      if(all(rownames(out) %in% orn)) { out <- out[orn,] }
      return(out)
    }
  }
}



#' Convert a snpStats SnpMatrix object to a dataframe
#' 
#' Converts a snpStats::SnpMatrix object to a dataframe where coding becomes 0,1,2,NA,
#' which represents genotypes as the number of copies of the reference allele.
#' @param SnpMat a snpStats::SnpMatrix object
#' @param matrix logical, whether to convert to a normal matrix (TRUE), or a data.frame (default, FALSE)
#' @export
#' @return data.frame object with genotype coding where 0,1,2 are the number of copies
#' of the reference allele (by default the letter latest in the alphabet), and NA is
#' missing data.
#' @seealso data.frame.to.SnpMatrix
#' @examples
#' samp <- matrix(sample(c(0:3),18,replace=T),
#'   ncol=3,dimnames=list(paste0("ID0",1:6),c("rs123","rs232","rs433")))
#' samp # note that 0's for SnpMatrix objects are missing values, and other
#' # values are the number of copies of the reference allele plus 1. The reference
#' # allele is by default the letter that is latest in the alphabet, e.g for C/G,
#' # G would be the reference, for T/A, T would be the reference
#' test.mat <- new("SnpMatrix", samp)
#' test.mat # preview does not show data, so use as.data.frame(test.mat)
#' is(as.data.frame(test.mat)[[1]]) # note that using 'as.data.frame' the type is 'raw'
#' SnpMatrix.to.data.frame(test.mat) 
#' # ^ now missing cells are NA, and genotypes are coded as number of copies of reference
SnpMatrix.to.data.frame <- function(SnpMat,matrix=FALSE) {
  #if(is(SnpMat)[1]=="snp.matrix") { SnpMat <- as(SnpMat,"SnpMatrix") }
  if(is(SnpMat)[1]!="SnpMatrix") { stop("SnpMat must be a SnpMatrix object") }
  cov.data <- as.data.frame(SnpMat)
  for(jj in 1:ncol(cov.data)) { 
    nuxt <- as.numeric(cov.data[,jj])-1
    nuxt[nuxt<0] <- NA
    cov.data[,jj] <- nuxt
    # assign(colnames(cov.data)[jj], nuxt)
  }
  if(matrix) { cov.data <- as.matrix(cov.data) }
  return(cov.data)
}


#' Convert a data.frame to a snpStats SnpMatrix object
#' 
#' Converts a dataframe to a snpStats::SnpMatrix object where the object contains
#' genotypes coded as number of copies of the reference allele: 0,1,2, and missing=NA.
#' This is an alternative to using new("SnpMatrix",data.frame()). Using 'new' the required 
#' format for the 'data.frame' argument is not as intuitive, as NA's are not allowed, and
#'  2 copies of the reference allele must be coded as 3, and 1 copy as 2, 0 copies as 1.
#' Note that this function will also accept data.frames/matrices coded in that way, and
#' will detect the coding automatically.
#' @param X a data.frame or matrix object containing allele codes 0,1,2 with missing=NA
#' @export
#' @seealso SnpMatrix.to.data.frame
#' @return SnpMatrix object
#' @examples
#' test.frame <- matrix(sample(c(0:2),18,replace=T),
#'   ncol=3,dimnames=list(paste0("ID0",1:6),c("rs123","rs232","rs433")))
#' test.frame[2,2] <- NA # set one genotype to missing
#' test.frame 
#' test.mat <- data.frame.to.SnpMatrix(test.frame) 
#' snp.mat <- new("SnpMatrix",test.frame) # note that does not handle the NAs/0s nicely
#' all.equal(print(as.data.frame(snp.mat)),print(as.data.frame(test.mat))) # shows offset by 1
#' test.mat
#' row.summary(test.mat) # call rate analysis by sample using snpStats function
#' col.summary(test.mat) # analysis by SNP
data.frame.to.SnpMatrix <- function(X){
  if(is.data.frame(X)) {
    if(any(sapply(lapply(X,is),"[",1) %in% c("character","factor"))) {
      for(cc in 1:ncol(X)) {
        X[[cc]] <- as.numeric(X[[cc]])
      }
    }
  } else { 
    if(!is.matrix(X)) { warning("X should be a matrix or data.frame, ",
                                "conversion is likely to fail if X is not sufficiently matrix-like")}
  }
  mxx <- maxna(X)
  if(mxx>3) { warning("Dataframe does not appear to contain allele codes") }
  X <- round(X)
  if(mxx==3) { X <- X-1 ; X[X<0] <- NA }
  NN <- as.matrix(X)
  #NN <- round(NN)
  SS <- as(NN,"SnpMatrix")
  return(SS)
}

#' Determine major or minor allele status for a set of SNPs
#' 
#' For a snpStats object or data.frame containing values 0,1,2, NA representing genotypes AA, AB, 
#' BB and no-call. Determines whether the reference allele is the major or minor allele 
#' Where the homozygous genotype coded as highest value = reference, e.g, if AA=0, AB=1, BB=2, 
#' then B is considered the reference here, and by the snpStats package. Combines this with
#' frequencies of the alleles to evaluate whether 'BB' is major or minor. Note that default
#' behaviour for a SnpMatrix is to code alleles alphabetically, so usually the reference allele
#' is the letter later in the alphabet, e.g, it is never an 'A' allele.
#' @param X a SnpMatrix object (see snpStats), or a data.frame coded by reference allele copies,
#' 0,1,2, with missing as NA
#' @param checks logical, whether to perform additional checks for valid values during conversion,
#' setting FALSE will give a slight increase in speed, but is not recommended.
#' @param tag.mono logical, whether to append a prefix of 'mono' to the major/minor factor code
#' for monomorphic SNPs (minor allele frequency = zero)
#' @param tag.neutral logical, whether to call SNPs with exactly equal A and B allele frequencies
#' 'neutral'; if this option is false, the default is to call neutral SNPs 'minor'.
#' @return returns a factor vector of the same length as the number of SNPs (columns) in the
#' SnpMatrix object/data.frame, indicating for each SNP whether the reference 'B' allele is the 
#' major or minor allele. If the SnpMatrix was creating using a snpStats import function the 
#' reference allele should be the nucleotide letter that is latest in the alphabet. So, for 
#' instance if a SNP is either T/A or A/T, then the reference (B) allele will be 'T'. This 
#' function indicates whether the 'B' allele is the major or minor allele (the major  allele
#' has the greatest frequency). This function can also code 'neutral' if both alleles have 
#' equal frequency when 'tag.neutral' is TRUE, and can add the prefix 'mono.' when 
#' 'tag.mono' is TRUE and one allele has 100% frequency, i.e, is monomorphic.
#' @seealso caseway
#' @export
#' @examples
#' cn <- c("rs123","rs232","rs433","rs234")
#' rn <- paste0("ID0",1:6)
#' samp <- rbind(c(0,2,0,2),c(0,0,0,2),c(0,1,1,2),c(0,0,2,NA),c(0,2,1,1),c(0,2,2,0))
#' dimnames(samp) <- list(rn,cn)
#' test.mat <- data.frame.to.SnpMatrix(samp)
#' majmin(test.mat)
#' col.summary(test.mat) # show call rates, MAF, allele frequencies
#' samp
#' majmin(samp) # also works on a matrix/data.frame with same structure as a Snp.Matrix
#' majmin(test.mat,tag.neutral=TRUE) # rs433 has equal A & B frequencies, and can be tagged neutral
#' majmin(test.mat,tag.mono=TRUE) # rs123 has zero B allele frequency, and can be tagged 'mono'
#' samp[2,2] <- 99 # insert invalid value into matrix
#' majmin(samp,TRUE,TRUE,TRUE) # warning for invalid value
#' majmin(samp,checks=FALSE) # invalid value is converted to NA without warning
majmin <- function(X,checks=TRUE,tag.mono=FALSE,tag.neutral=FALSE) {
  ## workhorse internal function ##
  #print(is(X)[1])
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
    if(a0==a2) { if(tag.neutral) { type <- "neutral"  } else { type <- "minor" } }
    if( length(which(c(A0,A1,A2)==0))==2 ) { if(tag.mono) { type <- paste("mono",type,sep=".") } }
    return(type)
  }
  ## main code ##
  if(length(Dim(X))!=2) { 
    if(length(Dim(X))==1) { return(do.mm(as.numeric(X))) } else {
      warning("invalid object for major/minor allele testing"); return(NA)
    }
  }
  snpmat <- snp.mat <- F
  if(is(X)[1] %in% "SnpMatrix") { snpmat <- snp.mat <- T} else {
    tt1 <- table(round(as.numeric(X[,1])));  if("3" %in% names(tt1)) { snpmat <- T }
  }
  if(snp.mat) {  
    ii <-  col.summary(X)
    #print(ii)
    all.typz <- c("minor","major")[as.numeric(round(ii$RAF,3)!=round(ii$MAF,3))+1]
    if(tag.neutral) { all.typz[round(ii$P.AA,6)==round(ii$P.BB,6)] <- "neutral" }
    if(tag.mono) { all.typz[round(ii$MAF,6)==0] <- paste("mono",all.typz[round(ii$MAF,6)==0],sep=".") }
  } else {
    all.typz <- apply(X,2,do.mm,snpmat=snpmat)
  }
  fl <- c("major","minor")
  if(tag.mono) { fl <- c(fl,"mono.major","mono.minor") }
  if(tag.neutral) { fl <- c(fl,"neutral") }
  return(factor(all.typz,levels=fl))
}


#' Find the direction of GWAS effects between cases and controls
#' 
#' After conducting an case-control association analysis on a SnpMatrix, e.g, 
#' GWAS using for example snp.rhs.tests from the snpStats package, it is not
#' always trivial to determine the direction of the effects with respect
#' to the reference allele. This function calculates which way around the
#' data is for case vs control (pheno) and will indicate with respect to 
#' cases whether they have more copies of the reference allele, or less, or also 
#' can highlight whether the heterozygous is the affected genotype. This function
#' works on a SnpMatrix or dataframe coded 0,1,2,NA (autodetects which). Note that
#' using a SnpMatrix, with het.effects=FALSE can be much faster (5-100x) than 
#' using a data.frame and/or setting het.effects=TRUE.
#' @param X a SnpMatrix or data.frame containing genotypes for a set of samples
#' @param pheno integer, must be coded as controls=0 and cases=1, or alternatively
#' using controls=1, cases=2 will be automatically detected and recoded to 0,1.
#' @param checks logical, whether to perform additional checks for valid values 
#' during conversion, setting FALSE will give a slight increase in speed, but is 
#' not recommended.
#' @param long logical, whether to use longer more explicit text in the resulting
#' categories describing effect directions, or to use abbreviated codes. For instance,
#' if long==TRUE, then SNPs where cases have more of the reference allele will produce
#' a resulting factor of "cases have more 2, less 0", whereas if long==FALSE, the result
#' would be "CasesRef+".
#' @param het.effects logical, whether to allow for the possibility that the risk
#' effect is seen for the heterozygous genotype rather than either homozygous allele,
#' if het.effects is TRUE, this will add to two additional possible categories to
#' the output vector, for and increase or decrease in the het allele for cases. If
#' het.effects is FALSE, then the borderline result will default to one side,
#' depending on how the allele is coded in the SnpMatrix object. Note that such
#' instances are unlikely to matter as the most commonly used 1df tests should 
#' be non-significant. Whereas if a 2df test is used, het.effects should be set TRUE.
#' @export
#' @return a factor with the same length as the number of SNPs in X. If het.effects
#' is FALSE and long is FALSE, this will take values 'CasesRef+' indicating cases
#' had a higher frequency of the reference allele than controls, or 'CasesRef-', 
#' indicating that cases had a lower frequency of the reference allele than controls.
#' If 'het.effects'=TRUE, then two more categories are added, 'CasesHet+' and 'CasesHet-',
#' indicating that cases had more of the heterozygous allele, or less, respectively, 
#' compared to controls. Or if 'long' is TRUE, then these four categories are respectively
#' changed to: 'cases have more 2, less 0', 'cases have more 0, less 2', and for het,
#' 'cases have more 1, less 0,2', 'cases have less 1, more 0,2'.
#' Note that these categories are based on absolute counts, and do not necessarily
#' reflect statistical differences in frequency, these should just be used to interpret
#' the direction of your findings, not to perform the analysis.
#' @seealso majmin
#' @examples
#' cn <- c("rs123","rs232","rs433","rs234","rs222")
#' rn <- paste0("ID0",1:10)
#' samp <- rbind(c(0,2,0,2,0),c(0,0,0,2,0),c(0,1,1,2,0),c(0,0,2,NA,1),c(0,2,1,1,0),
#'               c(0,2,2,0,2),c(0,2,0,2,2),c(0,1,NA,2,1),c(0,1,0,0,2),c(0,2,2,2,2))
#' dimnames(samp) <- list(rn,cn)
#' test.mat <- data.frame.to.SnpMatrix(samp)
#' phenotype <- c(rep(0,5),rep(1,5))
#' cbind(samp,phenotype) # show example data as a matrix with phenotype
#' caseway(test.mat,phenotype,long=TRUE) # long version
#' caseway(samp,phenotype,het.effects=TRUE) # also works using a data.frame
#' # example conducting association analysis, then adding interpretation of results
#' result <- data.frame(p.value=p.value(snp.rhs.tests(formula=phenotype~1,snp.data=test.mat)))
#' result[["direction"]] <- caseway(test.mat,phenotype)
#' result[["referenceIs"]] <- majmin(test.mat)
#' result # note that only rs222 is significant
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
      warning("invalid object for case/control effect direction testing"); return(NA)
    }
  }
  snpmat <- F
  if(is(X)[1] %in% "SnpMatrix") { 
    snpmat <- T
    if(!het.effects) {
      SSTS <- snpStats::single.snp.tests(pheno, snp.data=X, score=T)
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





##' log sum function @author Claudia Giambartolomei - internal
logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}



## functions from Chris W - adapted

# internal
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



##' Calculate approximate Bayes factors from p values and MAF
##'
##' This is a function to calculate approximate Bayes factors from p
##' values and MAF - for reference see Wakefield, J (2009) Bayes
##' factors for genome-wide association studies: comparison with
##' p-values. 
##' @param p p value
##' @param maf minor allele frequency
##' @param n0 number of controls
##' @param n1 number of cases
##' @param scale0 by default, =n0
##' @param scale1 by default, =n1
##' @references Genetic Epidemiology 33: 79–86.
##' @return the appoximate bayes factor calculated based on
##' the p-value, minor allele frequency and case/control sample
##' sizes
##' @export
##' @author Chris Wallace
##' @examples
##' abf(.0001,.2,9500,6670)
abf <- function(p,maf, n0=1000, n1=1000, scale0=n0, scale1=n1) { 
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
  abf <- 2 * log(sqrt(VW/V) * exp( - z^2 * W / (2 * VW) ))
  return(abf)
}




# internal
# snpStats imputation only works if there are correlated SNPs with non-missing values
# that can be used to interpolate missing SNPs. If any correlated SNPs are missing
# 'impute.missing' will leave these blank. This function mops up the remainder
# by randomly inserting values consistent with the minor allele frequency of each SNP
# (Nick)
randomize.missing2 <- function(X,verbose=FALSE,n.cores=1) {
  miss1 <- function(x) { 
    x <- as.numeric(x)
    TX <- table(c(round(x),0,1,2,3))-c(1,1,1,1) # to force zero counts to be in the table
    naz <- which(x==0)
    if(length(naz)>0 & length(TX)>0) {
      x[naz] <- sample(as.numeric(names(TX)),size=length(naz),
                       replace=T,prob=as.numeric(TX))
    }
    return(as.raw(x))
  }
  miss2 <- function(sel,aa,ab,bb) { 
    x <- X@.Data[,sel]
    naz <- which(x==0)
    p.vec <- c(aa,ab,bb)
    if(any(is.na(p.vec))) {
      if(all(is.na(p.vec))) {
        x[naz] <- as.raw(rsnp(length(naz),cr=1))
        return(x)
      } else {
        stop(p.vec)
        #Pleft <- 1-sum(p.vec,na.rm=T)  
      }
    }
    x[naz] <- sample(as.raw(1:3),size=length(naz),replace=T,prob=p.vec)
    return(x)
  }
  # randomly generate replacements for missing values using current distribution for each column of X
  typ <- is(X)[1]
  if(!typ %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { stop("X must be a SnpMatrix object or similar") }
  cs <- col.summary(X)
  FF <- nrow(X)
  nmiss <- FF - cs$Calls
  select <- nmiss>0
  if(length(which(select))<1) { return(X) }
  if(verbose) { cat(sum(nmiss),"missing values replaced with random alleles\n") }
  if(length(which(select))==1) { X[,select] <- miss2(select,cs$P.AA[select],cs$P.AB[select],cs$P.BB[select]); return(X) }
  if(n.cores>1) {
    do.chunk <- function(SEL) {
      ii <- mapply(miss2,SEL,cs$P.AA[SEL],cs$P.AB[SEL],cs$P.BB[SEL])
      return(ii)
    }
    posz <- which(select)
    Ls <- length(posz)
    split.to <- min(c(round(Ls/100),n.cores),na.rm=T)
    #if(n.cores>4) { split.to <- split.to * 4 } # divide more if using multicores
    stepz <- round(seq(from=1,to=Ls+1,length.out=round((split.to+1))))
    if((tail(stepz,1)) != Ls+1) { stepz <- c(stepz,Ls+1) }
    split.to <- length(stepz)-1
    outlist <- sel.range <- vector("list",split.to)
    for (cc in 1:split.to) {
      c1 <- stepz[cc]; c2 <- stepz[cc+1]-1  # check this in FN!
      sel.range[[cc]] <- posz[c(c1:c2)]
    }
    #prv(split.to,stepz,select,Ls,sel.range)
    outlist <- mclapply(sel.range,do.chunk)
    for (cc in 1:split.to) {
      X@.Data[,sel.range[[cc]]] <- outlist[[cc]]
    }
  } else {
    ii <- mapply(miss2,which(select),cs$P.AA[select],cs$P.AB[select],cs$P.BB[select])
    X@.Data[,select] <- ii
  }
  return(X)
}


# internal
# snpStats imputation only works if there are correlated SNPs with non-missing values
# that can be used to interpolate missing SNPs. If any correlated SNPs are missing
# 'impute.missing' will leave these blank. This function mops up the remainder
# by randomly inserting values consistent with the minor allele frequency of each SNP
# (Nick)
randomize.missing <- function(X,verbose=FALSE) {
  miss1 <- function(x) { 
    TX <- table(c(round(x),0,1,2))-c(1,1,1) # to force zero counts to be in the table
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
  if(verbose) { cat(sum(nmiss),"missing values replaced with random alleles\n") }
  if(length(which(select))==1) { X[,select] <- miss1(X[,select]); return(X) }
  X[,select] <- apply(X[,select],2,miss1)
  return(X)
}


#' Replace missing values in a SnpMatrix object with imputed values
#'
#' This function is a wrapper for the snpStats snp.imputation() and impute.snps() functions.
#' It allows a full imputation with one simple command, and facilitates stratified imputation for
#' subsets of samples, and the parameter 'by' allows subsetting by SNPs, which speeds up the
#' imputation for large sets by running it in smaller batches (as large ones are really slow).
#' The standard use of the snpStats imputation functions will still leave some NA's behind,
#' whereas the option 'random' ensures each missing value is replaced, even if just by an
#' allele chosen at random (using existing frequencies).
# (chris)
#' @param X SnpMatrix object with missing data that you wish to impute
#' @param strata factor, imputation can be done independently for different sample subsets,
#' use a vector here of the same length as the number of samples (nrow(X)) to code this
#' @param by integer, a parameter that is passed to impute.missing() to determine
#' what sized subsets to impute by (smaller subsets makes imputation faster)
#' @param random logical, when imputation is performed using the core function, impute.snps(),
#' usually there are some SNPs that remain missing, due to missingness in the most correlated
#' SNPs used to impute them. Setting this parameter TRUE, will replace this remainder with
#' random values (using the allele frequencies for the non missing data). Setting to FALSE
#' will leave these values as missing in which case you cannot rely on this function returning
#' a fully complete matrix.
#' @param data.frame logical, if TRUE then return the result as a data.frame instead of a
#' SnpMatrix object
#' @param verbose logical, if TRUE, more information about what is being done in the 
#' imputation will be provided. Additionally, if TRUE and the dataset being imputed is 
#' large enough to take a long time, then a progress bar will be displayed using
#' 'loop.tracker()'. If FALSE, no progress bar will be used, regardless of the size of dataset.
#' @param ... further arguments to 'impute.snps()' from snpStats, which is the core function
#' that performs the main imputation.
#' @return Returns a SnpMatrix with no missing values (as any originally present are imputed
#' or randomized), or if random=FALSE, it may contain some missing values. If data.frame=TRUE
#' then the returned object will be a data.frame rather than a SnpMatrix.
#' @export
#' @author Chris Wallace and Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' test.mat <- rSnpMatrix(nsnp=5,nsamp=10)
#' print(SnpMatrix.to.data.frame(test.mat))
#' out.mat <- impute.missing(test.mat)
#' print(SnpMatrix.to.data.frame(out.mat))
#' test.mat2 <- rSnpMatrix(nsnp=200,nsamp=100,call.rate=.98)
#' head(col.summary(test.mat2))
#' out.mat2 <- impute.missing(test.mat2,by=50)
#' cs <- col.summary(out.mat2)
#' head(cs)
#' # 'random' is true, so there should be no remaining missing values
#' paste(out.of(length(which(cs$Call.rate==1)),nrow(cs)),"SNPs imputed to 0% missing\n")
#' out.mat3 <- impute.missing(test.mat2,random=FALSE)
#' cs <- col.summary(out.mat3)
#' # random was FALSE, so some of these should still have missing values:
#' paste(out.of(length(which(cs$Call.rate==1)),nrow(cs)),"SNPs imputed to 0% missing\n")
#' # for instance, these
#' print(head(cs[which(cs$Call.rate<1),]))
#' # complicated example using both sample and SNP subsets for imputation
#' out.mat4 <- impute.missing(test.mat2, by=100, strata = c(rep(1,25),rep(2,50),rep(3,25)), random=FALSE)
impute.missing <- function (X,  strata = NULL, by=NULL, random=TRUE, data.frame = FALSE,
                            verbose=TRUE,...) {
  if(!is(X)[1]=="SnpMatrix") { stop("X must be a SnpMatrix object") }
  if(Dim(X)[1] < 2 | Dim(X)[2] <2) { stop("X must be at least 2 x 2") }
  N <- as(X, "numeric"); 
  done <- FALSE
  mem.max <- 0.02  #GB
  bp <- 1:ncol(X) # parameter needed for snp.imputation() ?
  if(!is.null(by)) {
    if(!is.numeric(by)) { warning("by must be numeric, ignoring invalid parameter"); by <- NULL }
    by <- round(by) # force integer
    if(by<3) { warning("by is too small to be useful, ignoring parameter"); by <- NULL }
    if((by*1.25)>ncol(X)) { by <- NULL } # ignore 'by' if num SNPs is not that large
  }
  if(any(Dim(X)!=Dim(N))) { stop("'as numeric' caused loss of dimensionality for unknown reason") }
  if (!is.null(strata)) {
    if(length(strata)==nrow(X)) {
      ## Partition by sample subsets, then re-send the subsets to this function for processing
      strata <- as.factor(strata)
      if (length(levels(strata)) > 20) 
        stop("too many levels in strata\n")
      for (i in levels(strata)) {
        if(verbose) { cat("\nImputing sample subset '", i, "' (strata)\n",sep="") }
        wh <- which(strata == i)
        N <- as.data.frame(N)
        N[wh, ] <- impute.missing(X[wh, , drop = FALSE], 
                                   data.frame = TRUE, by=by,...)
      }
      done <- TRUE
    } else {
      #print(length(strata));print(nrow(X))
      warning("'strata' should have length==nrow(X), all samples will be imputed together")
    }
  }
  if (is.null(strata) & is.numeric(by)) {
    ## Partition by SNP subsets, then re-send the subsets to this function for processing
    strata <- as.factor(make.split(ncol(X),by=by))
    # added by nick
    for (i in levels(strata)) {
      wh <- which(strata == i)
      if(verbose) { cat(" imputing SNP subset ", minna(wh), "-",maxna(wh),"\n",sep="") }
      anything <- impute.missing(X[,wh , drop = FALSE], data.frame = TRUE, ...)
      N <- as.data.frame(N)
      N[, wh] <- anything
    }
    done <- TRUE
  } else {
    if(!done) {
      # The core imputation, done using the snpStats functions 'snp.imputation' and 'impute.snps' #
      csumm <- col.summary(X)
      use <- csumm[, "Certain.calls"] == 1
      X2 <- X[, use]
      bp <- bp[use]
      imp <- (csumm[, "Call.rate"] < 1 & !is.na(csumm[, "Call.rate"]))[use]
      #if(sumna(imp)==175) { stst <- strata; byby <- by; prv(stst,byby)}
      if(verbose) { cat(sumna(imp), " SNP",if(sumna(imp)!=1) {"s"} else { "" }," in the set had missing values\n",sep="")  }
      ll <- length(which(imp)); dd <- 1 # mx <- max(which(imp),na.rm=T); 
      if(estimate.memory(X2)<mem.max) { verbose <- FALSE }
      for (i in which(imp)) {
        if(verbose) {
          loop.tracker(dd,ll); dd <- dd + 1 # track this loop as count, these pars not used elsewhere
        }
        supres <- capture.output(rule <- snp.imputation(X2[, -i, drop = FALSE], X2[, 
                                                   i, drop = FALSE], bp[-i], bp[i]))
        if (is.null(rule@.Data[[1]])) 
          next
        imp <- impute.snps(rules = rule, snps = X2[, rule@.Data[[1]]$snps, drop = FALSE], ...) 
        wh.na <- which(is.na(N[, i]))
        #gtWWW <- X2
        #prv(gtWWW)
        if(length(wh.na)>0) { 
          N[wh.na, colnames(X2)[i]] <- imp[wh.na]
        }
      }
    }
  }
  if(random) { N <- randomize.missing(N,verbose=verbose) } # replace any remaining missing values
  if(data.frame) {
    return(as.data.frame(N))
  }
  else {
    #print(Dim(N))
    #print(Dim(X@snps))
    N <- round(N)
    N[is.na(N)] <- -1 # replace missing with -1 [will become zero in next line]
    if(is(N)[1]=="list" | is(N)[1]=="data.frame") { N <- as.matrix(N) } # sometimes it becomes a list for some reason???
    return(new("SnpMatrix", data = (as.raw(N+1)),
                      nrow = nrow(N), ncol = ncol(N), dimnames = dimnames(N)))
  }
}




# internal
#' Remove trailing letter from non-unique rs-ids
#'
#' @examples
#' snp.ids <- rsnpid(25)
#' snp.ids[1:2] <- paste0(snp.ids[1:2],"b")
#' snp.ids[19:20] <- paste0(snp.ids[19:20],"c")
#' snp.ids[6:7] <- paste0(snp.ids[6:7],"d")
#' snp.ids[11:12] <- paste0(snp.ids[11:12],"a")
#' snp.ids
#' rmv.trail(snp.ids)
rmv.trail <- function(rs.ids,suffix=c("b","c","d","a")) {
  if(!is.character(suffix)) { stop("suffix must a character vector") }
  ind <- NULL
  for (cc in 1:length(suffix)) {
    ind <- c(ind,grep(suffix[cc],rs.ids))
  }
  ind <- unique(ind)
  X <- rs.ids[ind]
  nX <- nchar(X)
  last.chars <- substr(X,nX,nX)
  sufz <- (last.chars %in% c("a","b","c","d"))
  X[sufz] <- substr(X[sufz],1,nX[sufz]-1)
  rs.ids[ind] <- X
  return(rs.ids)
}

# internal
#' Add trailing letter(s) to non-unique rs-ids
#' @examples
#' snp.ids <- rsnpid(15)
#' snp.ids
#' add.trail(snp.ids)
#' snp.ids <- snp.ids[sample(15,30,replace=TRUE)]
#' snp.ids
#' add.trail(snp.ids)
add.trail <- function(rs.ids,suffix=c("b","c","d","a")) {
  rs.ids <- rmv.trail(rs.ids)
  for (txt in suffix) {
    dupz <- duplicated(rs.ids)
    if(length(which(dupz))>0) {
      rs.ids[dupz] <- paste0(rmv.trail(rs.ids[dupz]),txt)
    }
  }
  dupz <- duplicated(rs.ids)
  if(length(which(dupz))>0) { warning("more than ",length(suffix),
                                      " duplications of at least 1 individial rs-id, suffixes exhausted, duplicates remain")
  }
  return(rs.ids)
}


# internal, function generate a random MAF, then random SNP
rsnp <- function(n,A.freq.fun=runif,cr=.95, A.freq=NA) { 
  if(is.na(A.freq)) {  m <- A.freq.fun(1) } else { m <- A.freq }
  x <- sample(0:3, n, replace=T, prob=c(1-cr,cr*(m^2),cr*(2*(1-m)*m),cr*((1-m)^2)))
  return(x) 
}
 

# internal
rsnpid <- function(n) { 
  id.len <- sample(c(3:8),n,replace=T,prob=c(0.01, 0.01, 0.01, 0.10, 0.50, 0.37))
  each.id <- function(l) { sapply(l,function(n) { paste(replicate(n,sample(1:9,1)),collapse="",sep="") }) }
  sufz <- each.id(id.len)
  ids <- paste0("rs",sufz)
  return(ids)
}

# internal
ldfun <- function(n) { x <- runif(n); r2 <- runif(n); x[x<.6 & r2>.1] <- x[x<.6 & r2>.1]+.4 ; return(x) }

# internal
rsampid <- function(n,pref="ID0") { paste0(pref,pad.left(1:n,"0")) }

# internal
snpify.cont <- function(x,call.rate=.95,A.freq.fun=runif) { 
  if(length(x)<2) { stop("x must be longer than 1 element") }
  if(length(x)<5) { warning("for good simulation x should be fairly large, e.g, at least 10, better >100") }
  n <- length(x)
  if(all(x %in% 0:3)) { return(x) } # these are already snp-coded
  rr <- rank(x)
  sim1 <- rsnp(n,A.freq.fun=A.freq.fun,cr=call.rate) # A.freq=NA
  x[rr] <- sort(sim1)
  return(x)
}

# internal
get.biggest <- function(r2.mat) {
  mm <- max(r2.mat,na.rm=T)
  coord <- which(r2.mat==mm,arr.ind=T)
  if(!is.null(dim(coord))){ coord <- coord[1,] }
  return(coord)
}

# internal
get.top.n <- function(mat,n=10) {
  cr <- cor(mat,use="pairwise.complete")^2
  diag(cr) <- 0
  nn <- NULL
  while(length(nn)<n) { 
    coord <- get.biggest(r2.mat=cr)
    nn <- unique(c(nn,coord))
    cr[coord[1],coord[2]] <- NA
    cr[coord[2],coord[1]] <- NA
    #cr[,coord[1]] <- NA; cr[coord[1],] <- NA
    #cr[,coord[2]] <- NA; cr[coord[2],] <- NA
  }
  nn <- nn[1:n]
  new.mat <- mat[,nn]
  return(new.mat)
}


#' Create a SNP matrix with simulated data
#' 
#' A simple function to simulation random SnpMatrix objects for
#' testing purposes. Does not produce data with an 'LD' structure
#' so is not very realistic. 
#' @param nsnp number of SNPs (columns) to simulate
#' @param nsamp number of samples (rows) to simulate
#' @param call.rate numeric, percentage of SNPs called, e.g, .95 leaves 5% missing data,
#' or 1.0 will produce a matrix with no missing data
#' @param A.freq.fun function, a function to randomly generate frequencies of the 'A'
#' allele, this must produce 'n' results 0 < x < 1 with argument n, the default is
#' for uniform generation.
#' @export
#' @return returns a SnpMatrix object with nsnp columns and nsamp rows, with 1-call.rate% 
#' missing data, with allele frequencies generated by 'A.freq.fun'.
#' @examples
#' newMat <- rSnpMatrix(call.rate=.92) # specify target of 8% missing data
#' print(as.data.frame(newMat)) # to see actual values
rSnpMatrix <- function(nsamp=10, nsnp=5, call.rate=.95, A.freq.fun=runif) {
  dummy.vec <- rep(nsamp,nsnp)
  call.rate <- force.percentage(call.rate)
  exp.fac <- 10
  ld <- FALSE
  warn.mem <- 0.5 # 0.5GB
  if(!ld) {
    mat <- sapply(dummy.vec,rsnp,A.freq.fun=A.freq.fun,cr=call.rate)
  } else {
    # this didn't really work #
    if(estimate.memory(exp.fac*nsamp*nsnp)>warn.mem) { warning("when ld=T, large simulations can take a long time") }
    mat <- sapply(rep(nsamp,nsnp*exp.fac),rsnp,A.freq.fun=A.freq.fun,cr=call.rate)
    mat <- get.top.n(mat,nsnp)    
  }
  rownames(mat) <- rsampid(nsamp); colnames(mat) <- rsnpid(nsnp)
  #sm <- new("SnpMatrix",as.raw(mat),dimnames=list(rsampid(nsamp),rsnpid(nsnp)))
  sm <- new("SnpMatrix", data = as.raw(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  return(sm)
}





#' Extend an interval or SNP by distance in centimorgans (recombination distance)
#' 
#' It is straightforward to extend a genomic interval or position by a number of basepairs, or
#' a percentage, but extending by recombination units of centimorgans is more involved, requiring
#' annotation lookup. This function streamlines this process.
#' This function makes use of recombination rate hapmap reference files to calculate 
#' recombination distances for genome locations, in centimorgans. For a given position
#' (or vector), a window can be returned of a given extension on either side of the position,
#' for instance, 1 centimorgan to the left, and to the right of a SNP, giving a 2 centimorgan
#' range as a result.
#' @param ranges optional GRanges or RangedData object describing positions for which we want to
#' generate windows, removing the need to enter chr, start and end
#' @param chr character, an optional vector of chromosomes to combine with 'start' and 'end'
#'  to describe positions for which to generate recombination windows
#' @param start integer, an vector of start points for chromosome ranges
#' @param end integer, an vector of end points for chromosome ranges
#' @param window numeric, number of centimorgans to extend the window either side of
#' the range or location (can be a fraction)
#' @param bp.ext numeric, optional number of base-pairs to extend the window by in addition
#' to the centimorgan extension
#' @param rec.map recombination map object (list of 22 data.frames) generated using 
#' 'get.recombination.map()'; if you are performing many of these operations, loading this 
#' object into your workspace and passing it on to this function will save loading it each 
#' time, and provide a speed advantage. Only use an object generated by get.recombination.map(),
#'  as otherwise the results will almost certainly be meaningless.
#' @param info logical, whether to display the derived window size and number of hapmap SNPs within
#' the window for each window derived
#' @seealso get.recombination.map, get.nearby.snp.lists, exp.window.nsnp
#' @author Chris Wallace and Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @export
#' @examples
#' # not run, as initial download of the recombination map takes nearly a minute #
#' recwindow(chr=11,start=10000000,end=10000000,window=1,bp.ext=10000)
#' rd <- RangedData(ranges=IRanges(start=c(1.5,10.1)*10^7, end=c(1.55,10.1)*10^7),space=c(2,10))
#' rd # show original data
#' recwindow(rd) # now extended by the interval
#' recwindow(as(rd,"GRanges"),info=FALSE) # also works for GRanges
recwindow <- function(ranges=NULL,chr=NA,start=NA,end=start,window=0.1,bp.ext=0, rec.map=NULL, info=TRUE) {
  #www <- window ; ccc <- chr; sss <- start; prv(ccc,sss,www,bp.ext)
  if(!is.numeric(bp.ext)) { warning("bp.ext must be numeric, setting to zero"); bp.ext <- 0 }
  if(!is.numeric(window)) { warning("window must be numeric, setting to 0.1 centimorgans"); window <- 0.1 }
  if(!all(is.na(chr))) { if(any(!paste(chr) %in% paste(1:22))) { 
    stop("this function only works for autosomes 1-22 [e.g, no X,Y or formatting like 'chr2', etc]") } }
  chr <- as.numeric(chr)
  typ <- is(ranges)[1]
  if(typ %in% c("RangedData","GRanges")) { 
    if(typ=="GRanges") { ranges <- as(ranges,"RangedData") }
    ranges <- toGenomeOrder2(ranges,strict=T)
    ss <- start(ranges); ee <- end(ranges); cc <- chr2(ranges)
    out <- recwindow(chr=cc,start=ss,end=ee,window=window,bp.ext=bp.ext,info=info,rec.map=rec.map)
    if(length(out)==2) { out <- as.matrix(out); dim(out) <- c(1,2) } 
    outData <- RangedData(ranges=IRanges(start=out[,1],end=out[,2],names=rownames(ranges)),space=cc)
    outData <- toGenomeOrder2(outData,strict=TRUE)
    if(ncol(ranges)>0) {
      for (zz in 1:ncol(ranges)) { 
        if(typ=="GRanges") {
          outData[[colnames(ranges)[zz]]] <- mcols(ranges)[,colnames(ranges)[zz]]
        } else {
          outData[[colnames(ranges)[zz]]] <- ranges[[colnames(ranges)[zz]]]
        }
      }
    }
    if(is(ranges)[1]=="GRanges") { outData <- as(toGenomeOrder2(outData,strict=T),"GRanges") }
    return(outData)
  } else {
    if(all(!is.na(chr)) & all(!is.na(start)) & all(!is.na(end))) {
      if(length(chr)==length(start) & length(start)==length(end)) {
        if(length(chr)>1) {
          # run for a vector
          out <- matrix(ncol=2,nrow=length(chr)); colnames(out) <- c("start","end")
          for (dd in 1:length(chr)) {
            out[dd,] <- recwindow(chr=chr[dd],start=start[dd],end=end[dd],
                                  window=window,bp.ext=bp.ext,info=info,rec.map=rec.map)
          }
          return(out)
        } else {
          ## continue as normal, just a single coordinate/range to process
        }
      } else {
        stop("invalid input, start, end and chr need to be the same length")
      }
    } else {
      stop("invalid input, either use a RangedData object, or else chr, start and end")
    }
  }
  
  #rate.fn <- sprintf("/dunwich/scratch/chrisw/HapMap/rates_rel22/genetic_map_chr%s_b36.txt.gz",chr)
  #print(rate.fn)
  #rates <- read.table(gzfile(rate.fn),header=TRUE)
  if(is.list(rec.map)) { if(length(rec.map)==22) { rates <- rec.map[[chr]] } } else {
    if(is.null(rec.map)) { rates <- get.recombination.map()[[chr]] } else {
      if(is.character(rec.map)) { rates <- get.recombination.map(dir=rec.map)[[chr]] } else {
        stop("invalid value for rec.map entered")
      }
    }
  }
  cm.st <- rates[which.min(abs(rates$position-start)),3]
  cm.en <- rates[which.min(abs(rates$position-end)),3]
  
  mx <- max(window,1)
  kk <- rates[which.min(abs(rates[,3]-(cm.st-window))) : which.min(abs(rates[,3]-(cm.en+window))),]
  if(info) { cat("n hapmap snps in window =",nrow(kk),"\n") }
  from <- min(kk[,1])
  to <- max(kk[,1])
  ##
  if(info) {
    cat("new window size is\nleft: ",(start-from+bp.ext)/1000,"kb\tright: ",
      (to-end+bp.ext)/1000,"kb\ttotal: ",(to-from+(2*bp.ext))/1000,"kb\n",sep="")
  }
  if(info & bp.ext>0) { cat("in addition to cM distance, window was extended by",
                     bp.ext,"base pairs on either side\n")} 
  from <- max(c(0,(from-bp.ext)))
  to <- min(c((to+bp.ext),get.chr.lens()[chr][1]),na.rm=T)
  return(c(from,to))
}


#' Meta-analysis using odds ratio and standard error from 2 datasets
#' 
#' This function calculates meta analysis odds ratios, standard errors and p-values
#' using results from a table containing odds ratio and standard error data for analyses
#' of 2 different datasets (typically logistic regression, but other analyses can be
#' incorporated if an odds-ratio and SE can be derived, for instance one analysis might
#' be a case control logistic regression GWAS and the other a family TDT analysis).
#' @param X A data.frame with column names which should be entered in the parameters:
#' OR1, OR2, SE1, SE2, and optionally N1, N2. 
#' @param OR1 The column name of X containing odds ratios from the first analysis
#' @param OR2 Same as OR1 above but pertaining to the second analysis
#' @param SE1 The column name of X containing standard errors from the first analysis
#' @param SE2 Same as SE1 above but pertaining to the second analysis
#' @param N1 Only required if method="sample.size". Either the column name in X with the 
#' number of samples in the first analysis, of a vector of the same, or if N's is the
#'  same for all rows, a scalar' value can be entered
#' for each.
#' @param N2 Only required if method="sample.size". Same as N1 above but pertaining to analysis 2
#' @param method character, can be either 'beta', 'z.score' or 'sample.size', and upper/lower
#' case does not matter. 'Beta' is the default and will calculate meta-analysis weights using
#' the inverse variance method (based on standard errors), and will calculate the p-values
#' based on the weighted beta coefficients of the two analyses. 'Z.score' also uses inverse variance
#' but calculates p-values based on the weighted Z scores of the two analyses. 'Sample.size' uses
#' the sqrt of the sample sizes to weight the meta analysis and uses Z scores to calculate p values
#' like 'Z.score' does.#' 
#' @return The object returned should have the same number of rows and rownames as the data.frame
#'  X but columns are the meta analysis stastistics, namely:
#'   OR.meta, beta.meta, se.meta, z.meta, p.meta, which will contain the meta
#' analysis odds-ratio, beta-coefficient, standard error, z-score, and p-values respectively
#' for each row of X.
#' @export
#' @examples
#' X <- data.frame(OR_CC=c(1.8,1.15),OR_Fam=c(1.33,0.95),SE_CC=c(0.02,0.12),SE_Fam=c(0.07,0.5))
#' rownames(X) <- c("rs689","rs23444")
#' X
#' meta.me(X)
#' X <- data.frame(OR_CC=c(1.8,1.15),OR_CC2=c(1.33,0.95),
#'  SE_CC=c(0.02,0.12),SE_CC2=c(0.02,0.05),
#'  n1=c(5988,5844),n2=c(1907,1774))
#' # even with roughly the same number of samples the standard error will determine the influence of
#' # each analysis on the overall odds ratio, note here that the second SE for dataset goes
#' # from 0.5 to 0.05 and as a result the estimate of the odds ratio goes from 1.137 to 0.977,
#' i.e, from very close to OR1, changing to very close to OR2.
#' meta.me(X,OR2="OR_CC2",SE2="SE_CC2") 
#' # sample size and z-score methods give similar (but distinct) results
#' meta.me(X,OR2="OR_CC2",SE2="SE_CC2",N1="n1",N2="n2",method="sample.size") 
#' meta.me(X,OR2="OR_CC2",SE2="SE_CC2",N1="n1",N2="n2",method="z.score")  # N's will be ignored
meta.me <- function(X,OR1="OR_CC",OR2="OR_Fam",SE1="SE_CC",SE2="SE_Fam",
                    N1=NA,N2=NA,method=c("beta","z.score","sample.size")) {
  #N1=18856,N2=7638
  validz <- c("beta","z.score","sample.size")
  method <- tolower(method[1])
  if(!method %in% validz) { method <- validz[1]; warning("invalid method entered, using 'beta' method") }
  if(!is(X)[1]=="data.frame") { stop("X must be a data.frame") }
  cnx <- colnames(X)
  if(is.null(rownames(X))) { rownames(X) <- paste(1:nrow(X)) }
  if(!all(c(OR1,OR2,SE1,SE2) %in% cnx)) { stop("X must contain column names specified by OR1,OR2,SE1,SE2") } 
  ok <- FALSE
  if(method=="sample.size") {
    if(is.numeric(N1) & is.numeric(N2)) { if(length(N1)!=1 | length(N2)!=1) { 
      stop("N1, N2 should either be scalar integers, or column names containing N's") } else {
        N.coln <- FALSE
      } }
    if(is.character(N1) & is.character(N2)) { 
      if(!all(c(OR1,OR2,SE1,SE2) %in% cnx)) { 
        stop("N1,N2 must contain either be scalar integers or column names with N's") } else {
          N.coln <- TRUE
        }
    }
  }
  OR.CC <- X[,OR1]
  beta.CC  <- log(X[,OR1])
  se.CC <- X[,SE1]
  OR.family <- X[,OR2]
  beta.family  <- log(X[,OR2])
  se.family <- X[,SE2]
  z.CC <- beta.CC/se.CC
  z.family <- beta.family/se.family
  
  inv.CC <- 1 / (se.CC^2)
  inv.family <- 1 / (se.family^2)
  var.meta <- 1 / (inv.CC+inv.family)
  weight.CC <- inv.CC * var.meta
  weight.family <- inv.family * var.meta
  se.meta <- round(sqrt(var.meta), digit=3)
  if(method=="sample.size") {
    if(N.coln) {
      famN <- X[,N2]
      ccN <- X[,N1]
    } else {
      famN <- N2 # 3819*2  #3509*2   #  3819*2   #  10796
      ccN <- N1 # 6683+12173 # including CBR, or for UVA analyses use instead: 9416+6670
    }  
    WeightFam = sqrt(famN)/(sqrt(famN)+sqrt(ccN))
    WeightCC <- 1-WeightFam
    # beta calculated the same way for sample.size method using sample size weights
    beta.meta <- round((WeightCC * beta.CC) + (WeightFam * beta.family),digit=3) # beta based
    z.meta <- round((WeightCC * z.CC) + (WeightFam * z.family),digit=6) # n-based, z-based
  } else {
    # beta calculated the same way for z.score method and beta method using inverse variance
    beta.meta <- round((weight.CC * beta.CC) + (weight.family * beta.family),digit=3) # beta based
    if(method=="z.score") {
      z.meta <- round((weight.CC * z.CC) + (weight.family * z.family),digit=6) # z-based
    } else {
      # default (beta) method
      z.meta <- beta.meta/se.meta
    }
  }  
  OR.meta <- exp(beta.meta)
  p.meta <- 2*pnorm(-abs(z.meta))
  out <- (cbind(OR.meta,beta.meta,se.meta,z.meta,p.meta)) 
  colnames(out) <- c("OR.meta","beta.meta","se.meta","z.meta","p.meta") 
  rownames(out) <- rownames(X)
  return(out)
}



#' Obtain nearby SNP-lists within a recombination window
#' 
#' For a snp.id (or list), extend a window around that chromosome location in recombination
#' units (centimorgans) and return the list of SNPs from the current ChipInfo object that
#' lie in this window. This is a way of extracting SNPs in linkage disequilibrium with an 
#' index SNP, that could also be plausible causal candidates.
#' @param snpid.list character, list of snp-ids (e.g, rs-id or chip id) to obtain lists for.
#' SNPs must all be from the same chromosome - if ranges for SNPs spanning multiple ranges
#' are desired, you must use multiple calls. A warning will be given if SNPs from the same
#' karyotype band are entered as index SNPs, as in a typical GWAS analysis only one SNP would
#' be used like this from each region, ignore the warning if this is not the case for your
#' application.
#' @param cM numeric, the number of centimorgans to extend the window either side of each SNP
#' @param bp.ext numeric, optional number of base-pairs to extend the window by in addition
#' to the centimorgan extension
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to use The 
#' default when build is NULL is to use the build from the current ChipInfo annotation.
#' Note that build37 is slower as the recombination reference is in build 36, so conv.36.37() 
#' and conv.37.36() need to be used to retrieve the appropriate ranges for build 37.
#' @param excl.snps character, a list of rs-id or chip-ids of SNPs to exclude from the list
#' returned, as, for instance, they may have failed quality control such as call-rate.
#' @param name.by.bands give labels to each sublist returned by the karotype/cytoband name, 
#' but faster not to do this
#' @return Returns a list of vectors of snp-ids falling within the window(s) specified and not
#' in 'excl.snps'. Each snp in 'snpid.list' will correspond to an element in the list returned.
#' If name.by.bands is TRUE, then these list elements will each be named using the local
#' karyotype/cytoband location
#' @export
#' @seealso snps.in.range, get.recombination.map, recwindow, conv.37.36, conv.36.37, exp.window.nsnp
#' @examples
#' result <- get.nearby.snp.lists("rs900569")
#' get.nearby.snp.lists("rs900569",cM=0.2,excl.snps=result[[1]]) # get SNPs within 0.1-0.2cM
#' # note that the same query can return a different set with build 36 versus 37
#' get.nearby.snp.lists(c("rs689","rs4909944"),cM=0.001,build=36,name.by.bands=FALSE)
#' get.nearby.snp.lists(c("rs689","rs4909944"),cM=0.001,build=37,name.by.bands=FALSE) 
get.nearby.snp.lists <- function(snpid.list,cM=0.1,bp.ext=0,build=NULL,excl.snps=NULL,name.by.bands=TRUE) {
  #if(!exists("all.support")) { print(load("all.support.RData")) }
  all.support <- get.support()
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  snpic.list <- rs.to.id(snpid.list)
  cyto <- get.cyto(dir=getwd(),GRanges=FALSE); cyto[["gene"]] <- rownames(cyto)
  #which.snps <- match(snpid.list,mcols(all.support)$rs.id)
  #if(any(is.na(which.snps))) { stop(paste("NAs in dbSNP match:",paste(snpid.list[is.na(which.snps)],collapse=","))) }
  snps.locs <- Pos(snpid.list)
  snps.chrs <- Chr(snpid.list)
  if(build=="hg18") {
    snps.locs36 <- snps.locs;  snps.locs37 <- conv.36.37(chr=snps.chrs,pos=snps.locs)[,"start"] 
  } else { 
    snps.locs37 <- snps.locs;  snps.locs36 <- conv.37.36(chr=snps.chrs,pos=snps.locs)[,"start"]
  }
  next.chr <- unique(Chr(snpid.list)); if(length(next.chr)>1) { stop("enter snpids from only 1 chromosome at a time!") }
  if(any(snps.locs!=sort(snps.locs))) { 
    warning("snp-ids not in position order, rearrangement is preferred but will attempt to continue")
    sort.back <- match(snps.locs,sort(snps.locs))
  } else { sort.back <- 1:length(snps.locs) }
  ddz <- snpic.list[duplicated(snpic.list)]
  if(length(ddz)>0) { warning("dup SNPs:",ddz,"\n") }
  snp.rd <- RangedData(ranges=IRanges(start=snps.locs,end=snps.locs,names=snpic.list),
                       space=rep(next.chr,length(snps.locs)))
  snp.rd <- toGenomeOrder(snp.rd,strict=T) # think it autosorts anyway, but just in case
  if(name.by.bands) {
    #snp.rd <- annot.cnv(snp.rd,gs=cyto,quiet=TRUE); colnames(snp.rd) <- "band"
    bands <- Band.pos(ranges=snp.rd) #   snp.rd$band
  }
  #prv(next.chr,cM,bp.ext)
  ## recwindow uses build36 only, so convert back afterwards
  nxt.window <- lapply(snps.locs36, function(X,...) { recwindow(start=X,...) },
                       chr=next.chr,window=cM,bp.ext=bp.ext,info=FALSE)
  if(build=="hg18") {
    st.window <- sapply(nxt.window, "[",1)
    en.window <- sapply(nxt.window, "[",2)
  } else {
    #prv(next.chr,nxt.window)
    nncc <- rep(next.chr,length(nxt.window))
    st.window <- conv.36.37(chr=nncc,pos=sapply(nxt.window, "[",1))[,"start"]
    en.window <- conv.36.37(chr=nncc,pos=sapply(nxt.window, "[",2))[,"start"]
  }
  pozz <- start(all.support)
  n.snps <- vector("list",length(st.window))
  for(cc in 1:length(st.window)) {
    n.snps[[cc]] <- which(chr(all.support)==next.chr &
                            pozz>=st.window[cc] & 
                            pozz<=en.window[cc] &
                            (!rownames(all.support) %in% excl.snps) &
                            (!mcols(all.support)$rs.id %in% excl.snps) 
    )
  }
  grp.labs <- lapply(n.snps,function(X) { rownames(all.support)[X] })
  if(name.by.bands) {
    if(length(unique(bands))!=length(bands)) { warning("these bands are not unique ==> ",paste(bands[duplicated(bands)],collapse=",")) }
    grpz <- 1:length(bands)
    names(grp.labs) <- paste(grpz,bands,sep=":")
  }
  grp.labs <- grp.labs[sort.back]
  return(grp.labs)
}



# specific to me - but not quite finished with it yet!
#  SEE "metaboChipAlleleCodesFreqsFromPlink.RData"
# already copied this function to specific, so can delete when done
#' @examples
#' hwe.fn <- "snpdataout.hwe"
#' HH <- plink.allele.coding(hwe.fn)
#' head(HH)
plink.allele.coding <- function(hwe.fn) {
  hh <- read.table(hwe.fn,header=T)
  FF <- strsplit(paste(hh[[6]]),split="/",fixed=T)
  hh[["F1"]] <- as.numeric(sapply(FF,"[",1))
  hh[["F2"]] <- as.numeric(sapply(FF,"[",2))
  hh[["F3"]] <- as.numeric(sapply(FF,"[",3))
  hh[["Min"]] <- with(hh,((2*F1)+F2)/(2*(F1+F2+F3)))
  hh[["Maj"]] <- with(hh,((2*F3)+F2)/(2*(F1+F2+F3)))
  hh[["flip"]] <- apply(cbind(hh$A1,hh$A2),1,function(x) { all(x!=sort(x)) })
  hh[["RAF"]] <- hh$Min
  hh$RAF[!hh$flip] <- hh$Maj[!hh$flip]
  return(hh)
}






# iFunctions
# internal
# allows an sapply style function to only work on valid values
clean.fn <- function(x,fail=NA,fn=function(x) { x }) {
  if(!is.null(x)) { 
    x <- x[!is.na(x)]; x <- x[(x!="")]; return(fn(x)) 
  } else {  return(fail) } 
}




#' Convert a chr:pos1-pos2 vector to a matrix
#' 
#' Takes standard text positions, such as what you might see on the UCSC genome browser, such as 
#' chr1:10,000,234-11,000,567 for a range, or chrX:234,432 for a SNP, and converts to 
#' with cols: chr, start, end.
#' @param text character vector, format like chr:pos1-pos2
#' @export
#' @return a matrix of the same length as 'ranges' with columns chr, start and end, and
#' rownames will be the same as the original text vector.
#' @seealso Ranges.to.txt
#' @examples
#' txt <- Ranges.to.txt(rranges())
#' convert.textpos.to.data(txt)
convert.textpos.to.data <- function(text) {
  do.one <- function(X) {
    chr.pos <- strsplit(X,":",fixed=T)[[1]]  
    chr.txt <- gsub("chr","",chr.pos[[1]],ignore.case=T)
    chr <- chr.txt #allow for X , as.integer(chr.txt); if(is.na(chr)) { print(chr.txt) }
    pos.txt <- strsplit(chr.pos[[2]],"-",fixed=T)[[1]]
    pos12 <- as.integer(pos.txt)
    out <- c(chr,pos12[1],pos12[2])
    if(is.na(out[3])) { out[3] <- out[2] }
    names(out) <- c("chr","start","end")
    return(out)
  }
  return(t(sapply(text,do.one)))
}


#' Convert GRanges/RangedData to chr:pos1-pos2 vector
#' 
#' Takes a RangedData or GRanged object from some annotation lookup functions and converts to standard text
#' positions, such as what you might see on the UCSC genome browser, such as 
#' chr1:10,000,234-11,000,567 for a range, or chrX:234,432 for a SNP. Useful for printing
#' messages, concatenating positions to a single vector, or creating queries for databastes.
#' @param ranges A RangedData or GRanges object
#' @export
#' @return a text vector of the same length as 'ranges' with notation as described above
#' representing each position in the 'ranges' object
#' @seealso convert.textpos.to.data
#' @examples
#' Ranges.to.txt(rranges())
Ranges.to.txt <- function(ranges) {
  if(!is(ranges)[1] %in% c("RangedData","GRanges")) { stop("Not a GRanges or RangedData object") }
  text.out.a <- paste0("chr",chr2(ranges),":",format(start(ranges),scientific=F,trim=T))
  text.out.b <- paste0("-",format(end(ranges),scientific=F,trim=T))
  text.out.b[start(ranges)==end(ranges)] <- ""
  text.out <- paste0(text.out.a,text.out.b)
  return(text.out)
}




#' Select ranges only within the 22 autosomes in a ranged data object
#' 
#' Select only data from autosomes from a GRanges/RangedData object.
#' Will exclude X,Y, mitochondrial chromosome rows, and can automatically
#' detect whether chromosomes are coded as 'chr1' or just '1', etc.
#' @param ranges A RangedData or GRanges object
#' @export
#' @return an object of the same format as the input (ranges), except
#' with non-autosomal ranges removed.
#' @examples
#' rand.ranges <- rranges(chr.range=20:26)
#' rand.ranges # should include some non-autosomes
#' select.autosomes(rand.ranges) # only autosomes remain
select.autosomes <- function(ranges,deselect=FALSE) {
  #must.use.package("genoset",bioC=T)
  typ <- is(ranges)[1]
  if(!typ %in% c("RangedData","GRanges")) { 
    warning("not a RangedData or GRanges object"); return(ranges) 
  }
  if(length(unique(chr2(ranges))) < length(levels(chr2(ranges)))) {
    # this fixes the problem when a subset of a ranges object with less 
    #  chromosomes still has empty chr slots from previous object
    if(typ=="RangedData") { ranges <- ranges[as.numeric(unique(chr2(ranges)))] }
  } #else { cat("ok\n") }
  Chrz <- (rownames(chrInfo2(ranges)))
  chrz <- tolower(paste(Chrz))
  if(length(grep("chr",chrz))>0) {
    select1 <- which(chrz %in% paste("chr",1:22,sep=""))
    select2 <- which(chrz %in% paste(1:22))
    if(length(select2)>length(select1)) { select <- select2 } else { select <- select1 }
  } else {
    select <- which(chrz %in% paste(1:22))
  }
  if(deselect) {
    ok.chrs <- Chrz[!select]
  } else {
    ok.chrs <- Chrz[select]
  }
  if(typ=="RangedData") {
    return(ranges[ok.chrs])
  } else {
    return(ranges[chr2(ranges) %in% ok.chrs,])
  }
}





# internal
# Remove that pesky 'elementmetadata.' prefix from column names that have been converted from GRanges
emd.rmv <- function(X, rmv.genome=TRUE) {
  if(has.method("mcols",X)) {
    colnames(mcols(X)) <- gsub("elementMetadata.","",colnames(mcols(X)))
    ii <- which(colnames(mcols(X))=="genome")
    if(length(ii)>0) {
      gn <- mcols(X)[,ii[1]]
      if(length(unique(gn))==1) {
        mcols(X) <- mcols(X)[,-ii[1]]
      }
    }
  } else {
    if(has.method("colnames",X)) {
      colnames(X) <- gsub("elementMetadata.","",colnames(X))
      ii <- which(colnames(X)=="genome")
      if(length(ii)>0) {
        gn <- X[,ii[1]]
        if(length(unique(gn))==1) {
          X <- X[,-ii[1]]
        }
      }
    } else {
      stop("X must have column names, expecting GRanges, RangedData or data.frame")
    }
  }
  return(X)
}



# internal# iFunctions
chrNames2 <- function(X) {
  X <- toGenomeOrder2(X)
  XX <- chrIndices2(X)
  return(rownames(XX))
}

# internal # iFunctions
# version of toGenomeOrder() that is guaranteed to work for IRanges or GRanges
toGenomeOrder2 <- function(X,...) {
  if(has.method("toGenomeOrder",X)) {
    return(toGenomeOrder(X))
  } else {
    alreadyThere <-("strand" %in% colnames(X))
    out <- toGenomeOrder(as(X,"GRanges"),strict=T)
    X <- as(out,"RangedData")
    if(("strand" %in% colnames(X)) & !alreadyThere) {
      X <- X[,-which(colnames(X) %in% "strand")]
    }
    return(X)
  }
}

#internal # iFunctions
# version of chrInfo() that is guaranteed to work for IRanges or GRanges
chrInfo2 <- function(X) {
  if(has.method("chrInfo",X)) {
    return(chrInfo(X))
  } else {
    out <- chrInfo(as(X,"GRanges"))
    return(out)
  }
}

#internal # iFunctions
# version of chrIndices() that is guaranteed to work for IRanges or GRanges
chrIndices2 <- function(X,...) {
  if(has.method("chrIndices",X)) {
    return(chrIndices(X,...))
  } else {
    out <- chrIndices(as(X,"GRanges"))
    return(out)
  }
}

#internal # iFunctions
# version of chr() that is guaranteed to work for IRanges or GRanges
chr2 <- function(X) {
  if(has.method("chr",X)) {
    return(chr(X))
  } else {
    if(is(X)[1]=="RangedData") {
      return(space(X))
    } else {
      if(is.null(X)) { warning("X was NULL, expecting RangedData/GRanges"); return(NULL) }
      warning("chr2() function applies only to RangedData objects, attempting to pass ",is(X)[1]," to chr()")
      return(chr(X))
    }
  }
}



#' Determine whether a function can be applied to a class/object
#' 
#' Wrapper for 'showMethods', allows easy testing whether a function
#' (can be specified as a string, or the actual function itself (FUN)) can be
#' applied to a specific object or class of objects (CLASS)
#' @param FUN the function to test, can be specified as a string, or the actual function itself
#' @param CLASS  a specific object or a class of objects specified by a string, e.g, "GRanges"
#' @export
#' @return returns logical (TRUE/FALSE), or if the function is not S4 will return an error,
#' although this could potentially be because the function's package has not been loaded.
#' @examples
#' has.method(chr,"GRanges")
#' has.method("chr","GRanges")
#' g.example <- rranges()
#' has.method(start, g.example)
#' has.method(duplicated, g.example) # no such method
#' ## not run # has.method("notAFunction","GRanges") # should return error
has.method <- function(FUN,CLASS) {
  if(!is.character(CLASS)) { CLASS <- class(CLASS) }
  if(!is.character(FUN) & !is.function(FUN)) { stop("FUN must be an R function, as a string or function") }
  test <- showMethods(FUN,classes=CLASS,printTo=F)
  if(length(grep("not an S4 generic function",test))>0) {
    stop(FUN," was not an S4 generic function or required package not loaded")
  }
  return(!(length(grep("No methods",test))>0))
}


#' Convert RangedData/GRanges to a data.frame
#' 
#' Convert a RangedData/GRanges object to a data.frame with columns
#' chr, start and end. Default is to only translate the chromosome and
#' position information, which is faster. Using 'include.cols'=TRUE
#' allows all the columns from 'ranged' to be taken across to the resulting
#' data.frame.
#' @param ranged A RangedData or GRanges object
#' @param include.cols logical, whether to also bring across non-positional
#' columns to the resulting data.frame
#' @param use.names logical, whether to keep the rownames from the
#' original object for the output. Only has an effect when include.cols=FALSE,
#' otherwise original rownames are always kept.
#' @export
#' @seealso data.frame.to.ranged, data.frame.to.granges
#' @return A data.frame with columns chr, start and end, and depending on
#' chosen parameters, the same rownames as the input, and optionally the
#' same additional columns.
#' @examples
#' rd <- rranges(9,GRanges=FALSE, fakeids=TRUE)
#' rd[["fakecol"]] <- sample(nrow(rd))
#' rd[["rs.id"]] <- paste0("rs",sample(10000,9))
#' ranged.to.data.frame(rd)
#' ranged.to.data.frame(rd,,FALSE)
#' ranged.to.data.frame(rd,TRUE) # keep all the columns
#' data.frame.to.granges(ranged.to.data.frame(rd,TRUE)) # inverse returns original
ranged.to.data.frame <- function(ranged,include.cols=FALSE,use.names=TRUE) {
  if(!include.cols) {
    u <- Ranges.to.txt(ranged)
    v <- convert.textpos.to.data(u)
    if(!is.null(rownames(ranged)) & nrow(ranged)==nrow(v) & use.names) { rownames(v) <- rownames(ranged) }
    return(v)
  } else {
    u <- as(ranged,"data.frame")
    cn <- tolower(colnames(u))
    if(is(ranged)[1]=="RangedData") {
      if("names" %in% cn) { 
        rownames(u) <- u[["names"]]
        u <- u[,-which(cn=="names")]
      } 
      if("space" %in% cn) { colnames(u)[which(cn=="space")] <- "chr" }
    } else {
      if(is(ranged)[1]=="GRanges") {
        if("seqnames" %in% cn) { colnames(u)[which(cn=="seqnames")] <- "chr" }
      } else {
        warning("'ranged' should be RangedData or GRanges, coercion could fail")
      }
    }
    return(u)
  }
}

#' Convert a data.frame with positional information to GRanges
#' 
#' Convert a data.frame containing chromosome and position information
#' to a GRanges object. Assumes the position information is contained in
#' columns named 'chr', 'start' and 'end' respectively (not case sensitive) 
#' although you can enter alternative column names for each as parameters. 
#' 'seqnames' will be automatically detected as an alternative to 'chr' if 
#' present. Column names that are default GRanges slot names such as 'seqnames',
#' 'ranges', 'strand', 'seqlevels', etc, will be removed during conversion, so
#' rename these if you want them to be translated into the resulting GRanges
#' objects' column metadata. If there is a column 'pos' but no columns 'start'
#'  and 'end' this will be detected automatically without needing to change
#'  the default parameters and start will equal end equals pos (ie., SNPs).
#' @param dat a data.frame with chromosome and position information 
#' @param ... additional arguments to data.frame.to.ranged(), namely:
#' ids, start, end, width, chr, exclude and build
#' @export
#' @seealso ranged.to.data.frame, data.frame.to.ranged
#' @return A RangedData or GRanges object. If 'dat' doesn't
#' use the default column names, specify these using parameters
#' ids, start, and end or width. Exclude will remove prevent any 
#' column names of 'dat' specified not to be translated to the 
#' returned GRanges object. 'build' specifies the 'genome'
#' slot of the resulting object. 'ids' allows specification of
#' a column to be converted to the rownames of the new object.
#' @examples
#' chr <- sample(1:22,10)
#' start <- end <- sample(1000000,10)
#' df1 <- cbind(chr,start,end)
#' data.frame.to.granges(df1) # basic conversion
#' width <- rep(0,10)
#' df2 <- cbind(chr,start,width)
#' data.frame.to.granges(df2,end=NULL,width="width") # define ranges with start and width
#' id.col <- paste0("ID",1:10)
#' rs.id <- paste0("rs",sample(10000,10))
#' df3 <- cbind(chr,start,end,id.col,rs.id)
#' data.frame.to.granges(df3) # additional columns kept
#' df4 <- cbind(chr,start,end,id.col,rs.id, ranges=1:10)
#' data.frame.to.granges(df4) # 'ranges' column excluded as illegal name
#' data.frame.to.granges(df4, exclude="rs.id") # manually exclude column
#' df5 <- cbind(chr,start,end,rs.id)
#' rownames(df5) <- paste0("ID",1:10)
#' data.frame.to.granges(df5) # rownames are kept
#' data.frame.to.granges(df4,ids="id.col") # use column of 'dat' for rownames
data.frame.to.granges <- function(dat,...) {
  return(data.frame.to.ranged(dat=dat,...,GRanges=TRUE))
}


#' Convert a data.frame with positional information to RangedData/GRanges
#' 
#' Convert a data.frame containing chromosome and position information
#' to a RangedData or GRanges object. Assumes the position information is contained in
#' columns named 'chr', 'start' and 'end' respectively (not case sensitive) 
#' although you can enter alternative column names for each as parameters. 
#' 'seqnames' will be automatically detected as an alternative to 'chr' if 
#' present. If there is a column 'pos' but no columns 'start' and 'end' this
#' will be detected automatically without needing to change the default parameters
#' and start will equal end equals pos (ie., SNPs). Column names that are default 
#' GRanges slot names such as 'seqnames', 'ranges', 'strand', 'seqlevels', etc, will
#' be removed during conversion, so rename these if you want them to be translated 
#' into the resulting object.
#' @param dat a data.frame with chromosome and position information 
#' @param ids character string, an optional column name containing ids which
#' will be used for rownames in the new object, as long as the ids are unique.
#' If not, this option is overridden and the ids will simply be a normal column
#' in the new object.
#' @param start character, the name of a column in the data.frame contain
#' the start point of each range. Not case sensitive. In the case of SNP
#' data, a column called 'pos' will also be automatically detected without
#' modifying 'start' or 'end', and will be used for both start and end.
#' @param end character, the name of a column in the data.frame containing the
#' end point of each range, can also use 'width' as an alternative specifier,
#' in which case 'end' should be set to NULL. Not case sensitive. In the case of SNP
#' data, a column called 'pos' will also be automatically detected without
#' modifying 'start' or 'end', and will be used for both start and end.
#' @param width the name of a column in the data.frame containing 'width' of
#' ranges, e.g, SNPs would be width=0. This is optional, with 'start' and 'end'
#' being the default way to specify an interval. If using 'width' you must
#' also set 'end' to NULL.  Not case sensitive.
#' @param chr character, the name of the column in the data.frame containing
#' chromosome values. The default is 'chr' but 'seqnames' will also be
#' detected automatically even when chr='chr'. Not case sensitive.
#' @param exclude character string, and column names from the data.frame to 
#' NOT include in the resulting S4 object.
#' @param build the ucsc build for the result object which will apply to the
#' 'universe' (RangedData) or 'genome' slot (GRanges) of the new object.
#' @param GRanges logical, whether the resulting object should be GRanges (TRUE),
#' or RangedData (FALSE)
#' @export
#' @seealso ranged.to.data.frame, data.frame.to.ranged
#' @return A RangedData or GRanges object. If 'dat' doesn't use the default 
#' column names 'chr', 'start'/'end' or 'pos', specify these using parameters 
#' 'ids', 'start', and 'end' or 'width'. Exclude will remove prevent any 
#' column names of 'dat' specified not to be translated to the returned GRanges
#' object. 'build' specifies the 'genome' slot of the resulting object. 'ids' 
#' allows specification of a column to be converted to the rownames of the new object.
#' @examples
#' chr <- sample(1:22,10)
#' start <- end <- sample(1000000,10)
#' df1 <- cbind(CHR=chr,Start=start,enD=end)
#' print(df1)
#' data.frame.to.granges(df1) # not case sensitive!
#' width <- rep(0,10)
#' df2 <- cbind(chr,start,width)
#' data.frame.to.granges(df2,end=NULL,width="width") # define ranges with start and width
#' id.col <- paste0("ID",1:10)
#' rs.id <- paste0("rs",sample(10000,10))
#' df3 <- cbind(chr,start,end,id.col,rs.id)
#' data.frame.to.granges(df3) # additional columns kept
#' df4 <- cbind(chr,start,end,id.col,rs.id, ranges=1:10)
#' data.frame.to.granges(df4) # 'ranges' column excluded as illegal name
#' data.frame.to.granges(df4, exclude="rs.id") # manually exclude column
#' df5 <- cbind(chr,start,end,rs.id)
#' rownames(df5) <- paste0("ID",1:10)
#' data.frame.to.granges(df5) # rownames are kept
#' data.frame.to.granges(df4,ids="id.col") # use column of 'dat' for rownames
data.frame.to.ranged <- function(dat, ids=NULL,start="start",end="end",width=NULL,
                                 chr="chr",exclude=NULL,build=NULL,GRanges=FALSE) 
{
  ## abandon longer names as they clash with function names
  st <- paste(start); en <- paste(end); ch <- paste(chr); wd <- paste(width)
  if((!chr %in% colnames(dat)) & ("seqnames" %in% colnames(dat)) & GRanges) { ch <- "seqnames" }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  must.use.package(c("genoset","IRanges"),T)
  g.illegal <- tolower(c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
                         "isCircular", "start", "end", "width", "element"))
  if(is.matrix(dat)) { dat <- as.data.frame(dat,stringsAsFactors=FALSE) }
  if(!is.data.frame(dat)) { stop("Error: not a dataframe")}
  key.nms <- c(ids,st,en,ch,wd)
  tries <- 0
  #print(key.nms); print(colnames(dat))
  while(!all(key.nms %in% colnames(dat))) { 
    colnames(dat) <- tolower(colnames(dat)); key.nms <- tolower(key.nms)
    st <- tolower(st); en <- tolower(en); ch <- tolower(ch); wd <- tolower(wd)
    if(tries>2) {
      if((tolower(st)=="pos" | tolower(en)=="pos") & !(tolower(st)=="pos" & tolower(en)=="pos")) {
        st <- en <- "pos"
      } else {
        if(tolower(st)=="start" & tolower(en)=="end") { st <- en <- "pos" }
      }
    }
    key.nms <- c(ids,st,en,ch,wd)
    tries <- tries+1
    if(tries > 3) { if(!all(c(st,en,ch) %in% colnames(dat))) {
      warning("chromosome and position columns not found") } ; break }
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
  ## not sure why here are adding 'chr' to X and Y?
  #this was here before? :  if(length(ch)>0) { ch1 <- gsub("Y","chrY",gsub("X","chrX",gsub("chr","",dat[[ch]],ignore.case=T))) } else { ch1 <- NULL }
  if(length(ch)>0) { ch1 <- gsub("chr","",dat[[ch]],ignore.case=T) } else { ch1 <- NULL }
  if(length(st)>0) { st1 <- as.numeric(dat[[st]]) } else { st1 <- NULL }
  if(length(en)>0) { en1 <- as.numeric(dat[[en]]) } else { en1 <- NULL }
  if(length(wd)>0) { en1 <- st1+as.numeric(dat[[wd]]) } # { en1 <- st1+dat[[wd]] }
  #print(length(st1)); print(length(en1)); print(length(id)); print(length(ch1))
  outData <- GRanges(ranges=IRanges(start=st1,end=en1,names=id),seqnames=ch1); genome(outData) <- build[1]
  #outData <- RangedData(ranges=IRanges(start=st1,end=en1,names=id),space=ch1,universe=build[1])
  ###  ###  ###  outData <- toGenomeOrder2(outData,strict=T)
  # note when adding data subsequently that 'RangedData' sorts by genome order, so need
  # to re-sort any new data before adding.
  if(is.null(rownames(outData))) { rownames(outData) <- paste(1:nrow(outData)) }
  reorder <- match(rownames(outData),id)
  more.cols <- colnames(dat)[!colnames(dat) %in% key.nms]
  more.cols <- more.cols[!more.cols %in% exclude]
  if(is(outData)[1]=="GRanges") { more.cols <- more.cols[!more.cols %in% g.illegal] }
  if(length(more.cols)>0) {
    for (cc in 1:length(more.cols)) {
      u <- dat[[more.cols[cc]]][reorder]; #prv(u)
      if(is(outData)[1]=="GRanges") {
        mcols(outData)[[more.cols[cc]]] <- u
      } else {
        outData[[more.cols[cc]]] <- u
      }
    }
  }
  if(GRanges) {
    return(as(outData,"GRanges"))
  } else {
    cncn <- colnames(mcols(outData))
    outData <- as(outData,"RangedData")
    if(any(cncn %in% "strand")) {
      outData <- outData[,-which(cncn=="strand")]
    }
    #outData <- toGenomeOrder2(outData,strict=T)
    return(outData)
  }
}


#' Select chromosome subset of GRanges or RangedData object
#' 
#' One of the main differences between RangedData and GRanges is the way
#' of selecting the subset for a chromosome. RangedData just uses [n] where
#' 'n' is the chromosome name or number. Whereas GRanges, does not have a
#' method like this, so need to select using [chr(X)==chr.num,]
#' This wrapper allows selection of a chromosome or chromosomes regardless of
#' whether the object is RangedData or GRanges type.
#' @param X A GRanges or RangedData object
#' @param chr Vector, the chromosome(s) (number(s) or name(s)) to select
#' @export
#' @return returns an object of the same type as X, with only the chromosome
#' subset specified.
#' @examples
#' some.ranges <- rranges(100,chr.range=1:10)
#' chr.sel(some.ranges,6)
#' more.ranges <- rranges(10, chr.range=21:25)
#' chr.sel(more.ranges,1:22) # gives warning
#' select.autosomes(more.ranges)
chr.sel <- function(X,chr) {
  typ <- is(X)[1]
  if(!typ %in% c("RangedData","GRanges","ChipInfo")) { stop("not a ChipInfo, GRanges or RangedData object") }
  if(!(is.character(chr) | is.numeric(chr))) { stop("chr must be character or numeric type") }
  if(is.numeric(chr)) { if(!all(chr %in% 1:99)) { 
    stop("illegal chromosome index, valid range 1-99 [although 1-28 typical for human]") } }
  if(typ=="RangedData") { return(X[chr])}
  all.chr <- chr2(X)
  if(!all(chr %in% unique(all.chr))) { 
    if(!any(chr %in% unique(all.chr))) { 
      warning("none of the specified chromosome indices were present in the GRanges object, returning NULL")
      return(NULL)
    } else { 
      warning("some of the specified chromosome indices were not present in the GRanges object") 
    }
  }
  return(X[all.chr %in% chr,])
}


#' Simulate a GRanges or RangedData object
#' 
#' For testing purposes, this function will generate a S4 ranged object
#' based on the human genome. The default is to produce ranges selected
#' from chromosomes, with probability of a position in each chromosome
#' equal to the length of that chromosome versus the whole genome. The
#' maximum position allocated within each chromosome will be within
#' the length bounds of that chromosome. You can specify SNPs (ie., start
#' =end), but the default is for random ranges. You can alter the UCSC
#' build to base the chromosome lengths on, and you can specify whether
#' chromosomes should appear as chr1,chr2,... versus 1,2,..
#' @param n integer, number of rows to simulate
#' @param SNP logical, whether to simulate SNPs (width 1, when SNPs=TRUE)
#'  or just ranges (when SNP=FALSE)
#' @param chr.range integer vector of values from 1 to 26, to specify which
#' chromosomes to include in the simulated object. 23-26 are X,Y,XY,MT 
#' respectively.
#' @param chr.pref logical, if TRUE chromosomes will be coded as chr1,chr2,...,
#' versus 1,2,.. when chr.pref=FALSE
#' @param order logical, if TRUE the object returned will be in genomic order,
#' otherwise the order will be randomized
#' @param equal.prob logical, when FALSE (default), random positions will be
#' selected on chromosomes chosen randomly according to the their length (i.e,
#' assuming every point on the genome has equal probability of being chosen.
#' If equal.prob=TRUE, then chromosomes will be selected with equal probability,
#' so you could expect just as many MT (mitochondrial) entries as Chr1 entries.
#' @param GRanges logical, if TRUE the returned object will be GRanges format,
#' or if FALSE, then RangedData format
#' @param build character, to specify the UCSC version to use, which has a small
#' effect on the chromosome lengths. Use either "hg18" or "hg19". Will also 
#' accept build number, e.g, 36 or 37.
#' @param fakeids logical, whether to add rownames with random IDs (TRUE) or 
#' leave rownames blank (FALSE). If SNP=TRUE, then ids will be fake rs-ids.
#' @export
#' @return returns a ranged object (GRanges or RangedData) containing data
#' for 'n' simulated genomic ranges, such as SNPs or CNVs across chromosomes in
#' 'chr.range', using UCSC 'build'.
#' @examples
#' rranges()
#' rr <- rranges(SNP=TRUE,chr.pref=TRUE,fakeids=TRUE)
#' width(rr) # note all have width 1
#' rr
#' tt <- table(chr(rranges(1000)))
#' print(tt/sum(tt)) # shows frequencies at which the chr's were sampled
#' tt <- table(chr(rranges(1000,equal.prob=TRUE)))
#' print(tt/sum(tt)) # shows frequencies at which the chr's were sampled
rranges <- function(n=10,SNP=FALSE,chr.range=1:26,chr.pref=FALSE,order=TRUE,equal.prob=FALSE,
                    GRanges=TRUE,build=NULL, fakeids=FALSE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(!is.numeric(chr.range)) { stop("chr.range must be a numeric integer vector rangeing from 1 to 26") }
  chr.range <- unique(chr.range[chr.range<=26 & chr.range>=1]) 
  cL <- get.chr.lens(mito=TRUE,build=build)[c(1:24,24,25)]
  if(equal.prob) {
    cP <- rep(1/length(chr.range),length(chr.range))
  } else {
    cP <- cL[chr.range]/sum(cL[chr.range]) # probabilities of a location being in each chromosome
  }
  if(!is.numeric(n)) { stop("'n' must be a numeric integer vector, representing the number of rows to simulate") }
  nn <- round(force.scalar(n,min=1,max=10^9,default=10))
  chrs <- sort(sample(chr.range,size=nn,replace=TRUE,prob=cP))
  cnts <- table(chrs)
  chr.lengths <- get.chr.lens(mito=TRUE)[c(1:24,24,25)]
  dubb <- if(SNP) { 1 } else { 2 }
  randoms <- starts <- ends <- vector("list",length(chr.range))
  for (cc in 1:length(cnts)){
    ii <- match(names(cnts[cc]),paste(chr.range))
    if(!is.na(ii)) {
      randoms[[cc]] <- oo <- sort(sample(chr.lengths[chr.range[ii]],size=cnts[cc]*dubb,replace=TRUE))
    } else { randoms[[cc]] <- NA }
    #prv(oo)
    if(!SNP) { 
      kk <- length(randoms[[cc]])/2
      #print(kk)
      starts[[cc]] <- randoms[[cc]][1+(2*(0:(kk-1)))]
      ends[[cc]] <- randoms[[cc]][2*(1:kk)]
    } 
  }
  if(SNP) { starts <- ends <- randoms }
  starts <- unlist(starts); ends <- unlist(ends)
  chrs <- paste0(if(chr.pref) { "chr" } else { "" },chrs)
  chrs <- gsub("23","X",chrs);  chrs <- gsub("24","Y",chrs)
  chrs <- gsub("25","XY",chrs);  chrs <- gsub("26","MT",chrs)
  gg <- GRanges(ranges=IRanges(start=starts,end=ends),seqnames=chrs)
  if(!order) {
    gg <- gg[order(rnorm(nrow(gg))),]
  } else {
    gg <- toGenomeOrder(gg,strict=TRUE)
  }
  if(!GRanges) { gg <- as(gg,"RangedData"); gg <- gg[,-1] }
  if(n==0) { return(gg[-1,])} # returns empty ranges if that's what you really want
  if(fakeids) {
    if(SNP) { rownames(gg) <- rsnpid(nrow(gg)) } else { rownames(gg)  <- rsampid(nrow(gg)) }
  }
  return(gg)
}



#' Extract chromosome numbers from GRanges/RangedData 
#' 
#' Sometimes chromosomes are codeds as 1:22, sometimes there is also X,Y, etc, sometimes it's 
#' chr1, ch2, etc. This function extracts the set of chromosome labels used by a ranged object 
#' (ie, GRanges or RangedData) and converts the labels to numbers in a consistent way, so
#' 1:22, X, Y, XT, MT ==> 1:26, and optionally you can output the conversion table of codes to
#' numbers, then input this table for future conversions to ensure consistency.
#' @param ranged GRanges or RangedData object
#' @param warn logical, whether to display a warning when non autosomes are converted to numbers
#' @param table.out logical, whether to return a lookup table of how names matched to integers
#' @param table.in data.frame/matrix, col 1 is the raw text names, col 2 is the integer that should be assigned,
#'  col 3 is the cleaned text (of col 1) with 'chr' removed. the required form is outputted by this function if
#'  you set 'table.out=TRUE', so the idea is that to standardize coding amongst several RangedData objects you
#'  can save the table each time and ensure future coding is consistent with this. Note that chromosomes 1-22, X,
#'  Y, XY, and MT are always allocated the same integer, so table is only useful where there are extra NT, COX, HLA
#'  regions, etc.
#' @return a set of integers of length equal to the number of unique chromosomes in the ranged data.
#' @export
#' @examples
#' gg <- rranges(1000)
#' chrNames(gg); chrNums(gg)
#' gg <- rranges(1000,chr.pref=TRUE) # example where chromosomes are chr1, chr2, ...
#' chrNames(gg); chrNums(gg)
#' lookup <- chrNums(gg,table.out=TRUE)
#' lookup
#' gg2 <- rranges(10)
#' chrNums(gg2,table.in=lookup) # make chromosome numbers using same table as above
chrNums <- function(ranged,warn=FALSE,table.out=FALSE,table.in=NULL) {
  #must.use.package("genoset",bioC=T)
  typ <- is(ranged)[1]
  if(!typ %in% c("RangedData","GRanges")) { warning("not a GRanges or RangedData object"); return(NULL) }
  lookup <- c("X","Y","XY","MT")
  txt1 <- chrNames2(ranged)
  txt <- gsub("chr","",txt1,fixed=T)
  nums <- suppressWarnings(as.numeric(txt))
  num.na <- length(nums[is.na(nums)])
  if(num.na>0) { 
    if(warn) { warning(paste("chromosome numbers requested for non-autosomes, will assign numbers >=23 to letters",
                             paste(txt[is.na(nums)],collapse=","))) }
    aux.ind <- match(txt,lookup)
    nums[!is.na(aux.ind)] <- 22+aux.ind[!is.na(aux.ind)]
    unmatched <- txt[is.na(nums)]
    if(!is.null(table.in)) {
      if((all(table.in[,1] %in% unmatched)) | (all(unmatched %in% table.in[,1]))) {
        if(all(unmatched %in% table.in[,1])) {
          out <- table.in[,2][match(unmatched,table.in[,1])]
          nums[is.na(nums)] <- as.numeric(out)
        } else {
          out <- table.in[,2][match(unmatched,table.in[,1])]
          nums[is.na(nums)][!is.na(out)] <- as.numeric(out)[!is.na(out)]
          st.num <- max(c(22+length(lookup),as.numeric(table.in[,2])),na.rm=T)+1
          nums[is.na(nums)][is.na(out)] <- st.num:(st.num+length(nums[is.na(nums)][is.na(out)])-1)
        }
      } else {
        out <- table.in[,2][match(unmatched,table.in[,1])]
        nums[is.na(nums)][!is.na(out)] <- as.numeric(out)[!is.na(out)]
        st.num <- max(c(22+length(lookup),as.numeric(table.in[,2])),na.rm=T)+1
        nums[is.na(nums)][is.na(out)] <- st.num:(st.num+length(nums[is.na(nums)][is.na(out)])-1)
      }
    } else {
      nums[is.na(nums)] <- 27:(27+length(nums[is.na(nums)])-1)
    }
  }
  if(table.out) {
    out <- cbind(txt1,nums,txt)
    return(out)
  } else {
    return(sortna(as.numeric(nums)))
  }
}



#' Expand genomic locations to the ranges covering the 'n' closest SNPs
#' 
#' Sometimes for chip data we want to create windows around some locus, and
#' fixed distance [see flank()], recombination distance [see recwindow()] or a number of SNPs 
#' might be used. This function allows expansion of regions according to a set number of SNPs.
#' The result gives two regions for each row of a GRanges or RangedData object describing
#' the start and end of the left flanking 'nsnp' region, and right flanking 'nsnp' region
#' respectively.
#' @param ranged a GRanges or RangedData object describing the locations for
#' which we want to find regions encompassing 'nsnps' closest SNPs.
#' @param snp.info An object of type: ChipInfo, RangedData or GRanges, describing the set of SNPs
#' you are using (e.g, chip annotation). If left as null the ChipInfo object from get.support() 
#' with default options() will be used
#' @param nsnp Number of nearest SNPs to return for each location
#' @param add.chr logical, whether to add a chromosome column for the output object
#' @seealso nearest.snp, get.support, recwindow
#' @export
#' @return Two regions for each row of a the 'ranged' object describing
#' the start and end of the left flanking 'nsnp' region, and right flanking 'nsnp' region
#' respectively. If 'ranged' has rownames these should stay in the same order in the resulting
#' object. Chromosome will be the final column if you set add.chr=TRUE.
#' @examples
#' rngs <- rranges()
#' # not run - slow ~5 seconds # exp.window.nsnp(rngs)
#' # not run - slow ~5 seconds # exp.window.nsnp(rngs,add.chr=TRUE)
exp.window.nsnp <- function(ranged,snp.info=NULL,nsnp=10, add.chr=FALSE) {
  if(is.null(snp.info)) { snp.info <- get.support() }
  if(!is(snp.info)[1] %in% c("ChipInfo","RangedData","GRanges")) { 
    stop("snp.info must be of type: ChipInfo, RangedData or GRanges") }
  snp.info <- toGenomeOrder2(snp.info,strict=TRUE); rw.cnt <- 1
  all.fl <- matrix(ncol=4+(as.numeric(add.chr)), nrow=0)
  for(cc in chrNums(ranged)) {
    nxt.nm <- rownames(snp.info[paste(cc)]); pos <- start(snp.info[paste(cc)])
    rng <- chr.sel(ranged,paste(cc)) # ranged[paste(cc)]
    st.en.snp <- range.snp(ranged=rng,snp.info=snp.info)
    fl <- matrix(ncol=4, nrow=nrow(st.en.snp))
    fl[,2] <- start(rng); fl[,3] <- end(rng);
    for(dd in 1:nrow(st.en.snp)) {
      x1 <- pos[max(1,match(st.en.snp[dd,1],nxt.nm)-nsnp)]
      x2 <- pos[min(length(nxt.nm),match(st.en.snp[dd,2],nxt.nm)+nsnp)]
      #print(x1); print(x2); print(length(x1)); print(length(x2))
      fl[dd,1] <- x1[1]
      fl[dd,4] <- x2[1]
    }
    if(add.chr) { 
      chrz <- cc
      fl <- cbind(fl,chrz)
    }
    all.fl <- rbind(all.fl,fl)
  }
  fl <- (all.fl)
  fl[fl[,1]>fl[,2],1] <- fl[fl[,1]>fl[,2],2]
  fl[fl[,3]>fl[,4],1] <- fl[fl[,3]>fl[,4],4]
  cN <- c("left.start","left.end","right.start","right.end")
  if(add.chr) { cN <- c(cN,"chr") }
  colnames(fl) <- cN
  if(nrow(fl)==nrow(ranged) & !is.null(rownames(ranged))) { rownames(fl) <- rownames(ranged) }
  return(fl)
}




#' Find closest SNPs to the ends of ranges
#' 
#' For given genome ranges (GRanges/RangedData) will try to find the closest snps to the end of the ranges.
#' @param ranged A GRanges or RangedData object specifying the range(s) you wish to find SNPs near the
#' ends of. Alternatively leave this parameter as NULL and specify ranges using chr, pos
#' @param snp.info ChipInfo/GRanges/Ranged data object describing the SNPs relevant to your query, e.g, 
#' SNPs on the chip you are using. If left NULL, the SNP set used will be that retrieved by get.support()
#' which will depend on your options() settings, see ?get.support for more info
#' @param chr optional alternative to 'ranged' input, use in conjunction with 'pos' to specify the ranges
#' to find the SNPs near the ends of.
#' @param pos matrix with 2 columns for start, end positions, or a single column if all ranges are SNPs.
#' An optional alternative to 'ranged' input, use in conjunction with 'chr' to specify the ranges
#' to find the SNPs near the ends of.
#' @param nearest will preferably find an exact match but if nearest=T, will fall-back on nearest match, 
#' even if slightly outside the range.
#' @export
#' @return a list of SNP-ids (rownames of 'snp.info') fulfilling the criteria, the output vector (character)
#' should be the same length as the number of ranges entered.
#' @examples
#' end.snp(chr=c(1:3),pos=cbind(c(100000,200000,300000),c(30000000,4000000,10000000)))
#' end.snp(rranges())
end.snp <- function(ranged=NULL,snp.info=NULL,chr=NULL,pos=NULL,nearest=T) {
  return(start.snp(ranged=ranged,snp.info=snp.info,chr=chr,pos=pos,start=F,end=T,nearest=nearest))
}


#' Find closest SNPs to the starts and ends of ranges
#' 
#' For given genome ranges (GRanges/RangedData) will try to find the closest snps to the starts and ends
#' of the ranges.
#' @param ranged A GRanges or RangedData object specifying the range(s) you wish to find SNPs near the
#' starts/ends of. Alternatively leave this parameter as NULL and specify ranges using chr, pos
#' @param snp.info ChipInfo/GRanges/Ranged data object describing the SNPs relevant to your query, e.g, 
#' SNPs on the chip you are using. If left NULL, the SNP set used will be that retrieved by get.support()
#' which will depend on your options() settings, see ?get.support for more info
#' @param chr optional alternative to 'ranged' input, use in conjunction with 'pos' to specify the ranges
#' to find the SNPs near the starts/ends of.
#' @param pos matrix with 2 columns for start, end positions, or a single column if all ranges are SNPs.
#' An optional alternative to 'ranged' input, use in conjunction with 'chr' to specify the ranges
#' to find the SNPs near the starts/ends of.
#' @param nearest will preferably find an exact match but if nearest=T, will fall-back on nearest match, 
#' even if slightly outside the range.
#' @export
#' @return a list of SNP-ids (rownames of 'snp.info') fulfilling the criteria, the output will be a matrix
#' which should have the same number of rows as the number of ranges entered.
#' @examples
#' range.snp(chr=c(1:3),pos=cbind(c(100000,200000,300000),c(30000000,4000000,10000000)))
#' range.snp(rranges())
range.snp <- function(ranged=NULL,snp.info=NULL,chr=NULL,pos=NULL,nearest=T) {
  return(start.snp(ranged=ranged,snp.info=snp.info,chr=chr,pos=pos,start=T,end=T,nearest=nearest))
}


#' Find closest SNPs to the starts of ranges
#' 
#' For given genome ranges (GRanges/RangedData) will try to find the closest snps to the starts 
#' of the ranges.
#' @param ranged A GRanges or RangedData object specifying the range(s) you wish to find SNPs near the
#' starts of. Alternatively leave this parameter as NULL and specify ranges using chr, pos
#' @param snp.info ChipInfo/GRanges/Ranged data object describing the SNPs relevant to your query, e.g, 
#' SNPs on the chip you are using. If left NULL, the SNP set used will be that retrieved by get.support()
#' which will depend on your options() settings, see ?get.support for more info
#' @param chr optional alternative to 'ranged' input, use in conjunction with 'pos' to specify the ranges
#' to find the SNPs near the starts of.
#' @param pos matrix with 2 columns for start, end positions, or a single column if all ranges are SNPs.
#' An optional alternative to 'ranged' input, use in conjunction with 'chr' to specify the ranges
#' to find the SNPs near the starts of.
#' @param start logical whether to return the SNP nearest the range starts
#' @param end logical whether to return the SNP nearest the range ends
#' @param nearest will preferably find an exact match but if nearest=T, will fall-back on nearest match, 
#' even if slightly outside the range.
#' @export
#' @return a list of SNP-ids (rownames of 'snp.info') fulfilling the criteria, the output will be a vector
#' which will have the same length as the input. Unless start=TRUE and end=TRUE, then will return a matrix
#' which should have the same number of rows as the number of ranges entered. Note that end.snp() is 
#' equivalent to using this function when end=TRUE and start=FALSE, and range.snp() is the same as setting
#' start=TRUE and end=TRUE.
#' @examples
#' start.snp(chr=c(1:3),pos=cbind(c(100000,200000,300000),c(30000000,4000000,10000000)))
#' start.snp(rranges())
start.snp <- function(ranged=NULL,snp.info=NULL,chr=NULL,pos=NULL,start=T,end=F,nearest=T) {
  # will preferably find an exact match but if nearest=T, will fall-back on nearest match
  #must.use.package("genoset",T)
  nmz <- NULL
  if(!is(ranged)[1] %in% c("RangedData","GRanges")) {
    if(!is.null(chr) & !is.null(pos)) {
      if(is.null(dim(pos))) { st <- pos[1]; en <- pos[2] } else {
        st <- pos[,1]; en <- pos[,2]
      }
      if(length(st)>length(chr)) { chr <- rep(chr[1],length(st)) } else { chr <- chr[1:length(st)] }
    } else {
      stop("if not using 'ranged' input, then chr and pos must be valid")
    }
  } else {
    st <- start(ranged); en <- end(ranged); chr <- chr2(ranged)
    nmz <- rownames(ranged)
  }
  if(is.null(snp.info)) { snp.info <- get.support() }  # load default chip
  if(!is(snp.info)[1] %in% c("ChipInfo","RangedData","GRanges")) {
    stop("snp.info must be of type ChipInfo, RangedData or GRanges")
  } else {
    if(is.null(rownames(snp.info))) {  rownames(snp.info) <- paste(1:nrow(snp.info)) }
  }
  st.snps <- en.snps <- character(length(chr)) ; prch <- 0
  for (cc in 1:length(chr)) {
    if(chr[cc]!=prch) { 
      ref <- chr.sel(snp.info,paste(chr[cc])) # snp.info[paste(chr[cc])]
      st.ref <- start(ref); rnref <- rownames(ref)
      if(is.null(ref)) { stop(paste("snp.info did not contain chr",chr[cc])) }
    }
    #exact
    if(start) { 
      ind <- match(st[cc],st.ref)
      if(any(is.na(ind)) & nearest) {
        difs <- abs(st[cc]-st.ref)
        ind <- which(difs==min(difs,na.rm=T))[1]
      }
      if(length(ind)==0) { st.snps[cc] <- NA } else {
        st.snps[cc] <- rnref[ind]
      }
    }
    if(end){
      ind2 <- match(en[cc],st.ref)
      if(any(is.na(ind2)) & nearest) {
        difs <- abs(en[cc]-st.ref)
        ind2 <- which(difs==min(difs,na.rm=T))[1]
      }
      if(length(ind2)==0) { en.snps[cc] <- NA } else {
        en.snps[cc] <- rnref[ind2]
      }
    }
    prch <- chr[cc]
  }
  if(start & !end) {
    return(st.snps)
  }
  if(!start & end) {
    return(en.snps)
  }
  #otherwise looks like want both
  out <- cbind(st.snps,en.snps)
  if(!is.null(nmz)) { if(length(nmz)==nrow(out)) { rownames(out) <- nmz } }
  return(out)
}


#' Force a valid genomic range, given the inputted coordinates
#'
#' Enter a pair of genomic locations representing a range for a given chromosome and this
#' function will ensure that no position is less than 1 or greater than the relevant chromosome
#'  lengths. Anything below will be coerced to 1, and anything above to the chromosome length.
#' @param Pos must be numeric, length 2, e.g, c(20321,30123)
#' @param Chr chromosome label
#' @param snp.info optional object to take boundaries from, the maxima and minima for each
#' chromosome within this object will take the place of the chromsome lengths / 1.
#' @build ucsc build, only need to enter if this differs from getOption("ucsc")
#' @dir directory to use for download of chromosome lengths (only if you wish to
#' keep the chromosome length file)
#' @export
#' @examples
#' pss <- ps <- c(345035,345035); ch <- 1
#' force.chr.pos(ps,ch)
#' pss[1] <- 0
#' force.chr.pos(pss,ch) # won't allow zero
#' pss[1] <- -1
#' force.chr.pos(pss,ch) # won't allow negative
#' pss[1] <- 645035012
#' force.chr.pos(pss,ch) # won't allow pos > chromosome length
force.chr.pos <- function(Pos,Chr,snp.info=NULL,build=NULL,dir=NULL) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  # convert any non autosomes to numbers:
  Chr <- paste(Chr)
  Chr[grep("c6",Chr,ignore.case=T)] <- 6  # prevent issues with c6_COX, c6_QBL  
  Chr[grep("X",Chr,ignore.case=T)] <- 23
  Chr[grep("Y",Chr,ignore.case=T)] <- 24
  Chr[grep("M",Chr,ignore.case=T)] <- 25
  Chr[grep("NT",Chr,ignore.case=T)] <- 26  # prevent issues with NT_11387, etc
  Chr <- as.numeric(Chr)
  if(any(!paste(Chr) %in% paste(c(1:26)))) { stop("invalid chromosome(s) entered") }
  if(any(paste(Chr) == paste(26))) { warning("'NT' chromosome(s) entered, not supported, NAs produced") }
  if(length(Pos)==2 & is.numeric(Pos)) {
    if(is(snp.info)[1]!="RangedData" & is(snp.info)[1]!="GRanges") { 
      maxln <- get.chr.lens(dir=dir,mito=T,autosomes=FALSE,build=build)[Chr] 
    } else { 
      maxln <- end(tail(snp.info[paste(Chr)],1)) # force start and end to be within 1:chr.len
    }
    mbs <- min(max(1,Pos[1]),(maxln-1)); mbe <- min(max(2,Pos[2]),maxln)
    return(c(mbs,mbe))
  } else {
    Pos <- NA; warning("Pos needs to be numeric length 2, min, max")
  }
  return(round(Pos))
}




#' Retreive GO terms from biomart for a given gene list
#' 
#' Gene-ontology terms (GO-terms) are commonly used for testing for simple functional
#' enrichment for pathways, etc. This function can retrieve biological function, 
#' cellular component, or molecular description, depending on the parameters chosen.
#' @param gene.list a list of gene, use HGNC names, like COMT, HLA-C, CTLA4, etc.
#' @param bio logical, whether to return biological process GO terms
#' @param cel logical, whether to return cellular component GO terms
#' @param mol logical, whether to return molecular function GO terms
#' @return data.frame containing the gene name in the first column, chromosome in the
#' second column, and the GO terms in the third column, where one gene has multiple
#' GO terms, this will produce multiple rows, so there will usually be more rows
#' than genes entered. The data.frame can have 3,4 or 5 columns depending on
#' how many GO terms are selected.
#' @export
#' @examples
#' get.GO.for.genes(c("CTLA4","PTPN2","PTPN22")) # biological terms (default)
#' get.GO.for.genes(c("CTLA4","PTPN2","PTPN22"),cel=TRUE) # add cellular GO terms
get.GO.for.genes <- function(gene.list,bio=T,cel=F,mol=F) {
  must.use.package(c("biomaRt","genoset","gage"),T)
  ens <- useMart("ENSEMBL_MART_ENSEMBL",
                 dataset="hsapiens_gene_ensembl",
                 host="may2009.archive.ensembl.org",
                 path="/biomart/martservice",
                 archive=FALSE)
  ens <- useDataset("hsapiens_gene_ensembl",mart=ens)
  data(egSymb)
  base.attr <- c("hgnc_symbol", "chromosome_name")
  if(bio) { base.attr <- c(base.attr,"go_biological_process_description") }
  if(cel) { base.attr <- c(base.attr,"go_cellular_component_description") }
  if(mol) { base.attr <- c(base.attr,"go_molecular_function_description") }
#  dat <- getBM(attributes = c("hgnc_symbol", "chromosome_name",
#                              "start_position", "end_position", "band"), filters = "hgnc_symbol",
#               values = egSymb[,2], mart = ens)
  results <- getBM(attributes = base.attr, filters = "hgnc_symbol",
                   values = c(gene.list), mart = ens)
  return(results)
}




#' Select all ranges lying within a chromosome window
#' 
#' Input a ranged object (ie., GRanges or RangedData) and this function will
#' return the subset from chromosome 'chr' and within the base-pair range specified
#' by 'pos', in units of 'unit'. By default ranges with ANY overlap are returned, but
#' it can be specified that it must be full overlap. Duplicates can be removed.
#' @param ranged GRanges or RangedData object
#' @param chr a chromosome, e.g, 1,2,3,...,22,X,Y,XY,MT or however chromosomes are 
#' annotated in 'ranged'
#' @param pos a numeric range (length 2), with a start (minima) and end (maxima), specifying
#' the window on the chromosome to select ranges from, base-pair units are specified by 'unit'.
#' @param full.overlap logical, the default is to return objects with ANY overlap with the window,
#' whereas setting this as TRUE, will only return those that fully overlap
#' @param unit the unit of base-pairs that 'pos' is using, eg, "b", "kb", "mb", "gb"
#' @param rmv.dup logical, whether to remove duplicate ranges from the return result. The default
#' is not to remove duplicates.
#' @return an object of the same type as 'ranged', but only containing the rows that
#' were within the specified bounds
#' @export
#' @examples
#' iG <- get.immunog.locs()[2,] # select the 2nd iG region
#' ciG <- chr(iG)  #  get the chromosome
#' posiG <- c(start(iG),end(iG) # get the region start and end
#' rr <- rranges(10000) # create a large random GRanges object
#' in.window(rr,chr=ciG,pos=posiG) # set with ANY overlap of iG
#' in.window(rr,chr=ciG,pos=posiG,TRUE) # set with FULL overlap of iG
#' in.window(rr,chr=6,pos=c(25,35),unit="mb") # look between 25 - 35 MB on chr6 [ie, MHC]
in.window <- function(ranged,chr,pos,full.overlap=F, unit=c("b","kb","mb","gb"), rmv.dup=FALSE) {
  if(length(pos)>2 | !is.numeric(pos)) { warning("pos should be a start and end numeric range"); return(NULL) }
  if(length(pos)==1) { pos <- rep(pos,2) }
  all.chrz <- unique(chrNames2(ranged))
  if(length(chr)>1 | (!chr %in% all.chrz)) { warning("chr must be a value in",comma(all.chrz)); return(NULL) }
  typ <- is(ranged)[1]
  if(!any(typ %in% c("RangedData","IRanges","GRanges","RangesList"))) { 
    warning("'ranged' should be a RangedData type or similar"); return(NULL) }
  pos <- pos*make.divisor(unit,"unit")
  #unit <- tolower(unit[1]) ; mult <- switch(unit,b=0,kb=3,mb=6,gb=9); pos <- pos*10^mult
  # get set of genes in a position range for a chromosome
  chr.genez <- chr.sel(ranged,paste(chr)) 
  if(full.overlap) {
    ranged <- chr.genez[which(start(chr.genez)>min(pos) & end(chr.genez)<max(pos)),]
  } else {
    # any overlap
    ranged <- chr.genez[which(end(chr.genez)>min(pos) & start(chr.genez)<max(pos)),]
  }
  if(rmv.dup) {
    # remove duplicate genes/exons
    ranged <- ranged[!(duplicated(start(ranged)) & duplicated(end(ranged)) & duplicated(width(ranged))),]
  }
  return(ranged)
}



#' Plot genes to annotate figures with genomic axes
#' 
#' Quite often it is helpful to visualize genomic locations in the context of Genes
#' in the same region. This function makes it simple to overlay genes on plots
#' where the x-axis is chromosomal location.
#' @param chr chromosome number/name that the plot-range lies on
#' @param pos position range on the chromosome 'chr' that the plot lies on, in base-pairs
#' @param scl character, the scale that the x axis uses, ie, "b","kb","mb", or "gb", meaning
#' base-pairs, kilobases, megabases or gigabase-pairs.
#' @param y.ofs numeric, y-axis-offset, depending on what units are on your y-axis,
#' you may need to add an offset so that the gene annotation is drawn at an appropriate
#' level on the vertical axis, this value should be the centre of annotation
#' @param width depending on the range of your y-axis, you might want to expand or reduce
#' the vertical width of the gene annotation (in normal graph units)
#' @param txt logical, TRUE to include the names of genes on top of their representation
#' on the plot, or if FALSE, genes are drawn without labels.
#' @param chr.pos.offset if for some reason zero on the x-axis is not equal to 'zero' on
#' the chromsome, then this offset can correct the offset. For instance if you were using
#' a graph of the whole genome and you were plotting genes on chromosome 10, you would
#' set this offset to the combined lengths of chromosomes 1-9 to get the start point
#' in the correct place.
#' @param gs GRanges or RangedData object, this is annotation for the location of genes.
#' This will be retrieved using get.gene.annot() if 'gs' is NULL. THere may be several reasons
#' for passing an object directly to 'gs'; firstly speed, if making many calls then you won't
#' need to load the annotation every time; secondly, if you want to use an alternative annotation
#' you can create your own so long as it is a GRanges/RangedData object and contains a column
#' called 'gene' (which doesn't strictly have to contain gene labels, it could be any feature
#' you require, eg., transcript names, etc).
#' @param dir character, location to store file with the gene annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param box.col genes are drawn as boxes, this sets the colour of the boxes
#' @param txt.col this sets the colour of the label text (Gene names)
#' @param join.col for exons, or multipart genes, joins are made between the sections with
#' a central line, this sets the colour of that line.
#' @param ... further arguments to 'rect', the graphics function used to plot the 'genes'.
#' @export
#' @return Returns a data.frame, GRanges or RangedData object, depending on input parameters. Contained
#' will be HGNC gene labels, chromosome and start and end positions, other information depends on 
#' specific parameters documented above
#' @export
#' @examples
#' # EXAMPLE PLOT OF SOME SIMULATED SNPS on chr21-p11.1 #
#' loc <- c(9.9,10.2)
#' Band(chr=21,pos=loc*10^6)
#' rr <- in.window(rranges(50000),chr=21,pos=loc,unit="mb") # make some random MHC ranges
#' # create some SNPs and plot
#' rr3 <- rr; end(rr3) <- start(rr3) 
#' rownames(rr3) <- paste0("rs",sample(10^6,nrow(rr3)))
#' plot.ranges(rr3,col="blue",scl="mb",xlim=loc,xlab="Chr21 position (Mb)",ylab="")
#' # NOW add UCSC hg18 GENE annotation to the plot #
#' ## not run ## plot.gene.annot(chr=21,pos=c(9.95,10.1),scl="mb",y.ofs=1,build=36)
plot.gene.annot <- function(chr=1, pos=NA, scl=c("b","kb","mb","gb"), y.ofs=0, width=1, txt=T, chr.pos.offset=0,
                            gs=NULL, build=NULL, dir=NULL, box.col="green", txt.col="black", join.col="red", ...)
{
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  dir <- validate.dir.for(dir,"ano")
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(!is(gs)[1] %in% c("RangedData","GRanges")) { gs <- get.gene.annot(dir=dir,build=build,GRanges=FALSE) }
  if(is(gs)[1]=="GRanges") { gs <- as(gs,"RangedData") }
  if(!"gene" %in% colnames(gs)) { warning("didn't find 'gene' column in annotation") ; return(NULL) }
  Col <- c("green", "darkgreen")
  if(all(is.na(pos))) { pos <- c(1,Inf) } 
  # get set of genes in range of the graph section + remove duplicate genes/exons
  rng.genez <- in.window(gs,chr,pos,full.overlap=F, rmv.dup=T,unit=scl)
  if(nrow(rng.genez)<1) { warning("no genes found in range") ; return(NULL) }
  old.cc <- 1 ; tp <- 2
  # set vertical alignments for annotation
  old.auto <- F
  y.cent <- y.ofs
  y.bot <- y.ofs-(width/2)
  y.top <- y.ofs+(width/2)
  # text alignment
  tps <- y.bot + c(.18,.29,.46,.64,.82)[c(1,3,5)]*width
  # x position with scaling (e.g, Mb units = 10^6)
  #unit <- tolower(unit[1]) ; mult <- switch(unit,b=0,kb=3,mb=6,gb=9); pos <- pos*10^mult
  x.scl <- make.divisor(scl)
  cnrlo <- start(rng.genez)/x.scl+chr.pos.offset
  cnrhi <- end(rng.genez)/x.scl+chr.pos.offset
  gnnm <- (rng.genez$gene)
  n.genes <- length(cnrlo)
  txt.cex <- .75; if(n.genes>10) { txt.cex <- .5 } ; if(n.genes>100) { txt.cex <- .35 } # more genes = smaller labels
  #print(n.genes)
  for (cc in 1:n.genes) {
    # draw rectangle and label for each gene
    #prv(cnrlo[cc],y.top,cnrhi[cc], y.bot)
    rect(cnrlo[cc],y.top,cnrhi[cc], y.bot,border=box.col,...)
  }
  for (cc in 1:n.genes) {
    if (gnnm[old.cc]==gnnm[cc] & cc!=1)
    {
      link <- c(cnrhi[old.cc],cnrlo[cc])
      if(link[1]<link[2]) {
        lines(link,y=rep(y.cent,2),lwd=1,col=join.col,lty="dotted")
      }
    } else {
      if(txt) {
        if(cnrlo[cc] < min(pos/x.scl)) { 
          txt.x <- mean(c(min(pos/x.scl),min(cnrhi[cc],max(pos/x.scl))),na.rm=T)  
        } else { txt.x <- cnrlo[cc] }
        text(txt.x,tps[tp],gnnm[cc],col=txt.col,cex=txt.cex,pos=4,las=2,offset=0)
      }
    }
    old.cc <- cc  
    tp <- tp+1; if(tp==4) {tp <- 1}; #if(n.genes <=5) { tp <- 2 }
  }
  return(rng.genez)
}


#internal
make.divisor <- function(unit=c("b","kb","mb","gb"), par.name="scale (scl)") {
  valid.units <- c("k","m","g","b")
  unit <- tolower(unit[1]);
  unit <- substr(unit,1,1)
  if(!unit %in% valid.units) { warning("invalid entry to ",par.name," defaulting to base-pairs") ; unit <- "b" }
  divisor <- switch(unit,k=1000,m=10^6, g=10^9, b=1)
  return(divisor)
}

#internal
plotdf <- function(expr,fn="myTempPlot.pdf") {
  pdf(fn)
  { expr }
  dev.off()
  cat("wrote plot to",cat.path(getwd(),fn),"\n")
}



#' Plot the locations specified in a GRanges or RangedData object
#' 
#' GRanges and RangedData objects are used in bioconductor to store genomic locations and
#' ranges, such as transcripts, genes, CNVs and SNPs. This function allows simple
#' plotting of this data directly from the ranged object. SNPs will be plotted as dots 
#' and ranges as lines. Either can be plotted using vertical bars at the start/end of each
#' range. There are options for labelling and other graphical parameters.
#' @param ranged GRanges or RangedData object with genomic ranges. Should only contain
#' one chromosome, but if not, the first will be used
#' @param labels by default labels for each range are taken from the rownames of 'ranged',
#' but if you want to use another column in the ranged object, specify the column name
#' or number to use to label these ranges on the plot. Or else input a character
#' vector the same length as ranged for custom labels.
#' @param do.labs logical, whether or not to display these labels
#' @param skip.plot.new logical, whether to append to an existing plot (TRUE), or start
#' a new plot (FALSE --> default)
#' @param lty line type to use, see '?lines()' - not used for SNP data when v.lines=FALSE
#' @param alt.y alternative y-axis values (other than the default ordering from the input)
#' This can be a vector of length 1 or length(ranged), or else a column name in ranged to 
#' take the values from
#' @param v.lines TRUE will plot the ranges as pairs of vertical lines, occupying the full
#' vertical extent of the plot, whereas FALSE will plot the ranges as individual horizontal lines
#' @param ylim numeric, length 2, the y-axis limits for the plot, same a 'ylim' for ?plot()
#' @param xlim numeric, length 2, the x-axis limits for the plot, same a 'xlim' for ?plot(),
#' This shouldn't usually be needed as the automatic x-limits should work well,
#'  however is here in case fine tuning is required.
#' @param scl character, the scale that the x axis uses, ie, "b","kb","mb", or "gb", meaning
#' base-pairs, kilobases, megabases or gigabase-pairs.
#' @param col character, colour, same as 'col' argument for plot(), etc.
#' @param srt integer, text rotation in degrees (see par) for labels
#' @param pos integer, values of ‘1’, ‘2’, ‘3’ and ‘4’, respectively indicate positions below, 
#' to the left of, above and to the right of the specified coordinates. See 'pos' in graphics:text()
#' @param lwd line width, see '?lines()' - not used for SNP data when v.lines=FALSE
#' @param pch point type, see '?points()' - not used for ranged data
#' @param cex font/symbol size, see '?plot()' - passed to plot, points if using SNP data 
#' @param ... further arguments to 'plot', so long as skip.plot.new==FALSE.
#' @export
#' @return Plots the ranges specified in 'ranged' to the current plot, or to a new plot
#' @example
#' rr <- in.window(rranges(5000),chr=6,pos=c(28,32),unit="mb") # make some random MHC ranges
#' rownames(rr) <- paste0("range",1:length(rr))
#' # plot ranges vertically #
#' plot.ranges(rr,v.lines=TRUE)
#' # make some labels and plot as horizontal lines #
#' rr2 <- rr[1:5,]; mcols(rr2)[["GENE"]] <- c("CTLA9","HLA-Z","BS-1","FAKr","teST")
#' plot.ranges(rr2,label="GENE",scl="Mb",col="black",
#'             xlab="Chr6 position (megabases)",
#'             yaxt="n",ylab="",bty="n")
#' # create some SNPs and plot
#' rr3 <- rr; end(rr3) <- start(rr3) 
#' rownames(rr3) <- paste0("rs",sample(10^6,nrow(rr3)))
#' plot.ranges(rr3,col="blue",yaxt="n",ylab="",bty="n")
plot.ranges <- function(ranged,labels=NULL,do.labs=T,skip.plot.new=F,lty="solid", alt.y=NULL,
                        v.lines=FALSE,ylim=NULL,xlim=NULL,scl=c("b","Kb","Mb","Gb"),
                        col=NULL,srt=0,pos=4,pch=1,lwd=1,cex=1,...) {
  #if(is(ranges)[1]!="RangedData") { warning("need RangedData object") ; return(NULL) }
  chk <- chrNums(ranged)
  typ <- is(ranged)[1]
  if(!is.null(alt.y)) {
    if(is.numeric(alt.y)) {
      if(length(alt.y)==1 | length(alt.y)==length(ranged)) {
        yy <- alt.y
      } else {
        warning("alt.y ignored, must be same length as ranged, or else length 1"); alt.y <- NULL
      }
    } else {
      if(is.character(alt.y)) {
        if(typ=="GRanges") { 
          cn <- colnames(mcols(ranged)); df <- mcols(ranged)
        } else {
          cn <- colnames(ranged); df <- ranged
        }
        if(!alt.y %in% cn) { stop("alternative y.axis column name ",alt.y," not found in 'ranged'") }
        yy <- df[,alt.y]; rm(df)
      } else { 
        warning("invalid value for alt.y, ignoring"); alt.y <- NULL
      }
    }
  }
  if(!is.null(labels)) {
    labels <- paste(labels)
    if(is.character(labels)) {
      if(length(labels)==1 | length(labels)==length(ranged)) {
        if(length(labels)==1) {
          if(typ=="GRanges") { 
            cn <- colnames(mcols(ranged)); df <- mcols(ranged)
          } else {
            cn <- colnames(ranged); df <- ranged
          }
          if(!labels %in% cn) { stop("labels column name ",labels," not found in 'ranged'") }
          lab <- df[,labels]; rm(df)
        } else {
          lab <- labels
        }
      } else {
        warning("labels ignored, must be same length as ranged, or else length 1"); labels <- NULL
      }
    } else {
      warning("invalid value for labels, ignoring"); labels <- NULL
    }
  } else {
    lab <- rownames(ranged) 
  } 
  if(length(chk)>1) { 
    warning(length(chk)," chromosomes in 'ranged', only using the first, chr",chk[1]) 
    ranged <- chr.sel(ranged,1) 
  }
  if(all(width(ranged)<=1)) { theyAreSnps <- TRUE } else { theyAreSnps <- FALSE }
  scl <- make.divisor(scl)
  xl <- range(c(start(ranged),end(ranged)))
  xl <- xl + ((diff(xl)*0.1)*c(-1,1))
  xl <- xl/scl
  nr <- nrow(ranged); if(is.null(nr)) { nr <- length(ranged) }
  if(is.null(alt.y)) {
    yl <- c(0,(nr+2))
  } else {
    yl <- range(yy)
  }
  if(is.numeric(ylim) & length(ylim)==2) {
    ylim <- range(ylim)
    ydif <- diff(ylim)
    yl <- ylim
  }
  if(is.numeric(xlim) & length(xlim)==2) {
    xlim <- range(xlim)
    xdif <- diff(xlim)
    xl <- xlim
  }
  if(is.null(alt.y)) {
    YY <- seq(from=yl[1],to=yl[2],length.out=nr+2)[-1]
  } else {
    if(length(yy)==1) { YY <- rep(yy,length(nr)) } else { YY <- yy }
  }
  #print(YY)
  if(!is.null(col)) {
    if(length(col)==1) {
      col <- rep(col,times=nr) 
    } else {
      if(length(col)!=nr) { warning("col was not the same length as ranged, using first only"); col <- rep(col[1],nr) }
    }
  }
  if(is.null(col)) {
    if(nr>22) { colz <- rep("black",nr) } else { colz <- get.distinct.cols(nr) }
  } else { colz <- col[1:nr] }
  if(is.null(lab) & do.labs) { lab <- paste(1:nr) } # last resort
  if(!skip.plot.new) {
    position <- c(start(ranged[1,]),end(ranged[1,]))/scl
    Y <- YY[c(1,1)]
    #prv(position,Y)
    TY <- if(theyAreSnps) { "p" } else { "l" }
    if(v.lines) {
      plot(x=position, y=Y, xlim=xl, ylim=yl, type=TY, col="white", lty=lty, ...)
      abline(v=position,col=colz[1])
    } else {
      plot(x=position, y=Y, xlim=xl, ylim=yl, type=TY, col=colz[1], lty=lty, lwd=lwd, cex=cex, ...)
    }
    st <- 2
  } else {
    st <- 1
  }
  if(nr>1 | st==1) {
    for (cc in st:nr) {
      if(v.lines) {
        abline(v=c(start(ranged[cc,]),end(ranged[cc,]))/scl,col=colz[cc],lty=lty)
      } else {
        if(theyAreSnps) { 
          points(x=c(start(ranged[cc,]),end(ranged[cc,]))/scl,y=YY[c(cc,cc)],col=colz[cc], pch=pch, cex=cex)
        } else { 
          lines(x=c(start(ranged[cc,]),end(ranged[cc,]))/scl,y=YY[c(cc,cc)],col=colz[cc], lty=lty, lwd=lwd)
        }
      }
    }
  }
  if(do.labs) {
    for (cc in 1:nr) {
      if(v.lines) { YY <- rep(tail(YY,1),length(YY)) }
      V.scale <- (diff(head(YY,2))*0.5)
      if(length(V.scale)<1 | srt!=90) { V.scale <- 0 }
      text(x=start(ranged[cc,])/scl,y=YY[cc]+V.scale,labels=lab[cc],cex=0.6,pos=pos,offset=0,srt=srt)
    }
  }
}

#internal function
chrnums.to.txt <- function(X,do.x.y=TRUE) {
  cond <- paste(X) %in% paste(1:22)
  if(any(cond)) { X[cond] <-  paste0("chr",X[cond]) }
  if(do.x.y) {
    X <- gsub("X","chrX",X)
    X <- gsub("Y","chrY",X)
    X <- gsub("23","chrX",X)
    X <- gsub("24","chrY",X)
    X <- gsub("25","chrXY",X)
    X <- gsub("26","chrM",X)
    X <- gsub("chrXchrY","XY",X)
    X <- gsub("chrYchrX","YX",X) 
    X <- gsub("M","chrM",X)
    X <- gsub("XY","chrXY",X)
    X <- gsub("chrchr","chr",X)
  } else {
    X[X %in% paste(23:100)] <- paste0("chr",X[X %in% paste(23:100)])
  }
  return(X)
}

#internal function
chrnames.to.num <- function(X,keep.let=TRUE) {
  X <- tolower(X)
    if(!keep.let) {
      X <- gsub("chrM","26",X)
      X <- gsub("chrXY","25",X) 
      X <- gsub("chrY","24",X)
      X <- gsub("chrX","23",X)
    } else {
      X <- gsub("chrM","M",X)
      X <- gsub("chrXY","XY",X) 
      X <- gsub("chrY","Y",X)
      X <- gsub("chrX","X",X)
    }
  X <- gsub("chrchr","",X) 
  X <- gsub("chr","",X)
  X <- toupper(X)
  return(X)
}


#' Change the chromosome labels in a RangedData or GRanges object to string codes
#' 
#' @param ranged A GRanges or RangedData object
#' @param do.x.y logical, if TRUE then the usual numbers allocated to chromosomes, X,Y,XY, MT will
#' be allocated as 23,24,25,26 respectively. If false, these will just have 'chr' appended as a
#' prefix
#' @param keep logical, whether to keep additional metadata columns in the new object 
#' @seealso set.chr.to.numeric
#' @export
#' @return returns the 'ranged' object, but wherever a chromosome number was previously, a character
#' label, e.g, 'chr1', or 'X', will returned to replace the number, e.g, 1 or 23 respectively. 
#' If table.out is TRUE will return a list where the first element is the resulting object, and the second 
#' element is a table showing which numbers were converted to what label This table
#' can then be used for future conversions via the parameter 'table.in' to ensure consistency of
#' coding.
#' @examples
#' x <- rranges()
#' x
#' x <- set.chr.to.numeric(x) # make entirely numeric
#' x <- rranges(chr.range=20:26)
#' # next two will give warning about X, Y, etc
#' set.chr.to.char(x) # 23 = chrX, etc
#' set.chr.to.char(x,do.x.y=FALSE) # 23=chr23, etc
set.chr.to.char <- function(ranged,do.x.y=T,keep=T) {
  #must.use.package("genoset",bioC=T)
  typ <- is(ranged)[1]
  if(!typ %in% c("GRanges","RangedData")) { warning("not a GRanges or RangedData object"); return(NULL) }
  if(length(grep("chr",chrNames2(ranged)))<length(chrNames(ranged))) {
    ranged <- toGenomeOrder2(ranged,strict=TRUE)
    #prv(ranged)
    mychr2 <- mychr <- paste(chr2(ranged))
    RN <- rownames(ranged)
    #all.nams <- chrNames2(ranged)
    #all.nums <- chrNums(ranged,table.in=table.in)
    if(length(grep("23",paste(mychr2)))>0) { 
      warning("use of arbitrary chromosome numbers for non-autosomes (i.e, >=23)",
              "can lead to annotation issues, try to use labels, X, Y, MT, and XY where possible") }
    sar <- select.autosomes(ranged)
    if(nrow(sar)>0) {
      all.nums.t <- chrNums(sar,table.in=NULL,table.out=T) 
      all.nams <- all.nums.t[,1]
      all.nums <- all.nums.t[,2]
      #mychr2 <- all.nums.t[,2][match(mychr,all.nums.t[,1])]
      for (cc in 1:length(all.nams)) { mychr2[which(mychr==all.nams[cc])] <- paste("chr",all.nums[cc],sep="") }
    } else {
      # no autosomes
    }
    if(do.x.y) {
      mychr2 <- gsub("X","chrX",mychr2)
      mychr2 <- gsub("Y","chrY",mychr2)
      mychr2 <- gsub("23","chrX",mychr2)
      mychr2 <- gsub("24","chrY",mychr2)
      mychr2 <- gsub("25","chrXY",mychr2)
      mychr2 <- gsub("26","chrM",mychr2)
      mychr2 <- gsub("chrXchrY","XY",mychr2)
      mychr2 <- gsub("chrYchrX","YX",mychr2) 
      mychr2 <- gsub("MT","chrM",mychr2)
      mychr2 <- gsub("XY","chrXY",mychr2)
      mychr2 <- gsub("chrchr","chr",mychr2)
    } else {
      mychr2[mychr2 %in% paste(23:100)] <- paste0("chr",mychr2[mychr2 %in% paste(23:100)])
    }
    #print(tail(mychr2)); print((all.nums))
    #prv(mychr2)
    if(any(is.na(mychr2))) { prv(mychr2[which(is.na(mychr2))]) }
    if(is.null(RN) | length(RN)!=nrow(ranged)) { RN <- 1:nrow(ranged) } #make sure RN's are valid
    if(is(ranged)[1]=="GRanges") {
      all.chr <- chr2(ranged)
      out <- GRanges(ranges=IRanges(start=start(ranged),end=end(ranged),names=RN),seqnames=mychr2)
    } else {
      out <- RangedData(ranges=IRanges(start=start(ranged),end=end(ranged),names=RN),space=mychr2)
    }
    out <- toGenomeOrder2(out,strict=TRUE)
    # prv(out)
    ###return(out)
    # need to allow for different indexing of dataframe part for GRanges
    ncr <- switch(typ,RangedData=ncol(ranged),GRanges=ncol(mcols(ranged)))
    if(is.null(ncr)) { ncr <- 0 }
    if(ncr>0 & keep) {
      cn <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
      for(cc in 1:ncr) {
        if(typ=="GRanges") {
          mcols(out)[[paste(cn[cc])]] <- mcols(ranged)[[paste(cn[cc])]]
        } else {
          out[[paste(cn[cc])]] <- ranged[[paste(cn[cc])]]
        }
      }
    }
    return(out)
  } else {
    #cat("no change\n")
    return(ranged)      # change not needed
  }
}




#' Change the chromosome labels in a RangedData or GRanges object to numbers
#' 
#' @param ranged A GRanges or RangedData object
#' @param keep logical, whether to keep additional metadata columns in the new object 
#' @param table.in matrix/data.frame object, usually a result of a prior run of 
#' set.chr.to.numeric(table.out=TRUE), which shows for each label (column 1), what
#' chromosome number should correspond. A way of ensuring consistent coding in different
#' sets.
#' @param table.out logical, if FALSE, the output will just be the object with updated 
#' chromosome labels. If TRUE, then the output will be a list, where the first element
#' is the updated object and the second object is a table describing the coding
#' scheme used to convert from labels to numeric indices.
#' @seealso set.chr.to.char
#' @export
#' @return returns the 'ranged' object, but wherever a chromosome label was previously a character
#' label, e.g, 'chr1', or 'X', will return as a number, e.g, 1 or 23 respectively. If table.out
#' is TRUE will return a list where the first element is the resulting object, and the second 
#' element is a table showing which labels were converted to what number. This table
#' can then be used for future conversions via the parameter 'table.in' to ensure consistency of
#' coding.
#' @examples
#' char <- rranges(chr.pref=T)
#' char
#' set.chr.to.numeric(char)
#' # behaviour with X, Y, etc
#' char <- rranges(chr.range=c(20:26))
#' #' char
#' set.chr.to.numeric(char)
#' tab <- set.chr.to.numeric(char,table.out=TRUE)[[2]]
#' tab # codes used in conversion #
#' char <- rranges(chr.range=c(20:26))
#' set.chr.to.numeric(char, table.in=tab) # code using codes from 'tab'
set.chr.to.numeric <- function(ranged,keep=T,table.in=NULL,table.out=FALSE) {
  typ <- is(ranged)[1]
  if(!typ %in% c("GRanges","RangedData")) { warning("not a GRanges or RangedData object"); return(NULL) }
  if(table.out | suppressWarnings(any(is.na(as.numeric(paste(chr2(ranged))))))) {
    silly.name <- "adf89734t5b"
    ranged <- toGenomeOrder2(ranged,strict=T)
    if(typ=="GRanges") {
      mcols(ranged)[[silly.name]] <- paste(1:nrow(ranged))
    } else {
      ranged[[silly.name]] <- paste(1:nrow(ranged))
    }
    #prv(ranged)
    mychr2 <- mychr <- paste(chr2(ranged))
    #all.nams <- chrNames2(ranged)
    #all.nums <- chrNums(ranged,table.in=table.in)
    all.nums.t <- chrNums(ranged,table.in=table.in,table.out=T) 
    all.nams <- all.nums.t[,1]
    all.nums <- all.nums.t[,2]
    #mychr2 <- all.nums.t[,2][match(mychr,all.nums.t[,1])]
    for (cc in 1:length(all.nams)) { mychr2[which(mychr==all.nams[cc])] <- all.nums[cc] }
    #print(tail(mychr2)); print((all.nums))
    if(typ=="GRanges") {
      all.chr <- chr2(ranged)
      out <- GRanges(ranges=IRanges(start=start(ranged),end=end(ranged)),seqnames=mychr2,silly.name=mcols(ranged)[[silly.name]])
    } else {
      out <- RangedData(ranges=IRanges(start=start(ranged),end=end(ranged)),space=mychr2,silly.name=ranged[[silly.name]])
    }
    out <- toGenomeOrder2(out,strict=T)
    if(typ=="GRanges") {
      oo <- mcols(out)[["silly.name"]]
      rr <- mcols(ranged)[[silly.name]]
    } else {
      oo <- out[["silly.name"]]
      rr <- ranged[[silly.name]]
    }
    if(all(!is.na(oo))) {
      if(is.null(rownames(ranged))) { rownames(ranged) <- paste(1:nrow(ranged)) ; rmv.rn <- TRUE } else { rmv.rn <- FALSE }
      #iioo <- rownames(ranged)[match(oo,rr)]; print(iioo); print(is(iioo))
      rn <- narm(rownames(ranged)[match(oo,rr)])
      if(nrow(out)==length(rn) ) { rownames(out) <- rn } else { warning("rownames did not match number of rows") }
      if(rmv.rn) { rownames(out) <- NULL }
    } else {
      warning("index column was corrupted")
    }
    # prv(out)
    ncr <- switch(typ,RangedData=ncol(ranged),GRanges=ncol(mcols(ranged)))
    if(is.null(ncr)) { ncr <- 0 }
    if(ncr>0 & keep) {
      cn <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
      for(cc in 1:ncr) {
        if(typ=="GRanges") {
          mcols(out)[[paste(cn[cc])]] <- mcols(ranged)[[paste(cn[cc])]]
        } else {
          out[[paste(cn[cc])]] <- ranged[[paste(cn[cc])]]
        }
      }
    }
    cno <- switch(typ,RangedData=colnames(out),GRanges=colnames(mcols(out)))
    if(any(cno %in% "silly.name")) { out <- out[,-which(cno %in% "silly.name")] }
    cno <- switch(typ,RangedData=colnames(out),GRanges=colnames(mcols(out)))
    if(any(cno %in% silly.name)) { out <- out[,-which(cno %in% silly.name)] }
    if(table.out) {
      return(list(ranged=out,table.out=all.nums.t))
    } else {
      return(out)
    }
  } else {
    #cat("no change\n")
    cnr <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
    if(any(cnr %in% "silly.name")) { ranged <- ranged[,-which(cnr %in% "silly.name")] }
    return(ranged)      # change not needed
  }
}




#' Retrieve locations of Immunoglobin regions across the genome
#' 
#' Returns the locations of immunoglobin regions in the human genome, for a given build, as
#' a list by chromosome, text vector, or GRanges/RangedData object.
#' For instance, for CNV research, these regions are known to be highly structurally complex
#' and can lead to false positive CNV-calls, so are often excluded.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param text logical, whether to return locations as a text vector of the form: chrN:xxxx-xxxx
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @export
#' @return Returns a list, GRanges or RangedData object, depending on input parameters. Contained
#' will be immunoglobin chromosome, start and end positions.
#' @example
#' get.immunog.locs()
#' get.immunog.locs(bioC=FALSE)
#' get.immunog.locs(text=TRUE,build=37)
get.immunog.locs <- function(build=NULL,bioC=TRUE,text=FALSE,GRanges=TRUE) {
  nchr <- 22
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(build[1]=="hg19") {
    #hg19
    chr <- c(22,14,2,14)
    stz <- c(22385572,105994256,89156874,22090057)
    enz <- c(23265082,107281230,89630187,23021097)
  } else {
    # hg18
    chr <- c(22,14,2,14)
    stz <- c(20715572,105065301,88937989,21159897)
    enz <- c(21595082,106352275,89411302,22090937)
  }
  nmz <- c("ig_c22","ig_c14_a","ig_c2","ig_c14_b")
  reg.dat <- rep("immunoglobin",length(chr))
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=chr,
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=T)
    if(text) { outData <- Ranges.to.txt(outData) } else {
      if(GRanges) { outData <- as(outData,"GRanges") }
    }
  } else {
    outData <- vector("list",nchr); names(outData) <- paste("chr",1:nchr,sep="")
    for (cc in 1:nchr) {
      if(cc %in% chr) {
        outData[[cc]] <- list(start=stz[chr==cc],end=enz[chr==cc])
      }
    }
  }
  return(outData)
}


#' Return Centromere locations across the genome
#' 
#' Returns the locations of centromeres in the human genome, for a given build, as
#' a list by chromosome, text vector, or GRanges/RangedData object.
#' @param dir character, location to store file with the this annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame.
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @param text logical, whether to return locations as a text vector of the form: chrN:xxxx-xxxx
#' @param autosomes logical, if TRUE, only return results for autosomes, if FALSE, also include
#' X and Y.
#' @export
#' @return Returns a list, GRanges or RangedData object, depending on input parameters. Contained
#' will be centromere chromosome and start and end positions.
#' @example
#' get.centromere.locs()
#' get.centromere.locs(bioC=FALSE,autosomes=TRUE)
#' get.centromere.locs(text=TRUE)
get.centromere.locs <- function(dir=NULL,build=NULL,
                                bioC=TRUE,GRanges=TRUE,text=FALSE,autosomes=FALSE)
{
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  dir <- validate.dir.for(dir,c("ano"),warn=FALSE); success <- TRUE
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  local.file <- cat.path(dir$ano,"cyto")
  tt <- get.cyto(build=build,bioC=FALSE,dir=dir,GRanges=FALSE)
  chrn <- paste(1:22)
  if(!autosomes) { 
    chrn <- c(chrn,c("X","Y"))
  }
  nchr <- length(chrn)
  my.chr.range <- vector("list",nchr)
  names(my.chr.range) <- paste("chr",chrn,sep="")
  for (cc in 1:nchr) {
    just.centros <- tt[paste(tt[,5])=="acen",]
    just.chr <- just.centros[which(paste(just.centros[,1])==names(my.chr.range)[cc]),]
    my.chr.range[[cc]] <- list(start=min(just.chr[,2]), end=max(just.chr[,3]))
  }
  reg.dat <- rep("centromere",nchr)
  nmz <- paste(reg.dat,chrn,sep="_")
  stz <- sapply(my.chr.range,"[[",1)
  enz <- sapply(my.chr.range,"[[",2)
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=TRUE)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=gsub("chr","",chrn),
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=TRUE)
    if(text) { 
      outData <- Ranges.to.txt(outData) 
    } else {
      if(GRanges){
        outData <- as(outData,"GRanges")
      }
    }
  } else {
    outData <- my.chr.range 
  }
  return(outData)
}


#' Return Cytoband/Karyotype locations across the genome
#' 
#' Returns the locations of cytobands/karyotype-bands in the human genome, for a given build, as
#' a data.frame, or GRanges/RangedData object.
#' @param dir character, location to store file with the this annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @param refresh logical, whether to re-download the file if the existing file has become corrupted
#' @export
#' @return Returns a list, GRanges or RangedData object, depending on input parameters. Contained
#' will be centromere chromosome and start and end positions.
#' @example
#' get.centromere.locs()
#' get.centromere.locs(bioC=FALSE)
#' get.centromere.locs(text=TRUE)
get.cyto <- function(build=NULL,dir=NULL,bioC=TRUE,GRanges=TRUE,refresh=FALSE) {
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  local.file="cyto"
  if(is.null(dir)) {
    local.file <- cat.path("",local.file,suf=build,ext="tar.gz")
  } else { 
    local.file <- cat.path(dir,local.file,suf=build,ext="tar.gz")
  }
  if(!file.exists(local.file) | refresh) {
    golden.path <- paste("http://hgdownload.cse.ucsc.edu/goldenPath/",build,"/database/cytoBand.txt.gz",sep="")
    success <- tryCatch( download.file(url=golden.path,local.file,quiet=T),error=function(e) { F } )
    if(is.logical(success)) {
      if(!success) { warning("couldn't reach ucsc website! try sourcing cytoband data elsewhere"); return(NULL) } }
    tt <- reader(local.file,header=FALSE)
    if(is.null(dir)) { unlink(local.file) }
  } else {
    tt <- reader(local.file)
  }
  colnames(tt) <- c("chr","start","end","band","negpos")
  write.table(tt,file=local.file,col.names=T,row.names=F,sep="\t",quote=F)
  mychr <- gsub("chr","",tt$chr,fixed=T)
  fullbands <- paste(mychr,tt$band,sep="")
  if(bioC ) {
    st <- as.numeric(tt$start)
    en <- as.numeric(tt$end)
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=st,end=en,names=fullbands),space=mychr,
                          negpos=tt$negpos,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=T)
    if(GRanges) { outData <- as(outData,"GRanges") }
  } else {
    outData <- tt 
    if("band" %in% colnames(outData)) {
      ## make 'chr-band' rownames to be consistent with the RangedData object if bioC=T
      rownames(outData) <- fullbands
      #outData <- outData[,-which(colnames(outData) %in% "band")]
    }
  }
  return(outData)
}



#' Get HapMap recombination rates for hg18 (build 36)
#' 
#' Recombination rate files can be used to calculate recombination distances
#' for genome locations, in centimorgans. This function downloads these reference
#' files from the hapmap NCBI website. At the time of writing they were only 
#' availble for build 36. If using a more recent build I suggest using the
#' conversion function conv.37.36(), then recwindow(), then conv.36.37() to 
#' get recombination distances for other builds. If getOption("save.annot.in.current")
#' is <=0 then no files will be kept. Otherwise an object containing this mapping data
#' will be saved in the local directory if dir=NULL, or else in the directory specified.
#' Allowing this reference to be saved will greatly increase the speed of this function
#' for subsequent lookups
#' @param dir character, location to store binary file with the recombination maps for
#' chromosomes 1-22. If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param verbose logical, if the binary file is not already downloaded, when verbose
#' is TRUE, there will be some output to the console indicating the progress of the
#' download. If FALSE, all output is suppressed.
#' @param refresh logical, if you already have the binary file in the current directory,
#' this argument will let you re-download and re-generate this file, e.g, if the file
#' is modified or corrupted this will make a new one without having to manually delete it
#' @param compress logical, this argument is passed to 'save' and will result in a larger
#' binary file size, but quicker loading times, so 'FALSE' is recommended for faster retrieval.
#' @export
#' @return Returns a list object of length 22, containing the recombination map files
#' as 22 separate data.frame's.
#' @example
#' ## not run as it takes roughly 2 minutes to download and read-in ##
#' ## uncomment the following 3 lines to run:
#' ## rec.map <- get.recombination.map(getwd())
#' ## file.on.disk <- "rrates_genetic_map_chr_1_22_b36.RData"
#' ## if(file.exists(file.on.disk)) { unlink(file.on.disk) } # remove the downloaded file
get.recombination.map <- function(dir=NULL,verbose=TRUE,refresh=FALSE, compress=FALSE) {
  n.chr <- 22
  hap.dir <- "http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/latest/rates/"
  temp.dir <- "recombinationratesGF13fDR1er119"
  local.file <- "rrates_genetic_map_chr_1_22_b36.RData"
  if(!file.exists(temp.dir)) { dir.create(temp.dir) } 
  local.files=paste0(temp.dir,"/genetic_map_chr",1:n.chr,"_b36.txt")
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(!is.null(dir)) {
    local.files <- cat.path(dir,local.files,ext="txt")
    local.file <- cat.path(dir,local.file,ext="RData")
  }
  if(!file.exists(local.file) | refresh) {
    if(verbose) { cat("Downloading recombination data from: ",hap.dir,"\n") }
    hapmap.urls <- cat.path(dir=hap.dir,fn=basename(local.files))
    success <- TRUE
    for (cc in 1:n.chr) {
      #print(hapmap.urls[cc])
      success <- tryCatch( download.file(url=hapmap.urls[cc],local.files[cc],quiet=T),error=function(e) { F } )
      if(verbose) { loop.tracker(cc,n.chr*2) }
    }
    if(is.logical(success)) {
      if(!success) { warning("couldn't download at least one of the files from: ",hap.dir); return(NULL) } }
    map.files.list <- vector("list",n.chr)
    for (cc in 1:n.chr) {
      map.files.list[[cc]] <- read.table(local.files[cc],header=TRUE)
      if(is.data.frame(map.files.list[[cc]])) { 
        unlink(local.files[cc]) 
      } else { warning("downloaded map file was corrupt for chr",cc) }
      if(verbose) { loop.tracker(n.chr+cc,n.chr*2) }
    }
    if(file.exists(temp.dir)) { file.remove(temp.dir) }   # delete the temporary directory
  } else {
    map.files.list <- reader(local.file)
  }
  if(!is.null(dir)) { save(map.files.list,file=local.file,compress=compress) }
  if(length(map.files.list)!=n.chr) { stop("Unfortunately the object derived seems corrupted") }
  names(map.files.list) <- paste0("chr",1:n.chr)
  return(map.files.list)
}


#' Get exon names and locations from UCSC
#' 
#' Various R packages assist in downloading exonic information but often the input required is 
#' complex, or several lines of code are required to initiate, returning an object that
#' might require some manipulation to be useful. This function simplifies the job 
#' considerably, not necessarily requiring any arguments. The object returned can be
#' a standard data.frame or a bioconductor GRanges/RangedData object. The raw annotation
#' file downloaded will be kept in the working directory so that subsequent calls to
#' this function run very quickly, and also allow use offline.
#' @param dir character, location to store file with the gene annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param transcripts logical, if TRUE, return transcripts rather than exons
#' @param GRanges logical, if TRUE and bioC is also TRUE, then returned object will be GRanges, otherwise
#' it will be RangedData
#' @export
#' @return Returns a data.frame, GRanges or RangedData object, depending on input parameters. Contained
#' will be HGNC gene labels, chromosome, start and end positions, transcript id number and name
#' @example
#' ## not run as it takes too long to download for CRAN ##
#' ## get.exon.annot()
#' ## get.exon.annot(bioC=FALSE,build=37)
get.exon.annot <- function(dir=NULL,build=NULL,bioC=T, transcripts=FALSE, GRanges=TRUE) {
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  ## load exon annotation (store locally if not there already)
  from.scr <- T
  txt <- if(transcripts) { "trans" } else { "exon" }
  if(!is.null(dir)) {
    dir <- validate.dir.for(dir,"ano")
    ex.fn <- cat.path(dir$ano,pref=txt,"Annot",suf=build,ext="RData")
    if(file.exists(ex.fn)) {
      tS <- reader(ex.fn)
      if(transcripts) {
        if(is(tS)[1]=="data.frame") { from.scr <- F }
      } else {
        if(is(tS)[1]=="GRanges") { from.scr <- F }
      }
    }
  }
  must.use.package("GenomicFeatures",T)
  if(from.scr) {
    must.use.package("gage",T)
    # get transcripts from build table 'knownGene'
    success <- tryCatch(txdb <- makeTranscriptDbFromUCSC(genome=build,
                                                         tablename="knownGene")  ,error=function(e) { F } )
    if(is.logical(success)) { 
      if(!success) {
        warning("Couldn't reach build website! try again later or, \n",
                "if in europe/uk, there may still be a bug in rtracklayer; \n",
                "Installing the latest version of R and bioconductor\n",
                "and running biocLite('rtracklayer'), should fix this")
        return(NULL) }
    }
    if(!transcripts) {
      ex = exonsBy(txdb, by="gene")
      tS <- toGenomeOrder(unlist(ex),strict=T)
    } else {  
      tS = transcriptsBy(txdb, by="gene")
      data(egSymb)
      select <- match(names(tS),egSymb[,1])
      names(tS)[!is.na(select)] <- egSymb[,2][select[!is.na(select)]]
      tS <- as.data.frame(tS)
    }
  }
  if(exists("ex.fn")) { save(tS,file=ex.fn) }
  if(transcripts) {
    chrs <- paste(tS$seqnames); chrs <- gsub("chr","",chrs)
    tS$seqnames <- chrs
    if(all(colnames(tS)==c("element","seqnames","start","end","width","strand","tx_id","tx_name"))) {
      colnames(tS) <- c("gene","chr","start","end","width","strand","txid","txname")
    } else {
      cat(" unexpected colnames found using makeTranscriptDbFrombuild()\n")
      if(bioC) { cat(" therefore returning data.frame instead of RangedData object\n") ;
                 bioC <- F }
    }
    if(bioC) {
      tS <- RangedData(ranges=IRanges(start=tS$start,end=tS$end),
                       space=tS$chr,gene=tS$gene, strand=tS$strand,
                       txid=tS$txid, txname=tS$txname,universe=build)
      tS <- toGenomeOrder2(tS,strict=T)
      if(GRanges) { tS <- as(tS,"GRanges") }
    }
    return(tS)
  } else {
    rownames(tS) <- paste(1:nrow(tS))
    ei <- mcols(tS)[["exon_id"]]
    ei2 <- add.trail(ei,suffix=strsplit("abcdefghijklmnopqrstuvwxyz","")[[1]])
    mcols(tS)[["exon_name"]] <- ei2
    if(bioC) {
      if(GRanges) {
        return(tS)
      } else {
        return(as(tS,"RangedData"))
      }
    } else {
      return(ranged.to.data.frame(tS))
    }
  }
}

# iFunctions
# internal, tidy chromosome names using extra chromosomal annotation into rough chromosomes
tidy.extra.chr <- function(chr,select=FALSE) {
  # most relevant to hg18
  chr <- paste(chr)
  SEL_c6 <- grep("c6",chr,ignore.case=T)
  SEL_c5 <- grep("c5",chr,ignore.case=T)
  SEL_NT <- grep("NT",chr,ignore.case=T)
  # most relevant to hg19
  SEL_LRG <- grep("LRG",chr,ignore.case=T)
  SEL_HG <- grep("HG",chr,ignore.case=T)
  SEL_GL <- grep("GL",chr,ignore.case=T)
  SEL_HS <- grep("HSCHR",chr,ignore.case=T)
  if(select) {
    # create TRUE/FALSE as to whether list elements have weird chromosome codes
    all <- unique(c(SEL_c6,SEL_c5,SEL_NT,SEL_LRG,SEL_HG,SEL_GL,SEL_HS))
    return(!1:length(chr) %in% all)
  } else {
    # transform weird chromosomes into more palatable codes
    chr[SEL_c6] <- 6  # prevent issues with c6_COX, c6_QBL  
    chr[SEL_c5] <- 5  # prevent issues with c5_H2  
    chr[SEL_NT] <- "Z_NT"  # merge all NT regions to one label
    chr[SEL_LRG] <- "Z_LRG"  # merge all NT regions to one label
    chr[SEL_HG] <- "Z_HG"  # merge all NT regions to one label
    chr[SEL_GL] <- "Z_GL"  # merge all NT regions to one label
    X <- names(table(chr))
    X <- X[grep("HSCHR",X)]
    if(length(X)>0) {
      HSC <- gsub("_","",substr(gsub("HSCHR","",X),1,2))
      for(cc in 1:length(X)) {
        #cat("replacing ",X[cc]," with ",HSC[cc],"\n",sep="")
        chr[chr==X[cc]] <- HSC[cc]
      }
    }
    return(chr)
  }  
}


#' Get gene names and locations from biomart
#' 
#' Various R packages assist in downloading genomic information but often the input required is 
#' complex, or several lines of code are required to initiate, returning an object that
#' might require some manipulation to be useful. This function simplifies the job 
#' considerably, not necessarily requiring any arguments. The object returned can be
#' a standard data.frame or a bioconductor GRanges/RangedData object. The raw annotation
#' file downloaded will be kept in the working directory so that subsequent calls to
#' this function run very quickly, and also allow use offline.
#' @param dir character, location to store file with the gene annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param duplicate.report logical, whether to provide a report on the genes labels that are listed
#' in more than 1 row - this is because some genes span ranges with substantial gaps within them
#' @param one.to.one logical, as per above, some genes have duplicate entries, sometimes for simplicity
#' you want just one range per gene, if this parameter is set TRUE, one range per gene is enforced,
#' and only the widest range will be kept by default for each unique gene label
#' @param remap.extra logical, whether to remap chromosome annotation for alternative builds and
#' unconnected segments to the closest regular chromosome, e.g, mapping MHC mappings to chromosome 6
#' @param discard.extra logical, similar to above, but if TRUE, then any non-standard chromosome
#' genes will just be discarded
#' @param only.named logical, biomart annotation contains some gene segments without names, if TRUE, then
#' such will not be included in the returned object (note that this will happen also if one.to.one is TRUE)
#' @param ens.id logical, whether to include the ensembl id in the dataframe
#' @param refresh logical, if you already have the file in the current directory,
#' this argument will let you re-download and re-generate this file, e.g, if the file
#' is modified or corrupted this will make a new one without having to manually delete it
#' @param GRanges logical, if TRUE and bioC is also TRUE, then returned object will be GRanges, otherwise
#' it will be RangedData
#' @export
#' @return Returns a data.frame, GRanges or RangedData object, depending on input parameters. Contained
#' will be HGNC gene labels, chromosome and start and end positions, other information depends on 
#' specific parameters documented above
#' @example
#' ## not run as it takes too long to download for CRAN ##
#' ## get.gene.annot()
#' ## get.gene.annot(bioC=FALSE,build=37)
get.gene.annot <- function(dir=NULL,build=NULL,bioC=TRUE,duplicate.report=FALSE,
                           one.to.one=FALSE,remap.extra=FALSE,discard.extra=TRUE,only.named=FALSE,
                           ens.id=FALSE,refresh=FALSE,GRanges=TRUE) {
  # faster than exon, but only contains whole gene ranges, not transcripts
  # allows report on duplicates as some might be confused as to why some genes
  # have more than one row in the listing (split across ranges usually)
  # run with dir as NULL to refresh changes in COX
  must.use.package(c("biomaRt","genoset","gage"),T)
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  from.scr <- T
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(!is.null(dir)) {
    dir <- validate.dir.for(dir,"ano")
    utxt <- ""; if(one.to.one) { utxt <- "_unq" }
    if(ens.id) { utxt <- paste(utxt,"ens",sep="_") }
    gn.fn <- cat.path(dir$ano,"geneAnnot",pref=build,suf=utxt,ext="RData")
    if(file.exists(gn.fn) & !refresh) {
      dat <- get(paste(load(gn.fn)))
      from.scr <- F
    }
  }
  # colnames for output
  nm.list <- c("gene","chr","start","end","band")
  if(ens.id & bioC) { warning("ens.id=TRUE only has an effect when bioC=FALSE") }
  if(from.scr) {
    if(build=="hg18") {
      ens <- useMart("ENSEMBL_MART_ENSEMBL",
                     dataset="hsapiens_gene_ensembl",
                     host="may2009.archive.ensembl.org",
                     path="/biomart/martservice",
                     archive=FALSE)
    } else {
      ens <- useMart("ensembl")
    }
    ens <- useDataset("hsapiens_gene_ensembl",mart=ens)
  
    
    attr.list <- c("hgnc_symbol", "chromosome_name",
                   "start_position", "end_position", "band")
    if(ens.id) { attr.list <- c(attr.list,"ensembl_gene_id") }
    if(build=="hg18" & only.named & !ens.id) {
      data(egSymb)
      dat <- getBM(attributes = attr.list, filters = "hgnc_symbol",
                 values = egSymb[,2], mart = ens)
    } else {
      dat <- getBM(attributes = attr.list,  mart = ens)
    }
    if(exists("gn.fn")) { save(dat,file=gn.fn) }
  } 
  if(ens.id) { nm.list <- c(nm.list,"ens.id") }
  #return(dat)
  no.gene.names <- which(dat[[1]]=="")
  if(one.to.one | (only.named & !ens.id & length(no.gene.names)>0)) { dat <- dat[-no.gene.names,] }
  if(remap.extra) {
    dat$chromosome_name <- tidy.extra.chr(dat$chromosome_name)
    #dat$chromosome_name[grep("c6",dat$chromosome_name,ignore.case=T)] <- 6  # prevent issues with c6_COX, c6_QBL  
    #dat$chromosome_name[grep("c5",dat$chromosome_name,ignore.case=T)] <- 5  # prevent issues with c5_H2  
    #dat$chromosome_name[grep("NT",dat$chromosome_name,ignore.case=T)] <- "Z_NT"  # merge all NT regions to one label
  }
  if(discard.extra) {
    ## http://www.lrg-sequence.org/ ##
    # note that if remapping is already done, then these won't be discarded unless remapping failed
    tt <- tidy.extra.chr(dat$chromosome_name,select=TRUE)
    badz <- which(!tt)
    if(length(badz)>0) { dat <- dat[-badz,] } # remove LRG, GS, HG, NT, COX, etc annotation from set
  }
  if(bioC) {
    outData <- RangedData(ranges=IRanges(start=dat$start_position,end=dat$end_position),
                          space=dat$chromosome_name,gene=dat$hgnc_symbol, band=dat$band, universe=build)
    outData <- toGenomeOrder2(outData,strict=T)
    if(duplicate.report | one.to.one) {
      genez <- outData$gene
      dG <- which(duplicated(genez))
      if(!duplicate.report) {
        dup.genes <- genez[dG]
      } else {
        dup.genes <- gene.duplicate.report(outData)
      }
      stz <- start(outData); enz <- end(outData); wdz <- width(outData)
      to.del <- to.ch <- ch.st <- ch.en <- NULL
      n.dup <- length(dup.genes)
      if(one.to.one & n.dup>0) {
        #return(dup.genes)
        indz <- sapply(as.list(dup.genes),function(X) { which(genez %in% X) })
        # keep the range that is widest
        st.en <- lapply(indz,function(X) { c(stz[X][wdz[X]==max(wdz[X])][1],enz[X][wdz[X]==max(wdz[X])][1]) } )
        to.del <- dG
        to.ch <- sapply(indz,min,na.rm=T)
        ch.st <- sapply(st.en,"[",1) 
        ch.en <- sapply(st.en,"[",2) 
        start(outData)[to.ch] <- ch.st 
        end(outData)[to.ch] <- ch.en 
        outData <- outData[-to.del,]
        cat("kept widest ranges, merging",length(to.del)+length(to.ch),"duplicate gene labels to",length(to.ch),"labels\n")
      } 
      # but haven't done anything about them or removed them! 
    }
    if(GRanges) { outData <- as(outData, "GRanges") }
  } else {
    outData <- dat; colnames(outData) <- nm.list
  }
  return(outData)
}

#internal
gene.duplicate.report <- function(ga,full.listing=F,colname="gene",silent=FALSE) {
  # for a RangedData object, report on any multiple listings for the same gene
  if(is(ga)[1]!="RangedData") { warning("not a RangedData object") ; return(NULL) }
  if(colname=="gene") {
    if("gene" %in% tolower(colnames(ga)))
    { 
      gene.col <- (which(tolower(colnames(ga)) %in% c("gene","genes","geneid")))
    } else {
      gene.col <- 0
    }
  } else {
    if(colname %in% colnames(ga)) { 
      gene.col <- which(colnames(ga)==colname) 
    } else { 
      stop("colname not found in ga") 
    } 
  }
  if(length(gene.col)>0) { gene.col <- gene.col[1] } else { warning("no 'gene' column"); return(NULL) }
  colnames(ga)[gene.col] <- "gene" #force this colname
  duplicate.report <- T  ### when would this be FALSE???
  if(duplicate.report) {
    culprits <- unique(ga$gene[which(duplicated(ga$gene))])
    n.gene.multi.row <- length(culprits)
    culprit.ranges <- ga[ga$gene %in% culprits,]
    total.culprit.rows <- nrow(culprit.ranges)
    start.same.ct <- end.same.ct <- 0; which.ss <- NULL
    for (cc in 1:length(culprits)) { 
      mini <- (ga[ga$gene %in% culprits[cc],]) 
      if(full.listing) {
        cat(colname,":",culprits[cc],"# same start:",anyDuplicated(start(mini)),
            "# same end:",anyDuplicated(end(mini)),"\n") }
      start.same.ct <- start.same.ct+anyDuplicated(start(mini))
      end.same.ct <- end.same.ct+anyDuplicated(end(mini))
      if(anyDuplicated(start(mini)) | anyDuplicated(end(mini))) { which.ss <- c(which.ss,cc) }
    }
    if(!silent) {
      cat(" ",colname,"s with split ranges:\n"); print(culprits,quote=F); cat("\n")
      cat(" ",colname,"s with same start or end:\n"); print(culprits[which.ss],quote=F); cat("\n")
      cat(" total ",colname,"-segments with same start",start.same.ct,"; total with same end:",end.same.ct,"\n")
    }
  }
  return(culprits)
}


#' Derive Telomere locations across the genome
#' 
#' Returns the locations of telomeres in the human genome, for a given build, as
#' a list by chromosome, text vector, or GRanges/RangedData object.
#' @param dir character, location to store file with the this annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param kb The number of base pairs at the start and end of a chromosome that are defined as
#' belonging to the telomere can be a little arbitrary. This argument allows specification
#' of whatever threshold is required.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @param text logical, whether to return locations as a text vector of the form: chrN:xxxx-xxxx
#' @param autosomes logical, if TRUE, only return results for autosomes, if FALSE, also include
#' X and Y.
#' @param mito.zeros logical, Mitochondria have no telomeres (are circular) but for some purposes you
#' might want zero values in order to match with other annotation that includes all chromosomes and MT.
#' TRUE adds zeros for chrMT, and FALSE excludes chrMT.
#' @export
#' @return Returns a text vector, GRanges or RangedData object, depending on input parameters. Contained
#' will be telomere chromosome and start and end positions.
#' @example
#' get.telomere.locs()
#' get.telomere.locs(bioC=FALSE)
#' get.telomere.locs(text=TRUE)
get.telomere.locs <- function(dir=NULL,kb=10,build=NULL,bioC=TRUE,GRanges=TRUE,
                              text=FALSE,autosomes=FALSE,mito.zeros=FALSE)
{
  # the actual telomeres are typically about 10kb, but
  # for cnv-QC purposes want to exclude a larger region like 500kb
  # Mt have no telomeres, are circular, but for some purposes might want zero values in there
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  chr.lens <- get.chr.lens(dir=dir,build=build[1],autosomes=FALSE,mito=mito.zeros)
  n <- 1:22; if(!autosomes) { n <- c(n,"X","Y") } # Mt have no telomeres, are circular
  nchr <- length(n) #default
  if(mito.zeros) { n <- c(n,"M") }
  my.chr.range <- vector("list",length(n))
  names(my.chr.range) <- paste("chr",n,sep="")
  for (cc in 1:nchr) {
    one <- force.chr.pos(Pos=c(1,kb*1000),Chr=cc,build=build) # f..c..pos() makes sure is a valid range
    two <- force.chr.pos(Pos=chr.lens[cc]+c(-kb*1000,0),Chr=cc,build=build)
    my.chr.range[[cc]] <- list(start=c(one[1],two[1]),end=c(one[2],two[2]))
  }
  if(mito.zeros) {
    # add null values for the Mitochondrial chromosome
    cc <- cc+1; one <- c(1,1);  two <- chr.lens[cc]+c(0,0)
    my.chr.range[[cc]] <- list(start=c(one[1],two[1]),end=c(one[2],two[2]))
  }
  reg.dat <- rep("telomere",length(n)*2)
  chrz <- rep(n,each=2)
  nmz <- paste(reg.dat,chrz,rep(c("a","b"),times=length(n)),sep="_")
  stz <- as.vector(sapply(my.chr.range,"[[",1))
  enz <- as.vector(sapply(my.chr.range,"[[",2))
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=chrz,
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=T)
    if(text) { outData <- Ranges.to.txt(outData) } else { if(GRanges) { outData <- as(outData,"GRanges") } }
  } else {
    outData <- my.chr.range 
  }
  return(outData)
}



#' Get chromosome lengths from build database
#' 
#' Quick and easy way to retrieve human chromosome lengths. Can select from hg18/hg19 (ie, 
#'  build 36/37), or any future builds (hg20, etc) stored in the same location on the build website.
#'  Default is to return lengths for 22 autosomes, but can also retrieve X,Y 
#'  and Mitochondrial DNA lengths by 'autosomes=FALSE' or n=1:25. Even if not connected to 
#'  the internet can retrieve hard coded lengths for hg18/hg19.
#'  @param dir directory to retrieve/download the annotation from/to (defaults to current getwd())
#'  if dir is NULL then will automatically delete the annotation text file from the local directory
#'   after downloading
#'  @param build string, currently 'hg17','hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is getOption("ucsc"). Will also accept integers 17,18,19,35,36,37 as alternative arguments.
#'  @param autosomes logical, if TRUE, only load the lengths for the 22 autosomes, else load X,Y,[MT] as well
#'  @param len.fn optional file name to keep the lengths in
#'  @param mito logical, whether to include the length of the mitochondrial DNA (will not include unless autosomes is also FALSE)
#'  @param names logical, whether to name the chromosomes in the resulting vector
#'  @param delete.after logical, if TRUE then delete the text file that these lengths were downloaded to. 
#'  If FALSE, then the file will be kept, meaning future lookups will be faster, and available offline.
#'  @export
#'  @examples
#'  get.chr.lens(delete.after=TRUE) # delete.after simply deletes the downloaded txt file after reading
#'  get.chr.lens(build=35,autosomes=TRUE,delete.after=TRUE) # only for autosomes
#'  get.chr.lens(build="hg19",mito=TRUE,delete.after=TRUE) # include mitochondrial DNA length
get.chr.lens <- function(dir=NULL,build=NULL,autosomes=FALSE,len.fn="humanChrLens.txt",
                         mito=FALSE,names=FALSE, delete.after=FALSE, verbose=FALSE)
{
  # retrieve chromosome lengths from local annotation file, else download from build
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(dir)) { dir <- getwd() ; delete.after <- TRUE }
  if(is.null(build)) { build <- getOption("ucsc") }
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  chrlens.f <- cat.path(dir$ano,len.fn) # existing or future lengths file
  build <- ucsc.sanitizer(build)
  n <- 1:22; if(!autosomes) { n <- c(n,"X","Y","M") }
  hg18.backup <- c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,
                   146274826,140273252,135374737,134452384,132349534,114142980,106368585,
                   100338915,88827254,78774742,76117153,63811651,62435964,46944323,
                   49691432,154913754,57772954,16571)
  hg19.backup <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
                   146364022,141213431,135534747,135006516,133851895,115169878,107349540,
                   102531392,90354753,81195210,78077248,59128983,63025520,48129895,
                   51304566,155270560,59373566,16571)
  
  # backups for offline use
  if(build=="hg18") { offline.backup <- hg18.backup  } else { offline.backup <- hg19.backup }
  if(file.exists(chrlens.f))
  {
    # file seems to be in annotation directory already
    chrLens <- readLines(chrlens.f)
    if (length(chrLens)!=length(n))
    {
      #warning("Length of existing chromosome file didn't match expected:",length(n))
      notGot <- T
    } else {
      notGot <- F
      # we have the right length, but do we have the right version?
      if(build=="hg18" & length(which(chrLens %in% hg19.backup))>2) { notGot <- T }
      if(build=="hg19" & length(which(chrLens %in% hg18.backup))>2) { notGot <- T }
      names(chrLens) <- paste0("chr",n)
    }
  } else { notGot <- T }
  if (notGot | (!build %in% c("hg18","hg19"))) {
    #download from build
    if(verbose) { cat("attempting to download chromosome lengths from genome build ... ") }
    urL <- switch(build,
                  hg17="http://hgdownload.cse.ucsc.edu/goldenPath/hg17/database/chromInfo.txt.gz",
                  hg18="http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chromInfo.txt.gz",
                  hg19="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz",
                  hg20="http://hgdownload.cse.ucsc.edu/goldenPath/hg20/database/chromInfo.txt.gz")
    success <- T
    success <- tryCatch(download.file(urL, chrlens.f,quiet=T),error=function(e) { F } )
    if(!is.logical(success)) { success <- T }
    if(success) {
      if(verbose) {  cat("download successful\n") }
      chrL.f <- readLines(chrlens.f)
      len.lst <- strsplit(chrL.f,"\t")
      nmz <- sapply(len.lst,"[",1)
      lnz <- sapply(len.lst,"[",2)
      nnn <- paste0("chr",n)
      want.chr.names <- match(nnn,nmz)
      want.chr.names <- want.chr.names[!is.na(want.chr.names)]
      #print(want.chr.names)
      chrLens <- lnz[want.chr.names]
      names(chrLens) <- nnn
    } else {
      warning("couldn't reach build website, so have used offline versions of chr lengths")
      if(!build %in% c("hg18","hg19")) { warning("no offline version for build version:",build) }
      chrLens <- paste(offline.backup)[1:length(n)]; names(chrLens) <- paste0("chr",n)
      #print(n)
      delete.after <- F
    }
    if(length(dir)==1 & dir[1]=="" & delete.after) {
      unlink(chrlens.f)
    } else {
      writeLines(chrLens,con=chrlens.f) # save file for future use
    }
  }
  if(!mito & length(chrLens)>22) { chrLens <- chrLens[-grep("M",n)] }
  if(names) {
    return(chrLens)
  } else {
    return(as.numeric(chrLens))
  }
}


#Internal: Read in a plink formatted pedigree/family file
# This function will import a PLINK style
# ped file and return a data.frame object in the same form
read.ped.file <- function(fn,keepsix=TRUE) {
  rr1 <- reader(fn,header=TRUE)
  if(ncol(rr1)<6) { warning("invalid ped/fam file, should have at least 6 columns"); return(NULL) }
  if(any(colnames(rr1) %in% c("X0","X1","X2"))) { rr1 <- reader(fn,header=FALSE) }
  colnames(rr1) <- gsub("X.","",colnames(rr1))
  if(any(colnames(rr1)[1] %in% unique(rr1[,1]))) { rr1 <- reader(fn,header=FALSE) }
  colnames(rr1)[1:6] <- c("family","sample","father","mother","sex","phenotype")
  if(keepsix) { rr1 <- rr1[,1:6] }
  return(rr1)
}


#' Import a ped file to pedData format used by snpStats
#'
#' PLINK ped files (family information files) are not in the
#' same format as ped files required for use with snpStats, for
#' instance for the tdt.snp() function to conduct a transmission
#' disequilibrium test. This function will import a PLINK style
#' ped file and return a 'pedData' object in the correct form
#' for snpStats and other rpackages. The plink file format is:
#' column 1: family id, column 2: individual id: column 3:
#' father's ID or 0, column 4: mother's ID or 0, column 5: sex of subject,
#' column 7: phenotype of subject.
#' @param file character, the name of a valid PLINK 'ped'/'fam' file
#' @param correct.codes logical, if TRUE, then where coding seems
#' inconsistent, e.g, mother and father both the same, this will
#' try to fix this automatically
#' @param silent logical, when false, output will show what
#' operations and transformations are being done, when TRUE,
#' the function produces no interim text to the console
#' @export
#' @return return a pedData object (a data.frame)
#' @examples
#' # not run # fn <- "myfile.ped" # insert name of your own file
#' # not run # myPed <- read.pedData(fn); prv(myPed)
read.pedData <- function(file,correct.codes=TRUE,silent=FALSE) {
  rr <- read.ped.file(file,keepsix = TRUE)
  want <- c("familyid","individual","father","mother","sex","affected")
  have <- colnames(rr)
  if(!silent) { cat(paste("mapping column",have,"==>",want,"\n"),sep="") }
  if(ncol(rr)>6) { warning("expected 6 columns in pedfile, unexpected behaviour could result") }
  if(ncol(rr)<6) { stop("need at least 6 columns in pedfile to map to headings ",paste(want,collapse="")) }
  colnames(rr)[1:length(want)] <- want
  rr <- rr[order(rr$familyid),]
  rr[["member"]] <- unlist(tapply(rep(1,nrow(rr)),factor(rr$familyid),cumsum))
  rr <- rr[,c(2,1,7,3:6)]
  rr <- shift.rownames(rr,T)
  rr[["father"]] <- rr$member[match(rr$father,rownames(rr))]
  rr[["mother"]] <- rr$member[match(rr$mother,rownames(rr))]
  if(correct.codes) {
    badz <- with(rr,which(father==1 & mother==1))
    rr[badz,"mother"] <- 2
    badz <- with(rr,which(father==2 & mother==2))
    rr[badz,"father"] <- 1
  }
  return(rr)
}


# internal function
validate.dir.for <- function(dir,elements,warn=F) {
  # in case the 'dir' input list object is not the standardised list form, convert
  # allows flexible use of list or regular directory specifications in plumbCNV functions
  if(is.null(dir)) { cat("directory empty\n"); return(NULL) }
  if(!is.list(dir)) {
    if(warn) { cat(elements[cc],"'dir' object wasn't a list\n")}
    dir <- as.list(dir); names(dir)[1:length(dir)] <- elements[1:length(dir)] 
  }
  for (cc in 1:length(elements)) {
    if(!elements[cc] %in% names(dir)) { 
      dir[[paste(elements[cc])]] <- "" ;
      if(warn) { stop(paste("dir$",elements[cc]," was empty.. set to current\n",sep="")) } 
    }
  }
  return(dir)
}


#where x is a vector of gene labels, where multiple hits
# are delimited by sep, then this will convert to the form
# gene1, gene2, ..., genen, + length(genen>1) others
compact.gene.list <- function(x,n=3,sep=";",others=TRUE) {
  XX <- strsplit(x,sep,fixed=T)
  #prv(XX)
  XX <- lapply(XX,function(x) { x[x %in% c(""," ")] <- "unnamed gene"; if(length(x)==0) { x <- "unnamed gene" };return(x) })
 # prv(XX)
  lens <- sapply(XX,length)
  sel <- which(lens>n)
  all <- sapply(XX,function(x) { sel <- FALSE; if(length(x)>0) { sel <- 1:(min(length(x),n)) }; paste(x[sel],collapse=sep) })
  extrz <- lens[sel]-n
  all[sel] <- paste(all[sel],"+",extrz,if(others) { "others" } else { "" } )
  return(all)  
}
