

## USEFUL ONES !! HERE !!!


options(chip.info="/chiswick/data/ncooper/iChipData/ImmunoChip_ChipInfo.RData")
options(ucsc="hg18")
options(save.annot.in.current=1)


if(Sys.info()[["user"]]=="ncooper")
{
  source("~/github/iChip/firstscriptFunctions.R", echo=FALSE) # only needed for internal analyses
  source("~/github/plumbCNV/geneticsFunctions.R", echo=FALSE)
  source("~/github/iChip/ChipInfoClass.R", echo=FALSE)
  source('~/github/iChip/specificIFunctions.R', echo=FALSE)
}


require(snpStats)
require(reader)
require(genoset)

#


##file includes the generally useful functions: simple.date, out.of, randomize.missing
 # simple.date, out.of should now be in NCmisc


minna <- function(...) {
  min(...,na.rm=TRUE)
}
maxna <- function(...) {
  max(...,na.rm=TRUE)
}
meanna <- function(...) {
  mean(...,na.rm=TRUE)
}
medianna <- function(...) {
  median(...,na.rm=TRUE)
}
sdna <- function(...) {
  sd(...,na.rm=TRUE)
}
sumna <- function(...) {
  sum(...,na.rm=TRUE)
}
sortna <- function(...) {
  sort(..., na.last=TRUE)
}

# stats function convenience wrappers

#' Convert p-values to Z-scores
#' 
#' Simple conversion of two-tailed p-values to Z-scores
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
#' Simple conversion of Z-scores to two-tailed p-values
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
get.t1dbase.snps <- function(disease="T1D",build="hg19",show.codes=FALSE) {
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



#internal function to retrieve all.support from the global environment or else load from file
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


#internal function to propoerly sort chromosome labels as text
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
  if(is.null(build)) { build <- build(get.support()) }
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
  if(is.null(build)) { build <- build(get.support()) }
  build <- ucsc.sanitizer(build)
  char.lim <- 100
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.cyto(build=build,bioC=bioC,dir=dir)
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
  if(is.null(build)) { build <- build(get.support()) }
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
  if(is.null(build)) { build <- build(get.support()) }
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
  if(is.null(build)) { build <- build(get.support()) }
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
  if(is.null(build)) { build <- build(get.support()) }
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
  cyto <- get.cyto(build=build,bioC=TRUE,dir=dir)
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
  if(is.null(build)) { build <- build(get.support()) }
  must.use.package(c("biomaRt","genoset","gage"),T)
  build <- ucsc.sanitizer(build)
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


#' Retrieve 'n' closest SNP ids or positions near specified locus
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
  if(is.null(build)) { build <- build(get.support()) }
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



#internal
# standardize snp ids so they would always appear the same, all _,.;, etc replaced with _
# all names leading with a number preceeded with X. mainly 'make.names' standard R-conventions
clean.snp.ids <- function(snpid.list) {
  snpid.list <- make.names(snpid.list)
  snpid.list <- gsub(".","_",snpid.list,fixed=T)
  return(snpid.list)
}





#internal - to make valid GRanges from chr,pos or start/end
# build only affects 'universe' slot of RangedData, or 'genome' slot of GRanges

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
make.granges <- function(chr,pos=NULL,start=NULL,end=NULL,row.names=NULL,build="hg18") {
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
      extra <- data.frame(Chr=newchr,Pos=noopos)
      rownames(extra) <- ln
    } #else { cat("length is already the same\n") }
    Ind <- match(ranged.gr.37[["XMYINDEXX"]],orn)
    out <- data.frame(Chr=ranged.gr.37[["XMYCHRXX"]],Pos=start(ranged.gr.37),ind=Ind)
    rownames(out) <- RN
    if(length(orn)>length(RN)) {
      out <- out[,-3]
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
SnpMatrix.to.data.frame <- function(SnpMat) {
  #if(is(SnpMat)[1]=="snp.matrix") { SnpMat <- as(SnpMat,"SnpMatrix") }
  if(is(SnpMat)[1]!="SnpMatrix") { stop("SnpMat must be a SnpMatrix object") }
  cov.data <- as.data.frame(SnpMat)
  for(jj in 1:ncol(cov.data)) { 
    nuxt <- as.numeric(cov.data[,jj])-1
    nuxt[nuxt<0] <- NA
    cov.data[,jj] <- nuxt
    # assign(colnames(cov.data)[jj], nuxt)
  }
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


## for a snp dataset, determine whether reference allele is the major or minor allele 
# Where 2 copies coded as highest value = reference, e.g, if AA=0, AB=1, BB=2, then B is reference.
# Combines this with frequencies of the alleles to evaluate whether 'BB' is major or minor
# takes SnpMatrix or data.frame coded 0,1,2,NA as input
#'
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
#'
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
rSnpMatrix <- function(nsnp=5,nsamp=10,call.rate=.95,A.freq.fun=runif) {
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
#' @seealso get.recombination.map, get.nearby.snp.lists
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
#' @seealso snps.in.range, get.recombination.map, recwindow, conv.37.36, conv.36.37
#' @examples
#' result <- get.nearby.snp.lists("rs900569")
#' get.nearby.snp.lists("rs900569",cM=0.2,excl.snps=result[[1]]) # get SNPs within 0.1-0.2cM
#' # note that the same query can return a different set with build 36 versus 37
#' get.nearby.snp.lists(c("rs689","rs4909944"),cM=0.001,build=36,name.by.bands=FALSE)
#' get.nearby.snp.lists(c("rs689","rs4909944"),cM=0.001,build=37,name.by.bands=FALSE) 
get.nearby.snp.lists <- function(snpid.list,cM=0.1,bp.ext=0,build=NULL,excl.snps=NULL,name.by.bands=TRUE) {
  #if(!exists("all.support")) { print(load("all.support.RData")) }
  all.support <- get.support()
  if(is.null(build)) { build <- build(all.support) }
  build <- ucsc.sanitizer(build)
  snpic.list <- rs.to.id(snpid.list)
  cyto <- get.cyto(dir=getwd()); cyto[["gene"]] <- rownames(cyto)
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



