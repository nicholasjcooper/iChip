
if(Sys.info()[["user"]]=="ncooper")
{
  source("~/github/iChip/firstscriptFunctions.R") # only needed for internal analyses
  source("~/github/plumbCNV/geneticsFunctions.R")
  source("~/github/iChip/ChipInfoClass.R")
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



## USEFUL ONES !! HERE !!!


options(chip.info="/chiswick/data/ncooper/iChipData/ImmunoChip_ChipInfo.RData")
options(ucsc="hg18")
options(save.annot.in.current=1)


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
    if(is(all.support)[1]=="ChipInfo") { 
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
#' @return A character vector of SNP chip-ids, where the input was rs-ids, chip-ids or a mixture, any text
#' other than this will result in NA values being returned in the character vector output.
#' @export
#' @seealso id.to.rs, GENE.to.ENS, ENS.to.GENE
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' rs.to.id(c("rs689","rs9467354","rs61733845"))  # middle one has no chip id
rs.to.id <- function(rs.ids) {
  rs.ids <- clean.snp.ids(rs.ids)
  all.support <- get.support()
  if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found") }  ## load object: all.support [snp support for whole chip]
  idvec <- rownames(all.support)[match(rs.ids,mcols(all.support)$rs.id)]
  idvec2 <- rownames(all.support)[match(rs.ids,rownames(all.support))]
  idvec[is.na(idvec)] <- idvec2[is.na(idvec)]
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
Band.pos <- function(chr,pos=NA,start=NA,end=NA,ranges=NULL,build=NULL,dir=NULL,bioC=FALSE,one.to.one=TRUE) {
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
  if(max(Dim(pos))!=length(chr)) { stop("chr and pos must be of the same length") }
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
#' output that is sorted by genome order, regardless of the original order.
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
    ranges <- make.granges(chr=chr,pos=pos,row.names=ids,...)
  }
  # return(ranges)
  if(is(ranges)[1]=="RangedData") { ranges <- as(ranges, "GRanges") }
  if(is(ranges)[1] %in% c("RangedData","GRanges")) {
    wd <- width(ranges)
    if(all(wd==1)) { SNPs <- TRUE } else { SNPs <- FALSE }
    if(is.null(rownames(ranges))) { rownames(ranges) <- paste(1:nrow(ranges)) }    
    mcols(ranges)[["XMYINDEXX"]] <- orn <- rownames(ranges)
    mcols(ranges)[["XMYCHRXX"]] <- ocr <- chr2(ranges)
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
      noopos <- start(ranged[match(ln,ranged$XMYINDEXX),])
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
    }
    Ind <- match(ranged.gr.37[["XMYINDEXX"]],orn)
    out <- data.frame(Chr=ranged.gr.37[["XMYCHRXX"]],Pos=start(ranged.gr.37),ind=Ind)
    rownames(out) <- RN
    if(length(orn)>length(RN)) {
      out <- out[,-3]
      out <- rbind(out,extra)
      out <- out[orn,]
    }
    #return(out) 
  } else { warning("missing key columns for chr, snp-name")  }
  #print(outType)
  #return(out)
  #prv(out)
  ranged.rd <- toGenomeOrder2(data.frame.to.ranged(out))
  ranged.gr.37 <- as(ranged.rd,"GRanges")
  if(found.xy) {
    xy.ind <- match(xy.id,rownames(ranged.gr.37))
    lmis <- length(which(is.na(xy.ind)))
    if(lmis>0) { warning('liftOver function removed ",lmis," chrX/chrY ranges'); xy.ind <- narm(xy.ind) }
    if(!"XY" %in% seqlevels(ranged.gr.37)) { seqlevels(ranged.gr.37) <- c(seqlevels(ranged.gr.37),"XY") }
    seqnames(ranged.gr.37)[xy.ind] <- "XY"
  }
  if(outType=="GRanges") { 
    #return(ranged.gr.37)
    mind <- as.numeric(mcols(ranged.gr.37)[["ind"]])
    #prv(ranged.gr.37,mind)
    ranged.gr.37 <- ranged.gr.37[order(mind),]
    cn37 <- colnames(mcols(ranged.gr.37))
    if("ind" %in% cn37) { mcols(ranged.gr.37) <- mcols(ranged.gr.37)[,-which(cn37 %in% "ind")] }
    return(ranged.gr.37)
  } else {
    if(outType=="RangedData") {
      #if("ind" %in% colnames(ranged.gr.37)) { ranged.gr.37 <- ranged.gr.37[,-which(colnames(ranged.gr.37) %in% "ind")] }
      return(toGenomeOrder2(as(ranged.gr.37,"RangedData")))
    } else {
      mind <- as.numeric(mcols(ranged.gr.37)[["ind"]])
      #return(ranged.gr.37)
      out <- ranged.to.data.frame(ranged.gr.37,include.cols=FALSE,use.names=TRUE)
      out <- out[order(mind),]
      # if("ind" %in% colnames(out)) { out <- out[,-which(colnames(out) %in% "ind")] }
      #return(out)
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


#### HERE! #####
# remaining functions to test and document #
# randomize.missing
# impute.missing
# recwindow
# meta.me
# get.nearby.snp.lists
# plink.allele.coding

###

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




## HERE ##
# chris' function to get a centimorgan window from intervals
# vector input
# internal
## Given a region, add window cM either side and return.
## Everything is HapMap v2, release 22, build 36.
#' Extend an interval or SNP by distance in centimorgans (recombination distance)
#' 
#' It is straightforward to extend a genomic interval or position by a number of basepairs, or
#' a percentage, but extending by recombination units of centimorgans is more involved, requiring
#' annotation lookup. This function streamlines this process
recwindow <- function(ranged=NULL,chr=NA,start=NA,end=start,window=0.1,bp.ext=0, rec.map=NULL) {
  if(!is.numeric(bp.ext)) { warning("bp.ext must be numeric, setting to zero"); bp.ext <- 0 }
  if(!is.numeric(window)) { warning("window must be numeric, setting to 0.1 centimorgans"); window <- 0.1 }
  if(!is.na(chr)) { if(any(!paste(chr) %in% paste(1:22))) { 
    stop("this function only works for autosomes 1-22 [e.g, no X,Y or formatting like 'chr2', etc]") } }
  if(is(ranged)[1] %in% c("RangedData","GRanges")) { 
    if(is(ranged)[1]=="GRanges") { ranged <- as(ranged,"RangedData") }
    ranged <- toGenomeOrder2(ranged,strict=T)
    ss <- start(ranged); ee <- end(ranged); cc <- chr2(ranged)
   out <- recwindow(chr=cc,start=ss,end=ee,window=window,bp.ext=bp.ext)
    outData <- RangedData(ranges=IRanges(start=out[,1],end=out[,2],names=rownames(ranged)),space=cc)
    outData <- toGenomeOrder(outData,strict=TRUE)
    for (zz in 1:ncol(ranged)) { outData[[colnames(ranged)[zz]]] <- ranged[[colnames(ranged)[zz]]]  }
    if(is(ranged)[1]=="GRanges") { outData <- as(toGenomeOrder2(outData,strict=T),"GRanges") }
    return(outData)
  } else {
    if(all(!is.na(chr)) & all(!is.na(start)) & all(!is.na(end))) {
      if(length(chr)==length(start) & length(start)==length(end)) {
        if(length(chr)>1) {
          # run for a vector
          out <- matrix(ncol=2,nrow=length(chr)); colnames(out) <- c("start","end")
          for (dd in 1:length(chr)) {
            out[dd,] <- recwindow(chr=chr[dd],start=start[dd],end=end[dd],window=window,bp.ext=bp.ext)
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
  cat("n hapmap snps in window =",nrow(kk),"\n")
  from <- min(kk[,1])
  to <- max(kk[,1])
  ##
  cat("new window size is\nleft: ",(start-from+bp.ext)/1000,"kb\tright: ",(to-end+bp.ext)/1000,"kb\ttotal: ",(to-from+(2*bp.ext))/1000,"kb\n",sep="")
  if(bp.ext>0) { cat("in addition to cM distance, window was extended by",bp.ext,"base pairs on either side\n")} 
  from <- max(c(0,(from-bp.ext)))
  to <- min(c((to+bp.ext),get.chr.lens()[chr][1]),na.rm=T)
  return(c(from,to))
}



# to get meta analysis parameters from table containing case-control and family data beta, se values
meta.me <- function(X,OR1="OR_CC",OR2="OR_Fam",SE1="SE_CC",SE2="SE_Fam",N1=18856,N2=7638) {
  if(!is(X)[1]=="data.frame") { stop("X must be a data.frame") }
  cnx <- colnames(X)
  if(!all(c(OR1,OR2,SE1,SE2) %in% cnx)) { stop("X must contain column names specified by OR1,OR2,SE1,SE2") } 
  ok <- FALSE
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
  if(N.coln) {
    famN <- X[,N2]
    ccN <- X[,N1]
  } else {
    famN <- N2 # 3819*2  #3509*2   #  3819*2   #  10796
    ccN <- N1 # 6683+12173 # including CBR, or for UVA analyses use instead: 9416+6670
  }
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
  rownames(out) <- c("OR.meta","beta.meta","se.meta","z.meta","p.meta") #rownames(X)
  return(out)
}



## for a list of snp-ids from iChip, obtain the nearby SNP-lists within 0.1cm, etc
# name.by.bands labels each sublist by the band name, but faster not to do this
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
  ## recwindow uses build36 only, so convert back afterwards
  nxt.window <- lapply(snps.locs36, function(X,...) { recwindow(st=X,...) },chr=next.chr,window=cM,bp.ext=bp.ext)
  if(build=="hg18") {
    st.window <- sapply(nxt.window, "[",1)
    en.window <- sapply(nxt.window, "[",2)
  } else {
    st.window <- conv.36.37(chr=next.chr,pos=sapply(nxt.window, "[",1))[,"start"]
    en.window <- conv.36.37(chr=next.chr,pos=sapply(nxt.window, "[",2))[,"start"]
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




#' @examples
#' hwe.fn <- "snpdataout.hwe"
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
  hh$RAF[!hh$flip] <- hh$Maj
  return(hh)
}







## SPECIFIC TO ME ##



## imputation subscript ## 
# Nick
# used by conditional analysis and indistinguishable analyses scripts, does the imputation,
# saving imputed results as we go so don't need to keep recalculating the same SNPs
# will automatically proceed in the most efficient way possible
# for my ichip analysis, prv.file="allImputed.RData"
#' @param myData SnpMatrix object with missing data that you wish to impute
#' @param imp.file name to give the file where the result will be saved
#' @param smp.filt a list of sample (indexes or labels) to include
#' @param snp.filt a list of SNPs (indexes or labels) to include
#' @param prv.file character, if applicable the name of a previous imputation file, which
#' means you don't need to impute the same SNPs again if you've imputed them before
#' @param by integer, a parameter that is passed to impute.missing() to determine
#' what sized subsets to impute by (smaller subsets makes imputation faster)
#' @param ret.object logical, if TRUE then return the a SnpMatrix with missing values
#' filled in, if false, return the file name this object was written to (same as imp.file)
#' @return See ret.object. If ret.object=TRUE, then returns the same SnpMatrix with missing
#' values replaced.
#' @examples
#' test.mat <- rSnpMatrix(nsnp=5,nsamp=10)
#' ichip.imputation(test.mat, "myImp.RData")
#' ichip.imputation(test.mat, imp.file="myImp.RData",prv.file="myImp.RData")
#' test.mat2 <- rSnpMatrix(nsnp=500,nsamp=10)
#' ichip.imputation(test.mat2, imp.file="myImp.RData",prv.file="myImp.RData")
#' test.mat3 <- cbind(test.mat2[,1:200],rSnpMatrix(nsnp=500,nsamp=10))
#' ichip.imputation(test.mat3, imp.file="myImp.RData",prv.file="myImp.RData")
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



# make sure support file naming and symbol convention match internal conventions, see clean.snp.ids
#internal
clean.snp.support <- function(X) {
  if(is(X)[1]!="ChipInfo") { stop("this function only allows X to be a ChipInfo object") }
  mcols(X)[["dbSNP"]][is.na(X$dbSNP)] <- rownames(X)[is.na(mcols(X)[["dbSNP"]])]
  rownames(X) <- clean.snp.ids(rownames(X))
  mcols(X)[["dbSNP"]] <- clean.snp.ids(mcols(X)[["dbSNP"]])
  rownames(X) <- clean.snp.ids(rownames(X))
  return(X)
}


## FIX UVA XCHR name screw up ##
#internal
rmv.uva.dup.rows <- function(table1) {
  print(Dim(table1))
  bad.tab <- table1[which(is.na(id.to.rs(rownames(table1)))),]
  nr <- substr(rownames(bad.tab),1,nchar(rownames(bad.tab))-1)
  table1 <- table1[which(!rownames(table1) %in% paste(nr,"1",sep="")),] # remove duplicate rows
  print(Dim(table1))
  return(table1)
}


#internal
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


#internal
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




#internal OLD
# remove leading X from variable names (e.g, if original name started with a number and changed by make.names)
remove.X <- function(str,char="X") {
  bdz <- substr(str,1,1)
  str[bdz==char] <- substr(str,2,100000)[bdz==char]
  return(str)
}

#internal OLD
# add X to the start of any string which has a digit as the first character
add.x <- function(str) {
  bdz <- substr(str,1,1)
  numy <-(paste(bdz) %in% paste(c(0:9)))
  str[numy] <- paste("X",str[numy],sep="") 
  return(str)
}


# trivial and specific ... internal for some function somewhere?
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



##INTERNAL
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
### specific to ichip paper analysis
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


# specific to me
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
  p.m.reg <- names(p.meta[names(p.meta) %in% rs.to.id(tt$names[in.region])])
  p.m.out <- names(p.meta[names(p.meta) %in% rs.to.id(tt$names[out.region])])
  p.m.nov <- names(p.meta[names(p.meta) %in% rs.to.id(tt$names[potential.novel.region])])
  
  if(surround) {
    t1dsnps <- calibrate.cond.bonf(topsnplist,cm.window=0.2,bp.ext=0,build=37,qclist="snpsExcluded.txt",ret.snps=T)
  }
  t1region <- tt$gene %in% rownames(bonfs.filt)
  t1region.novel <- (tt$gene %in% rownames(bonfs.filt)) & (rs.to.id(tt$names) %in% rs.to.id(p.m.nov))
  t1regsnps <- rs.to.id(tt$names[t1region])
  t1regsnps.novel <- rs.to.id(tt$names[t1region.novel])
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
  p.meta.reg <- p.meta[names(p.meta) %in% rs.to.id(tt$names[in.region])]
  #p.meta.out <- p.meta[names(p.meta) %in% rs.to.id(tt$names[out.region])]
  p.meta.out <- p.meta[names(p.meta) %in% rs.to.id(p.m.out.not.t1)]
  do.the.meta.ones(p.meta.reg,"inside.dense")
  do.the.meta.ones(p.meta.out,"outside.dense")
  # QQ for t1d versus non dense regions #
  p.meta.reg <- p.meta[names(p.meta) %in% rs.to.id(t1regsnps)]
  p.meta.out <- p.meta[!names(p.meta) %in% rs.to.id(t1regsnps)]
  do.the.meta.ones(p.meta.reg,"inside.t1d")
  do.the.meta.ones(p.meta.out,"outside.t1d")
  
  return(c(lam1000.1,lam1000.2,lam1000.3,lam1000.4,lam1000.5))
}


#send cw: bonferonnis, numbers of snps


# function specific to getMetatable.R script
#  gets a list of equivalent SNPs to the current snp
get.equivs <- function(id,iden.list,na.fail=T) {
  id <- id.to.rs(id) 
  do.one <- function(id,iden.list,na.fail) {
    if(id %in% id.to.rs(unlist(iden.list))) {
      lnum <- which(sapply(iden.list,function(X) { id %in% id.to.rs(X) }))
      idens <- id.to.rs(iden.list[[lnum]])
      return(idens[-which(idens %in% id)])
    } else {
      if(na.fail) { return(NA) } else { return(id) }
    }
  }
  return(sapply(id,do.one,iden.list=iden.list,na.fail=na.fail))
}


# specific
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

#specific
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
      uva.lookup <- match(rs.to.id(excl.uva),rs.to.id(rownames(SNPQC)))
      print(length(which(is.na(uva.lookup))))
      excl.rules.uva <- mini.snp.qc(SNPQC[uva.lookup,],MAF=.005,CR=.95,HWE=4.89)
      num.excl.rules.uva <- length(excl.rules.uva)
      cat(out.of(length(excl.uva)-num.excl.rules.uva,length(excl.uva)),"UVA failing snps fail for other reasons\n")
    }
  }
  return(excl.rules)
}
####################



# what does this do?
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
      cat("Chr",chr2(ranged[ccc])[1],":\n")
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




data.frame.to.ranged2 <- function(dat,ids=NULL,start="start",end="end",width=NULL,
                                  chr="chr",exclude=NULL,build="hg18") 
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
  outData <- RangedData(ranges=IRanges(start=st1,end=en1,names=id),space=ch1,universe=build[1])
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
  kk2 <- which(!(rs.to.id(tt[,5]) %in% unique(rs.to.id(c(qc.excluded.snps,qc.cloud.fail)))))
  tt <- tt[kk2,]
  tt <- tt[,-3:-4]
  colnames(tt)[1:2] <- c("Chr","Pos")
  aabb <- AB(tt$names)
  forolly <- cbind(tt,aabb)
  colnames(forolly)[11:12] <- c("allele.A","allele.B")
  colnames(forolly)[10] <- "band"
  colnames(forolly)[3] <- "rsid"
  rownames(forolly) <- rs.to.id(forolly[,3])
  
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





# returns list of snps from old table, not in new
# prints summary of new vs old pvalues and odds ratios for those snps
# e.g, examine.no.longer.t1.snps(table1snpsfinaljan30,bonf.snps,table1a,identicals,T,T)
examine.no.longer.t1.snps <- function(table1snpsfinaljan30,bonf.snps,table1a,identicals,do.OR=FALSE,do.p=TRUE) {
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  excl <- reader(docs[2])
  uva.table <- reader(docs[1])
  why <- table1snpsfinaljan30[!id.to.rs(table1snpsfinaljan30) %in% id.to.rs(bonf.snps)]
  ct.fn <- "conditionalTests.csv"
  ct <- reader(ct.fn,stringsAsFactors=F)
  conditionals <- ct$TABLE1[ct$COND!=0]
  if(length(why)>0 & exists("conditionals")) {  why <- why[!id.to.rs(why) %in% id.to.rs(conditionals)] }
  if(length(why)>0) {   why <- why[!id.to.rs(why) %in% id.to.rs(identicals)] }
  if(length(why)==0) { return("none found a mystery")}  
  
  if(do.p) {
    # compare pvalues new vs old for SNPs in the old table 1, no longer in the new table 1
    ww <- cbind(table1a[rs.to.id(why),"P_Meta"],
                uva.table[rs.to.id(why),"P_Meta"],Chr(why),Pos(why))[order(Chr(why)),]
    xx <- cbind(table1a[rs.to.id(why),"P_CC"],
                uva.table[rs.to.id(why),"P_CC"])[order(Chr(why)),]
    zz<-cbind(xx,ww)
    colnames(zz) <- c("New_CC","Old_CC","New_Meta","Old_Meta","Chr","Pos")
    zz <- zz[,c(5,6,1,2,3,4)]
    rownames(zz) <- why[order(Chr(why))]
    print(zz,digits=4)
  }
  if(do.OR) {
    # compare odds ratios new vs old for SNPs in the old table 1, no longer in the new table 1
    ww2 <- cbind(table1a[rs.to.id(why),"OR_Fam"],
                 uva.table[rs.to.id(why),"OR_Fam"],Chr(why),Pos(why))[order(Chr(why)),]
    xx2 <- cbind(table1a[rs.to.id(why),"OR_CC"],
                 uva.table[rs.to.id(why),"OR_CC"])[order(Chr(why)),]
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




##### SPECIFIC TO ME ######




# calibrate the bonferroni threshold for conditional analyses
calibrate.cond.bonf <- function(snplist,cm.window=0.1,bp.ext=0,build=37,qclist="snpsExcluded.txt",ret.snps=FALSE) {
  all.snps.tested <- NULL
  qc.excluded.snps <- reader(qclist)
  qc.excluded.snps <- qc.excluded.snps[!rs.to.id(qc.excluded.snps) %in% rs.to.id(snplist)]
  
  for (cc in 1:22) {
    cat("chr",cc,"\n")
    snpid.list <- snplist[Chr(snplist) %in% cc]
    if(length(snpid.list)<1) { next }
    grp.labs <- get.nearby.snp.lists(snpid.list,cM=cm.window,bp.ext=bp.ext,build=37,excl.snps=qc.excluded.snps,name.by.bands=FALSE)
    all.snps.tested <- c(all.snps.tested,unlist(grp.labs))  
  }
  
  ast <- all.snps.tested[!rs.to.id(all.snps.tested) %in% rs.to.id(c(qc.excluded.snps,qc.cloud.fail))]
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




if(F) {
  all.support <- reader("/chiswick/data/ncooper/iChipData/ImmunoChip_ChipInfo.RData")
  
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
  p.m.reg <- p.meta[names(p.meta) %in% rs.to.id(tt$names[in.region])]
  p.m.out <- p.meta[names(p.meta) %in% rs.to.id(tt$names[out.region])]
  do.the.meta.ones(p.m.reg)
  do.the.meta.ones(p.m.out)
}





# wut?
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


