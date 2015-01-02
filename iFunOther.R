
#####################
## Other Functions ##
#####################


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


################## end other ##########################

