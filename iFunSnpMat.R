
#########################
## SnpMatrix Functions ##
#########################


#' Calculate Lambda inflation factors for SNP dataset
#' 
#' This function calculates SNP-wise or overall Lambda and Lambda_1000 statistics for inflation due
#' to population structure. It works on a SnpMatrix object or dataframe coded 0,1,2,NA (autodetects which).
#' @param X SnpMatrix or data.frame coded 0,1,2,NA containing genotypes
#' @param checks logical, whether to perform some sanity checks that will slightly slown down performance
#' @return output character vector, containing any of the following "all", "lambda","l1000", "both";
#' if snp.wise is false, the scalar value(s) specified by 'output', or otherwise a matrix
#' of parameters for each SNP, optionally limited to just the lambda column(s) by the value of 'output'.
#' @param snp.wise logical, if TRUE, return lambda statistics separately for each SNP (rather than an overall)
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
  if(is(X)[1] %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { snpmat <- T } else {
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
  if(is(SnpMat)[1]=="aSnpMatrix") { SnpMat <- SnpMat@.Data }
  if(is(SnpMat)[1]!="SnpMatrix") { stop("SnpMat must be a SnpMatrix/aSnpMatrix object") }
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
  if(is(X)[1] %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { snpmat <- snp.mat <- T} else {
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
  if(is(X)[1] %in% c("SnpMatrix","XSnpMatrix","aSnpMatrix","aXSnpMatrix")) { 
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



#' Multicore randomised replacement of missing genotypes
#' 
#' snpStats imputation only works if there are correlated SNPs with non-missing values
#' that can be used to interpolate missing SNPs. If any correlated SNPs are missing
#' 'impute.missing' will leave these blank. This function mops up the remainder
#' by randomly inserting values consistent with the minor allele frequency of each SNP.
#' This can be run using multiple cores to speed up progress for large matrices.
#' @param X a SnpMatrix object, presumably with missing genotypes
#' @param verbose logical, whether to report on progress, or suppress extraneous output
#' @param n.cores integer, maximum number of processing cores to use (limited by system)
#' @return returns a SnpMatrix object with all missing values replaced
#' @export
#' @examples
#' test <- rSnpMatrix(nsamp=1000,call.rate=.95) # create a random SnpMatrix
#' prv(test) 
#' col.summary(test) # shows some missing values
#' test2 <- randomize.missing(test,verbose=TRUE)
#' col.summary(test2) # shows none missing
randomize.missing <- function(X,verbose=FALSE,n.cores=1) {
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
  typ <- is(X)[1]
  if(!typ %in% c("SnpMatrix","aSnpMatrix")) { stop("X must be a SnpMatrix/aSnpMatrix object") }
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
  if(random) { N <- randomize.missing2(N,verbose=verbose) } # replace any remaining missing values
  if(data.frame) {
    return(as.data.frame(N))
  }
  else {
    #print(Dim(N))
    #print(Dim(X@snps))
    N <- round(N)
    N[is.na(N)] <- -1 # replace missing with -1 [will become zero in next line]
    if(is(N)[1]=="list" | is(N)[1]=="data.frame") { N <- as.matrix(N) } # sometimes it becomes a list for some reason???
    out.m <- (new("SnpMatrix", data = (as.raw(N+1)),
                      nrow = nrow(N), ncol = ncol(N), dimnames = dimnames(N)))
    if(typ=="SnpMatrix") { return(out.m) } else {
    	X@.Data <- out.m
    	return(X)
    }
  }
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




################## end snpmatrix ##########################
