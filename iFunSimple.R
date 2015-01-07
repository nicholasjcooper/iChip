

######################
## Simple Functions ##
######################

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





#' Convert a chr:pos1-pos2 vector to a matrix
#' 
#' Takes standard text positions, such as what you might see on the UCSC genome browser, such as 
#' chr1:10,000,234-11,000,567 for a range, or chrX:234,432 for a SNP, and converts to 
#' with cols: chr, start, end.
#' @param text character vector, format like chr:pos1-pos2
#' @export
#' @return a matrix of the same length as 'ranges' with columns chr, start and end, and
#' rownames will be the same as the original text vector.
#' @seealso ranges.to.txt
#' @examples
#' txt <- ranges.to.txt(rranges())
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



#' Make a compact version of gene annotation
#'
#' When adding gene annotation to genomic ranges, sometimes
#' there are many genes associated with a single feature, so
#' that compiling a table becomes awkward, if some rows contain
#' hundreds of genes. This function takes a character vector
#' of gene lists delimited by some separator and provides
#' a compact representation of the gene labels
#' 
#' @param x is a character vector of gene label listings, where multiple hits
#' are delimited by 'sep'
#' @param n number of genes to list before abbreviating
#' @param sep character, separator used to delimit genes in elements of x
#' @param others logical, TRUE to abbreviate with '+ # others' or FALSE to
#' append just the number of genes not listed.
#' @return a character vector with the form:
#' gene-1, gene-2, ..., gene-n, + length(gene-n) - n [others]
#' @export
#' @examples
#' my.genes <- c("ERAP1","HLA-C;CTLA4;IFIH","INS;MICA","AGAP1;APOE;DRDB1;FUT2;HCP5;BDNF;COMT")
#' compact.gene.list(my.genes)
#' compact.gene.list(my.genes,n=2,others=TRUE)
compact.gene.list <- function(x,n=3,sep=";",others=FALSE) {
  XX <- strsplit(x,sep,fixed=T)
  #prv(XX)
  XX <- lapply(XX,function(x) { x[x %in% c(""," ")] <- "unnamed gene"; if(length(x)==0) { x <- "unnamed gene" };return(x) })
 # prv(XX)
  lens <- sapply(XX,length)
  sel <- which(lens>n)
  all <- sapply(XX,function(x) { sel <- FALSE; if(length(x)>0) { sel <- 1:(min(length(x),n)) }; paste(x[sel],collapse=sep) })
  extrz <- lens[sel]-n
  if(others) { oth <- rep("others",length(extrz)); oth[extrz==1] <- "other" }
  all[sel] <- paste(all[sel],"+",extrz,if(others) { oth } else { "" } )
  return(all)  
}




################## end simple ##########################

