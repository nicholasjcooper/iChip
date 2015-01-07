###NAMESPACE ADDITIONS###
# Depends: R (>= 2.14), grDevices, graphics, stats, utils, reader, NCmisc
# Imports: genoset, IRanges, GenomicRanges, snpStats, Rcpp, GenomicFeatures, rtracklayer, biomaRt, methods, parallel
# Suggests: snpStats
# importClassesFrom("GenomicFeatures",TranscriptDb)
# importFrom("GenomicFeatures",makeTranscriptDbFromUCSC, exonsBy, transcriptsBy)
# importMethodsFrom("GenomicFeatures", exonsBy, transcriptsBy, as.list)
# importMethodsFrom("snpStats", effect.sign)
# importFrom("snpStats", row.summary, col.summary, read.pedfile, snp.imputation, impute.snps, single.snp.tests)
# importClassesFrom("rtracklayer", ChainFile)
# importClassesFrom("snpStats", SnpMatrix, XSnpMatrix, SingleSnpTests, SingleSnpTestsScore)
# importMethodsFrom("rtracklayer", liftOver, import.chain)
# importFrom("biomaRt", useMart, useDataset, getBM)
# importClassesFrom("biomaRt", Mart)
# importFrom(parallel, mclapply)
# import(grDevices, graphics, stats, utils, reader, genoset, NCmisc, Rcpp, methods, IRanges, GenomicRanges, genoset)
###END NAMESPACE###


# importNoClassesFrom("GenomicRanges", GRanges)
# importNoClassesFrom("IRanges", Rle)

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("humarray >= 1.0 is awesome.\n")
}

.onLoad <- function(libname, pkgname) {
  # library.dynam("humarray", pkgname, libname)
  options(chip.info="~/github/iChip/ImmunoChip_ChipInfo_New.RData") # if you can access the file, you won't need to change this path
  options(ucsc="hg19") # depends on which analysis, need to set something though
  options(save.annot.in.current=1)  # 1 = TRUE, stores annotation in current folder to speed up subsequent lookups
}




########################
## internal functions ##
########################

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


# internal
l10 <- function(x) { O <- log10(x); O[!is.finite(O)] <- NA; return(O) }
# internal
Cor <- function(...) { cor(...,use="pairwise.complete") }
# internal
pt2 <- function(q, df, log.p=FALSE) {  2*pt(-abs(q), df, log.p=log.p) }



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


#internal
# standardize snp ids so they would always appear the same, all _,.;, etc replaced with _
# all names leading with a number preceeded with X. mainly 'make.names' standard R-conventions
clean.snp.ids <- function(snpid.list) {
  snpid.list <- make.names(snpid.list)
  snpid.list <- gsub(".","_",snpid.list,fixed=T)
  return(snpid.list)
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



# internal
# snpStats imputation only works if there are correlated SNPs with non-missing values
# that can be used to interpolate missing SNPs. If any correlated SNPs are missing
# 'impute.missing' will leave these blank. This function mops up the remainder
# by randomly inserting values consistent with the minor allele frequency of each SNP
# (Nick)
randomize.missing2 <- function(X,verbose=FALSE) {
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



# iFunctions
# internal
# allows an sapply style function to only work on valid values
clean.fn <- function(x,fail=NA,fn=function(x) { x }) {
  if(!is.null(x)) { 
    x <- x[!is.na(x)]; x <- x[(x!="")]; return(fn(x)) 
  } else {  return(fail) } 
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
  if(nrow(X)==0) { return(character(0)) }
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



################# end internals #################
