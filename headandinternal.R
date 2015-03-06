###NAMESPACE ADDITIONS###
#' @importFrom BiocInstaller  biocVersion
#' @importFrom utils capture.output download.file read.table write.table head tail data
#' @importFrom stats family pnorm pt qnorm rchisq rnorm runif
#' @importFrom reader cat.path reader shift.rownames
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline lines points rect text plot
#' @importFrom methods as callNextMethod is new prototype representation setAs setClass setGeneric setMethod setValidity
#' @importClassesFrom GenomicRanges GIntervalTree GRanges GenomicRanges Seqinfo GenomicRangesORmissing
#' @importFrom GenomicRanges "seqinfo" "seqinfo<-" "seqnames"  "seqnames<-" "seqlevels"  "seqlevels<-"
#' @importFrom GenomicRanges GRanges seqlengths Seqinfo GRangesList
#' @importMethodsFrom GenomicRanges "genome<-" "names<-" genome length names seqlevels start end
#' @importMethodsFrom GenomicRanges width strand mcols  "mcols<-" show findOverlaps subsetByOverlaps
#' @importMethodsFrom GenomicRanges seqinfo  "seqinfo<-" seqnames  "seqnames<-" seqlevels  "seqlevels<-"
#' @importMethodsFrom GenomicRanges length names  "names<-" "dimnames<-" "["  "[<-"  "[["  "[[<-"  "$"  "$<-" cbind rbind
#' @importClassesFrom IRanges DataFrame RangedData Rle Vector Annotated
#' @importFrom IRanges "%over%" DataFrame IRanges  Rle  RangedData subjectHits queryHits showAsCell
#' @importMethodsFrom IRanges "colnames<-" "rownames<-" "universe<-" Rle subjectHits queryHits showAsCell
#' @importMethodsFrom IRanges as.data.frame as.list as.matrix cbind rbind colnames elementLengths
#' @importMethodsFrom IRanges end findOverlaps subsetByOverlaps gsub intersect is.unsorted lapply
#' @importMethodsFrom IRanges levels mean na.exclude nrow ncol order paste as.list head tail aggregate
#' @importMethodsFrom IRanges ranges rownames runLength runValue sapply space  flank  reduce resize
#' @importMethodsFrom IRanges start universe unlist Rle width  "start<-"  "width<-"  "end<-" ranges "ranges<-"
#' @importClassesFrom "GenomicFeatures" TranscriptDb
#' @importFrom "GenomicFeatures" makeTranscriptDbFromUCSC  exonsBy  transcriptsBy
#' @importMethodsFrom "GenomicFeatures"  exonsBy  transcriptsBy  as.list
#' @importClassesFrom "rtracklayer"  ChainFile
#' @importMethodsFrom "rtracklayer"  liftOver  import.chain
#' @importMethodsFrom "genoset"  chr  chrIndices  chrInfo  chrNames  genome  isGenomeOrder  locData  toGenomeOrder  universe
#' @importFrom "genoset"  chr  chrIndices  chrInfo  chrNames  chrOrder  
#' @importFrom "genoset" "genome"  "isGenomeOrder"  "locData"  "toGenomeOrder"  "universe" 
#' @importClassesFrom "genoset" RangedDataOrGenomicRanges
#' @importFrom "biomaRt"  useMart  useDataset  getBM
#' @importClassesFrom "biomaRt"  Mart
#' @importFrom parallel  mclapply
#' @import Rcpp BiocGenerics NCmisc
###END NAMESPACE###


# importNoClassesFrom("GenomicRanges", GRanges)
# importNoClassesFrom("IRanges", Rle, RangedData)
# doNotimportFrom AnnotationDbi head tail ncol as.list colnames get exists sample 
# doNotimportFrom(BiocGenerics,strand, "strand<-", colnames, cbind, rbind, unlist, order, rownames, ncol, as.vector, paste, as.data.frame)


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("humarray version 1.0.0\n")
}

.onLoad <- function(libname, pkgname) {
  # library.dynam("humarray", pkgname, libname)
  #~/github/iChip/ImmunoChip_ChipInfo_New.RData
  options(chip.info="") # if you can access the file, you won't need to change this path
  options(ucsc="hg19") # depends on which analysis, need to set something though
  #data("iChipRegionsB36", "egSymb", "ImmunoChipB37", "hg18ToHg19","hg38ToHg19","hg19ToHg18","hg19ToHg38",
  #     package=pkgname, envir=parent.env(environment()))
  options(save.annot.in.current=1)  # 1 = TRUE, stores annotation in current folder to speed up subsequent lookups
}



#require(GenomicRanges); require(IRanges); require(reader); require(genoset)


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




#' Manage flexible input for the build parameter
#' 
#' The genome annotation version for internals in this package should always be
#' of the form 'hgXX', where XX can be 15,16,17,18,19,38. However most functions
#' allow flexible entry of this parameter as a build number, e.g, 36,37,38, or as
#' 'build36', 'b36', etc. This function sanitizes various forms of input to the 
#' correct format for internal operations. 
#' @param build the input to be sanitized. 
#' @param allow.multiple logical, whether to force a single value, or allow a vector
#' of build strings as input
#' @param show.valid logical, if TRUE, show a list of supported values.
#' @return build string in the correct 'hgXX' format.
#' @export
#' @examples
#' ucsc.sanitizer(36)
#' ucsc.sanitizer("build38")
#' ucsc.sanitizer("b37")
#' ucsc.sanitizer(show.valid=TRUE)
ucsc.sanitizer <- function(build,allow.multiple=FALSE,show.valid=FALSE) {
  build.alt <- c("hg15","hg20","hg17","hg18","hg19","hg38",17,18,19,20,35,36,37,38,
                 "build35","build36","build37","build38","b35","b36","b37","b38")
  build.new <- c("hg15","hg38",rep(c("hg17","hg18","hg19","hg38"),times=5))
  if(show.valid) { return(cbind(valid=build.alt,mapsTo=build.new)) }
  build <- build.new[match(tolower(build),build.alt)]
  if(is.null(build)) { build <- "hg19"; warning("build was NULL (see getOption('ucsc')), set to hg19") }
  if(any(is.na(build))) { 
    warning("Illegal build parameter '",build[1],"', defaulting to hg19") 
    build[is.na(build)] <- "hg19" 
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




# internal
# Remove trailing letter from non-unique rs-ids
#
# #@examples
# snp.ids <- rsnpid(25)
# snp.ids[1:2] <- paste0(snp.ids[1:2],"b")
# snp.ids[19:20] <- paste0(snp.ids[19:20],"c")
# snp.ids[6:7] <- paste0(snp.ids[6:7],"d")
# snp.ids[11:12] <- paste0(snp.ids[11:12],"a")
# snp.ids
# rmv.trail(snp.ids)
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
# Add trailing letter(s) to non-unique rs-ids
# #@examples
# snp.ids <- rsnpid(15)
# snp.ids
# add.trail(snp.ids)
# snp.ids <- snp.ids[sample(15,30,replace=TRUE)]
# snp.ids
# add.trail(snp.ids)
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

# internal
rsampid <- function(n,pref="ID0") { paste0(pref,pad.left(1:n,"0")) }


# internal
rsnpid <- function(n) { 
  id.len <- sample(c(3:8),n,replace=T,prob=c(0.01, 0.01, 0.01, 0.10, 0.50, 0.37))
  each.id <- function(l) { sapply(l,function(n) { paste(replicate(n,sample(1:9,1)),collapse="",sep="") }) }
  sufz <- each.id(id.len)
  ids <- paste0("rs",sufz)
  return(ids)
}


# internal
# Remove that pesky 'elementmetadata.' prefix from column names that have been converted from GRanges
emd.rmv <- function(X, rmv.genome=TRUE) {
  requireNamespace("GenomicRanges"); requireNamespace("IRanges")
  if(has.method("mcols",X, where=environment(emd.rmv))) {
    colnames(mcols(X)) <- gsub("elementMetadata.","",colnames(mcols(X)))
    ii <- which(colnames(mcols(X))=="genome")
    if(length(ii)>0) {
      gn <- mcols(X)[,ii[1]]
      if(length(unique(gn))==1) {
        mcols(X) <- mcols(X)[,-ii[1]]
      }
    }
  } else {
    if(has.method("colnames",X, where=environment(emd.rmv))) {
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



chrOrder2 <- function (chr.names) {
  if(!is.character(chr.names)) { warning("expecting character() type for chr.names argument") }
  simple.names = gsub("^chr", "", chr.names)
  name.is.numeric = grepl("^[0-9]+$", simple.names, perl = T)
  numeric.names = chr.names[name.is.numeric][order(as.numeric(simple.names[name.is.numeric]))]
  non.numeric.names = chr.names[!name.is.numeric][order(chr.names[!name.is.numeric])]
  all.names = c(numeric.names, non.numeric.names)
  return(all.names)
}


# internal# iFunctions
chrNames2 <- function(X) {
  #requireNamespace("GenomicRanges"); requireNamespace("IRanges")
  if(nrow(X)==0) { return(character(0)) }
  out <- as.character(unique(seqnames(X)))
  return(out)
}


# internal from genoset
TGORD <- function (ds, strict = TRUE) {
  if (strict == TRUE) {
    if (!isTRUE(all.equal(chrOrder2(chrNames2(ds)), chrNames2(ds)))) {
      ds = ds[chrOrder2(chrNames2(ds))]
    }
  }
  row.order = order(as.integer(space(ds)), start(ds))
  if (is.unsorted(row.order)) {
    return(ds[row.order, , drop = FALSE])
  }
  else {
    return(ds)
  }
}


# internal from genoset
TGOGR <- function (ds, strict = TRUE) {
  if (strict == TRUE) {
    if (!isTRUE(all.equal(chrOrder(seqlevels(ds)), seqlevels(ds)))) {
      seqlevels(ds) = chrOrder(seqlevels(ds))
    }
  }
  row.order = order(as.integer(seqnames(ds)), start(ds))
  if (is.unsorted(row.order)) {
    ds = ds[row.order, , drop = FALSE]
  }
  return(ds)
}



# internal # iFunctions
# version of toGenomeOrder() that is guaranteed to work for IRanges or GRanges
toGenomeOrder2 <- function(X,...) {
  requireNamespace("GenomicRanges"); requireNamespace("IRanges"); requireNamespace("genoset")
  if(is(X)[1] %in% c("GRanges","RangedData","ChipInfo")) {
    if(is(X)[1]=="RangedData") {
      return(TGORD(X))
    } else {
      return(TGOGR(X))
    }
  } else {
    typ <- is(X)[1]
    if(!typ %in% c("GRanges","RangedData","ChipInfo")) { warning("unsupported type '",typ,"' for toGenomeOrder(), failure likely") }
    alreadyThere <-("strand" %in% colnames(X))
    out <- TGOGR(as(X,"GRanges"),strict=T) #genoset::
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
  requireNamespace("GenomicRanges"); requireNamespace("IRanges")
  if(is(X)[1] %in% c("GRanges","ChipInfo")) {
    return(genoset::chrInfo(X))
  } else {
    typ <- is(X)[1]
    if(!typ %in% c("GRanges","RangedData","ChipInfo")) { warning("unacceptable type '",typ,"' for chrInfo2(), failure likely") }
    out <- chrInfo(as(X,"GRanges"))
    return(out)
  }
}

#internal # iFunctions
# version of chrIndices() that is guaranteed to work for IRanges or GRanges
chrIndices2 <- function(X,...) {
  requireNamespace("GenomicRanges"); requireNamespace("IRanges")
  if(is(X)[1] %in% c("GRanges","ChipInfo")) {
    return(genoset::chrIndices(X,...))
  } else {
    typ <- is(X)[1]
    if(!typ %in% c("GRanges","RangedData","ChipInfo")) { warning("unacceptable type '",typ,"' for chrIndices2(), failure likely") }
    out <- chrIndices(as(X,"GRanges"))
    return(out)
  }
}

#internal # iFunctions
# version of chr() that is guaranteed to work for IRanges or GRanges
chr2 <- function(X) {
  requireNamespace("GenomicRanges"); requireNamespace("IRanges")
  if(is(X)[1] %in% c("GRanges","ChipInfo")) {
    return(genoset::chr(X))
  } else {
    if(is(X)[1]=="RangedData") {
      return(space(X))
    } else {
      if(is.null(X)) { warning("X was NULL, expecting RangedData/GRanges"); return(NULL) }
      warning("chr2() function applies only to RangedData objects, attempting to pass ",is(X)[1]," to chr()")
      return(genoset::chr(X))
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
