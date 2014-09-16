if(getwd()!= "/home/ncooper"){
  require(genoset)
}
# examples
# 
# snp.info <- ChipInfo(chr=all.support[,"Chr"],pos=all.support[,"Pos"],ids=rownames(all.support),chip="ImmunoChip",rs.id=all.support[,"dbSNP"],
#                      A1=all.support[,"A1"], A2=all.support[,"A2"])
# 
# gr.snp.info <- with(all.support,make.granges(chr=Chr,pos=Pos,row.names=rownames(all.support)))
# 
# snp.info <- ChipInfo(GRanges=gr.snp.info,chip="ImmunoChip",rs.id=all.support[,"dbSNP"],
#                      A1=all.support[,"A1"], A2=all.support[,"A2"])
# 
# snp.info[chr(snp.info)=="MT",] # look at the mitochondrial SNP
# QCcode(snp.info)[chr(snp.info)=="MT"] <- 1 # exclude it by changing the fail code
# snp.info[["MT"]] # revisit and see it now registers as 'fail'
# QCfail(snp.info)
# xx <- conv.36.37(chr=all.support[,"Chr"],pos=all.support[,"Pos"],ids=rownames(all.support))
# build(si36[["XY"]])
# si37 <- convTo37(snp.info)
# si36 <- convTo36(si37)
# snp.info <- ChipInfo(chr=all.support[,"Chr"],pos=all.support[,"Pos37"],ids=rownames(all.support),chip="ImmunoChip",rs.id=all.support[,"dbSNP"],
#                     A1=all.support[,"A1"], A2=all.support[,"A2"],build=37)


#' Class 'ChipInfo'
#' 
#' This class annotates a microarray SNP chip with data for each SNP including chromosome,
#' id, position, strand, 'rs' id, allele 1, allele 2 for each SNP of a microarray chip,
#' in either hg18 or hg19 (build 36/37) coordinates.
#' This package makes extension use of this class of annotation object for the working
#' microarray chip, e.g, default is ImmunoChip, but Metabochip is also built-in,
#' and you can also load your own annotation if using a different chip. The class
#' is basically a GRanges object, modified to always have columns for A1, A2 (alleles), 
#' rs.id, and a quality control flag. The default display is tidier than GRanges, it has
#' nice coersion to and frame data.frame and subsetting by chromosome using [[n]] has been
#' added, in addition to normal [i,j] indexing native to GRanges.
#' SLOTS
#'  seqnames, ranges, strand, seqinfo, chip, build, and elementMetadata, containing columns:
#'  A1, A2, QCcode and rs.id.
#' METHODS
#'  "[[", show, print, length, initialize
#'  build, chip, rs.id, A1, A2, QCcode, QCcode<-, QCpass, QCfail, convTo36, convTo37
#' COERCION
#'  can use 'as' to convert to and from: GRanges, RangedData, data.frame
setClass("ChipInfo",
         contains="GRanges",
         representation(
           seqnames="Rle",
           ranges="IRanges",
           strand="Rle",
           elementMetadata="DataFrame",
           seqinfo="Seqinfo",
           chip="character", # a single string
           build="character" # a single string
         ),
         prototype(
           seqnames=Rle(factor()),
           strand=Rle(strand()),
           seqinfo=Seqinfo(),
           ranges=IRanges(),
           chip=character(), 
           build=character(), 
           elementMetadata=DataFrame(A1=NULL, A2=NULL,  QCcode=integer(), rs.id=NULL)
         )
)

#' Length for ChipInfo
#' 
#' Simply returns the number of SNPs (internally using 'nrow')
#' @param x a ChipInfo object
#' @return integer
setMethod("length", "ChipInfo", function(x) length(x@ranges))

setGeneric("chip", function(x) standardGeneric("chip"))

#' Retrieve the Chip name for ChipInfo
#' 
#' Simply returns the name of the chip, e.g, 'ImmunoChip'
#' @param x a ChipInfo object
#' @return character
setMethod("chip", "ChipInfo", function(x) x@chip)
setGeneric("build", function(x) standardGeneric("build") )

#' Retrieve the build for ChipInfo
#' 
#' Returns the UCSC build of the chip object, e.g, "hg18" or "hg19"
#' @param x a ChipInfo object
#' @return character, "hg18" or "hg19"
setMethod("build", "ChipInfo", function(x) x@build)
setGeneric("rs.id", function(x,b=TRUE) standardGeneric("rs.id") )

#' Access rs-ids for ChipInfo
#' 
#' Returns the rs-ids for the chip object, e.g, "rs689", etc
#' Only if these are annotated internally, or else a vector of NAs
#' @param x a ChipInfo object
#' @return character vector of IDs (or NAs)
setMethod("rs.id", "ChipInfo", function(x,b=TRUE) { 
  u <- mcols(x) ;  
  if("rs.id" %in% colnames(u)) { 
    U <- u[,"rs.id"] 
    if(!b) { U <- gsub("b","",U) }
    return(U)
  } else { return(NULL) } 
})
setGeneric("A1", function(x) standardGeneric("A1") )

#' Access allele 1 for ChipInfo
#' 
#' Returns the letter for the first alleles for the chip object, 
#' e.g, "A","C","G","T", etc
#' Only if these are annotated internally, or else a vector of NAs
#' @param x a ChipInfo object
#' @return character vector of allele codes (or NAs)
setMethod("A1", "ChipInfo", function(x) { u <- mcols(x) ;  if("A1" %in% colnames(u)) { u[,"A1"] } else { NULL } })
setGeneric("A2", function(x) standardGeneric("A2") )

#' Access allele 2 for ChipInfo
#' 
#' Returns the letter for the second alleles for the chip object, 
#' e.g, "A","C","G","T", etc
#' Only if these are annotated internally, or else a vector of NAs
#' @param x a ChipInfo object
#' @return character vector of allele codes (or NAs)
setMethod("A2", "ChipInfo", function(x) { u <- mcols(x) ;  if("A2" %in% colnames(u)) { u[,"A2"] } else { NULL } })

setGeneric("QCcode", function(x) standardGeneric("QCcode") )

#' Access quality control pass or fail codes for ChipInfo
#' 
#' Returns the pass or fail codes for each SNP of the chip object, 
#' e.g, 0,1,..,n etc
#' Only if these are added manually, or else all will be 'pass' (=0)
#' @param x a ChipInfo object
#' @return integer vector of pass/fail codes
setMethod("QCcode", "ChipInfo", function(x) { u <- mcols(x) ;  if("QCcode" %in% colnames(u)) { u[,"QCcode"] } else { NULL } })
setGeneric("QCcode<-", function(x,values) standardGeneric("QCcode<-") )

#' Set quality control pass or fail codes for ChipInfo
#' 
#' Allows user to set the pass or fail codes for each SNP of the chip object, 
#' e.g, 0,1,..,n etc. 0 is always pass, >0 is always fail, but each integer
#' can be used to represent a different failure type, or for simplicity, stick
#' to 0 and 1, ie, just pass and fail.
#' @param x a ChipInfo object
#' @param values new pass/fail codes, e.g, 0,1,...,n
#' @return updates the object specified with new pass/fail codes for the 'QCcode' slot
setMethod("QCcode<-", "ChipInfo", function(x,values) {
  if(length(x)==length(values)) {
    if(is.numeric(values)) {
      mcols(x)[,"QCcode"] <- as.integer(values)
    } else {
      stop("only numeric values can be inserted into the QCcode column, 0=pass, higher integers are failure codes")
    }
  } else {
    stop("mismatching lengths, tried to insert ",length(values),"new values into ChipInfo with ",length(x)," rows")
  }
  return(x)
} )

setGeneric("QCpass", function(x) standardGeneric("QCpass") )
setGeneric("QCfail", function(x,type=NA) standardGeneric("QCfail") )


#' Filter ChipInfo to for only SNPs passing QC
#' 
#' Returns the subset of the ChipInfo object for which SNPs pass quality
#' control, according to the QCcodes() slot == 0.
#' @param x a ChipInfo object
#' @return ChipInfo object for which SNPs pass quality control
setMethod("QCpass", "ChipInfo", function(x) { 
  ii <- which(QCcode(x)==0)
  if(length(ii)>0) { return(x[ii,]) } else { warning("No SNPs passed QC"); return(NULL) } })

#' Filter ChipInfo to for only SNPs failing QC
#' 
#' Returns the subset of the ChipInfo object for which SNPs fail quality
#' control, according to the QCcodes() slot > 0.
#' @param x a ChipInfo object
#' @return ChipInfo object for which SNPs fail quality control
setMethod("QCfail", "ChipInfo", function(x,type=NA) { 
  ii <- which(QCcode(x)!=0)
  if(is.numeric(type)) { 
    if(type %in% 1:100) {
      ii <- which(QCcode(x)==type) 
    } else { 
      warning("type must be an integer between 1 and 100, returning all failures") 
    }
  }
  if(length(ii)>0) { return(x[ii,]) } else { warning("All SNPs passed QC"); return(NULL) } })


#' Subset ChipInfo by chromosome
#' 
#' Returns the subset of the ChipInfo object for which SNPs are on
#' the chromosome specified, by either number or character.
#' @param x a ChipInfo object
#' @param i a chromosome number or letter, i.e, one of seqlevels(x)
#' @return ChipInfo object for the subset of SNPs on chromosome i
setMethod("[[", "ChipInfo", function(x,i,j,...) { 
  dotArgs <- list(...)
  if (length(dotArgs) > 0)
    dotArgs <- dotArgs[names(dotArgs) != "exact"]
  if (!missing(j) || length(dotArgs) > 0)
    stop("invalid subsetting")
  if (missing(i))
    stop("subscript is missing")
  if (!is.character(i) && !is.numeric(i)) 
    stop("invalid subscript type")
  if (length(i) < 1L)
    stop("attempt to select less than one element")
  if (length(i) > 1L)
    stop("attempt to select more than one element")
  cn <- chrNames(x)
  if (is.numeric(i) && !is.na(i) && (i < 1L || i > length(cn)))
    stop("subscript out of bounds")
  # do the selection #
  if(i %in% paste(chr(x))) {
    out <- x[chr(x)==i,]
  } else {
    if(is.numeric(i)) {
      out <- x[match(chr(x),chrNames(x))==i,]
    } else {
      stop("unknown index")
    }
  }
  out@build <- x@build
  out@chip <- x@chip
  return(out)
} )


setGeneric("convTo37", function(x) standardGeneric("convTo37"))
          
#' Convert ChipInfo to build 37/hg19 coordinates
#' 
#' Returns the a ChipInfo object with positions updated to build
#' 37 coordinates, assuming that the existing object was in build 36,
#' or already in build 37 coordinates, and that the build() slot was
#' entered correctly. Ensure that the value of build(x) is correct before
#' running this function for conversion; for instance, if the coordinates 
#' are already build 37/hg19, but build(x)=="hg18" (incorrect value), then
#' these coordinates will be transformed in a relative manner rendering the
#' result meaningless.
#' @param x a ChipInfo object
#' @return ChipInfo object with the build updated to hg19 coordinates
#' @seealso convTo36
setMethod("convTo37", "ChipInfo", function(x) {
  if(ucsc.sanitizer(build(x))=="hg18") {
    u <- conv.36.37(ranges=as(x,"GRanges"))
    if(length(u)==length(x)) { 
      x@ranges <- u@ranges
      all.eq <- TRUE
      if(all.eq) { all.eq <- length(seqlevels(x))==length(seqlevels(u)) }
      if(all.eq) { all.eq <- all(sort(seqlevels(x))==sort(seqlevels(u))) }
      if(!all.eq) {
        #print(seqlevels(x)); print(seqlevels(u))
        #warning("conversion altered the chromosomes"); 
        seqlevels(x) <- c(seqlevels(x),seqlevels(u)[!seqlevels(u) %in% seqlevels(x)])  #x@seqinfo <- u@seqinfo 
      }
      xx <- as(x@seqnames,"character")
      uu <- as(u@seqnames,"character")
      if(any(xx!=uu)) { x@seqnames <- u@seqnames }
      x@build <- "hg19"
    } else { 
      stop("conversion to build37/hg19 failed, input had length ",length(x),"; output had length ",length(u)) 
    } 
  } else {
    if(ucsc.sanitizer(build(x))!="hg19") { 
      warning("input object was not tagged as hg18/build36 [@build], left unchanged") 
    } else {
      warning("object is already using hg19/build37, no change")
    }
  }
  return(x)
})

setGeneric("convTo36", function(x) standardGeneric("convTo36"))

#' Convert ChipInfo to build 36/hg18 coordinates
#' 
#' Returns the a ChipInfo object with positions updated to build
#' 36 coordinates, assuming that the existing object was in build 37,
#' or already in build 36 coordinates, and that the build() slot was
#' entered correctly. Ensure that the value of build(x) is correct before
#' running this function for conversion; for instance, if the coordinates 
#' are already build 36/hg18, but build(x)=="hg19" (incorrect value), then
#' these coordinates will be transformed in a relative manner rendering the
#' result meaningless.
#' @param x a ChipInfo object
#' @return ChipInfo object with the build updated to hg18 coordinates
#' @seealso convTo37
setMethod("convTo36", "ChipInfo", function(x) {
  if(ucsc.sanitizer(build(x))=="hg19") {
    u <- conv.37.36(ranges=as(x,"GRanges"))
    if(length(u)==length(x)) { 
      x@ranges <- u@ranges
      all.eq <- TRUE
      if(all.eq) { all.eq <- length(seqlevels(x))==length(seqlevels(u)) }
      if(all.eq) { all.eq <- any(sort(seqlevels(x))==sort(seqlevels(u))) }
      if(!all.eq) {
        #warning("conversion altered the chromosomes"); 
        seqlevels(x) <- c(seqlevels(x),seqlevels(u)[!seqlevels(u) %in% seqlevels(x)])  #x@seqinfo <- u@seqinfo 
      }
      xx <- as(x@seqnames,"character")
      uu <- as(u@seqnames,"character")
      if(any(xx!=uu)) { x@seqnames <- u@seqnames }
      x@build <- "hg18"
    } else { 
      stop("conversion to build36/hg18 failed") 
    } 
  } else {
    if(ucsc.sanitizer(build(x))!="hg18") { 
      warning("input object was not tagged as hg19/build37 [@build], left unchanged") 
    } else {
      warning("object is already using hg18/build36, no change")
    }
  }
  return(x)
})



#' Display a ChipInfo object
#' 
#' Returns a preview of a ChipInfo object to the console. This
#' is similar to a GRanges preview, but the seqlevels are hidden, the UCSC
#' build and chip name are displayed, start and end are merged to the virtual
#' label 'position' (as it's assume we are dealing with SNPs, not ranges), the strand
#' by default is hidden, and the integer codes for pass/fail in QCcodes() are 
#' displayed as 'pass' or 'fail', even though this is not how they are represented internally. 
#' @param object a ChipInfo object
#' @param up.to only SNPs at the start and end are generally displayed, however this
#' parameter specifies that when there are <= 'up.to' SNPs, then all SNPs will be displayed.
#' @param head.tail number of SNPs to display at start/end (only the head and tail are
#' shown as these objects are generally very large with >100K SNPs)
#' @param show.strand logical, by default the strand is hidden, particularly given that
#' the strand can vary between different datasets of the same chip. Setting to TRUE will
#' display the strand
#' @return displays a preview of the object on the console
setMethod("show", "ChipInfo", 
     function(object) { showChipInfo(object,up.to=10,head.tail=5,show.strand=FALSE) } )

#' Print a ChipInfo object to the console
#' 
#' See 'show' as the behaviour is very similar and ... are just arguments of 'show'.
#' The key difference with 'print' instead of 'show' is that by default the parameter
#' 'up.to' is set to 50, so that any ChipInfo object (or subset) of less than or equal
#' to 50 rows will be displayed in its entirety, rather than just the top/bottom 5 rows. 
#' @param object a ChipInfo object
#' @param ... further arguments to show()
#' @return displays a preview of the object on the console
setMethod("print", "ChipInfo", 
          function(x,...) { showChipInfo(x,...) } )


#' Constructor for ChipInfo annotation object
#' 
#' This class annotates a microarray SNP chip with data for each SNP including chromosome,
#' id, position, strand, 'rs' id, allele 1, allele 2 for each SNP of a microarray chip,
#' in either hg18 or hg19 (build 36/37) coordinates.
#' This package makes extension use of this class of annotation object for the working
#' microarray chip, e.g, default is ImmunoChip, but Metabochip is also built-in,
#' and you can also load your own annotation if using a different chip. The class
#' is basically a GRanges object, modified to always have columns for A1, A2 (alleles), 
#' rs.id, and a quality control flag. The default display is tidier than GRanges, it has
#' nice coersion to and frame data.frame and subsetting by chromosome using [[n]] has been
#' added, in addition to normal [i,j] indexing native to GRanges.
#' @param GRanges a GRanges object containing chromosome, start/end = position, and strand
#' information for the chip object to be created, also rownames should be used to code
#' the chip-ids for each SNP.
#' @param chr optional, alternative to using 'GRanges' to input SNP locations, enter here 
#' a vector of chromosome numbers/letters for each SNP. The recommended coding is: 
#' 1:22, X, Y, XY, MT
#' @param pos optional, vector of positions (integers), use in conjunction with 'chr' and
#'  'ids' as an alternative way to input SNP position information instead of GRanges.
#' @param ids optional, vector of SNP chip-ids, use in conjunction with 'chr' and
#'  'pos' as an alternative way to input SNP position information instead of GRanges.
#' @param chip character, name of the chip you are making this annotation for (only used
#' for labelling purposes)
#' @param build character, either "hg18" or "hg19". Will also accept build number, 36 or 37.
#' This indicates what coordinates the object is using, and will be taken into account by
#' conversion functions, and annotation lookup functions throughout this package.
#' @param rs.id 'rs' ids are standardized ids for SNPs, these usually differ from each chips'
#' own IDs for each snp. If you don't know these, or can't find them, they can be left blank,
#' but will render the functions 'rs.to.id()' and 'id.to.rs()' useless for this ChipInfo object.
#' @param A1 the first allele letter code for each SNP, e.g, usually "A","C","G", or "T", but
#' you can use any scheme you like. Can be left blank.
#' @param A2, as for A1, but for allele 2.
#' @param QCcode optional column to keep track of SNPs passing and failing QC. You can completely
#' ignore this column. It works based on integer codes, 0,1,2, you may wish to use simple 0 and 1,
#' for pass and fail respectively, or else 0 can be pass, and 1,2,... can indicate failure for 
#' different criteria. 0 will always be treated as a pass and anything else as a fail, so you
#' can code fails however you wish.
ChipInfo <- function(GRanges=NULL, chr=NULL, pos=NULL, ids=NULL, chip="unknown chip", build="",
                     rs.id=NULL, A1=NULL, A2=NULL, QCcode=NULL) {
  if(build!="") { build <- ucsc.sanitizer(build) }
  LL <- max(c(length(chr),length(GRanges)),na.rm=T)
  if(length(A1)!=LL | length(A2)!=LL) { A1 <- A2 <- rep(NA,times=LL) }
  if(length(rs.id)!=LL) { rs.id <- rep(NA,times=LL) } else { 
    if(any(duplicated(rs.id))) { rs.id <- add.trail(rs.id) } # appends letters to stop duplicates
  }
  if(length(QCcode)!=LL) { QCcode <- rep(0,LL) }
  if(is.null(GRanges)) {
    GRanges <- make.granges(chr=chr,pos=pos,row.names=ids)
  } else {
    if(is(GRanges)[1]!="GRanges") { GRanges <- as(GRanges,"GRanges") }
  }
  df <- DataFrame(A1=A1,A2=A2,QCcode=QCcode,rs.id=rs.id)
  #print(build)
  return(new("ChipInfo", seqnames=GRanges@seqnames, ranges=GRanges@ranges,  strand=GRanges@strand,
            elementMetadata=df, seqinfo=GRanges@seqinfo,
             chip=chip, build=build))
}


#' Initialize method for ChipInfo
#' 
#' Please use the 'ChipInfo' constructor
setMethod("initialize", "ChipInfo",
              function(.Object, ...){
          		  callNextMethod(.Object, ...)
          	  })


setAs("ChipInfo", "GRanges",
      function(from) { 
        #print(is(from)); print(from@seqnames)
        out <- GRanges(from@seqnames,ranges=from@ranges,strand=from@strand,
                       seqinfo=from@seqinfo,elementMetadata=from@elementMetadata,genome=build(from))
        return(out)
      }
)

setAs("ChipInfo", "RangedData",
      function(from) { 
        out <- as(as(from,"GRanges"),"RangedData")
        if("strand" %in% colnames(out)) { out <- out[,-which(colnames(out) %in% "strand")] }
        return(out)
      }
)

setAs("ChipInfo", "data.frame", function(from) { ranged.to.data.frame(as(from,"GRanges"),include.cols=TRUE) })

setAs("GRanges", "ChipInfo", 
      function(from) { 
        bb <- genome(from)
        if(all(is.na(bb))) { build <- "" } else {
          if(length(unique(bb))==1) { build <- ucsc.sanitizer(bb[1]) } else { build <- "" }  
        }
        cN <- colnames(mcols(from))
        ii <- match(toupper("A1"), toupper(cN))
        if(!is.na(ii)) { a1 <- mcols(from)[,ii] } else { a1 <- NULL }
        ii <- match(toupper("A2"), toupper(cN))
        if(!is.na(ii)) { a2 <- mcols(from)[,ii] } else { a2 <- NULL }
        ii <- match(toupper("rs.id"), toupper(cN))
        if(!is.na(ii)) { rss <- mcols(from)[,ii] } else { rss <- NULL }
        ii <- match(toupper("QCcode"), toupper(cN))
        if(!is.na(ii)) { qcc <- mcols(from)[,ii] } else { qcc <- NULL }
        ChipInfo("ChipInfo",GRanges=from,chip="unknown chip",build=build,rs.id=rss,A1=a1,A2=a2,QCcode=qcc)
      }
)

setAs("RangedData", "ChipInfo", function(from) { as(as(from,"GRanges"),"ChipInfo") } )

setAs("data.frame", "ChipInfo", 
      function(from) { 
        rr <- data.frame.to.granges(from,chr="seqnames") 
        return(as(as(rr,"GRanges"),"ChipInfo"))
      } 
)

setValidity("ChipInfo",
            function(object) {
              if (!is.character(chip(object)) || length(chip(object)) != 1 || is.na(chip(object))) {
                return("'chip' slot must be a single string") 
              }
              if (!is.character(build(object)) || length(build(object)) != 1 || is.na(build(object))) {
                return("'build' slot must be a single string") 
              } else {
                if(!build(object) %in% c("",ucsc.sanitizer(show.valid=T)[,1])) {
                  return("'build' must be a string, 36/37 or hg18/hg19") 
                }
              }
            }
)

#internal
showChipInfo <- function (x, margin = "", with.classinfo = FALSE, print.seqlengths = FALSE,...) 
{
  lx <- length(x)
  nc <- ncol(mcols(x))
  qc <- QCcode(x)
  bb <- build(x)
  if(bb=="") { bb <- "unknown" }
  if(length(qc)==lx) { QC <- rep("pass",lx); QC[qc>0] <- paste0("fail",QC[qc>0]) ; x$QCcode <- QC }
  cat("ChipInfo for ",chip(x)," with ", lx, " ", ifelse(lx == 1L, "SNP", 
                   "SNPs")," using ",bb," coordinates",":\n", sep = "")
  out <- makePrettyMatrixForCompactPrinting2(x, .makeNakedMatFromChipInfo,...)
  if (nrow(out) != 0L) 
    rownames(out) <- paste0(margin, rownames(out))
  print(out, quote = FALSE, right = TRUE)
}

#internal
.makeNakedMatFromChipInfo <- function (x,show.strand=TRUE) 
{
  lx <- length(x)
  nc <- ncol(mcols(x))
  if(!show.strand) {
    ans <- cbind(seqnames = as.character(seqnames(x)), ranges = showAsCell(ranges(x)))
  } else { 
    ans <- cbind(seqnames = as.character(seqnames(x)), ranges = showAsCell(ranges(x)),strand = as.character(strand(x)))
  }
  extraColumnNames <- GenomicRanges:::extraColumnSlotNames(x)
  if (length(extraColumnNames) > 0L) {
    ans <- do.call(cbind, c(list(ans), lapply(GenomicRanges:::extraColumnSlots(x), 
                                              showAsCell)))
  }
  if (nc > 0L) {
    tmp <- do.call(data.frame, c(lapply(mcols(x), showAsCell), 
                                 list(check.names = FALSE)))
    ans <- cbind(ans, `|` = rep.int("|", lx), as.matrix(tmp))
    if(all(colnames(ans)[1:2]==c("seqnames","ranges"))) { colnames(ans)[1:2] <- c("chr","pos") }
  }
  ans
}


#internal
makePrettyMatrixForCompactPrinting2 <- function (x, makeNakedMat.FUN,head.tail=6,up.to=50, show.strand=TRUE) 
{
  lx <- NROW(x)
  if(lx <= up.to) { head.tail <- up.to }
  nhead <- head.tail
  ntail <- head.tail
  if (lx < (nhead + ntail + 1L)) {
    ans <- makeNakedMat.FUN(x,show.strand=show.strand)
    ans_rownames <- .rownames3(names(x), lx)
  }
  else {
    top_idx <- 1:nhead
    if (nhead == 0) 
      top_idx <- 0
    bottom_idx = (lx - ntail + 1L):lx
    if (ntail == 0) 
      bottom_idx <- 0
    ans_top <- makeNakedMat.FUN(x[top_idx, , drop = FALSE],show.strand=show.strand)
    ans_bottom <- makeNakedMat.FUN(x[bottom_idx, , drop = FALSE],show.strand=show.strand)
    ans <- rbind(ans_top, matrix(rep.int("...", ncol(ans_top)), 
                                 nrow = 1L), ans_bottom)
    ans_rownames <- .rownames3(names(x), lx, top_idx, bottom_idx)
  }
  rownames(ans) <- format(ans_rownames, justify = "right")
  ans
}

#internal
.rownames3 <- function (names = NULL, len = NULL, tindex = NULL, bindex = NULL) 
{
  if (is.null(tindex) && is.null(bindex)) {
    if (len == 0L) 
      character(0)
    else if (is.null(names)) 
      paste0("[", seq_len(len), "]")
    else names
  }
  else {
    if (!is.null(names)) {
      c(names[tindex], "...", names[bindex])
    }
    else {
      s1 <- paste0("[", tindex, "]")
      s2 <- paste0("[", bindex, "]")
      if (all(tindex == 0)) 
        s1 <- character(0)
      if (all(bindex == 0)) 
        s2 <- character(0)
      c(s1, "...", s2)
    }
  }
}
