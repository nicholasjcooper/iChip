#if(getwd()!= "/home/ncooper"){
#  require(genoset)
#}

# examples for package
# 
# all.support <- chip.support()
# snp.info <- ChipInfo(chr=all.support[,"Chr"],pos=all.support[,"Pos"],ids=rownames(all.support),chip="ImmunoChip",rs.id=all.support[,"dbSNP"],
#                      A1=all.support[,"A1"], A2=all.support[,"A2"])
# 
# gr.snp.info <- with(all.support,makeGRanges(chr=Chr,pos=Pos,row.names=rownames(all.support)))
# 
# snp.info <- ChipInfo(GRanges=gr.snp.info,chip="ImmunoChip",rs.id=all.support[,"dbSNP"],
#                      A1=all.support[,"A1"], A2=all.support[,"A2"])
# 
# snp.info[chr(snp.info)=="MT",] # look at the mitochondrial SNP
# QCcode(snp.info)[chr(snp.info)=="MT"] <- 1 # exclude it by changing the fail code
# snp.info[["MT"]] # revisit and see it now registers as 'fail'
# QCfail(snp.info)
# xx <- conv.36.37(chr=all.support[,"Chr"],pos=all.support[,"Pos"],ids=rownames(all.support))
# ucsc(si36[["XY"]])
# si37 <- convTo37(snp.info)
# si36 <- convTo36(si37)
# snp.info <- ChipInfo(chr=all.support[,"Chr"],pos=all.support[,"Pos37"],ids=rownames(all.support),chip="ImmunoChip",rs.id=all.support[,"dbSNP"],
#                     A1=all.support[,"A1"], A2=all.support[,"A2"],build=37)


#' Class to represent SNP annotation for a microarray
#' 
#' This class annotates a microarray SNP chip with data for each SNP including chromosome,
#' id, position, strand, 'rs' id, allele 1, allele 2 for each SNP of a microarray chip,
#' in either hg18, hg19 or hg38 (build 36/37/38) coordinates.
#' This package makes extension use of this class of annotation object for the working
#' microarray chip, e.g, default is ImmunoChip, 
#' and you can also load your own annotation if using a different chip. The class
#' is basically a GRanges object, modified to always have columns for A1, A2 (alleles), 
#' rs.id, and a quality control flag. The default display is tidier than GRanges, it has
#' nice coersion to and frame data.frame and subsetting by chromosome using [[n]] has been
#' added, in addition to normal [i,j] indexing native to GRanges.
#' Note that with this package the first time annotation is used it might be slow, but
#' subsequent calls should be fast.
#' METHODS
#'  "[[", show, print, length, dim, rownames, initialize
#'  build, chip, rs.id, A1, A2, QCcode, QCcode<-, QCpass, QCfail
#'  convTo36, convTo37, convTo38
#' COERCION
#'  can use 'as' to convert to and from: GRanges, RangedData, data.frame
#'@section Fields: 
#'  \describe{
#'    \item{\code{seqnames}:}{Object of class \code{"Rle"}, containing chromosomes for each range, see GRanges.}
#'    \item{\code{ranges}:}{Object of class \code{"IRanges"}, containing genomic start and end, see GRanges.}
#'    \item{\code{strand}:}{Object of class \code{"Rle"}, containing plus or minus coding for forward or reverse strand, see GRanges.}
#'    \item{\code{seqinfo}:}{Object of class \code{"Seqinfo"}, containing chromosome listing, see GRanges.}
#'    \item{\code{chip}:}{Name, class \code{"character"}, containing user description of the chip, e.g, 'immunoChip'.}
#'    \item{\code{build}:}{Object of class \code{"character"}, annotation version, e.g, hg18, hg19, hg38, etc.}
#'    \item{\code{elementMetaData}:}{Object of class \code{"DataFrame"}, see GRanges, but with specific column names:
#'  A1, A2, QCcode and rs.id.}
#'  }
# @name ChipInfo-class
#' @aliases ChipInfo-method
#' @rdname ChipInfo-class
#' @exportClass ChipInfo
#' @author Nick Cooper
setClass("ChipInfo",
         contains="GRanges",
         slots=list(chip="character", 
                    build="character"),
         prototype=prototype(
              seqnames=Rle(factor()),
              ranges=IRanges::IRanges(),
              strand=Rle(GenomicRanges::strand()),
              elementMetadata=DataFrame(A1=NULL, A2=NULL,  QCcode=integer(), rs.id=NULL, chip.id=NULL),
              seqinfo=Seqinfo(),
              metadata=list(),
              chip=character(), 
              build=character()
            )
)



#' rownames method for GRanges objects
#' 
#' rownames: Returns the row names.
# @name rownames
#' @param x a GRanges object
#' @return rownames: Character vector of row names (e.g, SNP IDs).
#' @rdname GRanges-methods
#' @exportMethod rownames
setMethod("rownames", "GRanges", function(x) names(x@ranges))
setMethod("rownames<-", "GRanges", function(x,value) { names(x@ranges) <- value; return(x) })



# seqnames="Rle",
# ranges="IRanges",
# strand="Rle",
# elementMetadata="DataFrame",
# seqinfo="Seqinfo",
# metadata="list",
# ,
# prototype=prototype(
#   seqnames=IRanges::Rle(factor()),
#   ranges=IRanges::IRanges(),
#   strand=IRanges::Rle(GenomicRanges::strand()),
#   elementMetadata=IRanges::DataFrame(A1=NULL, A2=NULL,  QCcode=integer(), rs.id=NULL),
#   seqinfo=GenomicRanges::Seqinfo(),
#   chip=character(), 
#   build=character(),
# )

#' rownames method for ChipInfo objects
#' 
#' rownames: Returns the row names.
# @name rownames
#' @param x a ChipInfo object
#' @return rownames: Character vector of row names (SNP IDs).
#' @rdname ChipInfo-methods
#' @exportMethod rownames
setMethod("rownames", "ChipInfo", function(x) names(x@ranges))


#' dim method for ChipInfo objects
#' 
#' dim: Returns the dimension
# @name dim
#' @return dim: same as length
#' @rdname ChipInfo-methods
#' @exportMethod dim
setMethod("dim", "ChipInfo", function(x) dim(mcols(x)))


#' Length method for ChipInfo objects
#' 
#' length: Returns the number of rows
# @name length
#' @return length: integer, number of rows, same as inherited nrow()
#' @rdname ChipInfo-methods
#' @exportMethod length
setMethod("length", "ChipInfo", function(x) length(x@ranges))


#' Retrieve the Chip name for ChipInfo
#' 
#' Simply returns the name of the chip, e.g, 'ImmunoChip'
# @name chip
#' @param x a ChipInfo object
#' @return character string
#' @rdname chip-methods
#' @export
setGeneric("chip", function(x) standardGeneric("chip"))


#' @exportMethod chip
#' @rdname chip-methods
setMethod("chip", "ChipInfo", function(x) x@chip)


#' Retrieve the UCSC build for a ChipInfo object
#' 
#' Returns the UCSC build of the chip object, e.g, 'hg18', 'hg19', or 'hg38'
# @name ucsc
#' @param x a ChipInfo object
#' @return character, 'hg18', 'hg19', or 'hg38'
#' @rdname ucsc-methods
#' @export
setGeneric("ucsc", function(x) standardGeneric("ucsc") )

#' @rdname ucsc-methods
#' @exportMethod ucsc
setMethod("ucsc", "ChipInfo", function(x) x@build)


#' Access rs-ids for ChipInfo
#' 
#' Returns the rs-ids for the chip object, e.g, "rs689", etc
#' Only if these are annotated internally, or else a vector of NAs
# @name rs.id
#' @param x a ChipInfo object
#' @param b logical, whether to show 'b' suffixes on rs.ids which
#' are created in the background to allow duplicate ids to be uniquely
#' represented for lookup and reference purposes.
#' @return rs-ids: character vector of IDs (or NAs)
#' @rdname rs.id-methods
#' @export
setGeneric("rs.id", function(x,b=TRUE) standardGeneric("rs.id") )

#' @rdname rs.id-methods
#' @exportMethod rs.id
setMethod("rs.id", "ChipInfo", function(x,b=TRUE) { 
  u <- mcols(x) ;  
  if("rs.id" %in% colnames(u)) { 
    U <- u[,"rs.id"] 
    if(!b) { U <- gsub("b","",U) }
    return(U)
  } else { return(NULL) } 
})


#' Access chip-ids for ChipInfo
#' 
#' Returns the chip-ids for the chip object, e.g, "imm_1_898835", etc
#' Only if these are annotated internally, or else a vector of NAs
#' Note that the main purpose of this is because sometimes chip-ids
#' do not satisfy conditions to be an R column/row name, e.g, start
#' with a number, illegal characters, etc. So this allows certain
#' functions to return the actual chip names that would match the 
#' official manifest. These will largely be the same as the rownames,
#' but the rownames will always be valid R column names, converted from
#' the original using clean.snp.ids() [internal function]
# @name chip-id
#' @param x a ChipInfo object
#' @return chip ids: character vector of IDs (or NAs)
#' @rdname chipId-methods
#' @export
setGeneric("chipId", function(x) standardGeneric("chipId") )

#' @rdname chipId-methods
#' @exportMethod chipId
setMethod("chipId", "ChipInfo", function(x) { 
  u <- mcols(x) ;  
  if("chip.id" %in% colnames(u)) { 
    U <- u[,"chip.id"] 
    return(U)
  } else { return(NULL) } 
})

#' Access alleles for ChipInfo
#' 
#' A1/A2: Returns the letter for the A1/A2 alleles for the chip object, 
#' e.g, 'A','C','G','T', etc
#' Only if these are annotated internally, or else a vector of NAs
#' @param x a ChipInfo object
#' @rdname allele-methods
#' @family Alleles
#' @export
setGeneric("A1", function(x) standardGeneric("A1") )

#' @rdname allele-methods
#' @family Alleles
#' @exportMethod A1
setMethod("A1", "ChipInfo", function(x) { u <- mcols(x) ;  if("A1" %in% colnames(u)) { u[,"A1"] } else { NULL } })


#' @return character vector of allele codes (or NAs)
#' @rdname allele-methods
#' @family Alleles
#' @export
setGeneric("A2", function(x) standardGeneric("A2") )

#' @rdname allele-methods
#' @family Alleles
#' @exportMethod A2
setMethod("A2", "ChipInfo", function(x) { u <- mcols(x) ;  if("A2" %in% colnames(u)) { u[,"A2"] } else { NULL } })


#' Set quality control pass or fail codes for ChipInfo
#' 
#' A1<-/A2<-: Allows user to set the allele codes for each SNP of the chip object, 
#' e.g, A,C,G,T,K, etc. If you are using allele codes this is likely to 
#' necessary as each genotyping produces a different set of allele codes.
#' If using in conjunction with snpStats, remember that allele codes are
#' always flipped to be alphabetical, so the reference allele is the later
#' letter in the alphabet. Note, assignment to A2 needs to be done separately.
#' @param value new allele codes, e.g, A,C,G,T
#' @return A1<-: updates the ChipInfo object specified with new allele codes for the 'A1' slot
#' @rdname allele-methods
#' @family Alleles
#' @export
setGeneric("A1<-", function(x,value) standardGeneric("A1<-") )


#' @rdname allele-methods
#' @family Alleles
#' @exportMethod "A1<-"
setMethod("A1<-", "ChipInfo", function(x,value) {
  return(.updateAllele(x,value,"A1"))
} )




#' @return A2<-: updates the ChipInfo object specified with new allele codes for the 'A2' slot
#' @rdname allele-methods
#' @family Alleles
#' @export
setGeneric("A2<-", function(x,value) standardGeneric("A2<-") )


#' @rdname allele-methods
#' @family Alleles
#' @exportMethod "A2<-"
setMethod("A2<-", "ChipInfo", function(x,value) {
  return(.updateAllele(x,value,"A2"))
} )


#internal
.updateAllele <- function(x,value, allele="A1") {
  if(length(Dim(value))!=1) { stop("value must be a vector") }
  if(length(x)==length(value)) {
    if(is.character(value)) {
      mcols(x)[,allele] <- paste(value)
    } else {
      stop("only character values can be inserted into the ",allele," column, A,C,G,T, etc")
    }
  } else {
    stop("mismatching lengths, tried to insert ",length(value),"new allele codes into ChipInfo with ",length(x)," rows")
  }
  return(x)
}


#' Access quality control pass or fail codes for ChipInfo
#' 
#' Returns the pass or fail codes for each SNP of the chip object, 
#' e.g, 0,1,..,n etc
#' Only if these are added manually, or else all will be 'pass' (=0)
# @name QCcode
# why not? param x a ChipInfo object
#' @return integer vector of pass/fail codes
#' @rdname QC-methods
#' @family QC
#' @export
setGeneric("QCcode", function(x) standardGeneric("QCcode") )

#' @rdname QC-methods
#' @family QC
#' @exportMethod QCcode
setMethod("QCcode", "ChipInfo", function(x) { 
  u <- mcols(x) ;  if("QCcode" %in% colnames(u)) { u[,"QCcode"] } else { NULL } 
})


#' Set quality control pass or fail codes for ChipInfo
#' 
#' QCcode<-: Allows user to set the pass or fail codes for each SNP of the chip object, 
#' e.g, 0,1,..,n etc. 0 is always pass, >0 is always fail, but each integer
#' can be used to represent a different failure type, or for simplicity, stick
#' to 0 and 1, ie, just pass and fail.
#' @param x a ChipInfo object
#' @param value new pass/fail codes, e.g, 0,1,...,n
#' @return QCcode<-: updates the object specified with new pass/fail codes for the 'QCcode' slot
#' @rdname QC-methods
#' @family QC
#' @export
setGeneric("QCcode<-", function(x,value) standardGeneric("QCcode<-") )


#' @rdname QC-methods
#' @family QC
#' @exportMethod QCcode
setMethod("QCcode<-", "ChipInfo", function(x,value) {
  if(length(x)==length(value)) {
    if(is.numeric(value)) {
      mcols(x)[,"QCcode"] <- as.integer(value)
    } else {
      stop("only numeric values can be inserted into the QCcode column, 0=pass, higher integers are failure codes")
    }
  } else {
    stop("mismatching lengths, tried to insert ",length(value),"new values into ChipInfo with ",length(x)," rows")
  }
  return(x)
} )


#' Filter ChipInfo to for only SNPs passing QC
#' 
#' QCpass: Returns the subset of the ChipInfo object for which SNPs pass quality
#' control, according to the QCcodes() slot == 0.
# @name QCpass
# @param x a ChipInfo object
#' @return QCpass: ChipInfo object for which SNPs pass quality control
#' @rdname QC-methods
#' @family QC
#' @export
setGeneric("QCpass", function(x) standardGeneric("QCpass") )


#' Filter ChipInfo to for only SNPs failing QC
#' 
#' QCfail: Returns the subset of the ChipInfo object for which SNPs fail quality
#' control, according to the QCcodes() slot > 0.
# @name QCfail
# @param x a ChipInfo object
#' @param type integer between 1 and 100, failure type (user can assign own coding scheme)
#' @return QCfail: ChipInfo object for which SNPs fail quality control
#' @rdname QC-methods
#' @family QC
#' @export
setGeneric("QCfail", function(x,type=NA) standardGeneric("QCfail") )


#' @rdname QC-methods
#' @family QC
#' @exportMethod QCpass
setMethod("QCpass", "ChipInfo", function(x) { 
  ii <- which(QCcode(x)==0)
  if(length(ii)>0) { return(x[ii,]) } else { warning("No SNPs passed QC"); return(NULL) } })


#' @rdname QC-methods
#' @family QC
#' @exportMethod QCfail
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
#' @param j always leave missing, not applicable for this method.
#' @param ... further arguments - again there should not be any
#' @return ChipInfo object for the subset of SNPs on chromosome i
#' @rdname ChipInfo-subset
#' @exportMethod "[["
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


#' Convert ChipInfo between build 36/37/38 coordinates
#' 
#' Returns the a ChipInfo object with positions updated to build 36/37/38
#' coordinates, assuming that the build() slot was entered correctly. 
#' Ensure that the value of ucsc(x) is correct before
#' running this function for conversion; for instance, if the coordinates 
#' are already build 37/hg19, but ucsc(x)!="hg19" (incorrect value), then
#' these coordinates will be transformed in a relative manner rendering the
#' result meaningless.
# @name convTo37
#' @param x a ChipInfo object
#' @return convTo37: Returns a ChipInfo object with the build updated to hg19 coordinates
#' @family conversion
#' @rdname conv-methods
#' @export
setGeneric("convTo37", function(x) standardGeneric("convTo37"))
          
#' @family conversion
#' @aliases convTo37
#' @rdname conv-methods
#' @exportMethod convTo37
setMethod("convTo37", "ChipInfo", function(x) {
  if(ucsc.sanitizer(ucsc(x)) %in% c("hg18","hg38")) {
    if(ucsc.sanitizer(ucsc(x)) == c("hg18")) {
      u <- conv.36.37(ranges=as(x,"GRanges"))
    } else {
      u <- conv.38.37(ranges=as(x,"GRanges"))
    }
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
    if(ucsc.sanitizer(ucsc(x))!="hg19") { 
      warning("input object was not tagged as hg18/build36/hg38/build38 [@build], left unchanged") 
    } else {
      warning("object is already using hg19/build37, no change")
    }
  }
  return(x)
})



# @name convTo36
# @param x a ChipInfo object
#' @return convTo36: Returns a ChipInfo object with the build updated to hg18 coordinates
#' @family conversion
#' @rdname conv-methods
#' @export
setGeneric("convTo36", function(x) standardGeneric("convTo36"))


#' @family conversion
#' @aliases convTo36
#' @rdname conv-methods
#' @exportMethod convTo36
setMethod("convTo36", "ChipInfo", function(x) {
  if(ucsc.sanitizer(ucsc(x))=="hg19") {
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
    if(ucsc.sanitizer(ucsc(x))!="hg18") { 
      warning("input object was not tagged as hg19/build37 [@build], left unchanged") 
    } else {
      warning("object is already using hg18/build36, no change")
    }
  }
  return(x)
})



# @name convTo38
# @param x a ChipInfo object
#' @return convTo38: Returns a ChipInfo object with the build updated to hg38 coordinates
#' @family conversion
#' @rdname conv-methods
#' @export
setGeneric("convTo38", function(x) standardGeneric("convTo38"))


#' @family conversion
#' @aliases convTo38
#' @rdname conv-methods
#' @exportMethod convTo38
setMethod("convTo38", "ChipInfo", function(x) {
  if(ucsc.sanitizer(ucsc(x))=="hg18") { stop("can't convert 36 to 38; must convert from 36 to 37, then convert from 37 to 38") }
  if(ucsc.sanitizer(ucsc(x))=="hg19") {
    u <- conv.37.38(ranges=as(x,"GRanges"))
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
      x@build <- "hg38"
    } else { 
      stop("conversion to build38/hg38 failed") 
    } 
  } else {
    if(ucsc.sanitizer(ucsc(x))!="hg38") { 
      warning("input object was not tagged as hg19/build37 [@build], left unchanged") 
    } else {
      warning("object is already using hg18/build36, no change")
    }
  }
  return(x)
})


#' Display method for ChipInfo objects
#' 
#' show: Displays a preview of the object
# @name show
#' @param object a ChipInfo object
#' @return show: Displays a preview of the object
#' @rdname ChipInfo-methods
#' @exportMethod show
setMethod("show", "ChipInfo", 
          function(object) { showChipInfo(object,up.to=10,head.tail=5,show.strand=FALSE) } )


#' Print a ChipInfo object to the console
#' 
#' print: See 'show' as the behaviour is very similar and ... are just arguments of 'show'.
#' The key difference with 'print' instead of 'show' is that by default the parameter
#' 'up.to' is set to 50, so that any ChipInfo object (or subset) of less than or equal
#' to 50 rows will be displayed in its entirety, rather than just the top/bottom 5 rows. 
# @name print
# @param x a ChipInfo object
#' @param ... further arguments to showChipInfo()
#' @return print: Prints the object to terminal using 'showChipInfo()'.
#' @rdname ChipInfo-methods
#' @exportMethod print
setMethod("print", "ChipInfo", 
          function(x,...) { showChipInfo(x,...) } )


#' Constructor (wrapper) for ChipInfo annotation object
#' 
#' This class annotates a microarray SNP chip with data for each SNP including chromosome,
#' id, position, strand, 'rs' id, allele 1, allele 2 for each SNP of a microarray chip,
#' in either hg18, hg19 or hg38 (build 36/37/38) coordinates.
#' This package makes extension use of this class of annotation object for the working
#' microarray chip, e.g, default is ImmunoChip, but Metabochip is also built-in,
#' and you can also load your own annotation if using a different chip. The class
#' is basically a GRanges object, modified to always have columns for A1, A2 (alleles), 
#' rs.id, and a quality control flag. The default display is tidier than GRanges, it has
#' nice coersion to and frame data.frame and subsetting by chromosome using [[n]] has been
#' added, in addition to normal [i,j] indexing native to GRanges.
# @name ChipInfo
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
#' @param chip.id chip ids are the chip-specific ids for SNPs, these usually differ between chips'
#' even for the same snp. If you don't know these, or can't find them, they can be left blank,
#' but will render the function 'chip.id()' useless for this ChipInfo object. The main purpose
#' of this parameter is for when the real chip ids are not valid R row/column name strings, and
#' by using this column, some functions can return the real chip ids instead of the sanitized 
#' version
#' @param A1 the first allele letter code for each SNP, e.g, usually "A","C","G", or "T", but
#' you can use any scheme you like. Can be left blank.
#' @param A2, as for A1, but for allele 2.
#' @param QCcode optional column to keep track of SNPs passing and failing QC. You can completely
#' ignore this column. It works based on integer codes, 0,1,2, you may wish to use simple 0 and 1,
#' for pass and fail respectively, or else 0 can be pass, and 1,2,... can indicate failure for 
#' different criteria. 0 will always be treated as a pass and anything else as a fail, so you
#' can code fails however you wish.
#' @export
ChipInfo <- function(GRanges=NULL, chr=NULL, pos=NULL, ids=NULL, chip="unknown chip", build="",
                     rs.id=NULL, chip.id=NULL, A1=NULL, A2=NULL, QCcode=NULL) {
  if(build!="") { build <- ucsc.sanitizer(build) }
  LL <- max(c(length(chr),length(GRanges)),na.rm=T)
  if(length(A1)!=LL | length(A2)!=LL) { A1 <- A2 <- rep(NA,times=LL) }
  if(length(rs.id)!=LL) { rs.id <- rep(NA,times=LL) } else { 
    if(any(duplicated(rs.id))) { rs.id <- add.trail(rs.id) } # appends letters to stop duplicates
  }
  if(length(chip.id)!=LL) { chip.id <- rep(NA,times=LL) } else { 
    if(any(duplicated(narm(chip.id)))) { stop("chip.id shouldn't contain duplicates") } # appends letters to stop duplicates
  }
  if(length(QCcode)!=LL) { QCcode <- rep(0,LL) }
  if(is.null(GRanges)) {
    GRanges <- makeGRanges(chr=chr,pos=pos,row.names=ids)
  } else {
    if(is(GRanges)[1]!="GRanges") { GRanges <- as(GRanges,"GRanges") }
  }
  df <- DataFrame(A1=A1,A2=A2,QCcode=QCcode,rs.id=rs.id,chip.id=chip.id)
  #print(build)
  return(new("ChipInfo", seqnames=GRanges@seqnames, ranges=GRanges@ranges,  strand=GRanges@strand,
            elementMetadata=df, seqinfo=GRanges@seqinfo,
             chip=chip, build=build))
}


#' Initialize (constructor) method for ChipInfo
#' 
#' Use the 'ChipInfo()' wrapper to construct ChipInfo objects from scratch
#  @name ChipInfo
#' @param .Object An object generated from the ChipInfo class prototype,
#'  see methods:initialize
#' @param ... Additional arguments to initialize. None recommended.
#' @rdname ChipInfo-class
#' @exportMethod initialize
setMethod("initialize", "ChipInfo",
              function(.Object, ...){
          		  callNextMethod(.Object, ...)
          	  })


#' As("ChipInfo", "GRanges")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
setAs("ChipInfo", "GRanges",
      function(from) { 
        #print(is(from)); print(from@seqnames)
        out <- GRanges(from@seqnames,ranges=from@ranges,strand=from@strand,
                       seqinfo=from@seqinfo,elementMetadata=from@elementMetadata,genome=ucsc(from))
        return(out)
      }
)

#' As("ChipInfo", "GRanges")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
setAs("ChipInfo", "RangedData",
      function(from) { 
        out <- as(as(from,"GRanges"),"RangedData")
        if("strand" %in% colnames(out)) { out <- out[,-which(colnames(out) %in% "strand")] }
        colnames(out) <- gsub("elementMetadata.","",colnames(out),fixed=TRUE)
        colnames(out) <- gsub("elementMetadata","",colnames(out),fixed=TRUE)
        return(out)
      }
)


#' As("ChipInfo", "GRanges")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
setAs("ChipInfo", "data.frame", function(from) { ranged.to.data.frame(as(from,"GRanges"),include.cols=TRUE) })


#' As("GRanges", "ChipInfo")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
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
        ii <- match(toupper("chip.id"), toupper(cN))
        if(!is.na(ii)) { cid <- mcols(from)[,ii] } else { cid <- NULL }
        ii <- match(toupper("QCcode"), toupper(cN))
        if(!is.na(ii)) { qcc <- mcols(from)[,ii] } else { qcc <- NULL }
        ChipInfo("ChipInfo",GRanges=from,chip="unknown chip",build=build,rs.id=rss,chip.id=cid,A1=a1,A2=a2,QCcode=qcc)
      }
)

#' As("RangedData", "ChipInfo")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
setAs("RangedData", "ChipInfo", function(from) { as(as(from,"GRanges"),"ChipInfo") } )

#' As("data.frame", "ChipInfo")
#'
#' @name as
# @rdname ChipInfo-class
#' @export
setAs("data.frame", "ChipInfo", 
      function(from) { 
        rr <- df.to.GRanges(from,chr="seqnames") 
        return(as(as(rr,"GRanges"),"ChipInfo"))
      } 
)



# improved coersion functions for data.frames to RangedData/GRanges

#' As("data.frame", "RangedData")
#'
#' @name as
#' @export
setAs("data.frame", "RangedData", function(from) { return(df.to.ranged(from,GRanges=FALSE))  } )

#' As("data.frame", "GRanges")
#'
#' @name as
#' @export
setAs("data.frame", "GRanges", function(from) { return(df.to.ranged(from,GRanges=TRUE))  } )

#' As("RangedData", "data.frame")
#'
#' Note that for automatic conversion of a data.frame to RangedData/GRanges, a column named 'chr' 
#' or 'seqnames' in the data.frame is expected/required to make the conversion effectively. 
#' Otherwise use 'ranged.to.data.frame()'
#' @name as
#' @export
setAs("RangedData", "data.frame", function(from) { return(ranged.to.data.frame(from))  } )

#' As("GRanges", "data.frame")
#'
#' @name as
#' @export
setAs("GRanges", "data.frame", function(from) { return(ranged.to.data.frame(from))  } )



# No roxygen required i think?
setValidity("ChipInfo",
            function(object) {
              if (!is.character(chip(object)) || length(chip(object)) != 1 || is.na(chip(object))) {
                return("'chip' slot must be a single string") 
              }
              if (!is.character(ucsc(object)) || length(ucsc(object)) != 1 || is.na(ucsc(object))) {
                return("'build' slot must be a single string") 
              } else {
               # requireNamespace(humarray)
                if(!paste(ucsc(object)) %in% c("",paste(as.vector(ucsc.sanitizer(show.valid=T)[,1])))) {
                  return("'build' must be a string, 36/37/38 or hg18/hg19/hg38") 
                }
              }
            }
)

#' Display a ChipInfo object
#' 
#' Returns a preview of a ChipInfo object to the console. This
#' is similar to a GRanges preview, but the seqlevels are hidden, the UCSC
#' build and chip name are displayed, start and end are merged to the virtual
#' label 'position' (as it's assume we are dealing with SNPs, not ranges), the strand
#' by default is hidden, and the integer codes for pass/fail in QCcodes() are 
#' displayed as 'pass' or 'fail', even though this is not how they are represented internally.
#' This is called by the default 'show' method for ChipInfo objects. 
#' @param x a ChipInfo object
#' @param margin margin for display, usually ""
#' @param with.classinfo logical, whether to display class information
#' @param print.seqlengths logical, whether to display sequence lengths below
#' the main output listing (e.g, chromsomes). Usually tidier when this is FALSE.
#' @param ... hidden arguments including: 'head.tail'; number of SNPs to display 
#' at start/end (only the head and tail are shown as these objects are generally
#' very large with >100K SNPs); 'up.to'; only SNPs at the start and end are generally
#' displayed, however this parameter specifies that when there are <= 'up.to' SNPs,
#' then all SNPs will be displayed; 'show.strand'; logical, by default the strand is 
#' hidden, particularly given that the strand can vary between different datasets 
#' of the same chip. Setting to TRUE will display the strand.
#' @return print compact preview of the object to the standard output (terminal)
#' @seealso \code{\link{ChipInfo}}
#' @export
showChipInfo <- function (x, margin = "", with.classinfo = FALSE, print.seqlengths = FALSE,...) 
{
  lx <- length(x)
  nc <- ncol(mcols(x))
  qc <- QCcode(x)
  bb <- ucsc(x)
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
extraColumnSlots2 <- function(x) {
  sapply(extraColumnSlotNames2(x), slot, object = x, simplify = FALSE)
}

#internal
setGeneric("extraColumnSlotNames2",
           function(x) standardGeneric("extraColumnSlotNames2"))

setMethod("extraColumnSlotNames2", "ANY", function(x) character())


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
  extraColumnNames <- extraColumnSlotNames2(x)
  if (length(extraColumnNames) > 0L) {
    ans <- do.call(cbind, c(list(ans), lapply(extraColumnSlots2(x), 
                                              showAsCell)))
  }
  if (nc > 0L) {
    df <- mcols(x)
    if(tail(colnames(df),1)=="chip.id") { df <- df[,-ncol(df)] }
    tmp <- do.call(data.frame, c(lapply(df, showAsCell), 
                                 list(check.names = FALSE)))  ### hide chip.id here!
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



#' Select chromosome subset for ranged objects
#' 
#' Returns the object filtered for specific chromosomes for a ranged object
#' @param object a ChipInfo, GRanges or RangedData object
#' @param chr vector, string or numeric of which chromosome(s) to select
#' @return vector of chromosome values for each range/SNP
#' @rdname chrSel-methods
#' @export
setGeneric("chrSel", function(object,chr) standardGeneric("chrSel"))


#' Select chromosome subset for RangedData objects
#' 
#' Returns the object filtered for specific chromosomes for a RangedData object
#' @rdname chrSel-methods
#' @exportMethod chrSel
setMethod("chrSel", "RangedData", function(object,chr) {
  return(humarray::chrSelect(object,chr))
})


#' Select chromosome subset for GRanges objects
#' 
#' Returns the object filtered for specific chromosomes for a GRanges object
#' @rdname chrSel-methods
#' @exportMethod chrSel
setMethod("chrSel", "GRanges", function(object,chr) {
  return(humarray::chrSelect(object,chr))
})

#' Select chromosome subset for ChipInfo objects
#' 
#' Returns the object filtered for specific chromosomes for a GRanges object
#' @rdname chrSel-methods
#' @exportMethod chrSel
setMethod("chrSel", "ChipInfo", function(object,chr) {
  return(humarray::chrSelect(object,chr))
})



#' Chromosome method for RangedData objects
#' 
#' Return the list of chromosome values from a RangedData object
#' @param object RangedData object
#' @return vector of chromosome values for each range/SNP
#' @rdname chrm-methods
#' @export
setGeneric("chrm",function(object) standardGeneric("chrm"))

#' @rdname chrm-methods
#' @importMethodsFrom genoset  chr  
#' @exportMethod chrm
setMethod("chrm", "RangedData", function(object) {
  return(chr2(object))
})

#' @rdname chrm-methods
#' @importMethodsFrom genoset  chr  
#' @importFrom genoset  chr  
#' @exportMethod chrm
setMethod("chrm", "GRanges", function(object) {
  return(genoset::chr(object))
})

#' @rdname chrm-methods
#' @importMethodsFrom genoset  chr  
#' @importFrom genoset  chr  
#' @exportMethod chrm
setMethod("chrm", "ChipInfo", function(object) {
  return(genoset::chr(object))
})



#' @importMethodsFrom "genoset"  toGenomeOrder  
#' @exportMethod toGenomeOrder
setMethod("toGenomeOrder", "RangedData", function(ds) {
  return(TGORD(ds))
})



#' @importMethodsFrom "genoset"   chrIndices  
#' @exportMethod chrIndices
setMethod("chrIndices", "RangedData", function(object) {
  return(chrIndices2(object))
})



#' @importMethodsFrom "genoset"   chrInfo 
#' @exportMethod chrInfo
setMethod("chrInfo", "RangedData", function(object) {
  return(chrInfo2(object))
})




#' @importMethodsFrom "genoset"  chrNames  
#' @exportMethod chrNames
setMethod("chrNames", "RangedData", function(object) {
  return(chrNames2(object))
})


#' Plot method for GRanges objects
#' 
#' See plotRanges()
# @name plot
#' @param x a GRanges or RangedData object
#' @param y not used for plotRanges
#' @param ... further arguments, see plotRanges()
#' @rdname plot-methods
#' @aliases GRanges GRanges-method
#' @seealso \code{\link{plotRanges}}
#' @exportMethod plot
setMethod("plot", "GRanges", function(x,y,...) {
  plotRanges(ranged=x,...)
})


#' Plot method for RangedData objects
#' 
# @name plot
#' @rdname plot-methods
#' @aliases RangedData RangedData-method
#' @exportMethod plot
setMethod("plot", "RangedData", function(x,y,...) {
  plotRanges(ranged=x,...)
})
