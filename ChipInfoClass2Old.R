
setClass("ChipInfo",
         contains="GRanges",
         representation(
           seqnames="Rle",
           ranges="IRanges",
           strand="Rle",
           elementMetadata="DataFrame",
           seqinfo="Seqinfo",
           chip="character", # a single string
           build="character", # a single string
           rs.id="character", # an optional character vector of length N
           alt.pos="integer", # an optional integer vector of length N
           A1="character", # an optional character vector of length N
           A2="character", # an optional character vector of length N
           QCfail="integer" # an optional integer vector of length N
         ),
         prototype(
           seqnames=Rle(factor()),
           strand=Rle(strand()),
           chip=character(), 
           build=character(), 
           rs.id=NULL,
           alt.pos=NULL, 
           A1=NULL, 
           A2=NULL, 
           QCfail=integer() 
         )
)


setMethod("length", "ChipInfo", function(x) length(x@ranges))
setGeneric("chip", function(x) standardGeneric("chip"))
setMethod("chip", "ChipInfo", function(x) x@chip)
setGeneric("build", function(x) standardGeneric("build"))
setMethod("build", "ChipInfo", function(x) x@build)
setGeneric("rs.id", function(x) standardGeneric("rs.id"))
setMethod("rs.id", "ChipInfo", function(x) x@rs.id)
setGeneric("A1", function(x) standardGeneric("A1"))
setMethod("A1", "ChipInfo", function(x) x@A1)
setGeneric("A2", function(x) standardGeneric("A2"))
setMethod("A2", "ChipInfo", function(x) x@A2)
setGeneric("QCfail", function(x) standardGeneric("QCfail"))
setMethod("QCfail", "ChipInfo", function(x) x@QCfail)

setMethod("show", "ChipInfo", 
     function(object) { showChipInfo(object) } )

#snp.info <

#constructor
ChipInfo <- function(GRanges=NULL, chr=NULL, pos=NULL, ids=NULL, chip="unknown chip", build="hg18",
                     rs.id=NULL, A1=NULL, A2=NULL, QCfail=rep(0,nrow(GRanges))) {
  build <- ucsc.sanitizer(build)
  if(is.null(GRanges)) {
    GRanges <- make.chr.pos.granges(chr=chr,pos=pos)
  } else {
    if(is(GRanges)[1]!="GRanges") { GRanges <- as(GRanges,"GRanges") }
  }
  #return(GRanges) #); print(is(GRanges))
  if(build=="hg18") {
    alt.pos <- start(conv.36.37(GRanges))
  } else { 
    alt.pos <- start(conv.37.36(GRanges))
  }
  #return(GRanges@ranges)
  return(new("ChipInfo", seqnames=GRanges@seqnames, ranges=GRanges@ranges,  strand=GRanges@strand,
            elementMetadata=DataFrame(), seqinfo=GRanges@seqinfo,
             chip=chip, build=ucsc.sanitizer(build), rs.id=rs.id, 
             alt.pos=as.integer(alt.pos), A1=A1, A2=A2, QCfail=as.integer(QCfail)))
}


setMethod("initialize", "ChipInfo",
              function(.Object, ...){
          		  callNextMethod(.Object, ...)
          	  })

# setMethod("initialize",
#           "ChipInfo",
#           function(.Object, seqnames=Rle(factor()),
#                    strand=Rle(strand()),
#                    ranges=IRanges(),
#                    elementMetadata=DataFrame(),
#                    seqinfo=Seqinfo(),
#                    chip=character(), 
#                    build=character(), 
#                    rs.id=character(),
#                    alt.pos=integer(), 
#                    A1=character(), 
#                    A2=character(), 
#                    QCfail=integer()) {
#             .Object
#           })


setAs("ChipInfo", "GRanges",
      function(from) { 
        #print(is(from)); print(from@seqnames)
        out <- GRanges(from@seqnames,ranges=from@ranges,strand=from@strand,
                       seqinfo=from@seqinfo,elementMetadata=from@elementMetadata)
       # print(length(out)); print(head(rs.id(from)))
      #  if(length(rs.id(from))>0) { if(length(rs.id(from))==nrow(out)) { out[["rs.id"]] <- rs.id(from) } }
      #  if(length(A1(from))>0) { if(length(A1(from))==nrow(out)) { out[["A1"]] <- A1(from) } }
      #  if(length(A2(from))>0) { if(length(A2(from))==nrow(out)) { out[["A2"]] <- A2(from) } }
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

setAs("ChipInfo", "data.frame", function(from) { ranges.to.data.frame(as(from,"GRanges"),include.cols=TRUE) })

setAs("GRanges", "ChipInfo", 
      function(from) { 
        ChipInfo("ChipInfo",GRanges=from,chip="unknown chip",build="hg18",rs.id=NULL,A1=NULL,A2=NULL)
      }
)

setAs("RangedData", "ChipInfo", function(from) { as(as(from,"GRanges"),"ChipInfo") } )

setAs("data.frame", "ChipInfo", 
      function(from) { 
        rr <- data.frame.to.GRanges(from,chr="seqnames") 
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
                if(!build(object) %in% ucsc.sanitizer(show.valid=T)[,1]) {
                  return("'build' must be a string, 36,37 or hg18/hg19") 
                }
              }
              slot_lengths <- c(length(rs.id(object)),
                                length(A1(object)), length(A2(object)), length(QCfail(object)))
              #print(slot_lengths)
              if(any(!slot_lengths %in% c(0,1,length(object)))) {
                cat("ChipInfo had length",length(object),"\n")
                prv(slot_lengths)
                return("invalid lengths for at least one of the slots: 'rs.id', 'A1', 'A2' or 'QCfail'")
              } else { 
                if(length(A1(object))!=length(A2(object))) {
                  return("slots 'A1' and 'A2' should both be of the same length, with both present, or both empty")
                } else {  TRUE }
              }
            }
)




showChipInfo <- function (x, margin = "", with.classinfo = FALSE, print.seqlengths = FALSE) 
{
  lx <- length(x)
  nc <- ncol(mcols(x))
  rs <- rs.id(x);   a1 <- A1(x);   a2 <- A2(x);   qc <- QCfail(x)
  if(length(rs)==lx) { nc <- nc+1; RS <- rs } else { RS <- NULL }
  if(length(a1)==lx) { nc <- nc+1; A1 <- a1 } else { A1 <- NULL }
  if(length(a2)==lx) { nc <- nc+1; A2 <- a2 } else { A2 <- NULL }
  if(length(qc)==lx) { nc <- nc+1; QC <- rep("pass",lx); QC[qc>0] <- "fail" } else { QC <- NULL }
  cat("ChipInfo with ", lx, " ", ifelse(lx == 1L, "SNP", 
                   "SNPs"), " and ", nc, " metadata ", ifelse(nc == 1L, 
                                     "column", "columns"), ":\n", sep = "")
  out <- makePrettyMatrixForCompactPrinting2(x, .makeNakedMatFromChipInfo,nc=nc,a1=A1,a2=A2,rs=RS,qc=QC)
  if (nrow(out) != 0L) 
    rownames(out) <- paste0(margin, rownames(out))
  print(out, quote = FALSE, right = TRUE)
}


.makeNakedMatFromChipInfo <- function (x,idx=NA,nc=0,a1=NULL,a2=NULL,rs=NULL,qc=NULL) 
{
  lx <- length(x)
  ans <- cbind(seqnames = as.character(seqnames(x)), ranges = showAsCell(ranges(x)), 
               strand = as.character(strand(x)))
  ecn <- c(if(!is.null(a1)) { "A1" } else { NULL },if(!is.null(a2)) { "A2" } else { NULL },
           if(!is.null(rs)) { "rs.id" } else { NULL },if(!is.null(qc)) { "QCfail" } else { NULL })
  #extraColumnNames <- ecn
  extraColumnNames <- GenomicRanges:::extraColumnSlotNames(x)
  if (length(extraColumnNames) > 0L) {
    ans <- do.call(cbind, c(list(ans), lapply(GenomicRanges:::extraColumnSlots(x), 
                                              showAsCell)))
  }
  if (length(ecn) > 0 & is.numeric(idx)) {
    tmp <- do.call(data.frame, c(lapply(mcols(x), showAsCell), 
                                 list(check.names = FALSE)))
    codz <- substr(tolower(ecn),1,2)
    tmp <- as.data.frame(matrix(nrow=length(tmp),ncol=length(ecn)))
    colnames(tmp) <- ecn
    for (cc in 1:length(ecn)) {
      tmp[[cc]] <- get(codz[cc])[idx]
    }
    ans <- cbind(ans, `|` = rep.int("|", lx), as.matrix(tmp))
  }
  ans
}


makePrettyMatrixForCompactPrinting2 <- function (x, makeNakedMat.FUN,...) 
{
  lx <- NROW(x)
  nhead <- get_showHeadLines()
  ntail <- get_showTailLines()
  if (lx < (nhead + ntail + 1L)) {
    ans <- makeNakedMat.FUN(x,...)
    ans_rownames <- .rownames2(names(x), lx)
  }
  else {
    top_idx <- 1:nhead
    if (nhead == 0) 
      top_idx <- 0
    bottom_idx = (lx - ntail + 1L):lx
    if (ntail == 0) 
      bottom_idx <- 0
    ans_top <- makeNakedMat.FUN(x[top_idx, , drop = FALSE],idx=top_idx,...)
    ans_bottom <- makeNakedMat.FUN(x[bottom_idx, , drop = FALSE],idx=bottom_idx,...)
    ans <- rbind(ans_top, matrix(rep.int("...", ncol(ans_top)), 
                                 nrow = 1L), ans_bottom)
    ans_rownames <- .rownames2(names(x), lx, top_idx, bottom_idx)
  }
  rownames(ans) <- format(ans_rownames, justify = "right")
  ans
}


# 
# GenomicRanges:::showGenomicRanges
# function (x, margin = "", with.classinfo = FALSE, print.seqlengths = FALSE) 
# {
#   lx <- length(x)
#   nc <- ncol(mcols(x))
#   cat(class(x), " with ", lx, " ", ifelse(lx == 1L, "range", 
#                                           "ranges"), " and ", nc, " metadata ", ifelse(nc == 1L, 
#                                                                                        "column", "columns"), ":\n", sep = "")
#   out <- IRanges:::makePrettyMatrixForCompactPrinting(x, .makeNakedMatFromGenomicRanges)
#   if (with.classinfo) {
#     .COL2CLASS <- c(seqnames = "Rle", ranges = "IRanges", 
#                     strand = "Rle")
#     extraColumnNames <- extraColumnSlotNames(x)
#     .COL2CLASS <- c(.COL2CLASS, getSlots(class(x))[extraColumnNames])
#     classinfo <- makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
#     stopifnot(identical(colnames(classinfo), colnames(out)))
#     out <- rbind(classinfo, out)
#   }
#   if (nrow(out) != 0L) 
#     rownames(out) <- paste0(margin, rownames(out))
#   print(out, quote = FALSE, right = TRUE)
#   if (print.seqlengths) {
#     cat(margin, "---\n", sep = "")
#     showSeqlengths(x, margin = margin)
#   }
# }
