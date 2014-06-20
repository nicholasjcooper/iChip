
setClass("ChipInfo",
         contains="GRanges",
         representation(
           chip="character", # a single string
           build="character", # a single string
           rs.id="character", # an optional character vector of length N
           alt.pos="integer", # an optional integer vector of length N
           A1="character", # an optional character vector of length N
           A2="character", # an optional character vector of length N
           QCfail="integer" # an optional integer vector of length N
         ) )


setMethod("length", "ChipInfo", function(x) nrow(x@ranges))
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
     function(x) { GenomicRanges:::showGenomicRanges(x, with.classinfo = FALSE,print.seqlengths = FALSE)} )

#constructor
ChipInfo <- function(GRanges=NULL, chr=NULL, pos=NULL, chip="unknown chip", build="hg18",
                     rs.ids=NULL, A1=NULL, A2=NULL, QCfail=rep(0,nrow(GRanges))) {
  build <- ucsc.sanitizer(build)
  if(is.null(GRanges)) {
    GRanges <- make.chr.pos.granges(chr=chr,pos=pos)
  }
  return(GRanges) #); print(is(GRanges))
  if(build=="hg18") {
    alt.pos <- start(conv.36.37(GRanges))
  } else { 
    alt.pos <- start(conv.37.36(GRanges))
  }
  return(new("ChipInfo", chip=chip, build=ucsc.sanitizer(build),
      rs.ids=rs.ids, alt.pos=alt.pos, A1=A1, A2=A2, QCfail=QCfail))
}

setAs("ChipInfo", "GRanges",
      function(from) { 
        out <- from@ranges
        if(length(rs.id(from))==length(from)) { out[["rs.id"]] <- rs.id(from) }
        if(length(A1(from))==length(from)) { out[["A1"]] <- A1(from) }
        if(length(A2(from))==length(from)) { out[["A2"]] <- A2(from) }
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
        newChipInfo("ChipInfo",ranges=from@ranges,seqnames=from@seqnames,strand=from@strand,
                    chip="unknown chip",build="hg18",rs.id=NULL,A1=NULL,A2=NULL)
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
                if(!build %in% ucsc.sanitizer(show.valid=T)[,1]) {
                  return("'build' must be a string, 36,37 or hg18/hg19") 
                }
              }
              slot_lengths <- c(length(rs.ids(object)),
                                length(A1(object)), length(A2(object)), length(QCfail(object)))
              if(any(!slot_lengths %in% c(0,1,length(object)))) {
                return("invalid lengths for at least one of the slots: 'rs.ids', 'A1', 'A2' or 'QCfail'")
              } else { TRUE }
            }
)



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
