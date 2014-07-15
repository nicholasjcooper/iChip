require(genoset)

# examples
# 
# snp.info <- ChipInfo(chr=all.support[,"Chr"],pos=all.support[,"Pos"],ids=rownames(all.support),chip="ImmunoChip",rs.id=all.support[,"dbSNP"],
#                      A1=all.support[,"A1"], A2=all.support[,"A2"])
# 
# gr.snp.info <- with(all.support,make.chr.pos.granges(chr=Chr,pos=Pos,row.names=rownames(all.support)))
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
#snp.info <- ChipInfo(chr=all.support[,"Chr"],pos=all.support[,"Pos37"],ids=rownames(all.support),chip="ImmunoChip",rs.id=all.support[,"dbSNP"],
#                     A1=all.support[,"A1"], A2=all.support[,"A2"],build=37)


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
           elementMetadata=DataFrame(alt.pos=NULL,
               A1=NULL, A2=NULL,  QCcode=integer(), rs.id=NULL )
         )
)


setMethod("length", "ChipInfo", function(x) length(x@ranges))
setGeneric("chip", function(x) standardGeneric("chip"))
setMethod("chip", "ChipInfo", function(x) x@chip)
setGeneric("build", function(x) standardGeneric("build") )
setMethod("build", "ChipInfo", function(x) x@build)
setGeneric("rs.id", function(x) standardGeneric("rs.id") )
setMethod("rs.id", "ChipInfo", function(x) { u <- mcols(x) ;  if("rs.id" %in% colnames(u)) { u[,"rs.id"] } else { NULL } })
setGeneric("A1", function(x) standardGeneric("A1") )
setMethod("A1", "ChipInfo", function(x) { u <- mcols(x) ;  if("A1" %in% colnames(u)) { u[,"A1"] } else { NULL } })
setGeneric("A2", function(x) standardGeneric("A2") )
setMethod("A2", "ChipInfo", function(x) { u <- mcols(x) ;  if("A2" %in% colnames(u)) { u[,"A2"] } else { NULL } })

setGeneric("QCcode", function(x) standardGeneric("QCcode") )
setMethod("QCcode", "ChipInfo", function(x) { u <- mcols(x) ;  if("QCcode" %in% colnames(u)) { u[,"QCcode"] } else { NULL } })
setGeneric("QCcode<-", function(x,values) standardGeneric("QCcode<-") )

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
setMethod("QCpass", "ChipInfo", function(x) { 
  ii <- which(QCcode(x)==0)
  if(length(ii)>0) { return(x[ii,]) } else { warning("No SNPs passed QC"); return(NULL) } })
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
  if(i %in% paste(chr(snp.info))) {
    out <- snp.info[chr(snp.info)==i,]
  } else {
    if(is.numeric(i)) {
      out <- snp.info[match(chr(snp.info),chrNames(snp.info))==i,]
    } else {
      stop("unknown index")
    }
  }
  out@build <- x@build
  out@chip <- x@chip
  return(out)
} )


setGeneric("convTo37", function(x) standardGeneric("convTo37"))
          
setMethod("convTo37", "ChipInfo", function(x) {
  if(ucsc.sanitizer(build(x))=="hg18") {
    u <- conv.36.37(ranged=as(x,"GRanges"))
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
      stop("conversion to build37/hg19 failed") 
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

setMethod("convTo36", "ChipInfo", function(x) {
  if(ucsc.sanitizer(build(x))=="hg19") {
    u <- conv.37.36(ranged=as(x,"GRanges"))
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

setMethod("show", "ChipInfo", 
     function(object) { showChipInfo(object,up.to=10,head.tail=5,show.strand=FALSE) } )

setMethod("print", "ChipInfo", 
          function(x,...) { showChipInfo(x,...) } )


#constructor
ChipInfo <- function(GRanges=NULL, chr=NULL, pos=NULL, ids=NULL, chip="unknown chip", build="",
                     rs.id=NULL, A1=NULL, A2=NULL, QCcode=NULL) {
  if(build!="") { build <- ucsc.sanitizer(build) }
  LL <- max(length(chr),length(GRanges),na.rm=T)
  if(length(A1)!=LL | length(A2)!=LL) { A1 <- A2 <- rep(NA,times=LL) }
  if(length(rs.id)!=LL) { rs.id <- rep(NA,times=LL) }
  if(length(QCcode)!=LL) { QCcode <- rep(0,LL) }
  if(is.null(GRanges)) {
    GRanges <- make.chr.pos.granges(chr=chr,pos=pos,row.names=ids)
  } else {
    if(is(GRanges)[1]!="GRanges") { GRanges <- as(GRanges,"GRanges") }
  }
  df <- DataFrame(A1=A1,A2=A2,QCcode=QCcode,rs.id=rs.id)
  #print(build)
  return(new("ChipInfo", seqnames=GRanges@seqnames, ranges=GRanges@ranges,  strand=GRanges@strand,
            elementMetadata=df, seqinfo=GRanges@seqinfo,
             chip=chip, build=build))
}


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
        ChipInfo("ChipInfo",GRanges=from,chip="unknown chip",build=build,rs.id=NULL,A1=NULL,A2=NULL)
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
