
######################
## Ranged Functions ##
######################



#' Wrapper to construct GRanges object from chr,pos or chr,start,end
#' 
#' Slightly simplifies the creation of a GRanges object, allowing flexible input of
#' chr, pos, or chr,start,end, and specification of rownames and the 'genome' parameter
#' for specifying the build/coordinate type, e.g, hg18, build 37, etc. Designed for
#' a simplified GRanges object without metadata, and where the 'strand' data is of
#' no interest, so if strand/metadata is to be used, use the original GRanges() constructor.
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions for the GRanges object
#' @param pos integer/numeric, for SNPs, can enter positions just once in 'pos' instead of entering the same value
#' for start and end
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param start integer/numeric, specify the start position of ranges to encode in the new
#' GRanges object, alongside 'end' (do not use 'pos' if using start+end)
#' @param end integer/numeric, specify the end position of ranges to encode in the new
#' GRanges object, alongside 'start' (do not use 'pos' if using start+end)
#' @return Returns a GRanges object with the ranges, build and rownames specified. Rownames
#' will be 1:nrow if the 'row.names' parameter is empty. The strand information will default
#' to '+' for all entries, and the metadata will be empty (this function is only for creation
#' of a very basic GRanges object).
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.gene, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' g1 <- make.granges(chr=c(1,4,"X"),pos=c(132432,434342,232222))
#' g2 <- make.granges(chr=c(22,21,21),start=c(1,1,1),end=c(1000,10000,100000),row.names=c("1K","10K","100K"))
#' g1 ; g2
make.granges <- function(chr,pos=NULL,start=NULL,end=NULL,row.names=NULL,build=NULL) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is.null(start) & is.null(end) & !is.null(pos)) {
    dF <- cbind(paste(chr),round(as.numeric(pos)))
  } else {
    if(!is.null(start) & !is.null(end) & is.null(pos)) {
      pos <- cbind(round(as.numeric(start)),round(as.numeric(end)))
      dF <- cbind(paste(chr),pos) 
    } else {
      stop("must use either 'pos' or 'start' and 'end'")
    }
  }
  #prv(paste(row.names)); prv(dF)
  if(is.character(row.names)) { if(length(row.names)==nrow(dF)) { rownames(dF) <- row.names } else { warning("row.names had an incorrect length")} }
  if(!any(Dim(pos)==length(chr))) { stop("chr and pos must be of the same length") }
  if(length(Dim(pos))>1) { 
    if(!ncol(pos) %in% c(1,2)) { 
      stop("pos must be a vector of SNP locations, or a 2-column object with start and end coordinates") 
    } else {
      if(ncol(pos)==2) {
        if(any(pos[,2]<pos[,1])) { warning("end coordinates should be equal or greater than start coordinates") }
      }
    }
  }
  if(ncol(dF)==3) { 
    colnames(dF) <- c("chr","start","end") } else { colnames(dF) <- c("chr","pos") }
  #return(dF)
  ranged <- data.frame.to.granges(dF,start=colnames(dF)[2],end=tail(colnames(dF),1),build=build) 
  if(!any(rownames(ranged) %in% row.names)) {
    if(is.character(row.names)){ if(length(row.names)==nrow(ranged)) { rownames(ranged) <- row.names }}
  }
  return(ranged)
}




#' Convert from build 37 to build 36 SNP coordinates
#' 
#' Convert range or SNP coordinates between builds using a chain file. Depending on the chain file
#' this can do any conversion, but the default will use the hg19 to hg18 (37-->36) chain file
#' built into this package. The positions to convert can be entered using using chr, pos vectors,
#'  or a RangedData or GRanges object. This function is a wrapper for liftOver() from rtracklayer,
#' providing more control of input and output and 'defensive' preservation of order and length
#' of the output versus the input ranges/SNPs. 
#' @param chr character, an optional vector of chromosomes to combine with 'pos' to describe
#'  positions to convert from build hg19 to hg18
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' a ranges object if this is provided along with 'chr'
#' @param ranges optional GRanges or RangedData object describing positions for which conversion
#' should be performed. No need to enter chr, pos if using ranges
#' @param ids if the ranges have ids (e.g, SNP ids, CNV ids), then by including this parameter
#' when using chr, pos input, the output object will have these ids as rownames. For ranges input
#' these ids would already be in the rownames of the GRanges or RangedData object, so use of
#' this parameter should be unnecessary
#' @param ... additional arguments to make.granges(), so in other words, can use 'start' and
#' 'end' to specify ranges instead of 'pos'.
#' @return Returns positions converted from build 37 to 36. If using the 'ranges' parameter 
#' for position input, the object returned will be of the same format. If using chr and pos 
#' to input, then the object returned will be a data.frame with columns, chr and pos with 
#' rownames 'ids'. Output will be the same length as the input, which is not necessarily the
#'  case for liftOver() which does the core part of this conversion. Using vector or GRanges 
#'  input will give a resulting data.frame or GRanges object respectively that has the same
#'  order of rownames as the original input. Using RangedData will result in an output that
#'   is sorted by genome order, regardless of the original order.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso conv.36.37, convTo37, convTo36
#' @examples
#' gene.labs <- c("CTLA4","IL2RA","HLA-C")
#' pp <- Pos.gene(gene.labs,build=37)
#' gg <- GRanges(ranges=IRanges(start=pp$start,end=pp$end),seqnames=pp$chr)
#' conv.37.36(gg) # order of output is preserved   ### HERE!!! ###
#' rr <- as(gg,"RangedData")
#' conv.37.36(rr) # note the result is same as GRanges, but in genome order
conv.37.36 <- function(ranges=NULL,chr=NULL,pos=NULL,...,ids=NULL) {
  chain.file <- "~/github/iChip/hg19ToHg18.over.chain"
  return(conv.36.37(ranges=ranges,chr=chr,pos=pos,...,ids=ids,chain.file=chain.file))
}



#' Convert from build 36 to build 37 SNP coordinates
#' 
#' Convert range or SNP coordinates between builds using a chain file. Depending on the chain file
#' this can do any conversion, but the default will use the hg18 to hg19 (36-->37) chain file
#' built into this package. The positions to convert can be entered using using chr, pos vectors,
#'  or a RangedData or GRanges object. This function is a wrapper for liftOver() from rtracklayer,
#' providing more control of input and output and 'defensive' preservation of order and length
#' of the output versus the input ranges/SNPs. 
#' @param chr character, an optional vector of chromosomes to combine with 'pos' to describe
#'  positions to convert to an alternative build
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' a ranges object if this is provided along with 'chr'
#' @param ranges optional GRanges or RangedData object describing positions for which conversion
#' should be performed. No need to enter chr, pos if using ranges
#' @param ids if the ranges have ids (e.g, SNP ids, CNV ids), then by including this parameter
#' when using chr, pos input, the output object will have these ids as rownames. For ranges input
#' these ids would already be in the rownames of the GRanges or RangedData object, so use of
#' this parameter should be unnecessary
#' @param chain.file character, a file location for the liftOver chain file to use for the
#' conversion. If this argument is left blank the default UCSC file that converts from hg18
#' to hg19 will be used. However chain files for other conversions are available from
#' http://crossmap.sourceforge.net/, and you could also customize these or create your own.
#' So this function can be used for conversion between any in-out build combination, using
#' this argument, not just 36--<37.
#' @param ... additional arguments to make.granges(), so in other words, can use 'start' and
#' 'end' to specify ranges instead of 'pos'.
#' @return Returns positions converted from build 36 to 37 (or equivalent for alternative chain 
#' files). If using the 'ranges' parameter for position input, the object returned will be of
#' the same format. If using chr and pos to input, then the object returned will be a data.frame
#' with columns, chr and pos with rownames 'ids'. Output will be the same length as the input,
#' which is not necessarily the case for liftOver() which does the core part of this conversion.
#' Using vector or GRanges input will give a resulting data.frame or GRanges object respectively
#' that has the same order of rownames as the original input. Using RangedData will result in an
#' output that is sorted by genome order, regardless of the original order. If ranges has no
#' rownames, or if 'ids' is blank when using chr, pos, ids of the form rngXXXX will be generated
#' in order to preserve the original ordering of locations.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso conv.37.36, convTo37, convTo36
#' @references http://crossmap.sourceforge.net/
#' @examples
#' # various chain files downloadable from http://crossmap.sourceforge.net/ #
#' gene.labs <- c("CTLA4","IL2RA","HLA-C")
#' snp.ids <- c("rs3842724","rs9729550","rs1815606","rs114582555","rs1240708","rs6603785")
#' pp <- Pos(snp.ids); cc <- Chr(snp.ids)
#' conv.36.37(chr=cc,pos=pp,ids=snp.ids)
#' pp <- Pos(gene.labs)
#' gg <- GRanges(ranges=IRanges(start=pp$start,end=pp$end),seqnames=pp$chr)
#' conv.36.37(gg) # order of output is preserved
#' rr <- as(gg,"RangedData")
#' conv.36.37(rr) # note the result is same as GRanges, but in genome order
conv.36.37 <- function(ranges=NULL,chr=NULL,pos=NULL,...,ids=NULL,chain.file="~/github/iChip/hg18ToHg19.over.chain") {
  require(GenomicRanges); require(rtracklayer); require(genoset)
  if(!file.exists(chain.file)) { stop("couldn't find chain file: ",chain.file) }
  chn <- import.chain(chain.file)
  #toranged <- F
  outType <- is(ranges)[1]
     
  if(!is.null(chr) & (!is.null(pos) | all(c("start","end") %in% names(list(...))))) {
    if(is.null(pos) & length(chr)==1) {
      if(length(list(...)$start)>1) {
        warning("when using start/end, 'chr' must have the same length as 'start'") 
      }
    }
    if(is.null(ids)) { ids <- paste0("rng",1:(max(length(chr),length(pos)))) } 
    ranges <- make.granges(chr=chr,pos=pos,row.names=ids,...)
    orn <- ids
  } else {
    if(is.null(rownames(ranges))) { rownames(ranges) <- paste0("rng",1:nrow(ranges)) }    
    orn <- rownames(ranges)
  }
  # return(ranges)
  if(is(ranges)[1]=="RangedData") { ranges <- as(ranges, "GRanges") }
  if(is(ranges)[1] %in% c("RangedData","GRanges")) {
    wd <- width(ranges)
    if(all(wd==1)) { SNPs <- TRUE } else { SNPs <- FALSE }
    mcols(ranges)[["XMYINDEXX"]] <- rownames(ranges)
    mcols(ranges)[["XMYCHRXX"]] <- ocr <- chr2(ranges)
    #prv(orn,ocr)
    opos <- start(ranges)
    ranges <- set.chr.to.char(ranges)
    #print(head(ranged))
    ranged.gr <- ranges # as(ranges,"GRanges"); #toranged <- T
  } else {
    stop("input specified resulted in an invalid GRanges/RangedData 'ranged' object, type ",is(ranges)[1]) 
  } 
  # change CHR-XY to CHR-X prior to liftOver, then change back #
  xy.ind <- grep("XY",seqnames(ranged.gr))
  if(length(xy.ind)>0) {
    found.xy <- TRUE
    xy.id <- rownames(ranged.gr)[xy.ind]
    if(!"chrX" %in% seqlevels(ranged.gr)) { seqlevels(ranged.gr) <- c(seqlevels(ranged.gr),"chrX") }
    seqnames(ranged.gr)[xy.ind] <- "chrX"
  } else { found.xy <- FALSE }
  ranged.gr.37 <- liftOver(ranged.gr,chn)
  myfun <- function(x) { 
    data.frame(start=minna(start(x)),end=maxna(end(x))) 
  }
  if(!SNPs ) {
    new.coords.df <- do.call("rbind",lapply(ranged.gr.37,myfun))
    ranged.gr.37<-ranged.gr
    ranges(ranged.gr.37)<-with(new.coords.df,IRanges(start=start,end=end))
    #seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
    seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
    out <- ranged.gr.37
  } else {
    seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
    #seqnames(ranged.gr.37)<-gsub("chr","",seqnames(ranged.gr.37))
    out <- as(ranged.gr.37,"IRangesList")
    #seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
    #new.coords.df <- as.data.frame(ranged.gr.37)
  }
  # seqlevels(ranged.gr.37)<-gsub("chr","",seqlevels(ranged.gr.37))
  # out <- as(ranged.gr.37,"IRangesList")
  out <- as(out,"RangedData")
  #return(ranged.gr.37)
  #ranged.gr.37 <- set.chr.to.numeric(ranged.gr.37)
  #if(!toranged | T) { return(ranged.gr.37) }
  ranged.gr.37 <- out #toGenomeOrder2(out)
  #return(ranged.gr.37)
  if(all(c("XMYINDEXX","XMYCHRXX") %in% colnames(ranged.gr.37))) {
    RN <- ranged.gr.37[["XMYINDEXX"]]
    nr <- nrow(ranged.gr.37)
    MAXDISPLAY <- 50
    if(length(orn)>length(RN)) { 
      cat("conversion failed for",length(orn[!orn %in% RN]),"rows, original positions kept:\n") ;  
      failz <- orn[!orn %in% RN]
      cat(comma(head(failz,MAXDISPLAY))) 
      if(length(failz)>MAXDISPLAY) { cat(", ... and",length(failz)-MAXDISPLAY,"more\n")  } else { cat("\n") }
      ln <- orn[!orn %in% RN]
      #return(ranges)
      newchr <- gsub("chr","",chr2(ranges[match(ln,ranges$XMYINDEXX),]))
      noopos <- start(ranges[match(ln,ranges$XMYINDEXX),])
      hcc <- hard.coded.conv()
      ifNAthen0 <- function(X) { X[is.na(X)] <- 0; return(X) }
      h36 <- which(noopos %in% hcc$pos36 & newchr==ifNAthen0(hcc$chr[match(noopos,hcc$pos36)])); l36 <- length(h36)
      h37 <- which(noopos %in% hcc$pos37 & newchr==ifNAthen0(hcc$chr[match(noopos,hcc$pos37)])); l37 <- length(h37)
      if(l36>=l37 & l36>0) {
        noopos[h36] <- hcc$pos37[match(noopos[h36],hcc$pos36)]
        cat("found",l36,"of the missing SNP hg18-hg19 lookups in an internal table\n")
      } else {
        if(l36<l37 & l37>0) {
          noopos[h37] <- hcc$pos36[match(noopos[h37],hcc$pos37)]
          cat("found",l37,"of the missing SNP hg19-hg18 lookups in an internal table\n")
        } else {
          ## no matches to extras table
        }
      }
      extra <- data.frame(Chr=newchr,Start=noopos,End=noopos)
      rownames(extra) <- ln
    } #else { cat("length is already the same\n") }
    Ind <- match(ranged.gr.37[["XMYINDEXX"]],orn)
    out <- data.frame(Chr=ranged.gr.37[["XMYCHRXX"]],Start=start(ranged.gr.37),End=end(ranged.gr.37),ind=Ind)
    rownames(out) <- RN
    if(length(orn)>length(RN)) {
      #prv(out,extra)
      out <- out[,-4] # 4 is the 'ind' column
      out <- rbind(out,extra)
      out <- out[orn,]
    } #else { cat("length is now the same\n") }
    #return(out) 
  } else { warning("missing key columns for chr, snp-name")  }
  #print(outType)
  #return(out)
  #prv(out)
  ranged.rd <- toGenomeOrder2(data.frame.to.ranged(out))
  #print(colnames(ranged.rd))
  ranged.gr.37 <- as(ranged.rd,"GRanges")
  #print(colnames(ranged.gr.37))
  if(found.xy) {
    xy.ind <- match(xy.id,rownames(ranged.gr.37))
    lmis <- length(which(is.na(xy.ind)))
    if(lmis>0) { warning('liftOver function removed ",lmis," chrX/chrY ranges'); xy.ind <- narm(xy.ind) }
    if(!"XY" %in% seqlevels(ranged.gr.37)) { seqlevels(ranged.gr.37) <- c(seqlevels(ranged.gr.37),"XY") }
    seqnames(ranged.gr.37)[xy.ind] <- "XY"
  }
  if(outType=="GRanges") { 
    #return(ranged.gr.37)
    cn37 <- colnames(mcols(ranged.gr.37))
    if("ind" %in% cn37) { 
      mind <- as.numeric(mcols(ranged.gr.37)[["ind"]])
      #prv(ranged.gr.37,mind)
      ranged.gr.37 <- ranged.gr.37[order(mind),]
      mcols(ranged.gr.37) <- mcols(ranged.gr.37)[,-which(cn37 %in% "ind")] 
    } # else { warning("couldn't find index column, GRanges object not sorted in original order") }
    if(all(rownames(ranged.gr.37) %in% orn)) { ranged.gr.37 <- ranged.gr.37[orn,] }
    return(ranged.gr.37)
  } else {
    if(outType=="RangedData") {
      #if("ind" %in% colnames(ranged.gr.37)) { ranged.gr.37 <- ranged.gr.37[,-which(colnames(ranged.gr.37) %in% "ind")] }
      return(toGenomeOrder2(as(ranged.gr.37,"RangedData")))
    } else {
      #prv(ranged.gr.37)
      out <- ranged.to.data.frame(ranged.gr.37,include.cols=FALSE,use.names=TRUE)
      cn37 <- colnames(mcols(ranged.gr.37))
      if("ind" %in% cn37) { 
        mind <- as.numeric(mcols(ranged.gr.37)[["ind"]])
        out <- out[order(mind),,drop=FALSE]
        if("ind" %in% colnames(out)) { out <- out[,-which(colnames(out) %in% "ind")] }
      }
      if(is.null(dim(out))) { dim(out) <- c(length(out)/3,3) }
      if(all(rownames(out) %in% orn)) { out <- out[orn,] }
      return(out)
    }
  }
}




#' Extend an interval or SNP by distance in centimorgans (recombination distance)
#' 
#' It is straightforward to extend a genomic interval or position by a number of basepairs, or
#' a percentage, but extending by recombination units of centimorgans is more involved, requiring
#' annotation lookup. This function streamlines this process.
#' This function makes use of recombination rate hapmap reference files to calculate 
#' recombination distances for genome locations, in centimorgans. For a given position
#' (or vector), a window can be returned of a given extension on either side of the position,
#' for instance, 1 centimorgan to the left, and to the right of a SNP, giving a 2 centimorgan
#' range as a result.
#' @param ranges optional GRanges or RangedData object describing positions for which we want to
#' generate windows, removing the need to enter chr, start and end
#' @param chr character, an optional vector of chromosomes to combine with 'start' and 'end'
#'  to describe positions for which to generate recombination windows
#' @param start integer, an vector of start points for chromosome ranges
#' @param end integer, an vector of end points for chromosome ranges
#' @param window numeric, number of centimorgans to extend the window either side of
#' the range or location (can be a fraction)
#' @param bp.ext numeric, optional number of base-pairs to extend the window by in addition
#' to the centimorgan extension
#' @param rec.map recombination map object (list of 22 data.frames) generated using 
#' 'get.recombination.map()'; if you are performing many of these operations, loading this 
#' object into your workspace and passing it on to this function will save loading it each 
#' time, and provide a speed advantage. Only use an object generated by get.recombination.map(),
#'  as otherwise the results will almost certainly be meaningless.
#' @param info logical, whether to display the derived window size and number of hapmap SNPs within
#' the window for each window derived
#' @seealso get.recombination.map, get.nearby.snp.lists, exp.window.nsnp
#' @author Chris Wallace and Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @export
#' @examples
#' # not run, as initial download of the recombination map takes nearly a minute #
#' recwindow(chr=11,start=10000000,end=10000000,window=1,bp.ext=10000)
#' rd <- RangedData(ranges=IRanges(start=c(1.5,10.1)*10^7, end=c(1.55,10.1)*10^7),space=c(2,10))
#' rd # show original data
#' recwindow(rd) # now extended by the interval
#' recwindow(as(rd,"GRanges"),info=FALSE) # also works for GRanges
recwindow <- function(ranges=NULL,chr=NA,start=NA,end=start,window=0.1,bp.ext=0, rec.map=NULL, info=TRUE) {
  #www <- window ; ccc <- chr; sss <- start; prv(ccc,sss,www,bp.ext)
  if(!is.numeric(bp.ext)) { warning("bp.ext must be numeric, setting to zero"); bp.ext <- 0 }
  if(!is.numeric(window)) { warning("window must be numeric, setting to 0.1 centimorgans"); window <- 0.1 }
  if(!all(is.na(chr))) { if(any(!paste(chr) %in% paste(1:22))) { 
    stop("this function only works for autosomes 1-22 [e.g, no X,Y or formatting like 'chr2', etc]") } }
  chr <- as.numeric(chr)
  typ <- is(ranges)[1]
  if(typ %in% c("RangedData","GRanges")) { 
    if(typ=="GRanges") { ranges <- as(ranges,"RangedData") }
    ranges <- toGenomeOrder2(ranges,strict=T)
    ss <- start(ranges); ee <- end(ranges); cc <- chr2(ranges)
    out <- recwindow(chr=cc,start=ss,end=ee,window=window,bp.ext=bp.ext,info=info,rec.map=rec.map)
    if(length(out)==2) { out <- as.matrix(out); dim(out) <- c(1,2) } 
    outData <- RangedData(ranges=IRanges(start=out[,1],end=out[,2],names=rownames(ranges)),space=cc)
    outData <- toGenomeOrder2(outData,strict=TRUE)
    if(ncol(ranges)>0) {
      for (zz in 1:ncol(ranges)) { 
        if(typ=="GRanges") {
          outData[[colnames(ranges)[zz]]] <- mcols(ranges)[,colnames(ranges)[zz]]
        } else {
          outData[[colnames(ranges)[zz]]] <- ranges[[colnames(ranges)[zz]]]
        }
      }
    }
    if(is(ranges)[1]=="GRanges") { outData <- as(toGenomeOrder2(outData,strict=T),"GRanges") }
    return(outData)
  } else {
    if(all(!is.na(chr)) & all(!is.na(start)) & all(!is.na(end))) {
      if(length(chr)==length(start) & length(start)==length(end)) {
        if(length(chr)>1) {
          # run for a vector
          out <- matrix(ncol=2,nrow=length(chr)); colnames(out) <- c("start","end")
          for (dd in 1:length(chr)) {
            out[dd,] <- recwindow(chr=chr[dd],start=start[dd],end=end[dd],
                                  window=window,bp.ext=bp.ext,info=info,rec.map=rec.map)
          }
          return(out)
        } else {
          ## continue as normal, just a single coordinate/range to process
        }
      } else {
        stop("invalid input, start, end and chr need to be the same length")
      }
    } else {
      stop("invalid input, either use a RangedData object, or else chr, start and end")
    }
  }
  
  #rate.fn <- sprintf("/dunwich/scratch/chrisw/HapMap/rates_rel22/genetic_map_chr%s_b36.txt.gz",chr)
  #print(rate.fn)
  #rates <- read.table(gzfile(rate.fn),header=TRUE)
  if(is.list(rec.map)) { if(length(rec.map)==22) { rates <- rec.map[[chr]] } } else {
    if(is.null(rec.map)) { rates <- get.recombination.map()[[chr]] } else {
      if(is.character(rec.map)) { rates <- get.recombination.map(dir=rec.map)[[chr]] } else {
        stop("invalid value for rec.map entered")
      }
    }
  }
  cm.st <- rates[which.min(abs(rates$position-start)),3]
  cm.en <- rates[which.min(abs(rates$position-end)),3]
  
  mx <- max(window,1)
  kk <- rates[which.min(abs(rates[,3]-(cm.st-window))) : which.min(abs(rates[,3]-(cm.en+window))),]
  if(info) { cat("n hapmap snps in window =",nrow(kk),"\n") }
  from <- min(kk[,1])
  to <- max(kk[,1])
  ##
  if(info) {
    cat("new window size is\nleft: ",(start-from+bp.ext)/1000,"kb\tright: ",
      (to-end+bp.ext)/1000,"kb\ttotal: ",(to-from+(2*bp.ext))/1000,"kb\n",sep="")
  }
  if(info & bp.ext>0) { cat("in addition to cM distance, window was extended by",
                     bp.ext,"base pairs on either side\n")} 
  from <- max(c(0,(from-bp.ext)))
  to <- min(c((to+bp.ext),get.chr.lens()[chr][1]),na.rm=T)
  return(c(from,to))
}




#' Convert GRanges/RangedData to chr:pos1-pos2 vector
#' 
#' Takes a RangedData or GRanged object from some annotation lookup functions and converts to standard text
#' positions, such as what you might see on the UCSC genome browser, such as 
#' chr1:10,000,234-11,000,567 for a range, or chrX:234,432 for a SNP. Useful for printing
#' messages, concatenating positions to a single vector, or creating queries for databastes.
#' @param ranges A RangedData or GRanges object
#' @export
#' @return a text vector of the same length as 'ranges' with notation as described above
#' representing each position in the 'ranges' object
#' @seealso convert.textpos.to.data
#' @examples
#' Ranges.to.txt(rranges())
Ranges.to.txt <- function(ranges) {
  if(!is(ranges)[1] %in% c("RangedData","GRanges")) { stop("Not a GRanges or RangedData object") }
  text.out.a <- paste0("chr",chr2(ranges),":",format(start(ranges),scientific=F,trim=T))
  text.out.b <- paste0("-",format(end(ranges),scientific=F,trim=T))
  text.out.b[start(ranges)==end(ranges)] <- ""
  text.out <- paste0(text.out.a,text.out.b)
  return(text.out)
}






#' Select ranges only within the 22 autosomes in a ranged data object
#' 
#' Select only data from autosomes from a GRanges/RangedData object.
#' Will exclude X,Y, mitochondrial chromosome rows, and can automatically
#' detect whether chromosomes are coded as 'chr1' or just '1', etc.
#' @param ranges A RangedData or GRanges object
#' @export
#' @return an object of the same format as the input (ranges), except
#' with non-autosomal ranges removed.
#' @examples
#' rand.ranges <- rranges(chr.range=20:26)
#' rand.ranges # should include some non-autosomes
#' select.autosomes(rand.ranges) # only autosomes remain
select.autosomes <- function(ranges,deselect=FALSE) {
  #must.use.package("genoset",bioC=T)
  typ <- is(ranges)[1]
  if(!typ %in% c("RangedData","GRanges")) { 
    warning("not a RangedData or GRanges object"); return(ranges) 
  }
  if(length(unique(chr2(ranges))) < length(levels(chr2(ranges)))) {
    # this fixes the problem when a subset of a ranges object with less 
    #  chromosomes still has empty chr slots from previous object
    if(typ=="RangedData") { ranges <- ranges[as.numeric(unique(chr2(ranges)))] }
  } #else { cat("ok\n") }
  Chrz <- (rownames(chrInfo2(ranges)))
  chrz <- tolower(paste(Chrz))
  if(length(grep("chr",chrz))>0) {
    select1 <- which(chrz %in% paste("chr",1:22,sep=""))
    select2 <- which(chrz %in% paste(1:22))
    if(length(select2)>length(select1)) { select <- select2 } else { select <- select1 }
  } else {
    select <- which(chrz %in% paste(1:22))
  }
  if(deselect) {
    ok.chrs <- Chrz[!select]
  } else {
    ok.chrs <- Chrz[select]
  }
  if(typ=="RangedData") {
    return(ranges[ok.chrs])
  } else {
    return(ranges[chr2(ranges) %in% ok.chrs,])
  }
}






#' Convert RangedData/GRanges to a data.frame
#' 
#' Convert a RangedData/GRanges object to a data.frame with columns
#' chr, start and end. Default is to only translate the chromosome and
#' position information, which is faster. Using 'include.cols'=TRUE
#' allows all the columns from 'ranged' to be taken across to the resulting
#' data.frame.
#' @param ranged A RangedData or GRanges object
#' @param include.cols logical, whether to also bring across non-positional
#' columns to the resulting data.frame
#' @param use.names logical, whether to keep the rownames from the
#' original object for the output. Only has an effect when include.cols=FALSE,
#' otherwise original rownames are always kept.
#' @export
#' @seealso data.frame.to.ranged, data.frame.to.granges
#' @return A data.frame with columns chr, start and end, and depending on
#' chosen parameters, the same rownames as the input, and optionally the
#' same additional columns.
#' @examples
#' rd <- rranges(9,GRanges=FALSE, fakeids=TRUE)
#' rd[["fakecol"]] <- sample(nrow(rd))
#' rd[["rs.id"]] <- paste0("rs",sample(10000,9))
#' ranged.to.data.frame(rd)
#' ranged.to.data.frame(rd,,FALSE)
#' ranged.to.data.frame(rd,TRUE) # keep all the columns
#' data.frame.to.granges(ranged.to.data.frame(rd,TRUE)) # inverse returns original
ranged.to.data.frame <- function(ranged,include.cols=FALSE,use.names=TRUE) {
  if(!include.cols) {
    u <- Ranges.to.txt(ranged)
    v <- convert.textpos.to.data(u)
    if(!is.null(rownames(ranged)) & nrow(ranged)==nrow(v) & use.names) { rownames(v) <- rownames(ranged) }
    return(v)
  } else {
    u <- as(ranged,"data.frame")
    cn <- tolower(colnames(u))
    if(is(ranged)[1]=="RangedData") {
      if("names" %in% cn) { 
        rownames(u) <- u[["names"]]
        u <- u[,-which(cn=="names")]
      } 
      if("space" %in% cn) { colnames(u)[which(cn=="space")] <- "chr" }
    } else {
      if(is(ranged)[1]=="GRanges") {
        if("seqnames" %in% cn) { colnames(u)[which(cn=="seqnames")] <- "chr" }
      } else {
        warning("'ranged' should be RangedData or GRanges, coercion could fail")
      }
    }
    return(u)
  }
}

#' Convert a data.frame with positional information to GRanges
#' 
#' Convert a data.frame containing chromosome and position information
#' to a GRanges object. Assumes the position information is contained in
#' columns named 'chr', 'start' and 'end' respectively (not case sensitive) 
#' although you can enter alternative column names for each as parameters. 
#' 'seqnames' will be automatically detected as an alternative to 'chr' if 
#' present. Column names that are default GRanges slot names such as 'seqnames',
#' 'ranges', 'strand', 'seqlevels', etc, will be removed during conversion, so
#' rename these if you want them to be translated into the resulting GRanges
#' objects' column metadata. If there is a column 'pos' but no columns 'start'
#'  and 'end' this will be detected automatically without needing to change
#'  the default parameters and start will equal end equals pos (ie., SNPs).
#' @param dat a data.frame with chromosome and position information 
#' @param ... additional arguments to data.frame.to.ranged(), namely:
#' ids, start, end, width, chr, exclude and build
#' @export
#' @seealso ranged.to.data.frame, data.frame.to.ranged
#' @return A RangedData or GRanges object. If 'dat' doesn't
#' use the default column names, specify these using parameters
#' ids, start, and end or width. Exclude will remove prevent any 
#' column names of 'dat' specified not to be translated to the 
#' returned GRanges object. 'build' specifies the 'genome'
#' slot of the resulting object. 'ids' allows specification of
#' a column to be converted to the rownames of the new object.
#' @examples
#' chr <- sample(1:22,10)
#' start <- end <- sample(1000000,10)
#' df1 <- cbind(chr,start,end)
#' data.frame.to.granges(df1) # basic conversion
#' width <- rep(0,10)
#' df2 <- cbind(chr,start,width)
#' data.frame.to.granges(df2,end=NULL,width="width") # define ranges with start and width
#' id.col <- paste0("ID",1:10)
#' rs.id <- paste0("rs",sample(10000,10))
#' df3 <- cbind(chr,start,end,id.col,rs.id)
#' data.frame.to.granges(df3) # additional columns kept
#' df4 <- cbind(chr,start,end,id.col,rs.id, ranges=1:10)
#' data.frame.to.granges(df4) # 'ranges' column excluded as illegal name
#' data.frame.to.granges(df4, exclude="rs.id") # manually exclude column
#' df5 <- cbind(chr,start,end,rs.id)
#' rownames(df5) <- paste0("ID",1:10)
#' data.frame.to.granges(df5) # rownames are kept
#' data.frame.to.granges(df4,ids="id.col") # use column of 'dat' for rownames
data.frame.to.granges <- function(dat,...) {
  return(data.frame.to.ranged(dat=dat,...,GRanges=TRUE))
}


#' Convert a data.frame with positional information to RangedData/GRanges
#' 
#' Convert a data.frame containing chromosome and position information
#' to a RangedData or GRanges object. Assumes the position information is contained in
#' columns named 'chr', 'start' and 'end' respectively (not case sensitive) 
#' although you can enter alternative column names for each as parameters. 
#' 'seqnames' will be automatically detected as an alternative to 'chr' if 
#' present. If there is a column 'pos' but no columns 'start' and 'end' this
#' will be detected automatically without needing to change the default parameters
#' and start will equal end equals pos (ie., SNPs). Column names that are default 
#' GRanges slot names such as 'seqnames', 'ranges', 'strand', 'seqlevels', etc, will
#' be removed during conversion, so rename these if you want them to be translated 
#' into the resulting object.
#' @param dat a data.frame with chromosome and position information 
#' @param ids character string, an optional column name containing ids which
#' will be used for rownames in the new object, as long as the ids are unique.
#' If not, this option is overridden and the ids will simply be a normal column
#' in the new object.
#' @param start character, the name of a column in the data.frame contain
#' the start point of each range. Not case sensitive. In the case of SNP
#' data, a column called 'pos' will also be automatically detected without
#' modifying 'start' or 'end', and will be used for both start and end.
#' @param end character, the name of a column in the data.frame containing the
#' end point of each range, can also use 'width' as an alternative specifier,
#' in which case 'end' should be set to NULL. Not case sensitive. In the case of SNP
#' data, a column called 'pos' will also be automatically detected without
#' modifying 'start' or 'end', and will be used for both start and end.
#' @param width the name of a column in the data.frame containing 'width' of
#' ranges, e.g, SNPs would be width=0. This is optional, with 'start' and 'end'
#' being the default way to specify an interval. If using 'width' you must
#' also set 'end' to NULL.  Not case sensitive.
#' @param chr character, the name of the column in the data.frame containing
#' chromosome values. The default is 'chr' but 'seqnames' will also be
#' detected automatically even when chr='chr'. Not case sensitive.
#' @param exclude character string, and column names from the data.frame to 
#' NOT include in the resulting S4 object.
#' @param build the ucsc build for the result object which will apply to the
#' 'universe' (RangedData) or 'genome' slot (GRanges) of the new object.
#' @param GRanges logical, whether the resulting object should be GRanges (TRUE),
#' or RangedData (FALSE)
#' @export
#' @seealso ranged.to.data.frame, data.frame.to.ranged
#' @return A RangedData or GRanges object. If 'dat' doesn't use the default 
#' column names 'chr', 'start'/'end' or 'pos', specify these using parameters 
#' 'ids', 'start', and 'end' or 'width'. Exclude will remove prevent any 
#' column names of 'dat' specified not to be translated to the returned GRanges
#' object. 'build' specifies the 'genome' slot of the resulting object. 'ids' 
#' allows specification of a column to be converted to the rownames of the new object.
#' @examples
#' chr <- sample(1:22,10)
#' start <- end <- sample(1000000,10)
#' df1 <- cbind(CHR=chr,Start=start,enD=end)
#' print(df1)
#' data.frame.to.granges(df1) # not case sensitive!
#' width <- rep(0,10)
#' df2 <- cbind(chr,start,width)
#' data.frame.to.granges(df2,end=NULL,width="width") # define ranges with start and width
#' id.col <- paste0("ID",1:10)
#' rs.id <- paste0("rs",sample(10000,10))
#' df3 <- cbind(chr,start,end,id.col,rs.id)
#' data.frame.to.granges(df3) # additional columns kept
#' df4 <- cbind(chr,start,end,id.col,rs.id, ranges=1:10)
#' data.frame.to.granges(df4) # 'ranges' column excluded as illegal name
#' data.frame.to.granges(df4, exclude="rs.id") # manually exclude column
#' df5 <- cbind(chr,start,end,rs.id)
#' rownames(df5) <- paste0("ID",1:10)
#' data.frame.to.granges(df5) # rownames are kept
#' data.frame.to.granges(df4,ids="id.col") # use column of 'dat' for rownames
data.frame.to.ranged <- function(dat, ids=NULL,start="start",end="end",width=NULL,
                                 chr="chr",exclude=NULL,build=NULL,GRanges=FALSE) 
{
  ## abandon longer names as they clash with function names
  st <- paste(start); en <- paste(end); ch <- paste(chr); wd <- paste(width)
  if((!chr %in% colnames(dat)) & ("seqnames" %in% colnames(dat)) & GRanges) { ch <- "seqnames" }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  must.use.package(c("genoset","IRanges"),T)
  g.illegal <- tolower(c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
                         "isCircular", "start", "end", "width", "element"))
  if(is.matrix(dat)) { dat <- as.data.frame(dat,stringsAsFactors=FALSE) }
  if(!is.data.frame(dat)) { stop("Error: not a dataframe")}
  key.nms <- c(ids,st,en,ch,wd)
  tries <- 0
  #print(key.nms); print(colnames(dat))
  while(!all(key.nms %in% colnames(dat))) { 
    colnames(dat) <- tolower(colnames(dat)); key.nms <- tolower(key.nms)
    st <- tolower(st); en <- tolower(en); ch <- tolower(ch); wd <- tolower(wd)
    if(tries>2) {
      if((tolower(st)=="pos" | tolower(en)=="pos") & !(tolower(st)=="pos" & tolower(en)=="pos")) {
        st <- en <- "pos"
      } else {
        if(tolower(st)=="start" & tolower(en)=="end") { st <- en <- "pos" }
      }
    }
    key.nms <- c(ids,st,en,ch,wd)
    tries <- tries+1
    if(tries > 3) { if(!all(c(st,en,ch) %in% colnames(dat))) {
      warning("chromosome and position columns not found") } ; break }
  }
  if(!is.null(ids)) { 
    if(anyDuplicated(dat[[ids]])==0) { 
      id <- dat[[ids]] 
    } else { 
      key.nms <- key.nms[-match(ids,key.nms)] # allow non-unique ids as regular
      ids <- NULL
      warning("id must be unique to form rownames, will insert as a separate column") 
    }
  }
  if(is.null(ids)) { 
    if(!is.null(rownames(dat)) & all(rownames(dat)!=paste(1:nrow(dat)))) { 
      id <- rownames(dat)
    } else { 
      id <- paste(1:nrow(dat)) 
    }
  }
  ## not sure why here are adding 'chr' to X and Y?
  #this was here before? :  if(length(ch)>0) { ch1 <- gsub("Y","chrY",gsub("X","chrX",gsub("chr","",dat[[ch]],ignore.case=T))) } else { ch1 <- NULL }
  if(length(ch)>0) { ch1 <- gsub("chr","",dat[[ch]],ignore.case=T) } else { ch1 <- NULL }
  if(length(st)>0) { st1 <- as.numeric(dat[[st]]) } else { st1 <- NULL }
  if(length(en)>0) { en1 <- as.numeric(dat[[en]]) } else { en1 <- NULL }
  if(length(wd)>0) { en1 <- st1+as.numeric(dat[[wd]]) } # { en1 <- st1+dat[[wd]] }
  #print(length(st1)); print(length(en1)); print(length(id)); print(length(ch1))
  outData <- GRanges(ranges=IRanges(start=st1,end=en1,names=id),seqnames=ch1); genome(outData) <- build[1]
  #outData <- RangedData(ranges=IRanges(start=st1,end=en1,names=id),space=ch1,universe=build[1])
  ###  ###  ###  outData <- toGenomeOrder2(outData,strict=T)
  # note when adding data subsequently that 'RangedData' sorts by genome order, so need
  # to re-sort any new data before adding.
  if(is.null(rownames(outData))) { rownames(outData) <- paste(1:nrow(outData)) }
  reorder <- match(rownames(outData),id)
  more.cols <- colnames(dat)[!colnames(dat) %in% key.nms]
  more.cols <- more.cols[!more.cols %in% exclude]
  if(is(outData)[1]=="GRanges") { more.cols <- more.cols[!more.cols %in% g.illegal] }
  if(length(more.cols)>0) {
    for (cc in 1:length(more.cols)) {
      u <- dat[[more.cols[cc]]][reorder]; #prv(u)
      if(is(outData)[1]=="GRanges") {
        mcols(outData)[[more.cols[cc]]] <- u
      } else {
        outData[[more.cols[cc]]] <- u
      }
    }
  }
  if(GRanges) {
    return(as(outData,"GRanges"))
  } else {
    cncn <- colnames(mcols(outData))
    outData <- as(outData,"RangedData")
    if(any(cncn %in% "strand")) {
      outData <- outData[,-which(cncn=="strand")]
    }
    #outData <- toGenomeOrder2(outData,strict=T)
    return(outData)
  }
}



#' Select chromosome subset of GRanges or RangedData object
#' 
#' One of the main differences between RangedData and GRanges is the way
#' of selecting the subset for a chromosome. RangedData just uses [n] where
#' 'n' is the chromosome name or number. Whereas GRanges, does not have a
#' method like this, so need to select using [chr(X)==chr.num,]
#' This wrapper allows selection of a chromosome or chromosomes regardless of
#' whether the object is RangedData or GRanges type.
#' @param X A GRanges or RangedData object
#' @param chr Vector, the chromosome(s) (number(s) or name(s)) to select
#' @export
#' @return returns an object of the same type as X, with only the chromosome
#' subset specified.
#' @examples
#' some.ranges <- rranges(100,chr.range=1:10)
#' chr.sel(some.ranges,6)
#' more.ranges <- rranges(10, chr.range=21:25)
#' chr.sel(more.ranges,1:22) # gives warning
#' select.autosomes(more.ranges)
chr.sel <- function(X,chr) {
  typ <- is(X)[1]
  if(!typ %in% c("RangedData","GRanges","ChipInfo")) { stop("not a ChipInfo, GRanges or RangedData object") }
  if(nrow(X)==0) { warning("X has no ranges") ; return(X) }
  if(!(is.character(chr) | is.numeric(chr))) { stop("chr must be character or numeric type") }
  if(is.numeric(chr)) { if(!all(chr %in% 1:99)) { 
    stop("illegal chromosome index, valid range 1-99 [although 1-28 typical for human]") } }
  if(typ=="RangedData") { return(X[chr])}
  all.chr <- chr2(X)
  if(!all(chr %in% unique(all.chr))) { 
    if(!any(chr %in% unique(all.chr))) { 
      warning("none of the specified chromosome indices were present in the GRanges object, returning NULL")
      return(NULL)
    } else { 
      warning("some of the specified chromosome indices were not present in the GRanges object") 
    }
  }
  return(X[all.chr %in% chr,])
}


#' Simulate a GRanges or RangedData object
#' 
#' For testing purposes, this function will generate a S4 ranged object
#' based on the human genome. The default is to produce ranges selected
#' from chromosomes, with probability of a position in each chromosome
#' equal to the length of that chromosome versus the whole genome. The
#' maximum position allocated within each chromosome will be within
#' the length bounds of that chromosome. You can specify SNPs (ie., start
#' =end), but the default is for random ranges. You can alter the UCSC
#' build to base the chromosome lengths on, and you can specify whether
#' chromosomes should appear as chr1,chr2,... versus 1,2,..
#' @param n integer, number of rows to simulate
#' @param SNP logical, whether to simulate SNPs (width 1, when SNPs=TRUE)
#'  or just ranges (when SNP=FALSE)
#' @param chr.range integer vector of values from 1 to 26, to specify which
#' chromosomes to include in the simulated object. 23-26 are X,Y,XY,MT 
#' respectively.
#' @param chr.pref logical, if TRUE chromosomes will be coded as chr1,chr2,...,
#' versus 1,2,.. when chr.pref=FALSE
#' @param order logical, if TRUE the object returned will be in genomic order,
#' otherwise the order will be randomized
#' @param equal.prob logical, when FALSE (default), random positions will be
#' selected on chromosomes chosen randomly according to the their length (i.e,
#' assuming every point on the genome has equal probability of being chosen.
#' If equal.prob=TRUE, then chromosomes will be selected with equal probability,
#' so you could expect just as many MT (mitochondrial) entries as Chr1 entries.
#' @param GRanges logical, if TRUE the returned object will be GRanges format,
#' or if FALSE, then RangedData format
#' @param build character, to specify the UCSC version to use, which has a small
#' effect on the chromosome lengths. Use either "hg18" or "hg19". Will also 
#' accept build number, e.g, 36 or 37.
#' @param fakeids logical, whether to add rownames with random IDs (TRUE) or 
#' leave rownames blank (FALSE). If SNP=TRUE, then ids will be fake rs-ids.
#' @export
#' @return returns a ranged object (GRanges or RangedData) containing data
#' for 'n' simulated genomic ranges, such as SNPs or CNVs across chromosomes in
#' 'chr.range', using UCSC 'build'.
#' @examples
#' rranges()
#' rr <- rranges(SNP=TRUE,chr.pref=TRUE,fakeids=TRUE)
#' width(rr) # note all have width 1
#' rr
#' tt <- table(chr(rranges(1000)))
#' print(tt/sum(tt)) # shows frequencies at which the chr's were sampled
#' tt <- table(chr(rranges(1000,equal.prob=TRUE)))
#' print(tt/sum(tt)) # shows frequencies at which the chr's were sampled
rranges <- function(n=10,SNP=FALSE,chr.range=1:26,chr.pref=FALSE,order=TRUE,equal.prob=FALSE,
                    GRanges=TRUE,build=NULL, fakeids=FALSE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(!is.numeric(chr.range)) { stop("chr.range must be a numeric integer vector rangeing from 1 to 26") }
  chr.range <- unique(chr.range[chr.range<=26 & chr.range>=1]) 
  cL <- get.chr.lens(mito=TRUE,build=build)[c(1:24,24,25)]
  if(equal.prob) {
    cP <- rep(1/length(chr.range),length(chr.range))
  } else {
    cP <- cL[chr.range]/sum(cL[chr.range]) # probabilities of a location being in each chromosome
  }
  if(!is.numeric(n)) { stop("'n' must be a numeric integer vector, representing the number of rows to simulate") }
  nn <- round(force.scalar(n,min=1,max=10^9,default=10))
  chrs <- sort(sample(chr.range,size=nn,replace=TRUE,prob=cP))
  cnts <- table(chrs)
  chr.lengths <- get.chr.lens(mito=TRUE)[c(1:24,24,25)]
  dubb <- if(SNP) { 1 } else { 2 }
  randoms <- starts <- ends <- vector("list",length(chr.range))
  for (cc in 1:length(cnts)){
    ii <- match(names(cnts[cc]),paste(chr.range))
    if(!is.na(ii)) {
      randoms[[cc]] <- oo <- sort(sample(chr.lengths[chr.range[ii]],size=cnts[cc]*dubb,replace=TRUE))
    } else { randoms[[cc]] <- NA }
    #prv(oo)
    if(!SNP) { 
      kk <- length(randoms[[cc]])/2
      #print(kk)
      starts[[cc]] <- randoms[[cc]][1+(2*(0:(kk-1)))]
      ends[[cc]] <- randoms[[cc]][2*(1:kk)]
    } 
  }
  if(SNP) { starts <- ends <- randoms }
  starts <- unlist(starts); ends <- unlist(ends)
  chrs <- paste0(if(chr.pref) { "chr" } else { "" },chrs)
  chrs <- gsub("23","X",chrs);  chrs <- gsub("24","Y",chrs)
  chrs <- gsub("25","XY",chrs);  chrs <- gsub("26","MT",chrs)
  gg <- GRanges(ranges=IRanges(start=starts,end=ends),seqnames=chrs)
  if(!order) {
    gg <- gg[order(rnorm(nrow(gg))),]
  } else {
    gg <- toGenomeOrder(gg,strict=TRUE)
  }
  if(!GRanges) { gg <- as(gg,"RangedData"); gg <- gg[,-1] }
  if(n==0) { return(gg[-1,])} # returns empty ranges if that's what you really want
  if(fakeids) {
    if(SNP) { rownames(gg) <- rsnpid(nrow(gg)) } else { rownames(gg)  <- rsampid(nrow(gg)) }
  }
  return(gg)
}



#' Extract chromosome numbers from GRanges/RangedData 
#' 
#' Sometimes chromosomes are codeds as 1:22, sometimes there is also X,Y, etc, sometimes it's 
#' chr1, ch2, etc. This function extracts the set of chromosome labels used by a ranged object 
#' (ie, GRanges or RangedData) and converts the labels to numbers in a consistent way, so
#' 1:22, X, Y, XT, MT ==> 1:26, and optionally you can output the conversion table of codes to
#' numbers, then input this table for future conversions to ensure consistency.
#' @param ranged GRanges or RangedData object
#' @param warn logical, whether to display a warning when non autosomes are converted to numbers
#' @param table.out logical, whether to return a lookup table of how names matched to integers
#' @param table.in data.frame/matrix, col 1 is the raw text names, col 2 is the integer that should be assigned,
#'  col 3 is the cleaned text (of col 1) with 'chr' removed. the required form is outputted by this function if
#'  you set 'table.out=TRUE', so the idea is that to standardize coding amongst several RangedData objects you
#'  can save the table each time and ensure future coding is consistent with this. Note that chromosomes 1-22, X,
#'  Y, XY, and MT are always allocated the same integer, so table is only useful where there are extra NT, COX, HLA
#'  regions, etc.
#' @return a set of integers of length equal to the number of unique chromosomes in the ranged data.
#' @export
#' @examples
#' gg <- rranges(1000)
#' chrNames(gg); chrNums(gg)
#' gg <- rranges(1000,chr.pref=TRUE) # example where chromosomes are chr1, chr2, ...
#' chrNames(gg); chrNums(gg)
#' lookup <- chrNums(gg,table.out=TRUE)
#' lookup
#' gg2 <- rranges(10)
#' chrNums(gg2,table.in=lookup) # make chromosome numbers using same table as above
chrNums <- function(ranged,warn=FALSE,table.out=FALSE,table.in=NULL) {
  #must.use.package("genoset",bioC=T)
  typ <- is(ranged)[1]
  if(!typ %in% c("RangedData","GRanges")) { warning("not a GRanges or RangedData object"); return(NULL) }
  lookup <- c("X","Y","XY","MT")
  txt1 <- chrNames2(ranged)
  txt <- gsub("chr","",txt1,fixed=T)
  nums <- suppressWarnings(as.numeric(txt))
  num.na <- length(nums[is.na(nums)])
  if(num.na>0) { 
    if(warn) { warning(paste("chromosome numbers requested for non-autosomes, will assign numbers >=23 to letters",
                             paste(txt[is.na(nums)],collapse=","))) }
    aux.ind <- match(txt,lookup)
    nums[!is.na(aux.ind)] <- 22+aux.ind[!is.na(aux.ind)]
    unmatched <- txt[is.na(nums)]
    if(!is.null(table.in)) {
      if((all(table.in[,1] %in% unmatched)) | (all(unmatched %in% table.in[,1]))) {
        if(all(unmatched %in% table.in[,1])) {
          out <- table.in[,2][match(unmatched,table.in[,1])]
          nums[is.na(nums)] <- as.numeric(out)
        } else {
          out <- table.in[,2][match(unmatched,table.in[,1])]
          nums[is.na(nums)][!is.na(out)] <- as.numeric(out)[!is.na(out)]
          st.num <- max(c(22+length(lookup),as.numeric(table.in[,2])),na.rm=T)+1
          nums[is.na(nums)][is.na(out)] <- st.num:(st.num+length(nums[is.na(nums)][is.na(out)])-1)
        }
      } else {
        out <- table.in[,2][match(unmatched,table.in[,1])]
        nums[is.na(nums)][!is.na(out)] <- as.numeric(out)[!is.na(out)]
        st.num <- max(c(22+length(lookup),as.numeric(table.in[,2])),na.rm=T)+1
        nums[is.na(nums)][is.na(out)] <- st.num:(st.num+length(nums[is.na(nums)][is.na(out)])-1)
      }
    } else {
      nums[is.na(nums)] <- 27:(27+length(nums[is.na(nums)])-1)
    }
  }
  if(table.out) {
    out <- cbind(txt1,nums,txt)
    return(out)
  } else {
    return(sortna(as.numeric(nums)))
  }
}



#' Expand genomic locations to the ranges covering the 'n' closest SNPs
#' 
#' Sometimes for chip data we want to create windows around some locus, and
#' fixed distance [see flank()], recombination distance [see recwindow()] or a number of SNPs 
#' might be used. This function allows expansion of regions according to a set number of SNPs.
#' The result gives two regions for each row of a GRanges or RangedData object describing
#' the start and end of the left flanking 'nsnp' region, and right flanking 'nsnp' region
#' respectively.
#' @param ranged a GRanges or RangedData object describing the locations for
#' which we want to find regions encompassing 'nsnps' closest SNPs.
#' @param snp.info An object of type: ChipInfo, RangedData or GRanges, describing the set of SNPs
#' you are using (e.g, chip annotation). If left as null the ChipInfo object from chip.support() 
#' with default options() will be used
#' @param nsnp Number of nearest SNPs to return for each location
#' @param add.chr logical, whether to add a chromosome column for the output object
#' @seealso nearest.snp, chip.support, recwindow
#' @export
#' @return Two regions for each row of a the 'ranged' object describing
#' the start and end of the left flanking 'nsnp' region, and right flanking 'nsnp' region
#' respectively. If 'ranged' has rownames these should stay in the same order in the resulting
#' object. Chromosome will be the final column if you set add.chr=TRUE.
#' @examples
#' rngs <- rranges()
#' # not run - slow ~5 seconds # exp.window.nsnp(rngs)
#' # not run - slow ~5 seconds # exp.window.nsnp(rngs,add.chr=TRUE)
exp.window.nsnp <- function(ranged,snp.info=NULL,nsnp=10, add.chr=FALSE) {
  if(is.null(snp.info)) { snp.info <- chip.support() }
  if(!is(snp.info)[1] %in% c("ChipInfo","RangedData","GRanges")) { 
    stop("snp.info must be of type: ChipInfo, RangedData or GRanges") }
  snp.info <- toGenomeOrder2(snp.info,strict=TRUE); rw.cnt <- 1
  all.fl <- matrix(ncol=4+(as.numeric(add.chr)), nrow=0)
  for(cc in chrNums(ranged)) {
    nxt.nm <- rownames(snp.info[paste(cc)]); pos <- start(snp.info[paste(cc)])
    rng <- chr.sel(ranged,paste(cc)) # ranged[paste(cc)]
    st.en.snp <- range.snp(ranged=rng,snp.info=snp.info)
    fl <- matrix(ncol=4, nrow=nrow(st.en.snp))
    fl[,2] <- start(rng); fl[,3] <- end(rng);
    for(dd in 1:nrow(st.en.snp)) {
      x1 <- pos[max(1,match(st.en.snp[dd,1],nxt.nm)-nsnp)]
      x2 <- pos[min(length(nxt.nm),match(st.en.snp[dd,2],nxt.nm)+nsnp)]
      #print(x1); print(x2); print(length(x1)); print(length(x2))
      fl[dd,1] <- x1[1]
      fl[dd,4] <- x2[1]
    }
    if(add.chr) { 
      chrz <- cc
      fl <- cbind(fl,chrz)
    }
    all.fl <- rbind(all.fl,fl)
  }
  fl <- (all.fl)
  fl[fl[,1]>fl[,2],1] <- fl[fl[,1]>fl[,2],2]
  fl[fl[,3]>fl[,4],1] <- fl[fl[,3]>fl[,4],4]
  cN <- c("left.start","left.end","right.start","right.end")
  if(add.chr) { cN <- c(cN,"chr") }
  colnames(fl) <- cN
  if(nrow(fl)==nrow(ranged) & !is.null(rownames(ranged))) { rownames(fl) <- rownames(ranged) }
  return(fl)
}




#' Find closest SNPs to the ends of ranges
#' 
#' For given genome ranges (GRanges/RangedData) will try to find the closest snps to the end of the ranges.
#' @param ranged A GRanges or RangedData object specifying the range(s) you wish to find SNPs near the
#' ends of. Alternatively leave this parameter as NULL and specify ranges using chr, pos
#' @param snp.info ChipInfo/GRanges/Ranged data object describing the SNPs relevant to your query, e.g, 
#' SNPs on the chip you are using. If left NULL, the SNP set used will be that retrieved by chip.support()
#' which will depend on your options() settings, see ?chip.support for more info
#' @param chr optional alternative to 'ranged' input, use in conjunction with 'pos' to specify the ranges
#' to find the SNPs near the ends of.
#' @param pos matrix with 2 columns for start, end positions, or a single column if all ranges are SNPs.
#' An optional alternative to 'ranged' input, use in conjunction with 'chr' to specify the ranges
#' to find the SNPs near the ends of.
#' @param nearest will preferably find an exact match but if nearest=T, will fall-back on nearest match, 
#' even if slightly outside the range.
#' @export
#' @return a list of SNP-ids (rownames of 'snp.info') fulfilling the criteria, the output vector (character)
#' should be the same length as the number of ranges entered.
#' @examples
#' end.snp(chr=c(1:3),pos=cbind(c(100000,200000,300000),c(30000000,4000000,10000000)))
#' end.snp(rranges())
end.snp <- function(ranged=NULL,snp.info=NULL,chr=NULL,pos=NULL,nearest=T) {
  return(start.snp(ranged=ranged,snp.info=snp.info,chr=chr,pos=pos,start=F,end=T,nearest=nearest))
}


#' Find closest SNPs to the starts and ends of ranges
#' 
#' For given genome ranges (GRanges/RangedData) will try to find the closest snps to the starts and ends
#' of the ranges.
#' @param ranged A GRanges or RangedData object specifying the range(s) you wish to find SNPs near the
#' starts/ends of. Alternatively leave this parameter as NULL and specify ranges using chr, pos
#' @param snp.info ChipInfo/GRanges/Ranged data object describing the SNPs relevant to your query, e.g, 
#' SNPs on the chip you are using. If left NULL, the SNP set used will be that retrieved by chip.support()
#' which will depend on your options() settings, see ?chip.support for more info
#' @param chr optional alternative to 'ranged' input, use in conjunction with 'pos' to specify the ranges
#' to find the SNPs near the starts/ends of.
#' @param pos matrix with 2 columns for start, end positions, or a single column if all ranges are SNPs.
#' An optional alternative to 'ranged' input, use in conjunction with 'chr' to specify the ranges
#' to find the SNPs near the starts/ends of.
#' @param nearest will preferably find an exact match but if nearest=T, will fall-back on nearest match, 
#' even if slightly outside the range.
#' @export
#' @return a list of SNP-ids (rownames of 'snp.info') fulfilling the criteria, the output will be a matrix
#' which should have the same number of rows as the number of ranges entered.
#' @examples
#' range.snp(chr=c(1:3),pos=cbind(c(100000,200000,300000),c(30000000,4000000,10000000)))
#' range.snp(rranges())
range.snp <- function(ranged=NULL,snp.info=NULL,chr=NULL,pos=NULL,nearest=T) {
  return(start.snp(ranged=ranged,snp.info=snp.info,chr=chr,pos=pos,start=T,end=T,nearest=nearest))
}


#' Find closest SNPs to the starts of ranges
#' 
#' For given genome ranges (GRanges/RangedData) will try to find the closest snps to the starts 
#' of the ranges.
#' @param ranged A GRanges or RangedData object specifying the range(s) you wish to find SNPs near the
#' starts of. Alternatively leave this parameter as NULL and specify ranges using chr, pos
#' @param snp.info ChipInfo/GRanges/Ranged data object describing the SNPs relevant to your query, e.g, 
#' SNPs on the chip you are using. If left NULL, the SNP set used will be that retrieved by chip.support()
#' which will depend on your options() settings, see ?chip.support for more info
#' @param chr optional alternative to 'ranged' input, use in conjunction with 'pos' to specify the ranges
#' to find the SNPs near the starts of.
#' @param pos matrix with 2 columns for start, end positions, or a single column if all ranges are SNPs.
#' An optional alternative to 'ranged' input, use in conjunction with 'chr' to specify the ranges
#' to find the SNPs near the starts of.
#' @param start logical whether to return the SNP nearest the range starts
#' @param end logical whether to return the SNP nearest the range ends
#' @param nearest will preferably find an exact match but if nearest=T, will fall-back on nearest match, 
#' even if slightly outside the range.
#' @export
#' @return a list of SNP-ids (rownames of 'snp.info') fulfilling the criteria, the output will be a vector
#' which will have the same length as the input. Unless start=TRUE and end=TRUE, then will return a matrix
#' which should have the same number of rows as the number of ranges entered. Note that end.snp() is 
#' equivalent to using this function when end=TRUE and start=FALSE, and range.snp() is the same as setting
#' start=TRUE and end=TRUE.
#' @examples
#' start.snp(chr=c(1:3),pos=cbind(c(100000,200000,300000),c(30000000,4000000,10000000)))
#' start.snp(rranges())
start.snp <- function(ranged=NULL,snp.info=NULL,chr=NULL,pos=NULL,start=T,end=F,nearest=T) {
  # will preferably find an exact match but if nearest=T, will fall-back on nearest match
  #must.use.package("genoset",T)
  nmz <- NULL
  if(!is(ranged)[1] %in% c("RangedData","GRanges")) {
    if(!is.null(chr) & !is.null(pos)) {
      if(is.null(dim(pos))) { st <- pos[1]; en <- pos[2] } else {
        st <- pos[,1]; en <- pos[,2]
      }
      if(length(st)>length(chr)) { chr <- rep(chr[1],length(st)) } else { chr <- chr[1:length(st)] }
    } else {
      stop("if not using 'ranged' input, then chr and pos must be valid")
    }
  } else {
    st <- start(ranged); en <- end(ranged); chr <- chr2(ranged)
    nmz <- rownames(ranged)
  }
  if(is.null(snp.info)) { snp.info <- chip.support() }  # load default chip
  if(!is(snp.info)[1] %in% c("ChipInfo","RangedData","GRanges")) {
    stop("snp.info must be of type ChipInfo, RangedData or GRanges")
  } else {
    if(is.null(rownames(snp.info))) {  rownames(snp.info) <- paste(1:nrow(snp.info)) }
  }
  st.snps <- en.snps <- character(length(chr)) ; prch <- 0
  for (cc in 1:length(chr)) {
    if(chr[cc]!=prch) { 
      ref <- chr.sel(snp.info,paste(chr[cc])) # snp.info[paste(chr[cc])]
      st.ref <- start(ref); rnref <- rownames(ref)
      if(is.null(ref)) { stop(paste("snp.info did not contain chr",chr[cc])) }
    }
    #exact
    if(start) { 
      ind <- match(st[cc],st.ref)
      if(any(is.na(ind)) & nearest) {
        difs <- abs(st[cc]-st.ref)
        ind <- which(difs==min(difs,na.rm=T))[1]
      }
      if(length(ind)==0) { st.snps[cc] <- NA } else {
        st.snps[cc] <- rnref[ind]
      }
    }
    if(end){
      ind2 <- match(en[cc],st.ref)
      if(any(is.na(ind2)) & nearest) {
        difs <- abs(en[cc]-st.ref)
        ind2 <- which(difs==min(difs,na.rm=T))[1]
      }
      if(length(ind2)==0) { en.snps[cc] <- NA } else {
        en.snps[cc] <- rnref[ind2]
      }
    }
    prch <- chr[cc]
  }
  if(start & !end) {
    return(st.snps)
  }
  if(!start & end) {
    return(en.snps)
  }
  #otherwise looks like want both
  out <- cbind(st.snps,en.snps)
  if(!is.null(nmz)) { if(length(nmz)==nrow(out)) { rownames(out) <- nmz } }
  return(out)
}


#' Force a valid genomic range, given the inputted coordinates
#'
#' Enter a pair of genomic locations representing a range for a given chromosome and this
#' function will ensure that no position is less than 1 or greater than the relevant chromosome
#'  lengths. Anything below will be coerced to 1, and anything above to the chromosome length.
#' @param Pos must be numeric, length 2, e.g, c(20321,30123)
#' @param Chr chromosome label
#' @param snp.info optional object to take boundaries from, the maxima and minima for each
#' chromosome within this object will take the place of the chromsome lengths / 1.
#' @build ucsc build, only need to enter if this differs from getOption("ucsc")
#' @dir directory to use for download of chromosome lengths (only if you wish to
#' keep the chromosome length file)
#' @export
#' @examples
#' pss <- ps <- c(345035,345035); ch <- 1
#' force.chr.pos(ps,ch)
#' pss[1] <- 0
#' force.chr.pos(pss,ch) # won't allow zero
#' pss[1] <- -1
#' force.chr.pos(pss,ch) # won't allow negative
#' pss[1] <- 645035012
#' force.chr.pos(pss,ch) # won't allow pos > chromosome length
force.chr.pos <- function(Pos,Chr,snp.info=NULL,build=NULL,dir=NULL) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  # convert any non autosomes to numbers:
  Chr <- paste(Chr)
  Chr[grep("c6",Chr,ignore.case=T)] <- 6  # prevent issues with c6_COX, c6_QBL  
  Chr[grep("X",Chr,ignore.case=T)] <- 23
  Chr[grep("Y",Chr,ignore.case=T)] <- 24
  Chr[grep("M",Chr,ignore.case=T)] <- 25
  Chr[grep("NT",Chr,ignore.case=T)] <- 26  # prevent issues with NT_11387, etc
  Chr <- as.numeric(Chr)
  if(any(!paste(Chr) %in% paste(c(1:26)))) { stop("invalid chromosome(s) entered") }
  if(any(paste(Chr) == paste(26))) { warning("'NT' chromosome(s) entered, not supported, NAs produced") }
  if(length(Pos)==2 & is.numeric(Pos)) {
    if(is(snp.info)[1]!="RangedData" & is(snp.info)[1]!="GRanges") { 
      maxln <- get.chr.lens(dir=dir,mito=T,autosomes=FALSE,build=build)[Chr] 
    } else { 
      maxln <- end(tail(snp.info[paste(Chr)],1)) # force start and end to be within 1:chr.len
    }
    mbs <- min(max(1,Pos[1]),(maxln-1)); mbe <- min(max(2,Pos[2]),maxln)
    return(c(mbs,mbe))
  } else {
    Pos <- NA; warning("Pos needs to be numeric length 2, min, max")
  }
  return(round(Pos))
}







#' Select all ranges lying within a chromosome window
#' 
#' Input a ranged object (ie., GRanges or RangedData) and this function will
#' return the subset from chromosome 'chr' and within the base-pair range specified
#' by 'pos', in units of 'unit'. By default ranges with ANY overlap are returned, but
#' it can be specified that it must be full overlap. Duplicates can be removed.
#' @param ranged GRanges or RangedData object
#' @param chr a chromosome, e.g, 1,2,3,...,22,X,Y,XY,MT or however chromosomes are 
#' annotated in 'ranged'
#' @param pos a numeric range (length 2), with a start (minima) and end (maxima), specifying
#' the window on the chromosome to select ranges from, base-pair units are specified by 'unit'.
#' @param full.overlap logical, the default is to return objects with ANY overlap with the window,
#' whereas setting this as TRUE, will only return those that fully overlap
#' @param unit the unit of base-pairs that 'pos' is using, eg, "b", "kb", "mb", "gb"
#' @param rmv.dup logical, whether to remove duplicate ranges from the return result. The default
#' is not to remove duplicates.
#' @return an object of the same type as 'ranged', but only containing the rows that
#' were within the specified bounds
#' @export
#' @examples
#' iG <- get.immunog.locs()[2,] # select the 2nd iG region
#' ciG <- chr(iG)  #  get the chromosome
#' posiG <- c(start(iG),end(iG) # get the region start and end
#' rr <- rranges(10000) # create a large random GRanges object
#' in.window(rr,chr=ciG,pos=posiG) # set with ANY overlap of iG
#' in.window(rr,chr=ciG,pos=posiG,TRUE) # set with FULL overlap of iG
#' in.window(rr,chr=6,pos=c(25,35),unit="mb") # look between 25 - 35 MB on chr6 [ie, MHC]
in.window <- function(ranged,chr,pos,full.overlap=F, unit=c("b","kb","mb","gb"), rmv.dup=FALSE) {
  if(length(pos)>2 | !is.numeric(pos)) { warning("pos should be a start and end numeric range"); return(NULL) }
  if(length(pos)==1) { pos <- rep(pos,2) }
  all.chrz <- unique(chrNames2(ranged))
  if(length(chr)>1 | (!chr %in% all.chrz)) { warning("chr must be a value in",comma(all.chrz)); return(NULL) }
  typ <- is(ranged)[1]
  if(!any(typ %in% c("RangedData","IRanges","GRanges","RangesList"))) { 
    warning("'ranged' should be a RangedData type or similar"); return(NULL) }
  pos <- pos*make.divisor(unit,"unit")
  #unit <- tolower(unit[1]) ; mult <- switch(unit,b=0,kb=3,mb=6,gb=9); pos <- pos*10^mult
  # get set of genes in a position range for a chromosome
  chr.genez <- chr.sel(ranged,paste(chr)) 
  if(full.overlap) {
    ranged <- chr.genez[which(start(chr.genez)>min(pos) & end(chr.genez)<max(pos)),]
  } else {
    # any overlap
    ranged <- chr.genez[which(end(chr.genez)>min(pos) & start(chr.genez)<max(pos)),]
  }
  if(rmv.dup) {
    # remove duplicate genes/exons
    ranged <- ranged[!(duplicated(start(ranged)) & duplicated(end(ranged)) & duplicated(width(ranged))),]
  }
  return(ranged)
}



#' Plot the locations specified in a GRanges or RangedData object
#' 
#' GRanges and RangedData objects are used in bioconductor to store genomic locations and
#' ranges, such as transcripts, genes, CNVs and SNPs. This function allows simple
#' plotting of this data directly from the ranged object. SNPs will be plotted as dots 
#' and ranges as lines. Either can be plotted using vertical bars at the start/end of each
#' range. There are options for labelling and other graphical parameters.
#' @param ranged GRanges or RangedData object with genomic ranges. Should only contain
#' one chromosome, but if not, the first will be used
#' @param labels by default labels for each range are taken from the rownames of 'ranged',
#' but if you want to use another column in the ranged object, specify the column name
#' or number to use to label these ranges on the plot. Or else input a character
#' vector the same length as ranged for custom labels.
#' @param do.labs logical, whether or not to display these labels
#' @param skip.plot.new logical, whether to append to an existing plot (TRUE), or start
#' a new plot (FALSE --> default)
#' @param lty line type to use, see '?lines()' - not used for SNP data when v.lines=FALSE
#' @param alt.y alternative y-axis values (other than the default ordering from the input)
#' This can be a vector of length 1 or length(ranged), or else a column name in ranged to 
#' take the values from
#' @param v.lines TRUE will plot the ranges as pairs of vertical lines, occupying the full
#' vertical extent of the plot, whereas FALSE will plot the ranges as individual horizontal lines
#' @param ylim numeric, length 2, the y-axis limits for the plot, same a 'ylim' for ?plot()
#' @param xlim numeric, length 2, the x-axis limits for the plot, same a 'xlim' for ?plot(),
#' This shouldn't usually be needed as the automatic x-limits should work well,
#'  however is here in case fine tuning is required.
#' @param scl character, the scale that the x axis uses, ie, "b","kb","mb", or "gb", meaning
#' base-pairs, kilobases, megabases or gigabase-pairs.
#' @param col character, colour, same as 'col' argument for plot(), etc.
#' @param srt integer, text rotation in degrees (see par) for labels
#' @param pos integer, values of ‘1’, ‘2’, ‘3’ and ‘4’, respectively indicate positions below, 
#' to the left of, above and to the right of the specified coordinates. See 'pos' in graphics:text()
#' @param lwd line width, see '?lines()' - not used for SNP data when v.lines=FALSE
#' @param pch point type, see '?points()' - not used for ranged data
#' @param cex font/symbol size, see '?plot()' - passed to plot, points if using SNP data 
#' @param ... further arguments to 'plot', so long as skip.plot.new==FALSE.
#' @export
#' @return Plots the ranges specified in 'ranged' to the current plot, or to a new plot
#' @example
#' rr <- in.window(rranges(5000),chr=6,pos=c(28,32),unit="mb") # make some random MHC ranges
#' rownames(rr) <- paste0("range",1:length(rr))
#' # plot ranges vertically #
#' plot.ranges(rr,v.lines=TRUE)
#' # make some labels and plot as horizontal lines #
#' rr2 <- rr[1:5,]; mcols(rr2)[["GENE"]] <- c("CTLA9","HLA-Z","BS-1","FAKr","teST")
#' plot.ranges(rr2,label="GENE",scl="Mb",col="black",
#'             xlab="Chr6 position (megabases)",
#'             yaxt="n",ylab="",bty="n")
#' # create some SNPs and plot
#' rr3 <- rr; end(rr3) <- start(rr3) 
#' rownames(rr3) <- paste0("rs",sample(10^6,nrow(rr3)))
#' plot.ranges(rr3,col="blue",yaxt="n",ylab="",bty="n")
plot.ranges <- function(ranged,labels=NULL,do.labs=T,skip.plot.new=F,lty="solid", alt.y=NULL,
                        v.lines=FALSE,ylim=NULL,xlim=NULL,scl=c("b","Kb","Mb","Gb"),
                        col=NULL,srt=0,pos=4,pch=1,lwd=1,cex=1,...) {
  #if(is(ranges)[1]!="RangedData") { warning("need RangedData object") ; return(NULL) }
  chk <- chrNums(ranged)
  typ <- is(ranged)[1]
  if(!is.null(alt.y)) {
    if(is.numeric(alt.y)) {
      if(length(alt.y)==1 | length(alt.y)==length(ranged)) {
        yy <- alt.y
      } else {
        warning("alt.y ignored, must be same length as ranged, or else length 1"); alt.y <- NULL
      }
    } else {
      if(is.character(alt.y)) {
        if(typ=="GRanges") { 
          cn <- colnames(mcols(ranged)); df <- mcols(ranged)
        } else {
          cn <- colnames(ranged); df <- ranged
        }
        if(!alt.y %in% cn) { stop("alternative y.axis column name ",alt.y," not found in 'ranged'") }
        yy <- df[,alt.y]; rm(df)
      } else { 
        warning("invalid value for alt.y, ignoring"); alt.y <- NULL
      }
    }
  }
  if(!is.null(labels)) {
    labels <- paste(labels)
    if(is.character(labels)) {
      if(length(labels)==1 | length(labels)==length(ranged)) {
        if(length(labels)==1) {
          if(typ=="GRanges") { 
            cn <- colnames(mcols(ranged)); df <- mcols(ranged)
          } else {
            cn <- colnames(ranged); df <- ranged
          }
          if(!labels %in% cn) { stop("labels column name ",labels," not found in 'ranged'") }
          lab <- df[,labels]; rm(df)
        } else {
          lab <- labels
        }
      } else {
        warning("labels ignored, must be same length as ranged, or else length 1"); labels <- NULL
      }
    } else {
      warning("invalid value for labels, ignoring"); labels <- NULL
    }
  } else {
    lab <- rownames(ranged) 
  } 
  if(length(chk)>1) { 
    warning(length(chk)," chromosomes in 'ranged', only using the first, chr",chk[1]) 
    ranged <- chr.sel(ranged,1) 
  }
  if(all(width(ranged)<=1)) { theyAreSnps <- TRUE } else { theyAreSnps <- FALSE }
  scl <- make.divisor(scl)
  xl <- range(c(start(ranged),end(ranged)))
  xl <- xl + ((diff(xl)*0.1)*c(-1,1))
  xl <- xl/scl
  nr <- nrow(ranged); if(is.null(nr)) { nr <- length(ranged) }
  if(is.null(alt.y)) {
    yl <- c(0,(nr+2))
  } else {
    yl <- range(yy)
  }
  if(is.numeric(ylim) & length(ylim)==2) {
    ylim <- range(ylim)
    ydif <- diff(ylim)
    yl <- ylim
  }
  if(is.numeric(xlim) & length(xlim)==2) {
    xlim <- range(xlim)
    xdif <- diff(xlim)
    xl <- xlim
  }
  if(is.null(alt.y)) {
    YY <- seq(from=yl[1],to=yl[2],length.out=nr+2)[-1]
  } else {
    if(length(yy)==1) { YY <- rep(yy,length(nr)) } else { YY <- yy }
  }
  #print(YY)
  if(!is.null(col)) {
    if(length(col)==1) {
      col <- rep(col,times=nr) 
    } else {
      if(length(col)!=nr) { warning("col was not the same length as ranged, using first only"); col <- rep(col[1],nr) }
    }
  }
  if(is.null(col)) {
    if(nr>22) { colz <- rep("black",nr) } else { colz <- get.distinct.cols(nr) }
  } else { colz <- col[1:nr] }
  if(is.null(lab) & do.labs) { lab <- paste(1:nr) } # last resort
  if(!skip.plot.new) {
    position <- c(start(ranged[1,]),end(ranged[1,]))/scl
    Y <- YY[c(1,1)]
    #prv(position,Y)
    TY <- if(theyAreSnps) { "p" } else { "l" }
    if(v.lines) {
      plot(x=position, y=Y, xlim=xl, ylim=yl, type=TY, col="white", lty=lty, ...)
      abline(v=position,col=colz[1])
    } else {
      plot(x=position, y=Y, xlim=xl, ylim=yl, type=TY, col=colz[1], lty=lty, lwd=lwd, cex=cex, ...)
    }
    st <- 2
  } else {
    st <- 1
  }
  if(nr>1 | st==1) {
    for (cc in st:nr) {
      if(v.lines) {
        abline(v=c(start(ranged[cc,]),end(ranged[cc,]))/scl,col=colz[cc],lty=lty)
      } else {
        if(theyAreSnps) { 
          points(x=c(start(ranged[cc,]),end(ranged[cc,]))/scl,y=YY[c(cc,cc)],col=colz[cc], pch=pch, cex=cex)
        } else { 
          lines(x=c(start(ranged[cc,]),end(ranged[cc,]))/scl,y=YY[c(cc,cc)],col=colz[cc], lty=lty, lwd=lwd)
        }
      }
    }
  }
  if(do.labs) {
    for (cc in 1:nr) {
      if(v.lines) { YY <- rep(tail(YY,1),length(YY)) }
      V.scale <- (diff(head(YY,2))*0.5)
      if(length(V.scale)<1 | srt!=90) { V.scale <- 0 }
      text(x=start(ranged[cc,])/scl,y=YY[cc]+V.scale,labels=lab[cc],cex=0.6,pos=pos,offset=0,srt=srt)
    }
  }
}




#' Change the chromosome labels in a RangedData or GRanges object to string codes
#' 
#' @param ranged A GRanges or RangedData object
#' @param do.x.y logical, if TRUE then the usual numbers allocated to chromosomes, X,Y,XY, MT will
#' be allocated as 23,24,25,26 respectively. If false, these will just have 'chr' appended as a
#' prefix
#' @param keep logical, whether to keep additional metadata columns in the new object 
#' @seealso set.chr.to.numeric
#' @export
#' @return returns the 'ranged' object, but wherever a chromosome number was previously, a character
#' label, e.g, 'chr1', or 'X', will returned to replace the number, e.g, 1 or 23 respectively. 
#' If table.out is TRUE will return a list where the first element is the resulting object, and the second 
#' element is a table showing which numbers were converted to what label This table
#' can then be used for future conversions via the parameter 'table.in' to ensure consistency of
#' coding.
#' @examples
#' x <- rranges()
#' x
#' x <- set.chr.to.numeric(x) # make entirely numeric
#' x <- rranges(chr.range=20:26)
#' # next two will give warning about X, Y, etc
#' set.chr.to.char(x) # 23 = chrX, etc
#' set.chr.to.char(x,do.x.y=FALSE) # 23=chr23, etc
set.chr.to.char <- function(ranged,do.x.y=T,keep=T) {
  #must.use.package("genoset",bioC=T)
  typ <- is(ranged)[1]
  if(!typ %in% c("GRanges","RangedData")) { warning("not a GRanges or RangedData object"); return(NULL) }
  if(length(grep("chr",chrNames2(ranged)))<length(chrNames2(ranged))) {
    ranged <- toGenomeOrder2(ranged,strict=TRUE)
    #prv(ranged)
    mychr2 <- mychr <- paste(chr2(ranged))
    RN <- rownames(ranged)
    #all.nams <- chrNames2(ranged)
    #all.nums <- chrNums(ranged,table.in=table.in)
    if(length(grep("23",paste(mychr2)))>0) { 
      warning("use of arbitrary chromosome numbers for non-autosomes (i.e, >=23)",
              "can lead to annotation issues, try to use labels, X, Y, MT, and XY where possible") }
    sar <- select.autosomes(ranged)
    if(nrow(sar)>0) {
      all.nums.t <- chrNums(sar,table.in=NULL,table.out=T) 
      all.nams <- all.nums.t[,1]
      all.nums <- all.nums.t[,2]
      #mychr2 <- all.nums.t[,2][match(mychr,all.nums.t[,1])]
      for (cc in 1:length(all.nams)) { mychr2[which(mychr==all.nams[cc])] <- paste("chr",all.nums[cc],sep="") }
    } else {
      # no autosomes
    }
    if(do.x.y) {
      mychr2 <- gsub("X","chrX",mychr2)
      mychr2 <- gsub("Y","chrY",mychr2)
      mychr2 <- gsub("23","chrX",mychr2)
      mychr2 <- gsub("24","chrY",mychr2)
      mychr2 <- gsub("25","chrXY",mychr2)
      mychr2 <- gsub("26","chrM",mychr2)
      mychr2 <- gsub("chrXchrY","XY",mychr2)
      mychr2 <- gsub("chrYchrX","YX",mychr2) 
      mychr2 <- gsub("MT","chrM",mychr2)
      mychr2 <- gsub("XY","chrXY",mychr2)
      mychr2 <- gsub("chrchr","chr",mychr2)
    } else {
      mychr2[mychr2 %in% paste(23:100)] <- paste0("chr",mychr2[mychr2 %in% paste(23:100)])
    }
    #print(tail(mychr2)); print((all.nums))
    #prv(mychr2)
    if(any(is.na(mychr2))) { prv(mychr2[which(is.na(mychr2))]) }
    if(is.null(RN) | length(RN)!=nrow(ranged)) { RN <- 1:nrow(ranged) } #make sure RN's are valid
    if(is(ranged)[1]=="GRanges") {
      all.chr <- chr2(ranged)
      out <- GRanges(ranges=IRanges(start=start(ranged),end=end(ranged),names=RN),seqnames=mychr2)
    } else {
      out <- RangedData(ranges=IRanges(start=start(ranged),end=end(ranged),names=RN),space=mychr2)
    }
    out <- toGenomeOrder2(out,strict=TRUE)
    # prv(out)
    ###return(out)
    # need to allow for different indexing of dataframe part for GRanges
    ncr <- switch(typ,RangedData=ncol(ranged),GRanges=ncol(mcols(ranged)))
    if(is.null(ncr)) { ncr <- 0 }
    if(ncr>0 & keep) {
      cn <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
      for(cc in 1:ncr) {
        if(typ=="GRanges") {
          mcols(out)[[paste(cn[cc])]] <- mcols(ranged)[[paste(cn[cc])]]
        } else {
          out[[paste(cn[cc])]] <- ranged[[paste(cn[cc])]]
        }
      }
    }
    return(out)
  } else {
    #cat("no change\n")
    return(ranged)      # change not needed
  }
}




#' Change the chromosome labels in a RangedData or GRanges object to numbers
#' 
#' @param ranged A GRanges or RangedData object
#' @param keep logical, whether to keep additional metadata columns in the new object 
#' @param table.in matrix/data.frame object, usually a result of a prior run of 
#' set.chr.to.numeric(table.out=TRUE), which shows for each label (column 1), what
#' chromosome number should correspond. A way of ensuring consistent coding in different
#' sets.
#' @param table.out logical, if FALSE, the output will just be the object with updated 
#' chromosome labels. If TRUE, then the output will be a list, where the first element
#' is the updated object and the second object is a table describing the coding
#' scheme used to convert from labels to numeric indices.
#' @seealso set.chr.to.char
#' @export
#' @return returns the 'ranged' object, but wherever a chromosome label was previously a character
#' label, e.g, 'chr1', or 'X', will return as a number, e.g, 1 or 23 respectively. If table.out
#' is TRUE will return a list where the first element is the resulting object, and the second 
#' element is a table showing which labels were converted to what number. This table
#' can then be used for future conversions via the parameter 'table.in' to ensure consistency of
#' coding.
#' @examples
#' char <- rranges(chr.pref=T)
#' char
#' set.chr.to.numeric(char)
#' # behaviour with X, Y, etc
#' char <- rranges(chr.range=c(20:26))
#' #' char
#' set.chr.to.numeric(char)
#' tab <- set.chr.to.numeric(char,table.out=TRUE)[[2]]
#' tab # codes used in conversion #
#' char <- rranges(chr.range=c(20:26))
#' set.chr.to.numeric(char, table.in=tab) # code using codes from 'tab'
set.chr.to.numeric <- function(ranged,keep=T,table.in=NULL,table.out=FALSE) {
  typ <- is(ranged)[1]
  if(!typ %in% c("GRanges","RangedData")) { warning("not a GRanges or RangedData object"); return(NULL) }
  if(table.out | suppressWarnings(any(is.na(as.numeric(paste(chr2(ranged))))))) {
    silly.name <- "adf89734t5b"
    ranged <- toGenomeOrder2(ranged,strict=T)
    if(typ=="GRanges") {
      mcols(ranged)[[silly.name]] <- paste(1:nrow(ranged))
    } else {
      ranged[[silly.name]] <- paste(1:nrow(ranged))
    }
    #prv(ranged)
    mychr2 <- mychr <- paste(chr2(ranged))
    #all.nams <- chrNames2(ranged)
    #all.nums <- chrNums(ranged,table.in=table.in)
    all.nums.t <- chrNums(ranged,table.in=table.in,table.out=T) 
    all.nams <- all.nums.t[,1]
    all.nums <- all.nums.t[,2]
    #mychr2 <- all.nums.t[,2][match(mychr,all.nums.t[,1])]
    for (cc in 1:length(all.nams)) { mychr2[which(mychr==all.nams[cc])] <- all.nums[cc] }
    #print(tail(mychr2)); print((all.nums))
    if(typ=="GRanges") {
      all.chr <- chr2(ranged)
      out <- GRanges(ranges=IRanges(start=start(ranged),end=end(ranged)),seqnames=mychr2,silly.name=mcols(ranged)[[silly.name]])
    } else {
      out <- RangedData(ranges=IRanges(start=start(ranged),end=end(ranged)),space=mychr2,silly.name=ranged[[silly.name]])
    }
    out <- toGenomeOrder2(out,strict=T)
    if(typ=="GRanges") {
      oo <- mcols(out)[["silly.name"]]
      rr <- mcols(ranged)[[silly.name]]
    } else {
      oo <- out[["silly.name"]]
      rr <- ranged[[silly.name]]
    }
    if(all(!is.na(oo))) {
      if(is.null(rownames(ranged))) { rownames(ranged) <- paste(1:nrow(ranged)) ; rmv.rn <- TRUE } else { rmv.rn <- FALSE }
      #iioo <- rownames(ranged)[match(oo,rr)]; print(iioo); print(is(iioo))
      rn <- narm(rownames(ranged)[match(oo,rr)])
      if(nrow(out)==length(rn) ) { rownames(out) <- rn } else { warning("rownames did not match number of rows") }
      if(rmv.rn) { rownames(out) <- NULL }
    } else {
      warning("index column was corrupted")
    }
    # prv(out)
    ncr <- switch(typ,RangedData=ncol(ranged),GRanges=ncol(mcols(ranged)))
    if(is.null(ncr)) { ncr <- 0 }
    if(ncr>0 & keep) {
      cn <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
      for(cc in 1:ncr) {
        if(typ=="GRanges") {
          mcols(out)[[paste(cn[cc])]] <- mcols(ranged)[[paste(cn[cc])]]
        } else {
          out[[paste(cn[cc])]] <- ranged[[paste(cn[cc])]]
        }
      }
    }
    cno <- switch(typ,RangedData=colnames(out),GRanges=colnames(mcols(out)))
    if(any(cno %in% "silly.name")) { out <- out[,-which(cno %in% "silly.name")] }
    cno <- switch(typ,RangedData=colnames(out),GRanges=colnames(mcols(out)))
    if(any(cno %in% silly.name)) { out <- out[,-which(cno %in% silly.name)] }
    if(table.out) {
      return(list(ranged=out,table.out=all.nums.t))
    } else {
      return(out)
    }
  } else {
    #cat("no change\n")
    cnr <- switch(typ,RangedData=colnames(ranged),GRanges=colnames(mcols(ranged)))
    if(any(cnr %in% "silly.name")) { ranged <- ranged[,-which(cnr %in% "silly.name")] }
    return(ranged)      # change not needed
  }
}





#' Invert a ranged object
#' Select the empty space between ranges for the whole genome, for instance you may want
#' to overlap with everything NOT in a set of ranges.
#' @param X a ranged object, GRanges, RangedData or ChipInfo
#' @param inclusive logical, TRUE if the ends of ranges should be in the inverted object
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param pad.missing.autosomes logical, whether to add entire chromosomes to the inverted
#' range object when they are not contained within X
#' @return a ranged object of the same type as X, but with the inverse set of human genomic ranges selected
#' @export
#' @examples
#' X <- rranges()
#' invert.granges(X,inclusive=TRUE)
#' invert.granges(X)
#' invert.granges(X,pad.missing.autosomes=FALSE)
invert.granges <- function(X,inclusive=FALSE,build=NULL,pad.missing.autosomes=TRUE) {
  typ <- is(X)[1]
  if(!typ %in% c("GRanges","RangedData","ChipInfo")) { stop("invalid type for X; ",typ) }
  
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  X <- toGenomeOrder2(X)
  X <- set.chr.to.char(X)
  ch <- chrNames2(X)
  chrLs <- get.chr.lens(mito=T,names=T,build=build)
  chm <- ch
  chm[ch %in% c("chrXY","XY")] <- gsub("Y","",chm[ch %in% c("chrXY","XY")])
  ii <- match(chm,names(chrLs))
  if(any(is.na(ii))) { stop("contained chromosome name not in reference: ",comma(ch[is.na(ii)])) }
  chrL <- as.integer(chrLs[ii])
  #all.dat <- GRanges()
  if(length(ch)>0) { 
  	all.dat <- vector("list",length=length(ch)); names(all.dat) <- ch
  	offs <- if(inclusive) { 0 } else { 1 }
    for (cc in 1:length(ch)) {
      nxt.chr <- chr.sel(X,ch[cc])
      st <- as.integer(start(nxt.chr)); en <- as.integer(end(nxt.chr))
      new.st <- as.integer(c(1,en+offs))
      new.en <- as.integer(c(st-offs,chrL[cc]))
      #prv(new.st,new.en)
      if(any(new.en<new.st)) {
        ind <- (rep(which(new.en<new.st),each=5)+rep(c(-2,-1,0,1,2),length(which(new.en<new.st))))
        ind <- ind[ind %in% 1:length(new.st)] #; prv(ind)
        cat("Found illegal start/end in ",ch[cc],"\n")
        print(head(cbind(chr=(rep(ch[cc],length(new.st))),start=new.st,end=new.en)[ind,]))
        if(length(grep("19",ch[cc]))>0) { cat("you may be using incorrect value of 'build' (current is '",build,"')\n",sep="") }
      }
      all.dat[[cc]] <- make.granges(chr=rep(ch[cc],length(new.st)),start=new.st,end=new.en)
    }
  } else {
  	all.dat <- NULL
  	warning("X was empty")
  	if(!pad.missing.autosomes) { return(NULL) }
  }
  if(pad.missing.autosomes) {
    autoz <- paste0("chr",1:22)
    misn <- (!autoz %in% ch)
    if(any(misn)) { 
      mis.list <- vector("list",length(which(misn))); names(mis.list) <- autoz[misn]
      for(dd in 1:length(which(misn))) {
        mis.list[[dd]] <- make.granges(chr=autoz[which(misn)[dd]],st=1,end=chrLs[which(misn)[dd]])
      }
      all.dat <- c(all.dat,mis.list)
    }
  }
  myDat <- do.call("rbind",args=lapply(all.dat,as,"RangedData"))
  myDat <- toGenomeOrder2(myDat)
  myDat <- as(myDat,typ)
  return(myDat)
}


################## end ranged ##########################
