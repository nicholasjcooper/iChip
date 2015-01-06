
################################
## ChipInfo Support Functions ##
################################


#' Retrieve current ChipInfo annotation object
#' 
#' This function returns the current 'ChipInfo' annotation object, containing chromosome,
#' id, position, strand, 'rs' id, allele 1, allele 2 for each SNP of a microarray chip,
#' in either hg18 or hg19 (build 36/37) coordinates.
#' This package makes extension use of this class of annotation object for the working
#' microarray chip, e.g, default is ImmunoChip, but Metabochip is also built-in,
#' and you can also load your own annotation if using a different chip. The class
#' of the object used is 'ChipInfo' which is a GRanges object, modified to always
#' have columns for A1, A2 (alleles), rs.id, and a quality control flag. The
#' default display is tidier than GRanges, it has nice coersion to and frame data.frame
#' and indexing by chromosome using [[n]] has been added, in addition to normal [i,j]
#' indexing native to GRanges.
#' @param build character, either "hg18" or "hg19". Will also accept build number, 36 or 37.
#' @param refresh logical, TRUE to just load whatever object is already in memory (except
#' when first using a function in this package, there should be a ChipInfo object loaded),
#' or FALSE to reload from the original source. For instance you may wish to do this when 
#' you want to use a different chip, different build, or if the annotation has been modifed
#' via a manual correction).
#' @param alternate.file character, name of an alternative RData file containing a ChipInfo
#' object to use instead of the object found in getOption("chip.info"). Such files may
#' contain two objects, one for build 36 and one for build 37; when doing this, the object
#' names should contain the characters '36' and '37' respectively.
#' @return returns the current ChipInfo object [S4]
#' @seealso ChipInfo, build, rs.id, QCfail, convTo36, convTo37, A1, A2
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' chip.support() # shows the current ChipInfo object (default is 'ImmunoChip' build 36)
chip.support <- function(build=NULL,refresh=FALSE,alternate.file=NULL) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  ## NEED TO ADD SUPPORT HERE TO USE CORRECT OUT OF 36/37
  if(!exists("all.support",envir=globalenv())) { refresh <- T }  # change global environment to namespace of iChip package
  if(refresh) {
    use.options <- TRUE
    if(is.character(alternate.file)) {
      if(file.exists(alternate.file)) {
        use.options <- FALSE
      }
    }
    if(use.options) {
      file <- getOption("chip.info")
    } else {
      file <- alternate.file
    }
    all.support <- reader(file)
    if(length(all.support)>1) { typz <- sapply(lapply(all.support,is),"[",1) } else { typz <- is(all.support)[1] }
    if(all(typz=="ChipInfo")) { 
      assign("all.support",value=all.support,envir=globalenv())  # change global environment to namespace of iChip package
    } else {
      stop("object (all.support) in the file",file,
           "should have type ChipInfo, or else object all.support in the global environment has been modified")
    }
  }
  if(!exists("all.support",envir=globalenv())) { stop("ChipInfo data object 'all.support' not found") }  
  all.support <- get("all.support",envir=globalenv())  # change global environment to namespace of iChip package
  nnn <- (names(all.support))
  if(length(nnn)>1) { 
    # multiple objects in file, probably different builds
    any37 <- c(grep("37",nnn),grep("hg19",nnn))
    any36 <- c(grep("36",nnn),grep("hg18",nnn))
    if(build=="hg19" & length(any37)>0) { all.support <- all.support[[nnn[any37[1]]]]}
    if(build=="hg18" & length(any36)>0) { all.support <- all.support[[nnn[any36[1]]]]}
  }
  return(all.support)
}




#' Convert from chip ID labels to dbSNP rs-ids
#' 
#' Most SNPs will have an 'rs-id' from dbSNP/HapMap, and these are often the standard for reporting or
#' annotation lookup. These can differ from the IDs used on the chip. This functions looks at the current
#' snp support (ChipInfo object) and returns rs-ids in place of chip IDs.
#' @param ids character, meant to be a list of chip ids, but if rs-ids are present they will not be altered.
#' @return A character vector of SNP rs-ids, where the input was chip ids, rs-ids or a mixture, any text
#' other than this will result in NA values being returned in the character vector output.
#' @export
#' @seealso rs.to.id, GENE.to.ENS, ENS.to.GENE
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' id.to.rs(c("imm_11_2138800","rs9467354","vh_1_1108138")) # middle one is already a rs.id
id.to.rs <- function(ids) {
  ids <- clean.snp.ids(ids)
  all.support <- chip.support()
  if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found") }  ## load object: all.support [snp support for whole chip]
  rsvec <- mcols(all.support)$rs.id[match(ids,rownames(all.support))]
  rsvec2 <- mcols(all.support)$rs.id[match(ids,mcols(all.support)$rs.id)]
  rsvec[is.na(rsvec)] <- rsvec2[is.na(rsvec)]
  return(rsvec)
}


#' Convert from dbSNP rs-ids to chip ID labels
#' 
#' Most SNPs will have an 'rs-id' from dbSNP/HapMap, and these are often the standard for reporting or
#' annotation lookup. These can differ from the IDs used on the chip. This functions looks at the current
#' snp support (ChipInfo object) and looks up chip IDs based on rs-ids.
#' @param ids character, meant to be a list of rs-ids, but if chip-ids are present they will not be altered.
#' @param multi.list logical, some rs-ids could map to multiple chip ids. It is recommended that if that is
#' the case then a letter should be appended to duplicate rs-ids to make them unique in the ChipInfo object,
#' e.g, rs1234, rs1234b, rs1234c, etc. If multi.list is TRUE, then the id list will be returned as a list,
#' and any time an rs-id is entered without a letter suffix, all possible corresponding chip ids will be
#' listed.
#' @return A character vector of SNP chip-ids, where the input was rs-ids, chip-ids or a mixture, any text
#' other than this will result in NA values being returned in the character vector output. Or, if multi-list
#' is true, then returns a list instead, which takes more than 1 value where there are multiple chip-ids 
#' with the same rs-id; if there are no such rs-id duplicates the result will still be a list.
#' @export
#' @seealso id.to.rs, GENE.to.ENS, ENS.to.GENE
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @examples
#' rs.to.id(c("rs689","rs9467354","rs61733845"))  # middle one has no chip id
#' test.ids <- c("rs61733845","rs2227313","rs11577783","rs3748816","rs12131065","rs3790567","rs2270614")
#' rs.to.id(test.ids, multi.list=TRUE) # list with duplicates
rs.to.id <- function(rs.ids,multi.list=FALSE) {
  rs.ids <- clean.snp.ids(rs.ids)
  all.support <- chip.support()
  if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found") }  ## load object: all.support [snp support for whole chip]
  if(multi.list) {
    X0 <- rmv.trail(rs.ids)
    X1 <- paste0(X0,"b"); X2 <- paste0(X0,"c"); X3 <- paste0(X0,"d"); X4 <- paste0(X0,"a")
    X <- cbind(rs.to.id(X0),rs.to.id(X1),rs.to.id(X2),rs.to.id(X3),rs.to.id(X4))
    idvec <- apply(X,1,function(x) { unique(narm(x)) })
    if(!is.list(idvec)) { idvec <- as.list(idvec) }
      #warning("'multi.list' option was used, but no duplicate rs-ids found, so returning a vector, not a list") }
  } else {
    idvec <- rownames(all.support)[match(rs.ids,mcols(all.support)$rs.id)]
    idvec2 <- rownames(all.support)[match(rs.ids,rownames(all.support))]
    idvec[is.na(idvec)] <- idvec2[is.na(idvec)]
  }
  return(idvec)
}



#' Find chromosome for SNP ids, gene name or band
#' 
#' Allows retrieval of the chromosome associated with a SNP-id, HGNC gene label, karyotype band,
#' or vector of such ids. For SNPs the ids can be either chip ids, or rs-ids, but must be contained
#'  in the current annotation. Default behaviour is to assume 'id' are SNP ids, but if none are
#'  found in the SNP annotation, the id's will be passed to functions Pos.gene() and Pos.band() to
#'  see whether a result is found. This latter step will only happen if no SNP ids are retreived in
#'  the first instance, and if snps.only=TRUE, then genes and bands will not be searched and NA's 
#'  returned. If you are repeatedly searching for chromosomes for genes/bands, using the dedicated 
#'  Pos.gene and Pos.band functions would be slightly faster than relying on the fallback behaviour
#'  of the Chr() function.  See documentation for these functions for more information. The build
#'  used will be that in the current ChipInfo object.
#' @param ids character, a vector of rs-ids or chip-ids representing SNPs in the current ChipInfo
#'  annotation, or gene ids, or karyotype bands
#' @param dir character, only relevant when gene or band ids are entered, in this case 'dir' is the location
#' to download gene and cytoband information; if left as NULL, depending on the value of 
#' getOption("save.annot.in.current"), the annotation will either be saved in the working directory to 
#' speed-up subsequent lookups, or deleted after use.
#' @param snps.only logical, if TRUE, only search SNP ids, ignore the possibility of genes/cytobands.
#' @return A character vector of Chromosomes for each ids, with NA values where no result was found.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Pos
#' @examples
#' Chr(c("rs689","rs9467354","rs61733845"))
#' Chr("CTLA4")
#' Chr("13q21.31")
#' Chr(c("CTLA4","PTPN22"),snps.only=TRUE) # fails as these are genes
#' Chr(c("rs689","PTPN22","13q21.31")) # mixed input, will default to SNPs, as at least 1 was found
Chr <- function(ids,dir=NULL,snps.only=FALSE) {
  ic.chr <- function(ids) {
    ic.ids <- clean.snp.ids(ids)
    all.support <- chip.support()
    if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found")  }  ## load object: all.support [snp support for whole chip]
    outlist <- chr(all.support)[match(ic.ids,rownames(all.support))]
    return(outlist)
  }
  query <- rs.to.id(ids)
  if(!snps.only & all(is.na(query))) { 
    ## unless the 'snps.only' function is set, then if it looks like we have not been handed
    ## snp ids, then check for band ids or gene ids instead
    chr.ify <- function(X) { tolower(names(X)); return(gsub("chr","",X["chr"],ignore.case=TRUE)) }
    numpqs <- (length(grep("q",ids))+length(grep("p",ids)))
    if(numpqs==length(ids)) { try.band <- T } else { try.band <- F }
    if(try.band) {
      suppressWarnings(test <- Pos.band(ids,dir=dir))
      if(!is.null(test)) { return(chr.ify(test)) }
    }
    suppressWarnings(test <- Pos.gene(ids,dir=dir))
    if(!is.null(test)) { return(chr.ify(test)) } 
  } 
  ic <- ic.chr(query)
  return(ic)
}


#' Find the chromosome position for SNP ids, gene name or band
#' 
#' Allows retrieval of the the chromosome position associated with a SNP-id, HGNC gene label, 
#'  karyotype band, or vector of such ids. For SNPs the ids can be either chip ids, or rs-ids,
#'  but must be contained in the current annotation. Default behaviour is to assume 'id' are 
#'  SNP ids, but if none are found in the SNP annotation, the id's will be passed to functions
#'  Pos.gene() and Pos.band() to see whether a result is found. This latter step will only happen
#'  if no SNP ids are retreived in the first instance, and if snps.only=TRUE, then genes and bands
#'  will not be searched and NA's returned. If you are repeatedly searching for positions for 
#'  genes/bands, using the dedicated Pos.gene() and Pos.band() functions would be slightly faster
#'  than relying on the fallback behaviour of the Pos() function. Note that the position for
#'  genes and bands are not a single point, so the result will be a range with start and end, 
#'  see 'values' below. See documentation for these functions for more information.
#' @param ids character, a vector of rs-ids or chip-ids representing SNPs in the current ChipInfo annotation,
#'  or gene ids, or karyotype bands
#' @param dir character, only relevant when gene or band ids are entered, in this case 'dir' is the location
#' to download gene and cytoband information; if left as NULL, depending on the value of 
#' getOption("save.annot.in.current"), the annotation will either be saved in the working directory to 
#' speed-up subsequent lookups, or deleted after use.
#' @param snps.only logical, if TRUE, only search SNP ids, ignore the possibility of genes/cytobands.
#' @return When ids are SNP ids, returns a numeric vector of positions for each id, with NA values
#'  where no result was found. When ids are genes or karyotype bands, will return a data.frame with
#'  columns 'chr' [chromosome], 'start' [starting position of feature], 'end' [end position of feature], 
#'  and the band without the chromosome prefix, if ids are bands. Note that this function cannot
#'  retrieve multiple ranges for a single gene (e.g, OR2A1), which means you'd need to use Pos.gene().
#'  The coordinates used will be of version getOption(ucsc="hg18"), or ucsc(chip.support()), which
#'  should be equivalent.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr
#' @examples
#' Pos(c("rs689","rs9467354","rs61733845"))
#' Pos("CTLA4") # returns a range
#' Pos("13q21.31") # returns a range
#' Pos(c("CTLA4","PTPN22"),snps.only=TRUE) # fails as these are genes
#' Pos(c("rs689","PTPN22","13q21.31")) # mixed input, will default to SNPs, as at least 1 was found
Pos <- function(ids,dir=NULL,snps.only=FALSE) {
  all.support <- chip.support()
  ic.pos <- function(ic.ids) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { stop("ChipInfo data object 'all.support' not found") }  ## load object: all.support [snp support for whole chip]
    outlist <- start(all.support)[match(ic.ids,rownames(all.support))]
    return(outlist)
  }
  query <- rs.to.id(ids)
  if(!snps.only & all(is.na(query))) { 
    ## unless the 'snps.only' function is set, then if it looks like we have not been handed
    ## snp ids, then check for band ids or gene ids instead
    numpqs <- (length(grep("q",ids))+length(grep("p",ids)))
    if(numpqs==length(ids)) { try.band <- T } else { try.band <- F }
    if(try.band) {
      suppressWarnings(test <- Pos.band(ids,dir=dir))
      if(!is.null(test)) { return(test) }
    }
    suppressWarnings(test <- Pos.gene(ids,dir=dir))
    if(!is.null(test)) { return(test) } 
  } 
  ic <- ic.pos(query)
  return(ic)
}



#' Order rs-ids or ichip ids by chrosome and position
#' 
#' Simple function to sort a character list of SNP ids into genome order.
#' @param ids character, vector of SNP rs-ids or chip-ids, see rs.to.id()
#' @return the same vector 'ids', sorted by genome position
#' @seealso rs.to.id, id.to.rs, Chr, Pos
#' @export
#' @examples
#' snp.ids <- c("rs3842724","imm_11_2147527","rs689","rs9467354","rs61733845")
#' Chr(snp.ids) # shows each is on a different chromosome
#' Pos(snp.ids)
#' ids.by.pos(snp.ids)
#' Chr(ids.by.pos(snp.ids))
#' Pos(ids.by.pos(snp.ids))
ids.by.pos <- function(ids) {
  if(!is.character(ids)) { stop("ids must be a character vector") }
  pp <- Pos(ids)
  if(any(is.na(pp))) { stop("invalid id list, 'ids' must all be valid, with position information in the ChipInfo object") }
  ids <- ids[order(pp)]
  cc <- Chr(ids)
  ids <- ids[order.chr(cc)]
  return(ids)
}




#' Find the chromosome, start and end position for gene names
#' 
#' Allows retrieval of the the chromosome position associated with a HGNC gene label, 
#'  or vector of such labels. Note that the position returned for genes is not a 
#'  single point as for SNPs, so the result will be a chromosome, then a position range with
#'  start and end.
#' @param genes character, a vector of gene ids
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download gene annotation information to; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, else 
#' a data.frame
#' @param band logical, whether to include band/stripe in returned object
#' @param one.to.one logical, some genes have split ranges, TRUE merges these to give only 1 range 
#' per gene, NB: this is the default behaviour when using the more general Pos() function
#' @param remap.extra logical, if TRUE genes with chromosome annotation 'c6_cox' and 'c6_QBL' will
#'  be mapped to chromosome 6, and 'NT_xxxx' chromosome labels will all be mapped to 'Z_NT', etc
#' @param discard.extra logical, if TRUE then any gene hit with chromosome not in 1:22, X, Y, XY, MT, 
#' will be discarded.
#' @param warnings logical, whether to show warnings when some/all ids are not matched to the 
#' reference
#' @return Returns a data.frame with columns 'chr' [chromosome], 'start' [starting position of the
#'  gene],'end' [end position of the gene], or if bioC=TRUE, then returns a GRanges object with
#'  equivalent information, and if band=TRUE, then an extra column is added with band information
#'  If returning a data.frame, then it will be in the same order as 'genes'. If bioC=TRUE, then
#'  the result will be in genome order, regardless of the order of 'genes'.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.band, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' Pos.gene(c("CTLA4","PTPN22"))
#' Pos.gene("MICA",build=36)
#' Pos.gene("MICA",build=37)
#' Pos.gene(c("CTLA4","PTPN22"),bioC=TRUE,band=TRUE)
#' Pos.gene(c("CTLA4","OR2A1"),one.to.one=TRUE) # OR2A1 is split over two ranges
#' Pos.gene(c("CTLA4","OR2A1"),one.to.one=FALSE)
#' Pos.gene("RNU2-1",one.to.one=FALSE,bioC=T,remap.extra=FALSE) # RNU2-1 is split over several ranges
Pos.gene <- function(genes,build=NULL,dir=NULL,bioC=FALSE,band=FALSE,one.to.one=TRUE,
                     remap.extra=FALSE,discard.extra=TRUE,warnings=TRUE) {
  if(!is.character(genes)) { stop("'genes' must be a character vector of gene names") }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  char.lim <- 100
  ga <- get.gene.annot(dir=dir,build=build,bioC=bioC,one.to.one=one.to.one,
                       remap.extra=remap.extra,discard.extra=discard.extra,GRanges=TRUE)
  typ <- is(ga)[1]
  if(typ=="GRanges") {  mt <- match(genes,mcols(ga)$gene) } else { mt <- match(genes,ga$gene) }
  failz <- paste(genes[is.na(mt)],collapse=", "); if(nchar(failz)>char.lim) { failz <- paste(substr(failz,1,char.lim),",...",sep="") }
  if(length(mt)<1 | all(is.na(mt))) { 
    if(warnings) { warning("did not find any 'genes' features: ",failz) }; return(NULL) }
  if(any(is.na(mt))) { 
    cnt <- length(which(is.na(mt)))
    if(warnings) { warning("did not find the following ",cnt," 'genes' features: ",failz) }
  }
  if(!one.to.one) { 
    order.important <- FALSE  # not currently implemented
    if(order.important) {
     if(length(genes)>0) {
       mt <- NULL
       for(gg in 1:length(genes)) {
         if(typ=="GRanges") { 
           mt <- c(mt,which(mcols(ga)$gene %in% genes[gg]))
         } else {
           mt <- c(mt,which(ga$gene %in% genes[gg]))
         }
       }
     } else { warning("length of genes entered was zero?") }     
    } else {
      if(typ=="GRanges") { 
        mt <- which(mcols(ga)$gene %in% genes)
      } else {
        mt <- which(ga$gene %in% genes)
      }
    }
  }
  outlist <- ga[(narm(mt)),]
  if(typ=="GRanges") {
    if(!band) { mcols(outlist) <- mcols(outlist)[,-which(colnames(mcols(outlist)) %in% "band")]  }
    cnn <- colnames(mcols(outlist)); rnn <- mcols(outlist)$gene
  } else {
    if(!band) { outlist <- outlist[,-which(colnames(outlist) %in% "band")] }
    cnn <- colnames(outlist); rnn <- outlist$gene
  }
  if(one.to.one & ("gene" %in% cnn) & !anyDuplicated(genes)) {
    rownames(outlist) <- rnn
    if(all(genes %in% rownames(outlist))) {
      outlist <- outlist[match(genes,rnn),]
    }
    # outlist <- outlist[,-which(colnames(outlist) %in% "gene")]
  }
  if(bioC) { outlist <- as(outlist,"GRanges") }
  return(outlist)
}


#' Find the chromosome, start and end position for cytoband names
#' 
#' Allows retrieval of the the chromosome position of a karyotype/cytoband label, 
#'  or vector of such labels. Note that the position returned for bands is not a 
#'  single point as for SNPs, so the result will be a chromosome, then a position range with
#'  start and end, and lastly the band without the chromosome prefix
#' @param bands character, a vector of cytoband labels, chromosome[p/q]xx.xx ; 
#' e.g, 13q21.31, Yq11.221, 6p23, etc
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download cyto annotation information; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, else 
#' a data.frame
#' @return Returns a data.frame with columns 'chr' [chromosome], 'start' [starting position of the
#'  gene],'end' [end position of the gene] and 'band' [band without the chromosome prefix],
#'  or if bioC=TRUE, then returns a GRanges object with equivalent information.
#'  If returning a data.frame, then it will be in the same order as 'bands'. If bioC=TRUE, then
#'  the result will be in genome order, regardless of the order of 'bands'.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.gene, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' Pos.band("1p13.2")
#' Pos.band("Yq11.221",build=36)
#' Pos.band("Yq11.221",build=37)
#' Pos.band(c("13q21.31","1p13.2","2q33.2","6p23"),bioC=TRUE)
Pos.band <- function(bands,build=NULL,dir=NULL,bioC=FALSE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  char.lim <- 100
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.cyto(build=build,bioC=bioC,dir=dir,GRanges=FALSE)
  mt <- match(bands,rownames(ga))
  failz <- paste(bands[is.na(mt)],collapse=", "); if(nchar(failz)>char.lim) { failz <- paste(substr(failz,1,char.lim),",...",sep="") }
  msg <- ("format for bands is: chromosome[p/q]xx.xx ; e.g, 13q21.31, Yq11.221, 6p23, etc")
  if(length(mt)<1 | all(is.na(mt))) { 
    warning("did not find any 'bands' features: ",failz) ; warning(msg); return(NULL) }
  if(any(is.na(mt))) { 
    cat("format for bands is: chromosome[p/q]xx.xx ; e.g, 13q21.31, Xq27.1, 6p23, etc")
    cnt <- length(which(is.na(mt)))
    warning("did not find the following ",cnt," 'bands' features: ",failz) ; warning(msg)
  }
  outlist <- ga[sort(mt[!is.na(mt)]),]
  if(any(colnames(outlist) %in% "negpos")) { outlist <- outlist[,-which(colnames(outlist) %in% "negpos")] }
  if(all(bands %in% rownames(outlist)) & !bioC) {
    outlist <- outlist[bands,]
  }
  return(outlist)
}


#' Retrieve the cytoband(s) for snp ids, genes or locations
#' 
#' Allows retrieval of the the cytoband/karyotype label, based on multiple
#'  possible input featues, including SNP chip or rs-ids, HGNC gene labels, GRanges or
#'  RangedData object, chromosome and position vectors. The most robust way to use the
#'  function is to use the parameter names to imply the type of input, e.g, use the 'genes'
#'  parameter to input gene labels, the 'snps' parameter to enter SNP ids, etc. However,
#'  if you enter the first argument as a GRanges or RangedData object instead of using the
#'  'ranges' argument, this will be detected and automatically moved to the 'ranges' parameter.
#' @param genes character, an optional vector of gene ids, or RangedData/GRanges object
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions to retrieve the band from
#' @param ranges optional GRanges or RangedData object describing positions for which we want bands
#' @param snps optional SNP ids, e.g, chip ids or rs-ids, to retrieve the band they fall within
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download cyto annotation information; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param ... further arguments to Band.gene if entering gene names, or further arguments to 
#' Band.pos if entering ranges, or chr, pos/start/end
#' @return Returns a vector of bands, if any entries span more than one band, the bands will be
#' concatenated as character type, delimited by semicolons (;)
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.gene, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' Band(chr=1,pos=1234567) # using chr,pos vectors
#' rd <- RangedData(ranges=IRanges(start=87654321,end=87654321),space=1)
#' gr <- as(rd,"GRanges")
#' Band(rd)    # using RangedData, autodetects this parameter should be 'ranges' not 'genes'
#' Band(ranges=gr) # using GRanges
#' Band("SLC6A4")  # serotonin gene [5-HTT]
#' a.few.snps <- c("rs3842724","imm_11_2147527","rs9467354")
#' Band(a.few.snps) # using SNP ids in the 'genes' parameter (still works!)
#' Band(snps=a.few.snps) # using SNP ids with the dedicated 'snps' parameter is quicker
#' Band(chr="X",pos=8000000)
#' # Band() with longer ranges  #
#' Band(chr=12,start=40000000,end=50000000,build="hg19") # concatenates if range spans multiple bands
#' Band(chr=12,start=40000000,end=50000000,build="hg18") # one extra band in the older annotation
Band <- function(genes=NULL,chr=NULL,ranges=NULL,snps=NULL,build=NULL,dir=NULL,...) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(all(is.character(genes))) { 
    Band.gene(genes=genes,build=build,dir=dir,...,warnings=FALSE)
  } else {
    if(!is.null(snps)) {
      pp <- Pos(snps); ch <- Chr(snps)
      Band.pos(chr=ch,pos=pp,build=build,dir=dir)
    } else {
      if((is(genes)[1] %in% c("RangedData","GRanges")) & is.null(ranges) ) { ranges <- genes }
      Band.pos(chr=chr,ranges=ranges,build=build,dir=dir,...)
    }
  }
}


#' Retrieve the cytoband(s) for genes labels
#' 
#' Allows retrieval of the the cytoband/karyotype label for HGNC gene labels.
#' @param genes character, an optional vector of gene ids, or RangedData/GRanges object
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download cyto annotation information; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param append.chr logical, it is typical that the chromosome character preceeds cytoband labels,
#'  but if this parameter is set to FALSE, it will be left off.
#' @param data.frame logical, if data.frame is true, instead of returning a vector of full cytoband
#' labels, a data.frame will be returned.
#' @param warnings logical, if warnings=FALSE and SNP ids are entered instead of Gene labels,
#' then the function will automatically detect this and return the result of Band(snps='genes')
#' @return Returns a vector of bands, if any entries span more than one band, the bands will be
#' concatenated as character type, delimited by semicolons (;). If data.frame is true, instead of 
#' returning a vector of full cytoband labels, a data.frame will be returned with a 'chr'
#' [chromosome] column, 'band' cytoband label  without the chromosome prefix, and rownames 
#' equal to 'genes'
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.gene, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' a.few.snps <- c("rs3842724","imm_11_2147527","rs9467354")
#' Band.gene("HLA-C") # using chr,pos vectors
#' Band.gene(a.few.snps)  # fails with warning as these are SNPs, not genes
#' Band.gene(a.few.snps,warnings=FALSE) # with warnings=FALSE this recognises snps were entered and carries on
#' Band.gene("SLC6A4")  # serotonin gene [5-HTT]
#' Band.gene("SLC6A4",append.chr=F)
#' Band.gene("SLC6A4",data.frame=T)
Band.gene <- function(genes,build=NULL,dir=getwd(),append.chr=TRUE,data.frame=FALSE,warnings=TRUE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  pg <- Pos.gene(genes=genes,build=build,dir=dir,bioC=F,band=TRUE,one.to.one=TRUE,warnings=warnings)
  ## a bit contradictory what i've done here with 'warnings' :
  if(length(pg)<1 & !warnings) { 
    #if(warnings) { warning("perhaps parameter 'genes' should have been 'snps'?") } 
    return(Band(snps=genes)) 
  }
  if(data.frame) {
    if(all(c("start","end") %in% colnames(pg))) { pg <- pg[,-which(colnames(pg) %in% c("start","end"))] }
    if(all(c("gene") %in% colnames(pg))) { rownames(pg) <- pg[["gene"]] ; pg <- pg[,-which(colnames(pg) %in% c("gene"))] }
    out <- pg
  } else {
    if(append.chr) {
      out <- paste(pg[["chr"]],pg[["band"]],sep="")
    } else {
      out <- pg[["band"]]
    }
  }
  return(out)    
}



#' Find the gene(s) overlapping a chromosome location
#' 
#' Allows retrieval of genes intersected by a chromosome and position, which can be entered
#' using chr, pos/start/end vectors, or a RangedData or GRanges object
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions to retrieve the possible overlapping gene(s)
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' start or end if this is entered, and vice-versa
#' @param start integer, an optional vector of start points for chromosome ranges
#' @param end integer, an optional vector of end points for chromosome ranges
#' @param ranges optional GRanges or RangedData object describing positions for which we want genes,
#' removing the need to enter chr, pos, start or end
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download gene annotation information to; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, or
#' RangedData if 'ranges' is RangedData, else a data.frame
#' @param one.to.one logical, whether to concatenate multiple hits for the same range into one result,
#' or spread the result over multiple lines, one for each gene overlapped
#' @return Returns a set of genes separated by semicolons (if more than one) for each range entered.
#' If bioC=TRUE, returns the equivalent as a GRanges object, unless a RangedData object was used
#' for the ranges parameter, in which case a RangedData object would be returned. If one.to.one is
#' FALSE, then instead of concatenating multiple genes into one line per range, each is listed 
#' separately as a new row, with an index added to correspond to the original input order of ranges,
#' if bioC=TRUE; or just adds additional elements to the resulting vector if bioC=FALSE.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.band, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' Gene.pos(chr=6, start=31459636, end=31462760)
#' Gene.pos(chr=22, pos=3452345) # no gene here
#' Gene.pos(Chr("rs689"),Pos("rs689")) # combine with Chr() and Pos() to find gene(s) for a the rs689 SNP
#' Gene.pos(chr=1,start=114000000,end=115000000,build="hg19") # multiple genes in range
#' Gene.pos(chr=1,start=114000000,end=115000000,one.to.one=FALSE) # list separately
#' Gene.pos(Pos.gene(c("CTLA4","PTPN22"),bioC=TRUE)) # use the ranges object returned by Pos.gene()
Gene.pos <- function(chr=NA,pos=NA,start=NA,end=NA,ranges=NULL,build=NULL,dir=NULL,bioC=FALSE,one.to.one=TRUE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is(chr)[1] %in% c("RangedData","GRanges"))  { ranges <- chr } # in case first parameter is used
  typ <- is(ranges)[1]
  if(!typ %in% c("RangedData","GRanges")) {
    if(any(!is.na(pos))) { start <- pos; end <- pos }
    if(length(chr)==1 & length(start)>1) { chr <- rep(chr,times=length(start)) }
    if(length(chr)!=length(start)) { stop("chr vector must have same length as pos or start/end") }
    if(any(is.na(chr))) { stop("cannot have chr=NA") }
    Pos <- matrix(ncol=2,nrow=length(start))
    for (cc in 1:length(start)) {
      Pos[cc,] <- force.chr.pos(Pos=c(start[cc],end[cc]),Chr=chr[cc],dir=dir,build=build)
    }
    #if(any(tolower(substr(chr,1,3))!="chr")) { chr <- gsub("chr",chr,sep="") }
    #chr <- gsub("chrchr","chr",chr)
    testData <- RangedData(ranges=IRanges(start=Pos[,1],end=Pos[,2]),space=chr,index=1:length(chr),universe=build[1])
    testData <- toGenomeOrder2(testData,strict=T)
  } else {
    testData <- ranges # set.chr.to.char(ranges)
    if(typ=="GRanges") { testData <- as(ranges,"RangedData") }
    #chr <- chr(testData)
    if(!bioC) { warning("bioC was set false, but ranges argument in use so overriding") ; bioC <- T }
    if("index" %in% colnames(testData)) { 
      warning("'index' is a reserved column name for ranges objects passed to this function so will be replaced. Consider renaming this column if this is undesired") }
    testData[["index"]] <- 1:nrow(testData)
  }
  testData <- set.chr.to.numeric(testData,keep=T)
  #return(testData)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.gene.annot(build=build,dir=dir,GRanges=FALSE)
  #ga <- set.chr.to.numeric(ga,keep=F)
  newDataList <- vector("list",nrow(testData))
  overlaps <- findOverlaps(testData,ga)
  genez <- ga$gene[subjectHits(overlaps)]
  indexz <- queryHits(overlaps)
  #prv(indexz,genez)
  if(length(indexz)<1) { return(NA) }
  if(!one.to.one) {
    newData <- testData[queryHits(overlaps),]
    newData[["gnm.index"]] <- indexz
    newData[["gene"]] <- genez
  } else {
    out <- tapply(genez,factor(indexz),c,simplify=FALSE)
    out <- sapply(out,function(X) { X <- narm(unique(X)); X <- X[X!=""] ; paste(X,collapse=";") })
    newData <- testData
    newData[["gene"]] <- rep("intergenic",nrow(newData))
    if(!is.null(names(out))) {
      newData[["gene"]][as.numeric(names(out))] <- out
    } else {
      newData[["gene"]] <- out
    }
  }
  if(bioC) {
    if(typ!="RangedData") { newData <- as(newData,"GRanges") }
    return(newData)
  } else {
    OO <- newData[["gene"]][order(newData[["index"]])]
    OO <- narm(unique(OO)); OO <- OO[OO!=""]
    return(OO)
  }
  #if(!all(chr %in% chr2(ga))) { stop("invalid chromosome(s) entered") } # redundant i think
}





#' Find the cytoband(s) overlapping a chromosome location
#' 
#' Allows retrieval of cytobands/karyotypes intersected by a chromosome and position, which can be 
#' entered using chr, pos/start/end vectors, or a RangedData or GRanges object
#' @param chr character, an optional vector of chromosomes to combine with 'pos' or 'start'+'end'
#' (enter in ...) to describe positions to retrieve the possible overlapping cytoband(s)
#' @param pos integer, an optional vector of chromosome positions (for SNPs), no need to enter
#' start or end if this is entered, and vice-versa
#' @param start integer, an optional vector of start points for chromosome ranges
#' @param end integer, an optional vector of end points for chromosome ranges
#' @param ranges optional GRanges or RangedData object describing positions for which we want bands,
#' removing the need to enter chr, pos, start or end
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param dir character, 'dir' is the location to download gene annotation information to; if left as
#'  NULL, depending on the value of getOption("save.annot.in.current"), the annotation will either
#'  be saved in the working directory to speed-up subsequent lookups, or deleted after use.
#' @param bioC logical, if true then return position information as a GRanges object, or
#' RangedData if 'ranges' is RangedData, else a data.frame
#' @param one.to.one logical, whether to concatenate multiple hits for the same range into one result,
#' or spread the result over multiple lines, one for each cytoband overlapped
#' @return Returns a set of cytobands separated by semicolons (if more than one) for each range entered.
#' If bioC=TRUE, returns the equivalent as a GRanges object, unless a RangedData object was used
#' for the ranges parameter, in which case a RangedData object would be returned. If one.to.one is
#' FALSE, then instead of concatenating multiple cytobands into one line per range, each is listed 
#' separately as a new row, with an index added to correspond to the original input order of ranges,
#' if bioC=TRUE; or just adds additional elements to the resulting vector if bioC=FALSE.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.band, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' Band.pos(chr=6, start=31459636, end=31462760)
#' Band.pos(chr=22, pos=3452345) 
#' Band.pos(Chr("rs689"),Pos("rs689")) # combine with Chr() and Pos() to find the cytoband for a the rs689 SNP
#' Band.pos(chr=1,start=110000000,end=120000000,build="hg19") # multiple cytobands in range
#' Band.pos(chr=1,start=110000000,end=120000000,one.to.one=FALSE) # list separately
#' Band.pos(Pos.band(c("13q21.31","1p13.2"),bioC=TRUE)) # use the ranges object returned by Pos.band()
#' # note that three ranges are returned for each entry as the start/end overlap the adjacent ranges
Band.pos <- function(chr=NA,pos=NA,start=NA,end=NA,ranges=NULL,build=NULL,dir=NULL,bioC=FALSE,one.to.one=TRUE) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is(chr)[1] %in% c("RangedData","GRanges"))  { ranges <- chr } # in case first parameter is used
  typ <- is(ranges)[1]
  if(!typ %in% c("RangedData","GRanges")) {
    if(any(!is.na(pos))) { start <- pos; end <- pos }
    if(length(chr)==1 & length(start)>1) { chr <- rep(chr,times=length(start)) }
    if(length(chr)!=length(start)) { stop("chr vector must have same length as pos or start/end") }
    if(any(is.na(chr))) { stop("cannot have chr=NA") }
    Pos <- matrix(ncol=2,nrow=length(start))
    for (cc in 1:length(start)) {
      Pos[cc,] <- force.chr.pos(Pos=c(start[cc],end[cc]),Chr=chr[cc],dir=dir,build=build)
    }
    #if(any(tolower(substr(chr,1,3))!="chr")) { chr <- gsub("chr",chr,sep="") }
    #chr <- gsub("chrchr","chr",chr)
    testData <- RangedData(ranges=IRanges(start=Pos[,1],end=Pos[,2]),space=chr,index=1:length(chr),universe=build[1])
    testData <- toGenomeOrder2(testData,strict=T)
  } else {
    if(typ=="GRanges") { ranges <- as(ranges,"RangedData") }
    testData <- ranges # set.chr.to.char(ranges)
    if("index" %in% colnames(testData)) { warning("'index' is a reserved column name for ranges objects passed to this function so will be replaced. Consider renaming this column if this is undesired") }
    testData[["index"]] <- 1:nrow(testData)
  }
  testData <- set.chr.to.numeric(testData)
  #return(testData)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  cyto <- get.cyto(build=build,bioC=TRUE,dir=dir,GRanges=FALSE)
  cyto <- set.chr.to.numeric(cyto,keep=F)
  #if(any(colnames(cyto) %in% "negpos")) { cyto <- cyto[,-which(colnames(cyto) %in% "negpos")] }
  newDataList <- vector("list",nrow(testData))
  overlaps <- findOverlaps(testData,cyto)
  bandz <- rownames(cyto)[subjectHits(overlaps)]
  indexz <- queryHits(overlaps)
  if(length(indexz)<1) { return(NA) }
  if(!one.to.one) {
    newData <- testData[queryHits(overlaps),]
    newData[["gnm.index"]] <- indexz
    newData[["band"]] <- bandz
  } else {
    out <- tapply(bandz,factor(indexz),c,simplify=FALSE)
    if(one.to.one) { 
      out <- sapply(out,function(X) { paste(X,collapse=";") }) 
      newData <- testData
      #newData[["band"]] <- out 
      newData[["band"]] <- rep("",nrow(newData))
      if(!is.null(names(out))) {
        newData[["band"]][as.numeric(names(out))] <- out
      } else {
        newData[["band"]] <- out
      }
    }
  }
  if(bioC) {
    if(typ!="RangedData") { newData <- as(newData,"GRanges") }
    return(newData)
  } else {
    return(newData[["band"]][order(newData[["index"]])])
  }
  #if(!all(chr %in% chr2(ga))) { stop("invalid chromosome(s) entered") } # redundant i think
}



#' Returns the A and B allele for SNP ids
#' 
#' For a set of chip ids or rs ids, returns a two column matrix containing the A and B allele. 
#' For snpStats objects the default is that A,B are coded in alphabetical order, so A,C; A,T; 
#' C,T; C,G are possible A,B pairs. Allele codes are specific to each dataset, so you should
#' upload your allele codes into the current ChipInfo object to make the alleles produced by
#' this function meaningful.
#' @param ids character, a list of chip ids or rs-ids as contained in the current ChipInfo object
#' @return Returns a two column matrix containing the A and B allele.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso Chr, Pos, Pos.band, Band, Band.gene, Band.pos, Gene.pos
#' @examples
#' snp.ids <- c("rs3842724","rs9729550","rs1815606","rs114582555","rs1240708","rs6603785")
#' AB(snp.ids) 
AB <- function(ids) {
  all.support <- chip.support()
  ic.ab <- function(ic.ids) {
    ic.ids <- clean.snp.ids(ic.ids)
    if(!exists("all.support")) { all.support <- chip.support() }  ## load object: all.support [snp support for whole chip]
    outlist <- cbind(A1(all.support)[match(ic.ids,rownames(all.support))],A2(all.support)[match(ic.ids,rownames(all.support))])
    return(outlist)
  }
  out <- matrix(nrow=length(ids),ncol=2)
  ic.ab.id <- ic.ab(rs.to.id(ids))
  return(ic.ab.id)
}


#' Retrieve SNP ids or positions in specified range
#' 
#' @param chr integer, chromosome, should be a number from 1 to 25, where 23,24,25 are X,Y,MT
#' Alternatively chr can be a RangedData or GRanges object in which case SNP lists will be returned
#' in a list for each row of the ranges object.
#' @param start integer, genomic position to define the start of the range to look for SNPs,
#'  should be between 1 and the length of the chromosome 'chr'
#' @param end integer, genomic position to define the end of the range to look for SNPs,
#'  should be between 1 and the length of the chromosome 'chr', and >= start
#' @param ids logical, if TRUE will return snp ids (chip ids, for rs-ids, use id.to.rs on the output), 
#' or if FALSE will return the chromosome positions of the SNPs.
#' @param build integer whether to use build 36/37 parameters, 36/37 is preferred, but can enter
#' using any form recognised by ucsc.sanitizer()
#' @param dir string, specify the directory to store or look for annotation (defaults to current)
#' @export
#' @return Set of SNP ids (when ids=TRUE), or otherwise genomic positions within chromosome 'chr', that
#' fall within the genomic range described by the chr, start, and end parameters. Alternatively, if
#' chr is a RangedData or GRanges object then multiple SNP lists will be returned
#' in a list for each row of the ranges object.
#' @examples
#' snps.in.range(1,9000000,10000000)
#' snps.in.range(10,19000000,20000000,ids=T)
#' snps.in.range(10,19000000,20000000,ids=F) # return positions instead of rs-ids
snps.in.range <- function(chr, start=NA, end=start, ids=TRUE) { 
  # ids - whether to return ichip SNP ids or positions
  if(is(chr)[1]=="RangedData" | is(chr)[1]=="GRanges") {
    chrz <- chr2(chr); stz <- start(chr); enz <- end(chr)
    output <- vector("list",nrow(chr))
    for(cc in 1:nrow(chr)) {
      output[[cc]] <- snps.in.range(chrz[cc],stz[cc],enz[cc],ids=ids)
    }
    names(output) <- rownames(chr)
    return(output)
  }
  if(length(chr)>1) { warning("chr should be length 1, using only first entry"); chr <- chr[1] }
  if(length(start)>1) { warning("start should be length 1, using only first entry"); start <- start[1] }
  if(length(end)>1) { warning("end should be length 1, using only first entry"); end <- end[1] }
  if(start>end) { warning("start was higher than end, so switching") }
  the.range <- sort(c(start,end))
  all.support <- chip.support()
  #if(!exists("work.dir")) { if(is.null(dir)) { work.dir <- getwd() } else { work.dir <- dir } }
  if(!exists("all.support")) { all.support <- chip.support() }  ## load object: all.support [snp support for whole chip]
  all.chr <- chr(all.support)
  all.pos <- start(all.support)[all.chr %in% chr]
  if(length(all.pos)<1) { warning("no positions found for 'chr' specified"); return(NULL) }
  validz <- which(all.pos>=the.range[1] & all.pos<=the.range[2])
  if(ids) {
    out <- rownames(all.support)[all.chr %in% chr][validz]
  } else {
    out <- start(all.support)[(all.chr %in% chr)][validz]
  }
  return(out)
}


#' Retrieve the 'n' closest SNP ids or positions near specified locus
#' 
#' @param chr integer, chromosome, should be a number from 1 to 25, where 23,24,25 are X,Y,MT
#' @param pos integer, genomic position, should be between 1 and the length of the chromosome 'chr'
#' @param n integer, the number of nearest SNPs to seek, if there aren't enough in the annotation
#' then NAs will fill the gaps to force the return value length to equal 'n'
#' @param side character, can be 'either', 'left' or 'right' and specifies which side of the 'pos'
#' to look for nearest snps (where left is decreasing genomic position and right is increasing)
#' @param ids logical, if TRUE will return snp ids (chip ids, for rs-ids, use id.to.rs on the output), 
#' or if FALSE will return the chromosome positions of the SNPs.
#' @param limit integer, a limit on the maximum distance from the position 'pos' can be specified
#' @param build integer whether to use build 36/37 parameters, 36/37 is preferred, but can enter
#' using any form recognised by ucsc.sanitizer()
#' @param dir string, specify the directory to store or look for annotation (defaults to current)
#' @export
#' @seealso exp.window.nsnp, nearest.gene
#' @return Set of SNP ids (when ids=TRUE), or otherwise genomic positions within chromosome 'chr'.
#' If the number of SNPs on the chromosome or the bounds of the 'side' and 'limit' parameters
#' restrict the number returned to less than 'n' then the return value will be padded with NAs.
#' @examples
#' nearest.snp(1,159000000,n=10) # return ids
#' nearest.snp(1,159000000,n=10,build=37)
#' nearest.snp(1,159000000,n=10,build=36,ids=F) # return positions
#' nearest.snp(1,159000000,n=10,build=37,ids=F)
#' nearest.snp(6,25000000,n=10,build=37,ids=F,side="left")  # only SNPs to the left of the locus
#' nearest.snp(6,25000000,n=10,build=37,ids=F,side="right") # only SNPs to the right of the locus
nearest.snp <- function(chr, pos, n=1, side=c("either","left","right"),ids=TRUE,limit=NULL,build=NULL) { 
  # ids - whether to return ichip SNP ids or positions
  if(length(chr)>1) { warning("chr should be length 1, using only first entry"); chr <- chr[1] }
  if(length(pos)>1) { warning("pos should be length 1, using only first entry"); pos <- pos[1] }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  all.support <- chip.support(build=build)
  if(!exists("all.support")) { all.support <- chip.support() }  ## load object: all.support [snp support for whole chip]
  side <- tolower(side[1]); 
  if(!side %in% c("either","left","right")) {
    side <- "either"; warning("invalid side argument, defaulting to 'either'") }
  if(!is.null(limit)) { if(!is.numeric(limit)) { limit <- NULL; warning("invalid limit argument, defaulting to NULL") } }
  all.chr <- chr(all.support)
  all.pos <- start(all.support)[all.chr %in% chr]
  if(length(all.pos)<1) { warning("no positions found for 'chr' specified"); return(NULL) }
  difz <- pos-all.pos
  all.true <- difz==difz
  if(is.null(limit)) { lfilt <- all.true } else { lfilt <- abs(difz)<=limit }
  if(side=="left") { filt <- difz>0 & lfilt }
  if(side=="right") { filt <- difz<0 & lfilt }
  if(side=="either") { filt <- all.true & lfilt }
  Difz <- abs(difz[filt])
  if(length(Difz)<n)  { warning("fewer than ",n," positions found for 'chr' specified (within 'limit'), NAs returned") }
  indx <- order(Difz)[1:n]
  if(ids) {
    out <- rownames(all.support)[all.chr %in% chr][filt][indx]
  } else {
    out <- start(all.support)[(all.chr %in% chr)][filt][indx]
  }
  return(out)
}



#' Obtain nearby SNP-lists within a recombination window
#' 
#' For a snp.id (or list), extend a window around that chromosome location in recombination
#' units (centimorgans) and return the list of SNPs from the current ChipInfo object that
#' lie in this window. This is a way of extracting SNPs in linkage disequilibrium with an 
#' index SNP, that could also be plausible causal candidates.
#' @param snpid.list character, list of snp-ids (e.g, rs-id or chip id) to obtain lists for.
#' SNPs must all be from the same chromosome - if ranges for SNPs spanning multiple ranges
#' are desired, you must use multiple calls. A warning will be given if SNPs from the same
#' karyotype band are entered as index SNPs, as in a typical GWAS analysis only one SNP would
#' be used like this from each region, ignore the warning if this is not the case for your
#' application.
#' @param cM numeric, the number of centimorgans to extend the window either side of each SNP
#' @param bp.ext numeric, optional number of base-pairs to extend the window by in addition
#' to the centimorgan extension
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to use The 
#' default when build is NULL is to use the build from the current ChipInfo annotation.
#' Note that build37 is slower as the recombination reference is in build 36, so conv.36.37() 
#' and conv.37.36() need to be used to retrieve the appropriate ranges for build 37.
#' @param excl.snps character, a list of rs-id or chip-ids of SNPs to exclude from the list
#' returned, as, for instance, they may have failed quality control such as call-rate.
#' @param name.by.bands give labels to each sublist returned by the karotype/cytoband name, 
#' but faster not to do this
#' @return Returns a list of vectors of snp-ids falling within the window(s) specified and not
#' in 'excl.snps'. Each snp in 'snpid.list' will correspond to an element in the list returned.
#' If name.by.bands is TRUE, then these list elements will each be named using the local
#' karyotype/cytoband location
#' @export
#' @seealso snps.in.range, get.recombination.map, recwindow, conv.37.36, conv.36.37, exp.window.nsnp
#' @examples
#' result <- get.nearby.snp.lists("rs900569")
#' get.nearby.snp.lists("rs900569",cM=0.2,excl.snps=result[[1]]) # get SNPs within 0.1-0.2cM
#' # note that the same query can return a different set with build 36 versus 37
#' get.nearby.snp.lists(c("rs689","rs4909944"),cM=0.001,build=36,name.by.bands=FALSE)
#' get.nearby.snp.lists(c("rs689","rs4909944"),cM=0.001,build=37,name.by.bands=FALSE) 
get.nearby.snp.lists <- function(snpid.list,cM=0.1,bp.ext=0,build=NULL,excl.snps=NULL,name.by.bands=TRUE) {
  #if(!exists("all.support")) { print(load("all.support.RData")) }
  all.support <- chip.support()
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  snpic.list <- rs.to.id(snpid.list)
  cyto <- get.cyto(dir=getwd(),GRanges=FALSE); cyto[["gene"]] <- rownames(cyto)
  #which.snps <- match(snpid.list,mcols(all.support)$rs.id)
  #if(any(is.na(which.snps))) { stop(paste("NAs in dbSNP match:",paste(snpid.list[is.na(which.snps)],collapse=","))) }
  snps.locs <- Pos(snpid.list)
  snps.chrs <- Chr(snpid.list)
  if(build=="hg18") {
    snps.locs36 <- snps.locs;  snps.locs37 <- conv.36.37(chr=snps.chrs,pos=snps.locs)[,"start"] 
  } else { 
    snps.locs37 <- snps.locs;  snps.locs36 <- conv.37.36(chr=snps.chrs,pos=snps.locs)[,"start"]
  }
  next.chr <- unique(Chr(snpid.list)); if(length(next.chr)>1) { stop("enter snpids from only 1 chromosome at a time!") }
  if(any(snps.locs!=sort(snps.locs))) { 
    warning("snp-ids not in position order, rearrangement is preferred but will attempt to continue")
    sort.back <- match(snps.locs,sort(snps.locs))
  } else { sort.back <- 1:length(snps.locs) }
  ddz <- snpic.list[duplicated(snpic.list)]
  if(length(ddz)>0) { warning("dup SNPs:",ddz,"\n") }
  snp.rd <- RangedData(ranges=IRanges(start=snps.locs,end=snps.locs,names=snpic.list),
                       space=rep(next.chr,length(snps.locs)))
  snp.rd <- toGenomeOrder(snp.rd,strict=T) # think it autosorts anyway, but just in case
  if(name.by.bands) {
    #snp.rd <- annot.cnv(snp.rd,gs=cyto,quiet=TRUE); colnames(snp.rd) <- "band"
    bands <- Band.pos(ranges=snp.rd) #   snp.rd$band
  }
  #prv(next.chr,cM,bp.ext)
  ## recwindow uses build36 only, so convert back afterwards
  nxt.window <- lapply(snps.locs36, function(X,...) { recwindow(start=X,...) },
                       chr=next.chr,window=cM,bp.ext=bp.ext,info=FALSE)
  if(build=="hg18") {
    st.window <- sapply(nxt.window, "[",1)
    en.window <- sapply(nxt.window, "[",2)
  } else {
    #prv(next.chr,nxt.window)
    nncc <- rep(next.chr,length(nxt.window))
    st.window <- conv.36.37(chr=nncc,pos=sapply(nxt.window, "[",1))[,"start"]
    en.window <- conv.36.37(chr=nncc,pos=sapply(nxt.window, "[",2))[,"start"]
  }
  pozz <- start(all.support)
  n.snps <- vector("list",length(st.window))
  for(cc in 1:length(st.window)) {
    n.snps[[cc]] <- which(chr(all.support)==next.chr &
                            pozz>=st.window[cc] & 
                            pozz<=en.window[cc] &
                            (!rownames(all.support) %in% excl.snps) &
                            (!mcols(all.support)$rs.id %in% excl.snps) 
    )
  }
  grp.labs <- lapply(n.snps,function(X) { rownames(all.support)[X] })
  if(name.by.bands) {
    if(length(unique(bands))!=length(bands)) { warning("these bands are not unique ==> ",paste(bands[duplicated(bands)],collapse=",")) }
    grpz <- 1:length(bands)
    names(grp.labs) <- paste(grpz,bands,sep=":")
  }
  grp.labs <- grp.labs[sort.back]
  return(grp.labs)
}



################## end support ##########################

