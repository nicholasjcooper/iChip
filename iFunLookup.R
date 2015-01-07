
##################################
## General Annotation Functions ##
##################################


#' Download GWAS hits from t1dbase.org
#' 
#' Retrieve human disease top GWAS hits from t1dbase in either build hg18 or hg19 coords (b36/37).
#' 28 Diseases currently available
#' @param disease integer (1-28), or character (abbreviation), or full name of one of the listed
#' diseases. A full list of options can be obtained by setting show.codes=TRUE.
#' @param build character, either "hg18" or "hg19". Will also accept build number, 36 or 37.
#' @param show.codes logical, if set to TRUE, instead of looking up t1dbase, will simply return
#' a table of available diseases with their index numbers and abbreviations.
#' @return A character vector of SNP rs-ids
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @references PMID: 20937630
#' @examples
#' get.t1dbase.snps(disease="CEL",build=36) # get SNP ids for celiac disease in build-36/hg18
#' get.t1dbase.snps(disease="AS") # get SNP ids for Ankylosing Spondylitis in build-37/hg19
#' get.t1dbase.snps(show.codes=TRUE) # show codes/diseases available to download
#' get.t1dbase.snps(disease=27) # get SNP ids for Alopecia Areata
#' get.t1dbase.snps("Vitiligo")
get.t1dbase.snps <- function(disease="T1D",build=NULL,show.codes=FALSE) {
  disease.codes <- c("Type 1 Diabetes", "Crohns Disease","Rheumatoid Arthritis",
  "Systemic Scleroderma",  "Ulcerative Colitis","Inflammatory Bowel Disease",  "Multiple Sclerosis",
  "Bipolar Disorder",  "Diabetes Mellitus",  "Coronary Artery Disease",  "Hypertension",  "Celiac Disease",
  "Systemic Lupus Erythematosus",  "Ankylosing Spondylitis",  "Type 2 Diabetes",  "Sjogren Syndrome",
  "Graves' Disease",  "Juvenile Rheumatoid Arthritis",  "Vitiligo",  "Primary Biliary Cirrhosis",
  "Psoriasis",  "Idiopathic Membranous Nephropathy",  "Immunoglobulin A Deficiency",
  "Autoimmune Thyroid Disease",  "Juvenile Idiopathic Arthritis",  "Narcolepsy",  "Alopecia Areata",
  "Alzheimer's Disease")
  abbr <- c("T1D","CD","RA","SCL","UC","IBD","MS","BD","DM","CAD","HYP",
            "CEL","SLE","AS","T2D","SS","GD","JRA","VIT",
            "PBC","PSO","IMN","IGA","ATD","JIA","NAR","AA","AD")
  code.table <- cbind(Abbreviation=abbr,FullNames=disease.codes)
  if(show.codes) { cat("values for the 'disease' parameter can be specified by the following index numbers or abbreviations:\n")
                   print(code.table,quote=F) ; return() }
  disease <- disease[1]
  if(toupper(disease) %in% abbr) { 
    disN <- match(toupper(disease),toupper(abbr)) 
  } else {
    if(toupper(disease) %in% toupper(disease.codes)) {
      disN <- match(toupper(disease),toupper(disease.codes))
    } else {
      if(as.numeric(disease) %in% 1:length(disease.codes)) { 
        disN <- as.numeric(disease)
      } else {
        stop("Invalid input for 'disease', use show.codes=TRUE to see list of codes/abbreviations")
      }
    }
  }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(!build %in% c("hg18","hg19")) { stop("only hg18 and hg19 are supported for this function") }
  filenm <- cat.path(dir=getwd(),pref=tolower(abbr[disN]),"hits",suf=build,ext="tab")
  cat("attempting to download",abbr[disN],"hits from t1dbase\n")
  url36 <- paste("http://www.t1dbase.org/webservice/RegionDownloads/=/model/variantsTAB/species=Human&disease_id=",disN,"&type=assoc&build=GRCh36",sep="")
  url37 <- paste("http://www.t1dbase.org/webservice/RegionDownloads/=/model/variantsTAB/species=Human&disease_id=",disN,"&type=assoc&build=GRCh37",sep="")
  urL <- switch(build, hg18=url36,  hg19=url37)
  success <- T
  success <- tryCatch(download.file(urL ,filenm ,quiet=T),error=function(e) { F } )
  #prv(filenm,urL,success)
  if(!is.logical(success)) { success <- T }
  if(success) {
    t1dh <- readLines(filenm)
    firsts <- substr(t1dh,1,2)
    t1dh <- t1dh[firsts!="##"]
    len.lst <- strsplit(t1dh,"\t")
    rsids <- sapply(len.lst,"[",3)
    if(substr(rsids[1],1,2)!="rs") { rsids <- rsids[-1] }
    if(length(rsids)<1) { 
      cat("download successful but the list of hits for",disease.codes[disN],"was empty\n") 
    } else {
      cat("download successful for",disease.codes[disN],"\n")
    }
    return(unique(rsids))
  } else {
    stop("couldn't reach t1dbase website at: ",urL)
  }
}





#' Retreive GO terms from biomart for a given gene list
#' 
#' Gene-ontology terms (GO-terms) are commonly used for testing for simple functional
#' enrichment for pathways, etc. This function can retrieve biological function, 
#' cellular component, or molecular description, depending on the parameters chosen.
#' @param gene.list a list of gene, use HGNC names, like COMT, HLA-C, CTLA4, etc.
#' @param bio logical, whether to return biological process GO terms
#' @param cel logical, whether to return cellular component GO terms
#' @param mol logical, whether to return molecular function GO terms
#' @return data.frame containing the gene name in the first column, chromosome in the
#' second column, and the GO terms in the third column, where one gene has multiple
#' GO terms, this will produce multiple rows, so there will usually be more rows
#' than genes entered. The data.frame can have 3,4 or 5 columns depending on
#' how many GO terms are selected.
#' @export
#' @examples
#' get.GO.for.genes(c("CTLA4","PTPN2","PTPN22")) # biological terms (default)
#' get.GO.for.genes(c("CTLA4","PTPN2","PTPN22"),cel=TRUE) # add cellular GO terms
get.GO.for.genes <- function(gene.list,bio=T,cel=F,mol=F) {
  must.use.package(c("biomaRt","genoset","gage"),T)
  ens <- useMart("ENSEMBL_MART_ENSEMBL",
                 dataset="hsapiens_gene_ensembl",
                 host="may2009.archive.ensembl.org",
                 path="/biomart/martservice",
                 archive=FALSE)
  ens <- useDataset("hsapiens_gene_ensembl",mart=ens)
  data(egSymb)
  base.attr <- c("hgnc_symbol", "chromosome_name")
  if(bio) { base.attr <- c(base.attr,"go_biological_process_description") }
  if(cel) { base.attr <- c(base.attr,"go_cellular_component_description") }
  if(mol) { base.attr <- c(base.attr,"go_molecular_function_description") }
#  dat <- getBM(attributes = c("hgnc_symbol", "chromosome_name",
#                              "start_position", "end_position", "band"), filters = "hgnc_symbol",
#               values = egSymb[,2], mart = ens)
  results <- getBM(attributes = base.attr, filters = "hgnc_symbol",
                   values = c(gene.list), mart = ens)
  return(results)
}




#' Convert ensembl ids to HGNC gene ids 
#' 
#' Retrieve the gene IDs (HGNC) corresponding to a list of ensembl gene ids.
#' Note that this will not find all IDs found on ensembl.org, as it uses bioMart which
#' seems to be incomplete, but this only pertains to a small minority of genes, so this
#' function should have general utility for most applications. This is of course the case
#' at the time of writing - bioMart is likely to be updated at some point.
#' @param ens character, a list of ensembl gene ids, of the form ENSG00xxxxxxxxx
#' @param ... further arguments to get.gene.annot()
#' @param dir character, 'dir' is the location to download gene and cytoband information; if
#' left as NULL, depending on the value of getOption("save.annot.in.current"), the annotation
#' will either be saved in the working directory to speed-up subsequent lookups, or deleted 
#' after use.
#' @param build character, "hg18" or "hg19" (or 36/37) to show which reference to retrieve. The 
#' default when build is NULL is to use the build from the current ChipInfo annotation
#' @param name.dups logical, if TRUE then duplicates will have a suffix appended to force the
#' list to be unique (e.g, so it would be usable as rownames, or in a lookup table). Otherwise
#' duplicate entries will just appear in the list multiple times
#' @param name.missing logical, if TRUE then missing values will be named as MISSING_n (n=1
#'  to # of missing), ensuring a valid unique name if the results are to be used as rownames,
#' etc. If FALSE then these will be left as NA. 
#' @return Returns a vector of HGNC gene ids corresponding to the 'ens' ensembl ids entered,
#' any ids not found will be returned as MISSING_n (n=1 to # of missing), if name.missing=TRUE.
#' If name.missing is FALSE then missing will be set to NA. Similarly with 'name.dups', if
#' duplicates are found and name.dups is true, each will be appended with suffix _n; else
#' their names will be left as is.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso GENE.to.ENS, rs.to.id, id.to.rs
#' @examples
#' ENS.ids <- c("ENSG00000183214", "ENSG00000163599", "ENSG00000175354", "ENSG00000134460")
#' ENS.to.GENE(ENS.ids)
#' gene.ids <- c("HLA-B","IFIH1","fake_gene!","FUT2")
#' ENS.to.GENE(GENE.to.ENS(gene.ids)) # lookup fails for the fake id, gives warning
ENS.to.GENE <- function(ens,dir=NULL,build=NULL,name.dups=FALSE,name.missing=TRUE,...) {
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  must.use.package(c("biomaRt","genoset","gage"),T)
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.gene.annot(...,dir=dir,bioC=FALSE,ens.id=TRUE,GRanges=FALSE)
  ### now have the gene data with the ensembl ids ##
  #print(head(ga))
  indx <- match(ens,ga$ens.id)
  missin <- length(which(is.na(indx))); valid <- length(indx)-missin
  if(valid<1) { warning("did not find any ENSEMBL ids from 'ens' in the bioMart human gene reference"); return(NULL) }
  if(missin>0) { warning(out.of(missin,(valid+missin))," of 'ens' did not match any ENSEMBL ids in the bioMart human gene reference") }
  outData <- ga$gene[indx]
  if(name.missing & any(is.na(outData))) {
    outData[is.na(outData)] <- paste("MISSING",pad.left(1:length(which(is.na(outData))),"0"),sep="_")
  }
  if(any(duplicated(outData))) { 
    if(name.dups) { 
      cnt <- 2
      while(any(duplicated(outData))) { 
        if(cnt==2) {
          outData[duplicated(outData)] <- paste(outData[duplicated(outData)],cnt,sep="_")
        } else {
          outData[duplicated(outData)] <- gsub(paste("_",cnt-1,sep=""),paste("_",cnt,sep=""),outData[duplicated(outData)])
        }
        cnt <- cnt + 1
      }
    } else { 
      warning("duplicated gene names produced, select 'name.dups=TRUE' to append numbers to make these unique")
    }
  }
  return(outData)
}


#' Convert gene ids to ensembl ids
#' 
#' Retrieve the ensembl IDs corresponding to a list of common gene names (HGNC format).
#' @param genes character, gene labels, e.g, "APOE"
#' @param ... further arguments to get.gene.annot()
#' @param dir character, 'dir' is the location to download gene and cytoband information; if
#' left as NULL, depending on the value of getOption("save.annot.in.current"), the annotation
#' will either be saved in the working directory to speed-up subsequent lookups, or deleted 
#' after use.
#' @return Returns a vector of HGNC gene ids corresponding to the 'ens' ensembl ids entered,
#' any ids not found will be returned as NA.
#' @export
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
#' @seealso GENE.to.ENS, rs.to.id, id.to.rs
#' @examples
#' gene.ids <- c("MICA","PTPN2","IL2RA","APOE")
#' GENE.to.ENS(gene.ids)
GENE.to.ENS <- function(genes,dir=NULL,...) {
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  ga <- get.gene.annot(...,dir=dir,bioC=FALSE,ens.id=TRUE,GRanges=FALSE)
  ### now have the gene data with the ensembl ids ##
  indx <- match(genes,ga$gene)
  missin <- length(which(is.na(indx))); valid <- length(indx)-missin
  if(valid<1) { warning("did not find any gene ids from 'genes' in the bioMart human gene reference"); return(NULL) }
  if(missin>0) { warning("at least one of 'genes' did not match any gene ids in the bioMart human gene reference") }
  outData <- ga$ens.id[indx]
  return(outData)
}



#' Retrieve the 'n' closest GENE labels or positions near specified locus
#' 
#' @param chr integer, chromosome, should be a number from 1 to 25, where 23,24,25 are X,Y,MT
#' @param pos integer, genomic position, should be between 1 and the length of the chromosome 'chr'
#' @param n integer, the number of nearest GENEs to seek, if there aren't enough in the annotation
#' then NAs will fill the gaps to force the return value length to equal 'n'
#' @param side character, can be 'either', 'left' or 'right' and specifies which side of the 'pos'
#' to look for nearest genes (where left is decreasing genomic position and right is increasing)
#' @param ids logical, if TRUE will return GENE labels, 
#' or if FALSE will return the chromosome positions of the genes
#' @param limit integer, a limit on the maximum distance from the position 'pos' can be specified
#' @param build integer whether to use build 36/37 parameters, 36/37 is preferred, but can enter
#' using any form recognised by ucsc.sanitizer()
#' @param dir string, specify the directory to store or look for annotation (defaults to current)
#' @param ga RangedData object, e.g, result of get.gene.annot(); gene annotation to save download
#' time if repeatedly calling this function
#' @export
#' @seealso expand.nsnp, nearest.snp, get.gene.annot
#' @return Set of GENE ids (when ids=TRUE), or otherwise genomic positions within chromosome 'chr'.
#' If the number of gemes on the chromosome or the bounds of the 'side' and 'limit' parameters
#' restrict the number returned to less than 'n' then the return value will be padded with NAs.
#' @examples
#' nearest.gene(1,159000000,n=10) # return ids
#' nearest.gene(1,159000000,n=10,build=37)
#' nearest.gene(1,159000000,n=10,build=36,ids=F) # return positions
#' nearest.gene(1,159000000,n=10,build=37,ids=F)
#' nearest.gene(6,25000000,n=10,build=37,ids=F,side="left")  # only genes to the left of the locus
#' nearest.gene(6,25000000,n=10,build=37,ids=F,side="right") # only genes to the right of the locus
nearest.gene <- function(chr, pos, n=1, side=c("either","left","right"),ids=TRUE,limit=NULL,build=NULL, ga=NULL) { 
  # ids - whether to return ichip SNP ids or positions
  if(length(chr)>1) { warning("chr should be length 1, using only first entry"); chr <- chr[1] }
  if(length(pos)>1) { warning("pos should be length 1, using only first entry"); pos <- pos[1] }
  if(is.null(build)) { build <- getOption("ucsc") }
  chrom <- paste(chr)
  build <- ucsc.sanitizer(build)
  if(is(get.gene.annot)[1]=="RangedData") { 
    if(!"gene" %in% colnames(ga)) { ga <- NULL }
  } else { ga <- NULL }
  if(is.null(ga)) {  
    ga <- get.gene.annot(build=build,GRanges=F) 
    if(!exists("ga")) { stop("couldn't find gene annotation") }  ## load object: ga [gene database]
    ga <- ga[ga$gene!="",]
  }
  side <- tolower(side[1]); 
  if(!side %in% c("either","left","right")) {
    side <- "either"; warning("invalid side argument, defaulting to 'either'") }
  if(!is.null(limit)) { if(!is.numeric(limit)) { limit <- NULL; warning("invalid limit argument, defaulting to NULL") } }
  all.chr <- paste(chr2(ga))
  all.st <- start(ga)[all.chr %in% chrom]
  all.en <- end(ga)[all.chr %in% chrom]
  #prv(all.st,all.en)
  if(length(all.st)<1) { warning("no positions found for 'chr' specified"); return(NULL) }
  difzS <- pos-all.st
  difzE <- pos-all.en
  all.true <- difzS==difzS
  if(is.null(limit)) { lfilt <- all.true } else { lfilt <- abs(difzS)<=limit | abs(difzE)<=limit }
  within <- difzS>0 & difzE<0
  if(side=="left") { filt <- ( difzE>0 & lfilt ) | within }
  if(side=="right") { filt <- ( difzS<0 & lfilt ) | within }
  if(side=="either") { filt <- ( all.true & lfilt ) | within }
  #print(length(which(filt)))
  tab <- rbind(abs(difzS),abs(difzE))
  minz <- apply(tab,2,min,na.rm=T)
  Difz <- abs(minz[filt])
  if(length(Difz)<n)  { warning("fewer than ",n," genes found for 'chr' specified (within 'limit'), NAs returned") }
  indx <- order(Difz)[1:n]
 # prv(minz,Difz,filt,indx)
  subi <- ga[["gene"]][all.chr %in% chrom][filt]
  #prv(ga,subi,all.chr,chrom)
  if(ids) {
    out <- ga[["gene"]][all.chr %in% chrom][filt][indx]
  } else {
    out <- start(ga)[(all.chr %in% chrom)][filt][indx]
  }
  return(out)
}




#' Plot genes to annotate figures with genomic axes
#' 
#' Quite often it is helpful to visualize genomic locations in the context of Genes
#' in the same region. This function makes it simple to overlay genes on plots
#' where the x-axis is chromosomal location.
#' @param chr chromosome number/name that the plot-range lies on
#' @param pos position range on the chromosome 'chr' that the plot lies on, in base-pairs
#' @param scl character, the scale that the x axis uses, ie, "b","kb","mb", or "gb", meaning
#' base-pairs, kilobases, megabases or gigabase-pairs.
#' @param y.ofs numeric, y-axis-offset, depending on what units are on your y-axis,
#' you may need to add an offset so that the gene annotation is drawn at an appropriate
#' level on the vertical axis, this value should be the centre of annotation
#' @param width depending on the range of your y-axis, you might want to expand or reduce
#' the vertical width of the gene annotation (in normal graph units)
#' @param txt logical, TRUE to include the names of genes on top of their representation
#' on the plot, or if FALSE, genes are drawn without labels.
#' @param chr.pos.offset if for some reason zero on the x-axis is not equal to 'zero' on
#' the chromsome, then this offset can correct the offset. For instance if you were using
#' a graph of the whole genome and you were plotting genes on chromosome 10, you would
#' set this offset to the combined lengths of chromosomes 1-9 to get the start point
#' in the correct place.
#' @param gs GRanges or RangedData object, this is annotation for the location of genes.
#' This will be retrieved using get.gene.annot() if 'gs' is NULL. THere may be several reasons
#' for passing an object directly to 'gs'; firstly speed, if making many calls then you won't
#' need to load the annotation every time; secondly, if you want to use an alternative annotation
#' you can create your own so long as it is a GRanges/RangedData object and contains a column
#' called 'gene' (which doesn't strictly have to contain gene labels, it could be any feature
#' you require, eg., transcript names, etc).
#' @param dir character, location to store file with the gene annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param box.col genes are drawn as boxes, this sets the colour of the boxes
#' @param txt.col this sets the colour of the label text (Gene names)
#' @param join.col for exons, or multipart genes, joins are made between the sections with
#' a central line, this sets the colour of that line.
#' @param ... further arguments to 'rect', the graphics function used to plot the 'genes'.
#' @export
#' @return Returns a data.frame, GRanges or RangedData object, depending on input parameters. Contained
#' will be HGNC gene labels, chromosome and start and end positions, other information depends on 
#' specific parameters documented above
#' @export
#' @examples
#' # EXAMPLE PLOT OF SOME SIMULATED SNPS on chr21-p11.1 #
#' loc <- c(9.9,10.2)
#' Band(chr=21,pos=loc*10^6)
#' rr <- in.window(rranges(50000),chr=21,pos=loc,unit="mb") # make some random MHC ranges
#' # create some SNPs and plot
#' rr3 <- rr; end(rr3) <- start(rr3) 
#' rownames(rr3) <- paste0("rs",sample(10^6,nrow(rr3)))
#' plotRanges(rr3,col="blue",scl="mb",xlim=loc,xlab="Chr21 position (Mb)",ylab="")
#' # NOW add UCSC hg18 GENE annotation to the plot #
#' ## not run ## plotGeneAnnot(chr=21,pos=c(9.95,10.1),scl="mb",y.ofs=1,build=36)
plotGeneAnnot <- function(chr=1, pos=NA, scl=c("b","kb","mb","gb"), y.ofs=0, width=1, txt=T, chr.pos.offset=0,
                            gs=NULL, build=NULL, dir=NULL, box.col="green", txt.col="black", join.col="red", ...)
{
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  dir <- validate.dir.for(dir,"ano")
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(!is(gs)[1] %in% c("RangedData","GRanges")) { gs <- get.gene.annot(dir=dir,build=build,GRanges=FALSE) }
  if(is(gs)[1]=="GRanges") { gs <- as(gs,"RangedData") }
  if(!"gene" %in% colnames(gs)) { warning("didn't find 'gene' column in annotation") ; return(NULL) }
  Col <- c("green", "darkgreen")
  if(all(is.na(pos))) { pos <- c(1,Inf) } 
  # get set of genes in range of the graph section + remove duplicate genes/exons
  rng.genez <- in.window(gs,chr,pos,full.overlap=F, rmv.dup=T,unit=scl)
  if(nrow(rng.genez)<1) { warning("no genes found in range") ; return(NULL) }
  old.cc <- 1 ; tp <- 2
  # set vertical alignments for annotation
  old.auto <- F
  y.cent <- y.ofs
  y.bot <- y.ofs-(width/2)
  y.top <- y.ofs+(width/2)
  # text alignment
  tps <- y.bot + c(.18,.29,.46,.64,.82)[c(1,3,5)]*width
  # x position with scaling (e.g, Mb units = 10^6)
  #unit <- tolower(unit[1]) ; mult <- switch(unit,b=0,kb=3,mb=6,gb=9); pos <- pos*10^mult
  x.scl <- make.divisor(scl)
  cnrlo <- start(rng.genez)/x.scl+chr.pos.offset
  cnrhi <- end(rng.genez)/x.scl+chr.pos.offset
  gnnm <- (rng.genez$gene)
  n.genes <- length(cnrlo)
  txt.cex <- .75; if(n.genes>10) { txt.cex <- .5 } ; if(n.genes>100) { txt.cex <- .35 } # more genes = smaller labels
  #print(n.genes)
  for (cc in 1:n.genes) {
    # draw rectangle and label for each gene
    #prv(cnrlo[cc],y.top,cnrhi[cc], y.bot)
    rect(cnrlo[cc],y.top,cnrhi[cc], y.bot,border=box.col,...)
  }
  for (cc in 1:n.genes) {
    if (gnnm[old.cc]==gnnm[cc] & cc!=1)
    {
      link <- c(cnrhi[old.cc],cnrlo[cc])
      if(link[1]<link[2]) {
        lines(link,y=rep(y.cent,2),lwd=1,col=join.col,lty="dotted")
      }
    } else {
      if(txt) {
        if(cnrlo[cc] < min(pos/x.scl)) { 
          txt.x <- mean(c(min(pos/x.scl),min(cnrhi[cc],max(pos/x.scl))),na.rm=T)  
        } else { txt.x <- cnrlo[cc] }
        text(txt.x,tps[tp],gnnm[cc],col=txt.col,cex=txt.cex,pos=4,las=2,offset=0)
      }
    }
    old.cc <- cc  
    tp <- tp+1; if(tp==4) {tp <- 1}; #if(n.genes <=5) { tp <- 2 }
  }
  return(rng.genez)
}



#' Retrieve locations of Immunoglobin regions across the genome
#' 
#' Returns the locations of immunoglobin regions in the human genome, for a given build, as
#' a list by chromosome, text vector, or GRanges/RangedData object.
#' For instance, for CNV research, these regions are known to be highly structurally complex
#' and can lead to false positive CNV-calls, so are often excluded.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param text logical, whether to return locations as a text vector of the form: chrN:xxxx-xxxx
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @export
#' @return Returns a list, GRanges or RangedData object, depending on input parameters. Contained
#' will be immunoglobin chromosome, start and end positions.
#' @example
#' get.immunog.locs()
#' get.immunog.locs(bioC=FALSE)
#' get.immunog.locs(text=TRUE,build=37)
get.immunog.locs <- function(build=NULL,bioC=TRUE,text=FALSE,GRanges=TRUE) {
  nchr <- 22
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(build[1]=="hg19") {
    #hg19
    chr <- c(22,14,2,14)
    stz <- c(22385572,105994256,89156874,22090057)
    enz <- c(23265082,107281230,89630187,23021097)
  } else {
    # hg18
    chr <- c(22,14,2,14)
    stz <- c(20715572,105065301,88937989,21159897)
    enz <- c(21595082,106352275,89411302,22090937)
  }
  nmz <- c("ig_c22","ig_c14_a","ig_c2","ig_c14_b")
  reg.dat <- rep("immunoglobin",length(chr))
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=chr,
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=T)
    if(text) { outData <- ranges.to.txt(outData) } else {
      if(GRanges) { outData <- as(outData,"GRanges") }
    }
  } else {
    outData <- vector("list",nchr); names(outData) <- paste("chr",1:nchr,sep="")
    for (cc in 1:nchr) {
      if(cc %in% chr) {
        outData[[cc]] <- list(start=stz[chr==cc],end=enz[chr==cc])
      }
    }
  }
  return(outData)
}


#' Return Centromere locations across the genome
#' 
#' Returns the locations of centromeres in the human genome, for a given build, as
#' a list by chromosome, text vector, or GRanges/RangedData object.
#' @param dir character, location to store file with the this annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame.
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @param text logical, whether to return locations as a text vector of the form: chrN:xxxx-xxxx
#' @param autosomes logical, if TRUE, only return results for autosomes, if FALSE, also include
#' X and Y.
#' @export
#' @return Returns a list, GRanges or RangedData object, depending on input parameters. Contained
#' will be centromere chromosome and start and end positions.
#' @example
#' get.centromere.locs()
#' get.centromere.locs(bioC=FALSE,autosomes=TRUE)
#' get.centromere.locs(text=TRUE)
get.centromere.locs <- function(dir=NULL,build=NULL,
                                bioC=TRUE,GRanges=TRUE,text=FALSE,autosomes=FALSE)
{
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  dir <- validate.dir.for(dir,c("ano"),warn=FALSE); success <- TRUE
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  local.file <- cat.path(dir$ano,"cyto")
  tt <- get.cyto(build=build,bioC=FALSE,dir=dir,GRanges=FALSE)
  chrn <- paste(1:22)
  if(!autosomes) { 
    chrn <- c(chrn,c("X","Y"))
  }
  nchr <- length(chrn)
  my.chr.range <- vector("list",nchr)
  names(my.chr.range) <- paste("chr",chrn,sep="")
  for (cc in 1:nchr) {
    just.centros <- tt[paste(tt[,5])=="acen",]
    just.chr <- just.centros[which(paste(just.centros[,1])==names(my.chr.range)[cc]),]
    my.chr.range[[cc]] <- list(start=min(just.chr[,2]), end=max(just.chr[,3]))
  }
  reg.dat <- rep("centromere",nchr)
  nmz <- paste(reg.dat,chrn,sep="_")
  stz <- sapply(my.chr.range,"[[",1)
  enz <- sapply(my.chr.range,"[[",2)
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=TRUE)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=gsub("chr","",chrn),
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=TRUE)
    if(text) { 
      outData <- ranges.to.txt(outData) 
    } else {
      if(GRanges){
        outData <- as(outData,"GRanges")
      }
    }
  } else {
    outData <- my.chr.range 
  }
  return(outData)
}


#' Return Cytoband/Karyotype locations across the genome
#' 
#' Returns the locations of cytobands/karyotype-bands in the human genome, for a given build, as
#' a data.frame, or GRanges/RangedData object.
#' @param dir character, location to store file with the this annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @param refresh logical, whether to re-download the file if the existing file has become corrupted
#' @export
#' @return Returns a list, GRanges or RangedData object, depending on input parameters. Contained
#' will be centromere chromosome and start and end positions.
#' @example
#' get.centromere.locs()
#' get.centromere.locs(bioC=FALSE)
#' get.centromere.locs(text=TRUE)
get.cyto <- function(build=NULL,dir=NULL,bioC=TRUE,GRanges=TRUE,refresh=FALSE) {
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  local.file="cyto"
  if(is.null(dir)) {
    local.file <- cat.path("",local.file,suf=build,ext="tar.gz")
  } else { 
    local.file <- cat.path(dir,local.file,suf=build,ext="tar.gz")
  }
  if(!file.exists(local.file) | refresh) {
    golden.path <- paste("http://hgdownload.cse.ucsc.edu/goldenPath/",build,"/database/cytoBand.txt.gz",sep="")
    success <- tryCatch( download.file(url=golden.path,local.file,quiet=T),error=function(e) { F } )
    if(is.logical(success)) {
      if(!success) { warning("couldn't reach ucsc website! try sourcing cytoband data elsewhere"); return(NULL) } }
    tt <- reader(local.file,header=FALSE)
    if(is.null(dir)) { unlink(local.file) }
  } else {
    tt <- reader(local.file)
  }
  colnames(tt) <- c("chr","start","end","band","negpos")
  write.table(tt,file=local.file,col.names=T,row.names=F,sep="\t",quote=F)
  mychr <- gsub("chr","",tt$chr,fixed=T)
  fullbands <- paste(mychr,tt$band,sep="")
  if(bioC ) {
    st <- as.numeric(tt$start)
    en <- as.numeric(tt$end)
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=st,end=en,names=fullbands),space=mychr,
                          negpos=tt$negpos,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=T)
    if(GRanges) { outData <- as(outData,"GRanges") }
  } else {
    outData <- tt 
    if("band" %in% colnames(outData)) {
      ## make 'chr-band' rownames to be consistent with the RangedData object if bioC=T
      rownames(outData) <- fullbands
      #outData <- outData[,-which(colnames(outData) %in% "band")]
    }
  }
  return(outData)
}



#' Get HapMap recombination rates for hg18 (build 36)
#' 
#' Recombination rate files can be used to calculate recombination distances
#' for genome locations, in centimorgans. This function downloads these reference
#' files from the hapmap NCBI website. At the time of writing they were only 
#' availble for build 36. If using a more recent build I suggest using the
#' conversion function conv.37.36(), then recomWindow(), then conv.36.37() to 
#' get recombination distances for other builds. If getOption("save.annot.in.current")
#' is <=0 then no files will be kept. Otherwise an object containing this mapping data
#' will be saved in the local directory if dir=NULL, or else in the directory specified.
#' Allowing this reference to be saved will greatly increase the speed of this function
#' for subsequent lookups
#' @param dir character, location to store binary file with the recombination maps for
#' chromosomes 1-22. If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param verbose logical, if the binary file is not already downloaded, when verbose
#' is TRUE, there will be some output to the console indicating the progress of the
#' download. If FALSE, all output is suppressed.
#' @param refresh logical, if you already have the binary file in the current directory,
#' this argument will let you re-download and re-generate this file, e.g, if the file
#' is modified or corrupted this will make a new one without having to manually delete it
#' @param compress logical, this argument is passed to 'save' and will result in a larger
#' binary file size, but quicker loading times, so 'FALSE' is recommended for faster retrieval.
#' @export
#' @return Returns a list object of length 22, containing the recombination map files
#' as 22 separate data.frame's.
#' @example
#' ## not run as it takes roughly 2 minutes to download and read-in ##
#' ## uncomment the following 3 lines to run:
#' ## rec.map <- get.recombination.map(getwd())
#' ## file.on.disk <- "rrates_genetic_map_chr_1_22_b36.RData"
#' ## if(file.exists(file.on.disk)) { unlink(file.on.disk) } # remove the downloaded file
get.recombination.map <- function(dir=NULL,verbose=TRUE,refresh=FALSE, compress=FALSE) {
  n.chr <- 22
  hap.dir <- "http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/latest/rates/"
  temp.dir <- "recombinationratesGF13fDR1er119"
  local.file <- "rrates_genetic_map_chr_1_22_b36.RData"
  if(!file.exists(temp.dir)) { dir.create(temp.dir) } 
  local.files=paste0(temp.dir,"/genetic_map_chr",1:n.chr,"_b36.txt")
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(!is.null(dir)) {
    local.files <- cat.path(dir,local.files,ext="txt")
    local.file <- cat.path(dir,local.file,ext="RData")
  }
  if(!file.exists(local.file) | refresh) {
    if(verbose) { cat("Downloading recombination data from: ",hap.dir,"\n") }
    hapmap.urls <- cat.path(dir=hap.dir,fn=basename(local.files))
    success <- TRUE
    for (cc in 1:n.chr) {
      #print(hapmap.urls[cc])
      success <- tryCatch( download.file(url=hapmap.urls[cc],local.files[cc],quiet=T),error=function(e) { F } )
      if(verbose) { loop.tracker(cc,n.chr*2) }
    }
    if(is.logical(success)) {
      if(!success) { warning("couldn't download at least one of the files from: ",hap.dir); return(NULL) } }
    map.files.list <- vector("list",n.chr)
    for (cc in 1:n.chr) {
      map.files.list[[cc]] <- read.table(local.files[cc],header=TRUE)
      if(is.data.frame(map.files.list[[cc]])) { 
        unlink(local.files[cc]) 
      } else { warning("downloaded map file was corrupt for chr",cc) }
      if(verbose) { loop.tracker(n.chr+cc,n.chr*2) }
    }
    if(file.exists(temp.dir)) { file.remove(temp.dir) }   # delete the temporary directory
  } else {
    map.files.list <- reader(local.file)
  }
  if(!is.null(dir)) { save(map.files.list,file=local.file,compress=compress) }
  if(length(map.files.list)!=n.chr) { stop("Unfortunately the object derived seems corrupted") }
  names(map.files.list) <- paste0("chr",1:n.chr)
  return(map.files.list)
}


#' Get exon names and locations from UCSC
#' 
#' Various R packages assist in downloading exonic information but often the input required is 
#' complex, or several lines of code are required to initiate, returning an object that
#' might require some manipulation to be useful. This function simplifies the job 
#' considerably, not necessarily requiring any arguments. The object returned can be
#' a standard data.frame or a bioconductor GRanges/RangedData object. The raw annotation
#' file downloaded will be kept in the working directory so that subsequent calls to
#' this function run very quickly, and also allow use offline.
#' @param dir character, location to store file with the gene annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param transcripts logical, if TRUE, return transcripts rather than exons
#' @param GRanges logical, if TRUE and bioC is also TRUE, then returned object will be GRanges, otherwise
#' it will be RangedData
#' @export
#' @return Returns a data.frame, GRanges or RangedData object, depending on input parameters. Contained
#' will be HGNC gene labels, chromosome, start and end positions, transcript id number and name
#' @example
#' ## not run as it takes too long to download for CRAN ##
#' ## get.exon.annot()
#' ## get.exon.annot(bioC=FALSE,build=37)
get.exon.annot <- function(dir=NULL,build=NULL,bioC=T, transcripts=FALSE, GRanges=TRUE) {
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  ## load exon annotation (store locally if not there already)
  from.scr <- T
  txt <- if(transcripts) { "trans" } else { "exon" }
  if(!is.null(dir)) {
    dir <- validate.dir.for(dir,"ano")
    ex.fn <- cat.path(dir$ano,pref=txt,"Annot",suf=build,ext="RData")
    if(file.exists(ex.fn)) {
      tS <- reader(ex.fn)
      if(transcripts) {
        if(is(tS)[1]=="data.frame") { from.scr <- F }
      } else {
        if(is(tS)[1]=="GRanges") { from.scr <- F }
      }
    }
  }
  must.use.package("GenomicFeatures",T)
  if(from.scr) {
    must.use.package("gage",T)
    # get transcripts from build table 'knownGene'
    success <- tryCatch(txdb <- makeTranscriptDbFromUCSC(genome=build,
                                                         tablename="knownGene")  ,error=function(e) { F } )
    if(is.logical(success)) { 
      if(!success) {
        warning("Couldn't reach build website! try again later or, \n",
                "if in europe/uk, there may still be a bug in rtracklayer; \n",
                "Installing the latest version of R and bioconductor\n",
                "and running biocLite('rtracklayer'), should fix this")
        return(NULL) }
    }
    if(!transcripts) {
      ex = exonsBy(txdb, by="gene")
      tS <- toGenomeOrder(unlist(ex),strict=T)
    } else {  
      tS = transcriptsBy(txdb, by="gene")
      data(egSymb)
      select <- match(names(tS),egSymb[,1])
      names(tS)[!is.na(select)] <- egSymb[,2][select[!is.na(select)]]
      tS <- as.data.frame(tS)
    }
  }
  if(exists("ex.fn")) { save(tS,file=ex.fn) }
  if(transcripts) {
    chrs <- paste(tS$seqnames); chrs <- gsub("chr","",chrs)
    tS$seqnames <- chrs
    if(all(colnames(tS)==c("element","seqnames","start","end","width","strand","tx_id","tx_name"))) {
      colnames(tS) <- c("gene","chr","start","end","width","strand","txid","txname")
    } else {
      cat(" unexpected colnames found using makeTranscriptDbFrombuild()\n")
      if(bioC) { cat(" therefore returning data.frame instead of RangedData object\n") ;
                 bioC <- F }
    }
    if(bioC) {
      tS <- RangedData(ranges=IRanges(start=tS$start,end=tS$end),
                       space=tS$chr,gene=tS$gene, strand=tS$strand,
                       txid=tS$txid, txname=tS$txname,universe=build)
      tS <- toGenomeOrder2(tS,strict=T)
      if(GRanges) { tS <- as(tS,"GRanges") }
    }
    return(tS)
  } else {
    rownames(tS) <- paste(1:nrow(tS))
    ei <- mcols(tS)[["exon_id"]]
    ei2 <- add.trail(ei,suffix=strsplit("abcdefghijklmnopqrstuvwxyz","")[[1]])
    mcols(tS)[["exon_name"]] <- ei2
    if(bioC) {
      if(GRanges) {
        return(tS)
      } else {
        return(as(tS,"RangedData"))
      }
    } else {
      return(ranges.to.data.frame(tS))
    }
  }
}



#' Get human gene names and locations from biomart
#' 
#' Various R packages assist in downloading genomic information but often the input required is 
#' complex, or several lines of code are required to initiate, returning an object that
#' might require some manipulation to be useful. This function simplifies the job 
#' considerably, not necessarily requiring any arguments. The object returned can be
#' a standard data.frame or a bioconductor GRanges/RangedData object. The raw annotation
#' file downloaded will be kept in the working directory so that subsequent calls to
#' this function run very quickly, and also allow use offline.
#' @param dir character, location to store file with the gene annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param duplicate.report logical, whether to provide a report on the genes labels that are listed
#' in more than 1 row - this is because some genes span ranges with substantial gaps within them
#' @param one.to.one logical, as per above, some genes have duplicate entries, sometimes for simplicity
#' you want just one range per gene, if this parameter is set TRUE, one range per gene is enforced,
#' and only the widest range will be kept by default for each unique gene label
#' @param remap.extra logical, whether to remap chromosome annotation for alternative builds and
#' unconnected segments to the closest regular chromosome, e.g, mapping MHC mappings to chromosome 6
#' @param discard.extra logical, similar to above, but if TRUE, then any non-standard chromosome
#' genes will just be discarded
#' @param only.named logical, biomart annotation contains some gene segments without names, if TRUE, then
#' such will not be included in the returned object (note that this will happen also if one.to.one is TRUE)
#' @param ens.id logical, whether to include the ensembl id in the dataframe
#' @param refresh logical, if you already have the file in the current directory,
#' this argument will let you re-download and re-generate this file, e.g, if the file
#' is modified or corrupted this will make a new one without having to manually delete it
#' @param GRanges logical, if TRUE and bioC is also TRUE, then returned object will be GRanges, otherwise
#' it will be RangedData
#' @export
#' @return Returns a data.frame, GRanges or RangedData object, depending on input parameters. Contained
#' will be HGNC gene labels, chromosome and start and end positions, other information depends on 
#' specific parameters documented above
#' @example
#' ## not run as it takes too long to download for CRAN ##
#' ## get.gene.annot()
#' ## get.gene.annot(bioC=FALSE,build=37)
get.gene.annot <- function(dir=NULL,build=NULL,bioC=TRUE,duplicate.report=FALSE,
                           one.to.one=FALSE,remap.extra=FALSE,discard.extra=TRUE,only.named=FALSE,
                           ens.id=FALSE,refresh=FALSE,GRanges=TRUE) {
  # faster than exon, but only contains whole gene ranges, not transcripts
  # allows report on duplicates as some might be confused as to why some genes
  # have more than one row in the listing (split across ranges usually)
  # run with dir as NULL to refresh changes in COX
  must.use.package(c("biomaRt","genoset","gage"),T)
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  from.scr <- T
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(!is.null(dir)) {
    dir <- validate.dir.for(dir,"ano")
    utxt <- ""; if(one.to.one) { utxt <- "_unq" }
    if(ens.id) { utxt <- paste(utxt,"ens",sep="_") }
    gn.fn <- cat.path(dir$ano,"geneAnnot",pref=build,suf=utxt,ext="RData")
    if(file.exists(gn.fn) & !refresh) {
      dat <- get(paste(load(gn.fn)))
      from.scr <- F
    }
  }
  # colnames for output
  nm.list <- c("gene","chr","start","end","band")
  if(ens.id & bioC) { warning("ens.id=TRUE only has an effect when bioC=FALSE") }
  if(from.scr) {
    if(build=="hg18") {
      ens <- useMart("ENSEMBL_MART_ENSEMBL",
                     dataset="hsapiens_gene_ensembl",
                     host="may2009.archive.ensembl.org",
                     path="/biomart/martservice",
                     archive=FALSE)
    } else {
      ens <- useMart("ensembl")
    }
    ens <- useDataset("hsapiens_gene_ensembl",mart=ens)
  
    
    attr.list <- c("hgnc_symbol", "chromosome_name",
                   "start_position", "end_position", "band")
    if(ens.id) { attr.list <- c(attr.list,"ensembl_gene_id") }
    if(build=="hg18" & only.named & !ens.id) {
      data(egSymb)
      dat <- getBM(attributes = attr.list, filters = "hgnc_symbol",
                 values = egSymb[,2], mart = ens)
    } else {
      dat <- getBM(attributes = attr.list,  mart = ens)
    }
    if(exists("gn.fn")) { save(dat,file=gn.fn) }
  } 
  if(ens.id) { nm.list <- c(nm.list,"ens.id") }
  #return(dat)
  no.gene.names <- which(dat[[1]]=="")
  if(one.to.one | (only.named & !ens.id & length(no.gene.names)>0)) { dat <- dat[-no.gene.names,] }
  if(remap.extra) {
    dat$chromosome_name <- tidy.extra.chr(dat$chromosome_name)
    #dat$chromosome_name[grep("c6",dat$chromosome_name,ignore.case=T)] <- 6  # prevent issues with c6_COX, c6_QBL  
    #dat$chromosome_name[grep("c5",dat$chromosome_name,ignore.case=T)] <- 5  # prevent issues with c5_H2  
    #dat$chromosome_name[grep("NT",dat$chromosome_name,ignore.case=T)] <- "Z_NT"  # merge all NT regions to one label
  }
  if(discard.extra) {
    ## http://www.lrg-sequence.org/ ##
    # note that if remapping is already done, then these won't be discarded unless remapping failed
    tt <- tidy.extra.chr(dat$chromosome_name,select=TRUE)
    badz <- which(!tt)
    if(length(badz)>0) { dat <- dat[-badz,] } # remove LRG, GS, HG, NT, COX, etc annotation from set
  }
  if(bioC) {
    outData <- RangedData(ranges=IRanges(start=dat$start_position,end=dat$end_position),
                          space=dat$chromosome_name,gene=dat$hgnc_symbol, band=dat$band, universe=build)
    outData <- toGenomeOrder2(outData,strict=T)
    if(duplicate.report | one.to.one) {
      genez <- outData$gene
      dG <- which(duplicated(genez))
      if(!duplicate.report) {
        dup.genes <- genez[dG]
      } else {
        dup.genes <- gene.duplicate.report(outData)
      }
      stz <- start(outData); enz <- end(outData); wdz <- width(outData)
      to.del <- to.ch <- ch.st <- ch.en <- NULL
      n.dup <- length(dup.genes)
      if(one.to.one & n.dup>0) {
        #return(dup.genes)
        indz <- sapply(as.list(dup.genes),function(X) { which(genez %in% X) })
        # keep the range that is widest
        st.en <- lapply(indz,function(X) { c(stz[X][wdz[X]==max(wdz[X])][1],enz[X][wdz[X]==max(wdz[X])][1]) } )
        to.del <- dG
        to.ch <- sapply(indz,min,na.rm=T)
        ch.st <- sapply(st.en,"[",1) 
        ch.en <- sapply(st.en,"[",2) 
        start(outData)[to.ch] <- ch.st 
        end(outData)[to.ch] <- ch.en 
        outData <- outData[-to.del,]
        cat("kept widest ranges, merging",length(to.del)+length(to.ch),"duplicate gene labels to",length(to.ch),"labels\n")
      } 
      # but haven't done anything about them or removed them! 
    }
    if(GRanges) { outData <- as(outData, "GRanges") }
  } else {
    outData <- dat; colnames(outData) <- nm.list
  }
  return(outData)
}


#' Derive Telomere locations across the genome
#' 
#' Returns the locations of telomeres in the human genome, for a given build, as
#' a list by chromosome, text vector, or GRanges/RangedData object.
#' @param dir character, location to store file with the this annotation.
#' If NULL then getOption("save.annot.in.current")>=1 will result in
#' this file being stored in the current directory, or if <=0, then this file will not
#' be stored.
#' @param kb The number of base pairs at the start and end of a chromosome that are defined as
#' belonging to the telomere can be a little arbitrary. This argument allows specification
#' of whatever threshold is required.
#' @param build string, currently 'hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is build-36/hg-18. Will also accept integers 36,37 as alternative arguments.
#' @param bioC logical, whether to return the annotation as a ranged S4 object (GRanges or
#' RangedData), or as a data.frame
#' @param GRanges logical, whether to return a GRanges object, or FALSE to return RangedData
#' @param text logical, whether to return locations as a text vector of the form: chrN:xxxx-xxxx
#' @param autosomes logical, if TRUE, only return results for autosomes, if FALSE, also include
#' X and Y.
#' @param mito.zeros logical, Mitochondria have no telomeres (are circular) but for some purposes you
#' might want zero values in order to match with other annotation that includes all chromosomes and MT.
#' TRUE adds zeros for chrMT, and FALSE excludes chrMT.
#' @export
#' @return Returns a text vector, GRanges or RangedData object, depending on input parameters. Contained
#' will be telomere chromosome and start and end positions.
#' @example
#' get.telomere.locs()
#' get.telomere.locs(bioC=FALSE)
#' get.telomere.locs(text=TRUE)
get.telomere.locs <- function(dir=NULL,kb=10,build=NULL,bioC=TRUE,GRanges=TRUE,
                              text=FALSE,autosomes=FALSE,mito.zeros=FALSE)
{
  # the actual telomeres are typically about 10kb, but
  # for cnv-QC purposes want to exclude a larger region like 500kb
  # Mt have no telomeres, are circular, but for some purposes might want zero values in there
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  chr.lens <- get.chr.lens(dir=dir,build=build[1],autosomes=FALSE,mito=mito.zeros)
  n <- 1:22; if(!autosomes) { n <- c(n,"X","Y") } # Mt have no telomeres, are circular
  nchr <- length(n) #default
  if(mito.zeros) { n <- c(n,"M") }
  my.chr.range <- vector("list",length(n))
  names(my.chr.range) <- paste("chr",n,sep="")
  for (cc in 1:nchr) {
    one <- force.chr.pos(Pos=c(1,kb*1000),Chr=cc,build=build) # f..c..pos() makes sure is a valid range
    two <- force.chr.pos(Pos=chr.lens[cc]+c(-kb*1000,0),Chr=cc,build=build)
    my.chr.range[[cc]] <- list(start=c(one[1],two[1]),end=c(one[2],two[2]))
  }
  if(mito.zeros) {
    # add null values for the Mitochondrial chromosome
    cc <- cc+1; one <- c(1,1);  two <- chr.lens[cc]+c(0,0)
    my.chr.range[[cc]] <- list(start=c(one[1],two[1]),end=c(one[2],two[2]))
  }
  reg.dat <- rep("telomere",length(n)*2)
  chrz <- rep(n,each=2)
  nmz <- paste(reg.dat,chrz,rep(c("a","b"),times=length(n)),sep="_")
  stz <- as.vector(sapply(my.chr.range,"[[",1))
  enz <- as.vector(sapply(my.chr.range,"[[",2))
  if(bioC | text) {
    must.use.package(c("genoset","IRanges"),bioC=T)
    outData <- RangedData(ranges=IRanges(start=stz,end=enz,names=nmz),space=chrz,
                          reg=reg.dat,universe=build[1])
    outData <- toGenomeOrder2(outData,strict=T)
    if(text) { outData <- ranges.to.txt(outData) } else { if(GRanges) { outData <- as(outData,"GRanges") } }
  } else {
    outData <- my.chr.range 
  }
  return(outData)
}



#' Get chromosome lengths from build database
#' 
#' Quick and easy way to retrieve human chromosome lengths. Can select from hg18/hg19 (ie, 
#'  build 36/37), or any future builds (hg20, etc) stored in the same location on the build website.
#'  Default is to return lengths for 22 autosomes, but can also retrieve X,Y 
#'  and Mitochondrial DNA lengths by 'autosomes=FALSE' or n=1:25. Even if not connected to 
#'  the internet can retrieve hard coded lengths for hg18/hg19.
#'  @param dir directory to retrieve/download the annotation from/to (defaults to current getwd())
#'  if dir is NULL then will automatically delete the annotation text file from the local directory
#'   after downloading
#'  @param build string, currently 'hg17','hg18' or 'hg19' to specify which annotation version to use. 
#'  Default is getOption("ucsc"). Will also accept integers 17,18,19,35,36,37 as alternative arguments.
#'  @param autosomes logical, if TRUE, only load the lengths for the 22 autosomes, else load X,Y,[MT] as well
#'  @param len.fn optional file name to keep the lengths in
#'  @param mito logical, whether to include the length of the mitochondrial DNA (will not include unless autosomes is also FALSE)
#'  @param names logical, whether to name the chromosomes in the resulting vector
#'  @param delete.after logical, if TRUE then delete the text file that these lengths were downloaded to. 
#'  If FALSE, then the file will be kept, meaning future lookups will be faster, and available offline.
#'  @export
#'  @examples
#'  get.chr.lens(delete.after=TRUE) # delete.after simply deletes the downloaded txt file after reading
#'  get.chr.lens(build=35,autosomes=TRUE,delete.after=TRUE) # only for autosomes
#'  get.chr.lens(build="hg19",mito=TRUE,delete.after=TRUE) # include mitochondrial DNA length
get.chr.lens <- function(dir=NULL,build=NULL,autosomes=FALSE,len.fn="humanChrLens.txt",
                         mito=FALSE,names=FALSE, delete.after=FALSE, verbose=FALSE)
{
  # retrieve chromosome lengths from local annotation file, else download from build
  if(is.null(dir)) { if(any(getOption("save.annot.in.current")<1)) { dir <- NULL } else { dir <- getwd() } }
  if(is.null(dir)) { dir <- getwd() ; delete.after <- TRUE }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  dir <- validate.dir.for(dir,c("ano"),warn=F)
  chrlens.f <- cat.path(dir$ano,len.fn) # existing or future lengths file
  n <- 1:22; if(!autosomes) { n <- c(n,"X","Y","M") }
  hg18.backup <- c(247249719,242951149,199501827,191273063,180857866,170899992,158821424,
                   146274826,140273252,135374737,134452384,132349534,114142980,106368585,
                   100338915,88827254,78774742,76117153,63811651,62435964,46944323,
                   49691432,154913754,57772954,16571)
  hg19.backup <- c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,
                   146364022,141213431,135534747,135006516,133851895,115169878,107349540,
                   102531392,90354753,81195210,78077248,59128983,63025520,48129895,
                   51304566,155270560,59373566,16571)
  
  # backups for offline use
  if(build=="hg18") { offline.backup <- hg18.backup  } else { offline.backup <- hg19.backup }
  if(file.exists(chrlens.f))
  {
    # file seems to be in annotation directory already
    chrLens <- readLines(chrlens.f)
    if (length(chrLens)!=length(n))
    {
      #warning("Length of existing chromosome file didn't match expected:",length(n))
      notGot <- T
    } else {
      notGot <- F
      # we have the right length, but do we have the right version?
      if(build=="hg18" & length(which(chrLens %in% hg19.backup))>2) { notGot <- T }
      if(build=="hg19" & length(which(chrLens %in% hg18.backup))>2) { notGot <- T }
      names(chrLens) <- paste0("chr",n)
    }
  } else { notGot <- T }
  if (notGot | (!build %in% c("hg18","hg19"))) {
    #download from build
    if(verbose) { cat("attempting to download chromosome lengths from genome build ... ") }
    urL <- switch(build,
                  hg17="http://hgdownload.cse.ucsc.edu/goldenPath/hg17/database/chromInfo.txt.gz",
                  hg18="http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chromInfo.txt.gz",
                  hg19="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz",
                  hg20="http://hgdownload.cse.ucsc.edu/goldenPath/hg20/database/chromInfo.txt.gz")
    success <- T
    success <- tryCatch(download.file(urL, chrlens.f,quiet=T),error=function(e) { F } )
    if(!is.logical(success)) { success <- T }
    if(success) {
      if(verbose) {  cat("download successful\n") }
      chrL.f <- readLines(chrlens.f)
      len.lst <- strsplit(chrL.f,"\t")
      nmz <- sapply(len.lst,"[",1)
      lnz <- sapply(len.lst,"[",2)
      nnn <- paste0("chr",n)
      want.chr.names <- match(nnn,nmz)
      want.chr.names <- want.chr.names[!is.na(want.chr.names)]
      #print(want.chr.names)
      chrLens <- lnz[want.chr.names]
      names(chrLens) <- nnn
    } else {
      warning("couldn't reach build website, so have used offline versions of chr lengths")
      if(!build %in% c("hg18","hg19")) { warning("no offline version for build version:",build) }
      chrLens <- paste(offline.backup)[1:length(n)]; names(chrLens) <- paste0("chr",n)
      #print(n)
      delete.after <- F
    }
    if(length(dir)==1 & dir[1]=="" & delete.after) {
      unlink(chrlens.f)
    } else {
      writeLines(chrLens,con=chrlens.f) # save file for future use
    }
  }
  if(!mito & length(chrLens)>22) { chrLens <- chrLens[-grep("M",n)] }
  if(names) {
    return(chrLens)
  } else {
    return(as.numeric(chrLens))
  }
}



#' Obtain a listing of known T1D associated genomic regions
#'
#' This function uses a full list of ichip dense regions combined with a list of t1d
#' SNPs to get the t1d regions. For type 1 diabetes researchers.
#' @param dense.reg GRanges or RangedData object, only use if you need to provide for a
#' build other than 36 or 37 (hg18/hg19).
#' @param build e.g, 36/hg18 or 37/hg19, if left as NULL current getOption('ucsc') will
#' be used.
#' @param invert logical, set to TRUE if you wish to get the set of NON-T1D regions.
#' @return a GRanges object with the specified type 1 diabetes (or inverse) ranges
#' @export
#' @examples
#' # not run due to exceeding CRAN time limit
#' # t1d.reg <- get.t1d.regions()
#' # non.t1d <- get.t1d.regions(build=36,invert=TRUE)
get.t1d.regions <- function(dense.reg=NULL,build=NULL,invert=FALSE) {
  #source("~/github/iChip/iFunctions.R")
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  def.build <- getOption("ucsc")
  if(is.null(dense.reg)) { 
  	dense.reg <- reader("~/github/iChip/iChipFineMappingRegionsB36.RData") 
  	if(build!="hg18") {
     	if(build=="hg19") { 
     		dense.reg <- conv.36.37(dense.reg) 
        } else { 
        	stop("automatic loading for dense regions is only supported for build 36 or 37") 
        }
    }
  }  
  if(is(dense.reg)[1]=="GRanges") { dense.reg <- as(dense.reg,"RangedData") }
  if(!is(dense.reg)[1]=="RangedData") { stop("dense.reg must be RangedData or GRanges") }
  ichip.regions <- dense.reg
  rs.ids <- get.t1dbase.snps()
  locs <- Pos(rs.ids); chrs <- Chr(rs.ids)
  good <- !is.na(locs) & !is.na(chrs)
  locs <- locs[good]; chrs <- chrs[good]
  if(build!=def.build) {
  	if(!all(c(build,def.build) %in% c("hg18","hg19"))) { stop("if build is not equal to getOption('ucsc') then build and this getOption('ucsc') parameter must both be equivalent to either hg18 or hg19") }
  	if(build=="hg18") {
  		#prv(chrs,locs)
		locs <- conv.37.36(chr=chrs,pos=locs)[,"start"]
  	} else {
  		locs <- conv.36.37(chr=chrs,pos=locs)[,"start"]
  	}
  }
  t1dgr <- makeGRanges(chr=chrs,pos=locs,build=build)
  #prv(t1dgr,ichip.regions)
  t1d.regions <- suppressWarnings(subsetByOverlaps(set.chr.to.numeric(as(ichip.regions,"GRanges")),t1dgr))
#  t1dgr <- as(makeGRanges(chr=chrs,pos=locs,build=build),"RangedData")
#  t1d.regions <- find.overlaps(ichip.regions,ref=t1dgr,thresh=0.000000000001,ranges.out=TRUE)
  #prv(t1d.regions)
  if(invert) { t1d.regions <- invGRanges(t1d.regions,build=build) }
  return(t1d.regions)
}


# return subset of a ranged object that overlaps ichip dense mapped regions
get.t1d.subset <- function(X,t1d.only=TRUE,build=36,ichip.regions=NULL,T1D.regions=NULL,invert=FALSE) {
  source("~/github/iChip/iFunctions.R")
  if(is.null(ichip.regions)) {
    ichip.regions <- reader("~/github/iChip/iChipFineMappingRegionsB36.RData")
  }
  if(t1d.only) {
    if(is.null(T1D.regions)) {
      T1D.regions <- get.t1d.regions(ichip.regions,build=build,invert=invert)
    }
    T1D.regions <- as(T1D.regions,"RangedData")
    filt.sd <- find.overlaps(X,ref=T1D.regions,thresh=0.000000000000001,ranges.out=TRUE)
  } else {
    filt.sd <- find.overlaps(X,ref=ichip.regions,thresh=0.000000000000001,ranges.out=TRUE)
  }
  return(filt.sd)
}

#return subset of a ranged object that overlaps genes
get.genic.subset <- function(X,DB="gene",...) {
  filt.sd <- find.overlaps(X,db=DB,thresh=0.00000000000001,ranges.out=TRUE,...)
  return(filt.sd)
}


################## end annotation ##########################

