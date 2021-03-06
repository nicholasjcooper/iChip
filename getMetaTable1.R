## This file is used to generate table 1 by reading the UVA results, iChip regions,
# and finding the top meta results in each region, then tabulating


bonf <- .05/135836   #3.68   3.23*(10^-7)
source("~/github/iChip/iFunctions.R")
work.dir <- "/chiswick/data/ncooper/iChipData"
setwd(work.dir)
library(reader)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")
# get hard coded lists of SNPs passing, failing QC, appearing in table 1, etc #tab1snps
source('~/github/iChip/hardCodedSnpLists.R', echo=FALSE)

print(load("all.support.RData"))
new.table <- TRUE
table.only <- F  #F #TRUE
new.table.fn <- "nick.meta.table4.RData"
qc.excluded.snps <- reader("snpsExcluded.txt")
qc.excluded.snps <- qc.excluded.snps[!rs.to.ic(qc.excluded.snps) %in% rs.to.ic(table1snpsfinaljan30)]

# once this has been run once, can speed up by setting this to FALSE
if(!exists("TR") | (table.only) ) {
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  if(new.table) {
    table1 <- reader(new.table.fn)
  } else {
    table1 <- reader(docs[1])
  }
  ## FIX UVA XCHR name screw up ##
  table1 <- rmv.uva.dup.rows(table1)
  ##
  rownames(table1) <- clean.snp.ids(rownames(table1))
  table1 <- rmv.uva.dup.rows(table1)
  table1 <- table1[which(!is.na(ic.to.rs(rownames(table1)))),]
  excl <- reader(docs[2])
  prv(table1)
  prv(excl)
  
  t1d.prior.snps <- get.t1d.snps(ucsc="hg19")
  t1d.prior.snps <- ic.to.rs(narm(rs.to.ic(t1d.prior.snps)))
  
  excl.index <- which((rownames(table1) %in% clean.snp.ids(rownames(excl))) | is.na(table1[,"Position"]))
  if(length(excl.index)>0) { table1a <- table1[-excl.index,] } else { table1a <- table1 }
  table1a <- table1a[order(table1a[,"Position"]),]
  table1a <- table1a[order(table1a[,"CHR"]),]
  #print(table1a)
  #table1a <- table1a[order(table1a[,11]),]
  table1a <- rmv.uva.dup.rows(table1a)
  table1a <- table1a[!rownames(table1a) %in% c("rs1790121","rs12920271"),]  # shouldn't be in t1a
  
  if(table.only) { stop() }
  p.cols <- c(2,9,11,16,18,20)  #c(3,10,12,15,18)-1
  prv.large(table1a[,p.cols],rows=20,cols=7)
  poz37 <- as.numeric(table1a[,"Position"])
  ## poz36 <- lookup in all.support# iChip regions
  if(F) {
    print(load("/chiswick/data/ncooper/iChipData/dense.ic.regions.b36.RData"))
    ## gives  regions.gr
    cyto <- get.cyto(); cyto[["gene"]] <- rownames(cyto)
    cy.ranges <- toGenomeOrder(as(regions.gr,"RangedData"))
    cy.rangesX <- annot.cnv(cy.ranges,gs=cyto)  ### YES!!!! ###
    ichip.regions <- cy.rangesX
    save(ichip.regions,file="/chiswick/data/ncooper/iChipData/ichip.regions.RData")
  } else {
    print(load("/chiswick/data/ncooper/iChipData/ichip.regions.RData"))
    library(GenomicRanges)
    library(rtracklayer)
    #print(load("/home/oliver/IC/PAPER_PLOT/R/summary.stat.regions.gr.RData"))
    regions.rd <- ichip.regions
    colnames(regions.rd)[2] <- "Band"
    ch <- import.chain("/home/oliver/R/stuff/hg18ToHg19.over.chain")
    regions.gr <- as(regions.rd,"GRanges")
    regions.gr.37<-liftOver(regions.gr,ch)
    new.coords.df<-do.call("rbind",lapply(regions.gr.37,function(x) data.frame(start=min(start(x)),end=max(end(x)))))
    regions.gr.37<-regions.gr
    ranges(regions.gr.37)<-with(new.coords.df,IRanges(start=start,end=end))
    seqlevels(regions.gr.37)<-gsub("chr","",seqlevels(regions.gr.37))
    ichip.regions.37 <- toGenomeOrder(as(regions.gr.37,"RangedData"))
  }
  print(load("all.support.RData")) # 
  #snp.info <- get(paste(load("/chiswick/data/ncooper/immunochipRunTest/ANNOTATION/snpinfo.RData")))
  #lookup.ind <- match(clean.snp.ids(all.support$SNP),clean.snp.ids(rownames(snp.info)))
  poz <- poz37
  table.ranges <- RangedData(ranges=IRanges(start=poz,end=poz,names=ic.to.rs(rownames(table1a))),
                             space=table1a[,"CHR"],OR=table1a[,"OR_CC"],p.value=table1a[,"P_CC"],fam.OR=table1a[,"OR_Fam"],
                             fam.p.value=table1a[,"P_Fam"],meta.OR=table1a[,"OR_Meta"],meta.p.value=table1a[,"P_Meta"])
  
  cyto <- get.cyto(ucsc="hg19"); cyto[["gene"]] <- rownames(cyto)

  # rename any regions which are in the same band as another : ichip.regions.37 #
  oo <- which(duplicated(ichip.regions.37$Band))
  suffz <- rep("b",times=length(oo))
  suffz[duplicated(ichip.regions.37$Band[oo])] <- "c"
  ichip.regions.37$Band[oo] <- paste(ichip.regions.37$Band[oo],suffz,sep="_")
  #
  gs <- select.autosomes(toGenomeOrder(ichip.regions.37[,-1]))
  iiii <- (colnames(gs) %in% "Band"); if(length(iiii)>0) { colnames(gs)[iiii] <- "gene" }
  table.ranges1 <- annot.cnv(select.autosomes(toGenomeOrder(table.ranges)),gs=gs,ucsc="hg19")
  gs[["gene"]] <- paste("EXT",gs[["gene"]],sep="_")
  start(gs) <- pmax(1,(start(gs)-50000))
  end(gs) <- end(gs)+50000
  table.ranges2 <- annot.cnv(table.ranges,gs=gs,ucsc="hg19")
  gs2 <- cyto
  gs2[["gene"]] <- paste("OTHER",gs2[["gene"]],sep="_")
  table.ranges3 <- annot.cnv(table.ranges,gs=gs2,ucsc="hg19")
  table.ranges1[["gene"]][which(table.ranges1$gene=="")] <- table.ranges2[["gene"]][which(table.ranges1$gene=="")]
  table.ranges1[["gene"]][which(table.ranges1$gene=="")] <- table.ranges3[["gene"]][which(table.ranges1$gene=="")]
  
  table.ranges <- table.ranges1
  MHC <- RangedData(ranges=IRanges(start=25000000,end=35000000),space=6)
  remv <- queryHits(findOverlaps(table.ranges[6],MHC))
  #new6 <- table.ranges[6][-remv,]
  st6 <- (chrIndices(table.ranges)[6,1]-1)
  TR <- table.ranges[-(remv+st6),]
  if(length(queryHits(findOverlaps(TR[6],MHC)))) { cat("MHC not removed\n") } else { cat("MHC successfully removed\n") }
  
  
  # chr(cy.ranges) <- gsub("chr","",chr(cy.ranges))
  # cy.chr <- gsub("chr","",chr(cy.ranges))
  # annot.cnv(table.ranges,gs=cy.ranges)
  # table.rangesX <- annot.cnv(table.ranges,gs=cy.ranges)
  
  #topsnplist <- reader(fn="/chiswick/data/ncooper/iChipData/topsnplist.txt")
  #table.ranges <- annot.cnv(table.ranges,gs=gs)
}

if(T) {
  jj <- order(TR[["meta.p.value"]])
  tt <- as.data.frame(TR)[jj,]

#  kk2 <- all.support$SNP[match(topsnplist,all.support$dbSNP)]
  kk2 <- which(!(rs.to.ic(tt[,5]) %in% unique(rs.to.ic(c(qc.excluded.snps,qc.cloud.fail)))))
  tt <- tt[kk2,]
  save(tt,file="compiledTableAllResultsPassingQC.RData")
  kk <- which(as.numeric(tt[["meta.p.value"]])<(10^-5))
  tt <- tt[kk,]
  prv.large(tt[,-2:-4],rows=20,cols=6)
}


#11p15.5: equals are c("rs3842752","rs3842729","rs3842727","rs6357","rs4074905","rs4930043") but rs689 is better
#1p13.2: equals are c("rs6679677","rs2476601")  choose rs2476601 over the other [highlights() has hack for this]

if(T) {
  ## extract infuva-controls.txt  3	txt	1	2	5	6	12	1nplist.txt")
  #table.ranges <- annot.cnv(table.ranges,gs=gs)
  #do for each 'row' (regional summary)
  out.list <- tapply(tt$meta.p.value,tt$gene,highlights) # main stats
  check.list <- tapply(tt$meta.p.value,tt$gene,multihit.check)
  identicals <- NULL
  if(any(check.list)) {
    wh.list <- tapply(tt$meta.p.value,tt$gene,multihit.return)[which(check.list)]
    #prv(wh.list)
    iden.list <- list()
    cat("\nREPORT: SNPs with identical top p-values in regions\n")
    dubs <- tapply(tt$names,tt$gene,c)[which(check.list)]
    for(cc in 1:length(dubs)) { 
      idzz <- dubs[[cc]][wh.list[[cc]]]
      cat(names(wh.list)[cc],": ",idzz,"\n")
      iden.list[[cc]] <- idzz
      identicals <- c(identicals,idzz)
    } 
    cat("\n")
  }
  out.list$`1p13.2`[1] <- 2 ### cheap hack!
  out.snps <- tapply(tt$names,tt$gene,"[",1) #  top snp (1st because sorted)
  out.snps[["1p13.2"]] <- tapply(tt$names,tt$gene,"[",2)[["1p13.2"]]   ### cheap hack!
  grp.snps <- tapply(tt$names,tt$gene,c) # snp list
  # convert list to data.frame and format
  out.frame <- cbind(sapply(out.list,"[",1),sapply(out.list,"[",2),sapply(out.list,"[",3),sapply(out.list,"[",4))
  colnames(out.frame) <- c("whichSNP","best.p.value","region.hits","hits<bonferroni")
  out.frame <- cbind(out.frame,out.frame[,3]-out.frame[,4])
  colnames(out.frame)[5] <- c("hits>=bonferroni")
  out.frame[,1] <- as.character(out.frame[,1]) 
  out.frame[,1] <- out.snps
  ###
  top.snps <- out.frame[as.numeric(out.frame[,"best.p.value"])<bonf,"whichSNP"]
 # names(top.snps)[1] <- "ChrX"
  top.snps.dat <- out.frame[as.numeric(out.frame[,"best.p.value"])<bonf,]
 # save(top.snps, top.snps.dat, file="topSNPsNEWqc.RData")
  
  #print(load("topSNPsNEWqc.RData"))
  #print(load("finalMetaTopHits.RData"))
  
  # replace SNPs where P=0 and can't tell best
 # think 'b' is ok for rs34536443(b)
 ### top.snps <- gsub("imm_19_10324118","rs34536443",top.snps)
 ### top.snps.dat[,"whichSNP"] <- gsub("imm_19_10324118","rs34536443",top.snps.dat[,"whichSNP"])
  for (dd in 1:length(iden.list)) {
    next.lst <- iden.list[[dd]]
    first.snp <- next.lst[1]
    dubious <- function(x) { (unchecked(ic.to.rs(x)) | (ic.to.rs(x) %in% qc.cloud.fail)) }
    better.ns <- which(!dubious(next.lst[-1]))+1
    if(dubious(first.snp) & any(!dubious(next.lst[-1]))) {
      # replace with a better choice
      if(exists("table1snpsfinaljan30")) {
        repl <- ((next.lst[better.ns])[next.lst[better.ns] %in% table1snpsfinaljan30])[1]
      } else {
        repl <- (next.lst[better.ns])[1]
      }
      cat("replacing SNP: ",first.snp," with known/checked SNP: ",repl,"\n")
      top.snps <- gsub(first.snp,repl,top.snps)
      top.snps.dat[,"whichSNP"] <- gsub(first.snp,repl,top.snps.dat[,"whichSNP"])
    } else {
      if(exists("table1snpsfinaljan30")) {
        repl <- ((next.lst[better.ns])[next.lst[better.ns] %in% table1snpsfinaljan30])[1]
        if(!is.na(repl)) {
          cat("replacing SNP: ",first.snp," with SNP that was in the original table 1: ",repl,"\n")
          top.snps <- gsub(first.snp,repl,top.snps)
          top.snps.dat[,"whichSNP"] <- gsub(first.snp,repl,top.snps.dat[,"whichSNP"])
        } else { if((dubious(first.snp))) { cat("first snp",first.snp,"is unchecked but no good replacement\n") } }
      } else { if((dubious(first.snp))) { cat("first snp",first.snp,"is unchecked but no good replacement\n") } }
    }
  }
  top.snps <- gsub("rs3842752","rs3842727",top.snps)
  top.snps.dat[,"whichSNP"] <- gsub("rs3842752","rs3842727",top.snps.dat[,"whichSNP"])
 ## top.snps <- gsub("imm_16_11258712","rs193778",top.snps)
 ## top.snps.dat[,"whichSNP"] <- gsub("imm_16_11258712","rs193778",top.snps.dat[,"whichSNP"])
  top.snps.dat[,"whichSNP"] <- ic.to.rs(top.snps.dat[,"whichSNP"])  
  p10.5 <- as.numeric(top.snps.dat[,2])>=bonf
  pbonf <- as.numeric(top.snps.dat[,2])<bonf
  non.bonfs <- top.snps.dat[p10.5,]
  bonfs <- top.snps.dat[pbonf,]
  
  bands2 <- rownames(non.bonfs)
  b2 <- gsub("EXT_","",gsub("OTHER_","",bands2))
  b2[duplicated(b2)]
  bands1 <- rownames(bonfs)
  b1 <- gsub("EXT_","",gsub("OTHER_","",bands1))
  b1[duplicated(b1) | rev(duplicated(rev(b1)))]
  to.chuck <- c(bands1[duplicated(b1)])  # after inspection, these are those in the bonf table to remove
  duplicates <- bands1[duplicated(b1) | rev(duplicated(rev(b1)))] # show table of both members of dup pairs
  superceeded <- bands2[b2 %in% b1]  # these have a passing hit passing bonferroni on same arm
  if(length(superceeded)>0) {  non.bonfs.filt <- non.bonfs[-(match(superceeded,rownames(non.bonfs))),] } else { non.bonfs.filt <- non.bonfs }
  if(length(to.chuck)>0) { bonfs.filt <- bonfs[-(match(to.chuck,rownames(bonfs))),] }
  if(length(which(not.top.snps %in% bonfs.filt[,1]))>0) {
     bonfs.filt <- bonfs.filt[-(narm(match(not.top.snps,bonfs.filt[,1]))),] } # not.top.snps from 'hardcoded' file
  bonf.snps <- bonfs.filt[,1]
  non.bonf.snps <- non.bonfs.filt[,1]
  
  # manually pickup these overlaps
  #OTHER_12q14.1 delete as too close to 12q13.2
  #OTHER_12q24.21 delete as too close to 12q24.11;12q24.12;12q24.13
  #OTHER_12q24.13 delete as within 12q24.11;12q24.12;12q24.13
  

  # these haven't had cloud QC
  list.to.check <- bonf.snps[unchecked(bonf.snps)]
  # these failed cloud qc and need next best for replacement
  list.to.get.new <- paste(bonf.snps[bonf.snps %in% qc.cloud.fail])
  
  save(bonfs.filt,non.bonfs.filt,bonf.snps,non.bonf.snps,
     file=cat.path(work.dir,fn="finalMetaTopHits",suf=simple.date(time=F),ext="RData"))
  
  if(length(list.to.check)>0) {
   cat("please check:",paste(list.to.check,collapse=","),"\n as these haven't had cloud QC\n")
  } 
  if(length(list.to.get.new)>0) {
   cat("please replace:",paste(list.to.get.new,collapse=","),"\n as these failed cloud QC\n")
  }
}

prv.large(bonfs.filt[order(Chr(bonfs.filt[,"whichSNP"])),],cols=10,rows=60)

