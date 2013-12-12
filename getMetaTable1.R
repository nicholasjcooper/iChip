bonf <- 3.23*(10^-7)
source("~/github/iChip/iFunctions.R")
setwd("/chiswick/data/ncooper/iChipData")
library(reader)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")

## these sent to chris ollie 10DEC13
table1snpsfinal <- c("rs6679677","rs113010081","rs75793288","rs11954020","rs12212193","rs1538171","rs1574285","rs61839660","rs12416116","rs3842727","rs1701704","rs151234","rs12927355","rs8067378","rs16939895","rs74956615","rs6043409","rs917911","rs7511678","rs3024493","rs10865035","rs35667974","rs3087243","rs17630466","rs2611215","rs722988","rs653178","rs11069349","rs1350275","rs56994090","rs72727394","rs34593439","rs8056814","rs1615504","rs516246","rs80054410","rs4820830","rs229533")

qc.cloud.fail <- c("imm_19_10324118","rs1016431","rs1160544","rs12743005","rs12927355","rs151234","rs1788097","rs202535","rs229533","rs2793108","rs3024493","rs34934650","rs3826110","rs4656327","rs4759229","rs4935063","rs56994090","rs5763751","rs5763779","rs61955086","rs666619","rs6691977","rs6710297","rs6907898","rs72729320","rs8074179","rs9585056","rs34536443","rs36038753") #,"rs12212193") #,"rs61839660" last 3 new, 2:1, last just failed 99% call rate
ok <- c("rs3024493","rs61839660","rs56994090","rs12927355","rs151234","rs229533")
qc.cloud.fail <- qc.cloud.fail[!qc.cloud.fail %in% ok]
good.snps.checked <- c("rs6679677","rs2476601","rs11552449","rs2269240","rs12568515","rs7511678","rs3024505","rs2111485",
                       "rs1534422","rs2309837","rs984971","rs3087243","rs62625034","rs75793288","rs6827756","rs6840119","rs4505848",
                       "rs2611215","rs1445898","rs71624119","rs12212193","rs2326451","rs1010474","rs7753008","rs10272724","rs7790800",
                       "rs10758594","rs61839660","rs12416116","rs12722563","rs689","rs3842727","rs7122407","rs2651830","rs11215766",
                       "rs705704","rs4378452","rs917911","rs3176775","seq-t1d-12-56449572-T-C","rs3184504","rs11066320","rs9557217","rs6573857","rs174213",
                       "rs4900384","rs8022656","rs55791667","rs34593439","rs151233","rs9746695","rs12922409","rs12924112","rs4238595",
                       "rs6498184","rs8054218","rs2290438","rs757411","rs2435200","rs1893217","rs17207042","rs74956615","rs12978105",
                       "rs10408526","rs516246","rs6043409","rs2250261","rs9981624","rs1989870","rs4820830","rs229536",
                       "rs1701704","rs10865035","rs113010081","rs11954020","rs12927355","rs1350275","rs151234",
                       "rs1574285","rs1615504","rs16939895","rs17630466","rs3024493","rs653178","rs80054410",
                       "rs35667974","rs8067378","rs11069349","rs72727394","rs8056814","rs1538171","rs722988","rs56994090",
                       "rs229533","rs7960225","rs4767000")  #rs722988 singleton to check


if(T) {
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  table1 <- reader(docs[1])
  excl <- reader(docs[2])
  prv(table1)
  prv(excl)
  table1a <- table1[-which(rownames(table1) %in% rownames(excl)),]
  table1a <- table1a[order(table1a[,3]),]
  table1a <- table1a[order(table1a[,2]),]
  #table1a <- table1a[order(table1a[,11]),]
  prv.large(table1a[,c(3,10,12,15,18)-1],rows=100,cols=7)
  poz37 <- as.numeric(table1a[,3])
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
    c<-import.chain("/home/oliver/R/stuff/hg18ToHg19.over.chain")
    regions.gr <- as(regions.rd,"GRanges")
    regions.gr.37<-liftOver(regions.gr,c)
    new.coords.df<-do.call("rbind",lapply(regions.gr.37,function(x) data.frame(start=min(start(x)),end=max(end(x)))))
    regions.gr.37<-regions.gr
    ranges(regions.gr.37)<-with(new.coords.df,IRanges(start=start,end=end))
    seqlevels(regions.gr.37)<-gsub("chr","",seqlevels(regions.gr.37))
    ichip.regions.37 <- toGenomeOrder(as(regions.gr.37,"RangedData"))
  }
  print(load("all.support.RData")) # 
  #snp.info <- get(paste(load("/chiswick/data/ncooper/immunochipRunTest/ANNOTATION/snpinfo.RData")))
  #lookup.ind <- match(clean.snp.ids(all.support$SNP),clean.snp.ids(rownames(snp.info)))
  
  table.ranges <- RangedData(ranges=IRanges(start=poz,end=poz,names=table1a[,1]),
                             space=table1a[,2],OR=table1a[,9],p.value=table1a[,11],fam.OR=table1a[,12],
                             fam.p.value=table1a[,14],meta.OR=table1a[,15],meta.p.value=table1a[,17])
  
  cyto <- get.cyto(ucsc="hg19"); cyto[["gene"]] <- rownames(cyto)

  gs <- select.autosomes(toGenomeOrder(ichip.regions.37[,-1]))
  iiii <- (colnames(gs) %in% "Band"); if(length(ii)>0) { colnames(gs)[iiii] <- "gene" }
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
  
  topsnplist <- reader(fn="/chiswick/data/ncooper/iChipData/topsnplist.txt")
  #table.ranges <- annot.cnv(table.ranges,gs=gs)
  jj <- order(TR[["meta.p.value"]])
  tt <- as.data.frame(TR)[jj,]

#  kk2 <- all.support$SNP[match(topsnplist,all.support$dbSNP)]
  kk2 <- which(!tt[,5] %in% qc.cloud.fail)
  tt <- tt[kk2,]
  kk <- which(as.numeric(tt[["meta.p.value"]])<(10^-5))
  tt <- tt[kk,]
  prv.large(tt,rows=100,cols=6)
}



if(T) {
  ## extract infuva-controls.txt  3	txt	1	2	5	6	12	1nplist.txt")
  #table.ranges <- annot.cnv(table.ranges,gs=gs)
  #do for each 'row' (regional summary)
  out.list <- tapply(tt$meta.p.value,tt$gene,highlights) # main stats
  out.snps <- tapply(tt$names,tt$gene,"[",1) #  top snp (1st because sorted)
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
  top.snps <- gsub("imm_19_10324118","rs34536443",top.snps)
  top.snps.dat[,"whichSNP"] <- gsub("imm_19_10324118","rs34536443",top.snps.dat[,"whichSNP"])
  top.snps <- gsub("rs3842752","rs3842727",top.snps)
  top.snps.dat[,"whichSNP"] <- gsub("rs3842752","rs3842727",top.snps.dat[,"whichSNP"])
  
  p10.5 <- as.numeric(top.snps.dat[,2])>=(3.23*(10^-7))
  pbonf <- as.numeric(top.snps.dat[,2])<(3.23*(10^-7))
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
  non.bonfs.filt <- non.bonfs[-(match(superceeded,rownames(non.bonfs))),]
  bonfs.filt <- bonfs[-(match(to.chuck,rownames(bonfs))),]
  
  bonf.snps <- bonfs.filt[,1]
  non.bonf.snps <- non.bonfs.filt[,1]
  
  # manually pickup these overlaps
  #OTHER_12q14.1 delete as too close to 12q13.2
  #OTHER_12q24.21 delete as too close to 12q24.11;12q24.12;12q24.13
  #OTHER_12q24.13 delete as within 12q24.11;12q24.12;12q24.13
  
  # these haven't had cloud QC
  list.to.check <- paste(bonf.snps[(!bonf.snps %in% qc.cloud.fail) & (!bonf.snps %in% good.snps.checked)])
  # these failed cloud qc and need next best for replacement
  list.to.get.new <- paste(bonf.snps[bonf.snps %in% qc.cloud.fail])
  
  save(bonfs.filt,non.bonfs.filt,bonf.snps,non.bonf.snps,file="~/Downloads/finalMetaTopHits11DEC.RData")
  
  
}
