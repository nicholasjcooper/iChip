#working to generate ichip annotation (not used)



ichip.dir <- "/chiswick/data/ncooper/iChipData/"
impute.dir <- "/chiswick/data/ncooper/imputation/common/"
data.dir <- paste0(impute.dir,"ALL.integrated_phase1_SHAPEIT_16-06-14.nosing/")
temp.txt <- "nextCHR.txt"

source("~/github/imputer/imputeFunctions.R")


options(ucsc="hg19")

#update allele info in ichip object
if(F) {
  
  rm(all.support)
  ichip37 <- chip.support(build=37)
  ichip36 <- chip.support(build=36)
  
  
  #### generate fresh all.support ####
  print(load(cat.path(ichip.dir,"temp.ichip-data1.RData")))
  all.support <- annotated.snp.support
  
  for (cc in c(2:22,"X","MT")) {
    print(load(cat.path(ichip.dir,"temp.ichip-data",suf=cc,ext="RData")))
    all.support <- rbind(all.support,annotated.snp.support)
  }
  
  save(all.support,file=cat.path(ichip.dir,"all.support.refreshed.RData"))
  
  print(load(cat.path(ichip.dir,"all.support.refreshed.RData")))
  rownames(all.support) <- clean.snp.ids(rownames(all.support))
  all.support[["Pos37"]] <- conv.36.37(chr=all.support$Chr,pos=all.support$Pos,ids=rownames(all.support))[,"start"]
  
  all.support[["dbSNP"]][is.na(all.support[["dbSNP"]])] <- rownames(all.support)[is.na(all.support[["dbSNP"]])]
  dil.support <- all.support; rm(all.support)
  ####################################
  
  #### generate fresh all.rr ####
  chrnum <- 1
  next.rn <- cat.path(impute.dir,fn="CHR",suf=chrnum,ext="RData")
  print(load(next.rn))
  all.rr <- rr
  for (chrnum in c(2:22)) {
    next.rn <- cat.path(impute.dir,fn="CHR",suf=chrnum,ext="RData")
    print(load(next.rn))
    all.rr <- rbind(all.rr,rr)
  }
  rownames(all.rr) <- clean.snp.ids(rownames(all.rr))
  save(all.rr,file=cat.path(impute.dir,"all.rr.RData"))
  print(load(cat.path(impute.dir,"all.rr.RData")))
  ###############################
  mchip <- reader("/chiswick/data/ncooper/iChipData/ImChip_T1D_20130913/metamerged_ImmChip_CC_fam_MAF_09132013.txt")
  mchip <- mchip[,c(1,2,3,4,5,6,8)]
  rownames(mchip) <- clean.snp.ids(rownames(mchip))
  mchip2 <- cbind(mchip[,1:3],alphabetize.alleles(mchip[,4],mchip[,5], mchip[,7]))
  
  immunobase <- reader("/chiswick/data/ncooper/iChipData/iChipAnnot.txt")
  immunobase[["rs.id"]] <- rownames(immunobase); rownames(immunobase) <- immunobase[[1]]
  immunobase <- immunobase[,-1]; colnames(immunobase) <- c("Chr","Pos","Alleles","rs.id")
  
  # only two ichip dups by chr,pos
  #immunobase[dup.pairs(i.ch.ps),]
  #Chr       Pos Alleles      rs.id
  #imm_12_54649934  12  56363666     A/C  rs2069407    #18421
  #imm_12_54649933  12  56363666     A/C  rs3213122    #19148
  extrai <- immunobase[18421,]
  immunobase <- immunobase[-18421,]
  
  
  # we have: dil.support (dbSNP, Chr, Pos , Pos37),  [hg18, 19]  193997*11, index chip.id
  #          all.rr (id, position),      [hg19]             131757*7 , index rs.id .. position maybe better
  #          ichip36/37 (rs.id, chr(), start()) [hg18/19]  196524  , index chip.id
  #          mchip (SNP_ID, CHR, Position)  [hg19]  154939*7  , index chip.id
  #          immunobase (rs.id, Chr, Pos) [hg19]   175868*4   , index chip.id? varies.. pos only complete
  #
  # rules
  # if the ichip is rs-id, then the rs-id should be rs-id
  # if the rs-id is the chip id, and the chip id is not rs-id, and there is an immunobase-id, use that
  # if the rs-id is the chip id, and the chip id is not rs-id, and there is no immunobase-id, use the rr-id
  # if no rr-id, use all.support. if no all.support, use mchip.
  # for each tabulate the number of rs-ids amongst the 4 refs, all.support, rr, mchip and immunobase, and the % for most common one
  chipids <- rownames(ichip37)
  
  categs <- c("chipid.is.rs","chipid.in.immuno.with.pos.match","chipid.in.immuno.without.pos.match")
  
  new.rs <- character(length(chipids))
  a1 <- a2 <- character(length(chipids))
  ch <- chr(ichip37); ps <- pos(ichip37); rsid <- rs.id(ichip37)
  ch.ps <- paste(ch,ps,sep="_")
  
  #chip id is an rs id?
  is.rs <- (substr(rownames(ichip37),1,2)=="rs"); is.rs <- rep(F,length(new.rs))
  new.rs[is.rs] <- chipids[is.rs]
  i.ch.ps <- paste(immunobase[["Chr"]],immunobase[["Pos"]],sep="_")
  ib.lookup <- match(ch.ps,i.ch.ps)
  not.ib <- (is.na(ib.lookup))
  new.rs[!is.rs & !not.ib] <- immunobase$rs.id[ib.lookup[!is.rs & !not.ib]]
  still.missing <- (new.rs=="")
  
  d.ch.ps <- paste(dil.support[["Chr"]],dil.support[["Pos37"]],sep="_")
  indo <- match(ch.ps[still.missing],d.ch.ps)
  drs <- dil.support2[indo,"dbSNP"]
  wh <- substr(drs,1,2)=="rs"
  wh[is.na(wh)] <- FALSE
  new.rs[still.missing][wh] <- drs[indo][wh]
  new.rs[is.na(new.rs)] <- ""
  still.missing <- (new.rs=="")
  
  muse <- substr(mchip$SNP_ID,1,2)=="rs"
  mref <- paste(mchip[muse,"CHR"],mchip[muse,"Position"],sep="_")
  
  indd <- match(ch.ps[still.missing],mref)
  new.rs[still.missing][!is.na(indd)] <- mchip[muse,"SNP_ID"][indd[!is.na(indd)]]
  still.missing <- (new.rs=="")
  
  indy <- match(ps[still.missing],all.rr$position)
  new.rs[still.missing][!is.na(indy)] <- all.rr[,"id"][indy[!is.na(indy)]]
  still.missing <- (new.rs=="")
  
  save(new.rs,file="rsFromManySources2.RData")
  new.rs[new.rs==""] <- chipids[new.rs==""]
  
  new.rs[(grep("AMBIG",new.rs))] <- rs.id(ichip37[chipids[(grep("AMBIG",new.rs))],])
  save(new.rs,file="rsFromManySources3.RData") # but contains errrors  - original seems superior
  
  ## alleles
  dil.support <- dil.support[-which(is.na(dil.support[,4]) & is.na(dil.support[,5])),]
  #
  indx <- match(rownames(ichip37),rownames(dil.support))
  need.alleles <- is.na(indx)
  a1[!need.alleles] <- dil.support$allele.A[indx[!need.alleles]]
  a2[!need.alleles] <- dil.support$allele.B[indx[!need.alleles]]
  
  #######
  
  print(load("/chiswick/data/ncooper/iChipData/ImmunoChip_ChipInfo.RData"))
  table(paste0(A1(ichip37[[19]]),A2(ichip37[[19]])),exclude=NULL)
  
  
  bad.na1 <- which(is.na(A1(sup.for.impute)))
  bad.na2 <- which(is.na(A2(sup.for.impute)))
  any.bad <- unique(c(bad.na1,bad.na2))
  if(length(any.bad)>0) {
    cat("missing alleles\n")
    print(sup.for.impute[any.bad,])
  }
  equal.ones <- which(A1(sup.for.impute)==A2(sup.for.impute))
  if(length(equal.ones)>0) {
    cat("equal alleles\n")
    print(sup.for.impute[equal.ones,])
  }
  
  mcols(ichip36["X17_57823058",])[,"A2"] <- "T"
  mcols(ichip37["X17_57823058",])[,"A2"] <- "T"
  save(ichip36,ichip37,file="/chiswick/data/ncooper/iChipData/ImmunoChip_ChipInfo.RData")
  cs <- col.summary(gt.for.impute)
  cs[rs.to.id("rs689"),]
  mcols(ichip36["imm_11_2138800",])[,"A1"] <- "A"
  mcols(ichip36["imm_11_2138800",])[,"A2"] <- "T"
  mcols(ichip37["imm_11_2138800",])[,"A1"] <- "A"
  mcols(ichip37["imm_11_2138800",])[,"A2"] <- "T"
  save(ichip36,ichip37,file="/chiswick/data/ncooper/iChipData/ImmunoChip_ChipInfo.RData")
  
  # to do
  mcols(ichip36["imm_3_50851128",])[,"rs.id"] <- "rs1458385b"
  mcols(ichip37["imm_3_50851128",])[,"rs.id"] <- "rs1458385b"
  mcols(ichip36["imm_3_50851128",])[,"A2"] <- "C"
  mcols(ichip37["imm_3_50851128",])[,"A2"] <- "C"
  
  mcols(ichip36["imm_3_50896281",])[,"A1"] <- "C"
  mcols(ichip36["imm_3_50896281",])[,"A2"] <- "T"
  mcols(ichip37["imm_3_50896281",])[,"A1"] <- "C"
  mcols(ichip37["imm_3_50896281",])[,"A2"] <- "T"
  
  # to do
  mcols(ichip36["imm_2_234139989",])[,"rs.id"] <- "rs13387539b"
  mcols(ichip37["imm_2_234139989",])[,"rs.id"] <- "rs13387539b"
  mcols(ichip36["imm_2_234139989",])[,"A2"] <- "G"
  mcols(ichip37["imm_2_234139989",])[,"A2"] <- "G"
  mcols(ichip36["imm_2_234139989",])[,"A1"] <- "A"
  mcols(ichip37["imm_2_234139989",])[,"A1"] <- "A"
  
  # to do
  mcols(ichip36["imm_1_2486652",])[,"rs.id"] <- "rs2227313b"
  mcols(ichip37["imm_1_2486652",])[,"rs.id"] <- "rs2227313b"
  mcols(ichip36["imm_1_2486652",])[,"A1"] <- "C"
  mcols(ichip36["imm_1_2486652",])[,"A2"] <- "T"
  mcols(ichip37["imm_1_2486652",])[,"A1"] <- "C"
  mcols(ichip37["imm_1_2486652",])[,"A2"] <- "T"
  
  mcols(ichip36["imm_1_2476740",])[,"rs.id"] <- "rs11577783b"
  mcols(ichip37["imm_1_2476740",])[,"rs.id"] <- "rs11577783b"
  mcols(ichip36["imm_1_2476740",])[,"A2"] <- "G"
  mcols(ichip37["imm_1_2476740",])[,"A2"] <- "G"
  mcols(ichip36["imm_1_2476740",])[,"A1"] <- "A"
  mcols(ichip37["imm_1_2476740",])[,"A1"] <- "A"
  
  mcols(ichip36["rs3734837",])[,"A2"] <- "G"
  mcols(ichip37["rs3734837",])[,"A2"] <- "G"
  mcols(ichip36["rs3734837",])[,"A1"] <- "A"
  mcols(ichip37["rs3734837",])[,"A1"] <- "A"
  
  mcols(ichip36["rs3132662",])[,"A2"] <- "A"
  mcols(ichip37["rs3132662",])[,"A2"] <- "A"
  mcols(ichip36["rs3132662",])[,"A1"] <- "G"
  mcols(ichip37["rs3132662",])[,"A1"] <- "G"
  
  mcols(ichip36["rs9261424",])[,"A2"] <- "G"
  mcols(ichip37["rs9261424",])[,"A2"] <- "G"
  mcols(ichip36["rs9261424",])[,"A1"] <- "C"
  mcols(ichip37["rs9261424",])[,"A1"] <- "C"
  #rownames("imm_1_2478204"
  
  mcols(ichip36["imm_1_2480888",])[,"rs.id"] <- "rs7554260"
  mcols(ichip37["imm_1_2480888",])[,"rs.id"] <- "rs7554260"
  mcols(ichip36["imm_1_2485683",])[,"rs.id"] <- "imm_1_2485683"
  mcols(ichip37["imm_1_2485683",])[,"rs.id"] <- "imm_1_2485683"
  start(ichip36)[match("imm_1_2485683",rownames(ichip36))] <- 2485683
  start(ichip37)[match("imm_1_2485683",rownames(ichip37))] <- 2495826
  
  mcols(ichip36["imm_2_234149786",])[,"rs.id"] <- "imm_2_234149786"
  mcols(ichip37["imm_2_234149786",])[,"rs.id"] <- "imm_2_234149786"
  mcols(ichip36["imm_2_234138651",])[,"rs.id"] <- "imm_2_234138651"
  mcols(ichip37["imm_2_234138651",])[,"rs.id"] <- "imm_2_234138651"
  start(ichip36)[match("imm_2_234138651",rownames(ichip36))] <- 234138651
  start(ichip37)[match("imm_2_234138651",rownames(ichip37))] <- 234473912
  
  mcols(ichip36["imm_2_234144460",])[,"rs.id"] <- "imm_2_234144460"
  mcols(ichip37["imm_2_234144460",])[,"rs.id"] <- "imm_2_234144460"
  mcols(ichip36["imm_2_234147184",])[,"rs.id"] <- "rs62192778"
  mcols(ichip37["imm_2_234147184",])[,"rs.id"] <- "rs62192778"
  start(ichip36)[match("imm_2_234144460",rownames(ichip36))] <- 234144460
  start(ichip37)[match("imm_2_234144460",rownames(ichip37))] <- 234475415
  
  
  save(ichip36,ichip37,file="/chiswick/data/ncooper/iChipData/ImmunoChip_ChipInfo_New.RData")
  
  ## REPORT ON ALL DUPLICATED POSITIONS ##
  for(chrnum in c(1,2,5:13,16:19)) { show((sup[chr(sup)==chrnum,])[dup.pairs(start(sup[chr(sup)==chrnum])),]) }
  ## fix allele codes 
}