
source("~/github/iChip/iFunctions.R")
work.dir <- "/chiswick/data/ncooper/iChipData/"
setwd(work.dir)
library(reader)

#setwd("~/Downloads/")
setwd("/chiswick/data/ncooper/iChipData/")
#print(load("all.resultsSat_Nov_30_13.23.02_2013.RData"))
#print(load("all.resultsFri_Dec__6_19.47.18_2013.RData"))
#print(load("all.resultsTue_Dec__3_20.16.49_2013.RData"))
print(load("/chiswick/data/ncooper/iChipData/all.resultsThu_Dec_12_16.11.37_2013.RData"))

print(load("all.resultsTue_Dec_10_22.43.19_2013.RData"))
all.results[[1]] <- all.results[[1]][-2]
all.results[[14]] <- all.results[[14]][-1]
#print(load("topSNPs.RData"))
#print(load("finalMetaTopHits.RData"))
print(load("~/Downloads/finalMetaTopHits11DEC.RData"))
top.snps <- bonf.snps
top.snps.dat <- bonfs.filt


#print(load("all.resultsSat_Nov_30_13.23.02_2013.RData"))
top.snps <- gsub("imm_19_10324118","rs34536443",top.snps)
top.snps.dat[,"whichSNP"] <- gsub("imm_19_10324118","rs34536443",top.snps.dat[,"whichSNP"])
top.snps <- gsub("rs3842752","rs3842727",top.snps)
top.snps.dat[,"whichSNP"] <- gsub("rs3842752","rs3842727",top.snps.dat[,"whichSNP"])


tab1snps <- c("rs6679677","rs3024493","rs10865035","rs35667974","rs3087243","rs113010081","rs75793288","rs17630466","rs11954020","rs12212193","rs1538171","rs1574285","rs61839660","rs12416116","rs3842727","rs1701704","rs917911","rs653178","rs56994090","rs72727394","rs34593439","rs151234","rs12927355","rs8056814","rs8067378","rs16939895","rs1615504","rs74956615","rs516246","rs6043409","rs80054410","rs4820830","rs229533","rs2611215","rs722988","rs11069349")

qc.cloud.fail <- c("imm_19_10324118","rs1016431","rs1160544","rs12743005","rs12927355","rs151234","rs1788097","rs202535","rs229533","rs2793108","rs3024493","rs34934650","rs3826110","rs4656327","rs4759229","rs4935063","rs56994090","rs5763751","rs5763779","rs61955086","rs666619","rs6691977","rs6710297","rs6907898","rs72729320","rs8074179","rs9585056","rs34536443","rs36038753") #,"rs61839660" last 2 new

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

badz <- c("imm_19_10324118","rs1016431","rs1160544","rs12927355","rs151234","rs1788097","rs202535","rs229533","rs3024493","rs3826110","rs4759229","rs56994090","rs5763751","rs5763779","rs61839660","rs6691977","rs6907898","rs72729320","rs8074179","rs9585056")
goodz <- c("rs74956615","rs7790800","rs2309837","rs9746695","rs151233","rs17207042","rs2250261","rs229536","rs3024505","rs8054218","rs705704","rs8022656","rs1989870","rs4820830","rs12722496","rs7511678","rs2326451","rs55791667","rs757411","rs9557217")
ok <- c("rs3024493","rs61839660","rs56994090","rs12927355","rs151234","rs229533")
goodz <- goodz[!badz %in% ok]
badz <- badz[!badz %in% ok]


weird <- c("ccc_2_100113638_G_A","ccc_2_100118573_G_C","ccc_2_100126398_C_A",
"ccc_2_100127085_C_T","ccc_2_100174734_C_T","ccc_2_100175954_C_T",
"ccc_2_100176881_A_G","ccc_2_100177022_A_G","ccc_22_35871557_C_A",
"imm_12_54649933","rs1690500","rs181581",
"rs2304237","rs35397078_rs8078409","rs62212954_rs66733675",
"rs72687017","rs72687027","rs73966411",
"seq_NOVEL_11243","seq_NOVEL_11254","seq_NOVEL_11256",
"seq_NOVEL_3186","seq_NOVEL_7020","seq_NOVEL_99")

to.remove.and.recalc.t1 <- bonf.snps[bonf.snps %in% badz]

qc.cloud.fail <- qc.cloud.fail[!qc.cloud.fail %in% ok]

# these table 1's must fail
# mustfail <- tab1snps[which((tab1snps %in% qc.cloud.fail) & (!tab1snps %in% ok))]
# uu <- tab1snps[tab1snps %in% badz]
# tab1snps[tab1snps %in% uu] <- goodz[match(tab1snps[tab1snps %in% uu],badz)]
# top.snps[top.snps %in% badz]


allsnpN <- c(c(1,3,5,9,10,15,16,17,19,22,24,26,27,28,30,31,34,39),c(2,7,21,27,30,33),c(8,13,35,36,37,38),11,4)
# chr 11           chr 13           chr 22           chr 42           chr 43 
# "rs2476601"      "rs6691977"      "rs3087243"     "rs75793288"      "rs2611215" 
# chr 10.1         chr 10.2         chr 10.3         chr 12.1         chr 13 
# "rs61839660"     "rs12416116"       "rs722988"       "rs917911"      "rs9585056" 
# chr 14.2         chr 15.2         chr 16.1         chr 16.2         chr 16.4 
# "rs56994090"     "rs34593439"     "rs12927355"       "rs151234"     "rs12927355" 
# chr 17          chr 19.1         chr 22.2          chr 12            chr 3 
# "rs36038753"     "rs34536443"       "rs229533"      "rs3024505" "imm_3_46488935" 
# chr 12.3         chr 16.1        chr 16.4          chr 18.2         chr 4.1 
# "rs3184504"     "rs12927355"     "rs12927355"     "rs71178827"     "rs34046593" 
# chr 6.2          chr 19.2           chr 20           chr 21         chr 22.1 
# "rs2326451"       "rs485186"       "rs202535"      "rs9981624"      "rs5763779" 
# chr 5           chr 21 
# "rs11950831"      "rs2111485" 

sum <- lapply(all.results,function(Y) { lapply(Y,suck.bic,dif=3) } )
marg.liks <- lapply(all.results,function(Y) { lapply(Y,do.bic.marg,dif=3) } )
#what <- lapply(all.results,function(Y) { lapply(Y,do.bic,dif=3) } )
vs.max.liks <- lapply(all.results,function(Y) { lapply(Y,do.bic.max,dif=3) } )
clear.not.tops <- lapply(all.results,function(Y) { lapply(Y,clearly.suck,thresh=20) } )
names(vs.max.liks) <- paste("chr",1:22)
#summary(unlist(sum))


# what is  a list by Chr, then subregion, then MAR/EXP for marginal lik
#what <- lapply(all.results,function(Y) { lapply(Y,do.bic,dif=3) } )





print(load("all.support.RData"))

#topsnplist <- reader(sl,work.dir)
topsnplist <- tab1snps   #top.snps   #readLines(cat.path(fn=sl,dir=work.dir))
grp <- list()
for (cc in 1:22) { grp[[cc]] <- topsnplist[topsnplist %in% all.support$dbSNP[all.support$Chr==cc]] }
for (cc in 1:22) { names(vs.max.liks[[cc]]) <- grp[[cc]]  }

bf3 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-3]) })})
bf5 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-5.3]) })})
bf10 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-10]) })})
bf20 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-20]) })})
bf30 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-30]) })})
bf100 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-100]) })})
bf3 <- unlist(bf3); bf5 <- unlist(bf5); bf10 <- unlist(bf10)
bf20 <- unlist(bf20); bf30 <- unlist(bf30); bf100 <- unlist(bf100)

cbind(bf3,bf5,bf10,bf20,bf30,bf100)


## remove QC failing SNPs from list and recalculate MAX comparison log BF values
# for(cc in 1:22) {
#   if(length(vs.max.liks[[cc]])<1) { next }
#   for(dd in 1:length(vs.max.liks[[cc]])) {
#     poz <- 1
#     while(names(vs.max.liks[[cc]][[dd]])[poz] %in% qc.cloud.fail) {
#       nm <- names(vs.max.liks[[cc]][[dd]])[poz]
#       vs.max.liks[[cc]][[dd]][poz] <- NA
#       vs.max.liks[[cc]][[dd]] <- vs.max.liks[[cc]][[dd]] - vs.max.liks[[cc]][[dd]][poz+1]
#       cat("removed",nm,"from grp",dd,"chr",cc,"\n")
#       poz <- poz + 1
#     }
#   }
# }

backup <- vs.max.liks; vs.max.liks <- backup
## recalc BFs based on top SNP from meta, recalculate MAX comparison log BF values
for(cc in 1:22) {
  if(length(vs.max.liks[[cc]])<1) { next }
  for(dd in 1:length(vs.max.liks[[cc]])) {
      nm <- names(vs.max.liks[[cc]])[[dd]]
      mxmx <- vs.max.liks[[cc]][[dd]][nm] 
      vs.max.liks[[cc]][[dd]] <- vs.max.liks[[cc]][[dd]] + (vs.max.liks[[cc]][[dd]][1]-mxmx)
  }
}

save(vs.max.liks,file="IndistinguishableListDEC11.RData")
