
source("~/github/iChip/iFunctions.R")
work.dir <- "/chiswick/data/ncooper/iChipData/"
setwd(work.dir)
library(reader)

#setwd("~/Downloads/")
setwd("/chiswick/data/ncooper/iChipData/")
#print(load("all.resultsSat_Nov_30_13.23.02_2013.RData"))
#print(load("all.resultsFri_Dec__6_19.47.18_2013.RData"))
#print(load("all.resultsTue_Dec__3_20.16.49_2013.RData"))

# overall
#print(load("all.resultsFri_Jan_24_10.10.55_2014.RData"))
print(load("all.resultsFri_Feb__7_10.34.44_2014.RData"))

# conditional
#print(load("indcond.resultsFri_Jan_24_13.31.32_2014.RData"))
print(load("indcond.resultsFri_Feb__7_14.38.04_2014.RData"))

if(F) {
  # main table 1 results
  #print(load("all.resultsTue_Dec_10_22.43.19_2013.RData"))
  print(load("all.resultsThu_Dec_12_18.03.38_2013.RData"))
  all.results[[1]] <- all.results[[1]][-2]
  all.results[[14]] <- all.results[[14]][-1]
  
  #BIC for rs689 == 22949, = BF mx - 11474.5  for rs3842727 = 11577.47 ... so not even in ballpark
} else {
  # latest conditional bayes results
#  print(load("/chiswick/data/ncooper/iChipData/indcond.resultsFri_Jan_24_16.27.45_2014.RData"))
  print(load("indcond.resultsFri_Feb__7_14.38.04_2014.RData"))
  vs.max.liks <- lapply(all.results,do.bic.max,dif=3)
  ct.fn <- "conditionalTests.csv"
  ct <- reader(ct.fn,stringsAsFactors=F)
  ct[,1] <- ic.to.rs(ct[,1]); ct[,2] <- ic.to.rs(ct[,2])
  condlist <- apply(ct,1,function(X) { paste(c(narm(X[1:2]))) })
  # ^instead of print(load("conditionalSummaryObjects.RData"))
  vs.max.liks <- lapply(vs.max.liks,function(X) { 
         if(any(duplicated(names(X)))) { print("removed duplicates") }; X[!duplicated(names(X))] }) #remove duplicates
  # get grpz,cond,condlist,glm.results 
  #condlist <- list("rs35667974","rs61839660",c("rs61839660", "rs10795791"),
  #                 "rs3842727",c("rs3842727", "rs11043003"),c("rs3842727", "rs11043003", "rs6357"),
  #                 "rs12927355","rs8067378","rs16939895","rs74956615")
  condon <- sapply(condlist,function(x) paste(x,collapse="_"))
  names(vs.max.liks) <- condon

  # manually remove some duplicates
#  vs.max.liks[[8]] <- (vs.max.liks[[8]][-2])
#  vs.max.liks[[12]] <- (vs.max.liks[[12]][-2])
  
  bf3 <- lapply(vs.max.liks,function(X) { length(X[X>-3]) })
  bf5 <- lapply(vs.max.liks,function(X) { length(X[X>-5.29]) })
  bf10 <- lapply(vs.max.liks,function(X) { length(X[X>-10]) })
  bf20 <- lapply(vs.max.liks,function(X) { length(X[X>-20]) })
  bf30 <- lapply(vs.max.liks,function(X) { length(X[X>-30]) })
  bf100 <- lapply(vs.max.liks,function(X) { length(X[X>-100]) })
  cbind(bf3,bf5,bf10,bf20,bf30,bf100)

#  vs.max.liks <- tapply(vs.max.liks,factor(ct$CHR),c)

  to.move.to.prime <- vs.max.liks[which(ct$PRIME==1)]
  names(to.move.to.prime) <- ct$TABLE1[which(ct$PRIME==1)]
  to.keep.as.cond <- vs.max.liks[which(ct$COND!=0)]
  names(to.keep.as.cond) <- ct$TABLE1[which(ct$COND!=0)]
  st.fn <- "/chiswick/data/ncooper/iChipData/spectable1bf.csv"
  spec <- reader(st.fn)
  nms.for.prime <- ic.to.rs(rownames(spec)[spec$MULTI==0])
  nms.for.cond <- ic.to.rs(rownames(spec)[spec$MULTI==1])
  save(vs.max.liks,file="/chiswick/data/ncooper/iChipData/conditionalBFlists_rs689.RData")
}
#print(load("topSNPs.RData"))
#print(load("finalMetaTopHits.RData"))

# get hard coded lists of SNPs passing, failing QC, appearing in table 1, etc #tab1snps
source('~/github/iChip/hardCodedSnpLists.R', echo=FALSE)

## DONT NEED THESE NOW
#to.remove.and.recalc.t1 <- bonf.snps[bonf.snps %in% badz]
#qc.cloud.fail <- qc.cloud.fail[!qc.cloud.fail %in% ok]

# these table 1's must fail
# mustfail <- tab1snps[which((tab1snps %in% qc.cloud.fail) & (!tab1snps %in% ok))]
# uu <- tab1snps[tab1snps %in% badz]
# tab1snps[tab1snps %in% uu] <- goodz[match(tab1snps[tab1snps %in% uu],badz)]
# top.snps[top.snps %in% badz]


sum <- lapply(all.results,function(Y) { lapply(Y,suck.bic,dif=3) } )
marg.liks <- lapply(all.results,function(Y) { lapply(Y,do.bic.marg,dif=3) } )
#what <- lapply(all.results,function(Y) { lapply(Y,do.bic,dif=3) } )
vs.max.liks <- lapply(all.results,function(Y) { lapply(Y,do.bic.max,dif=3) } )
clear.not.tops <- lapply(all.results,function(Y) { lapply(Y,clearly.suck,thresh=20) } )
names(vs.max.liks) <- paste("chr",1:22,sep="")
#summary(unlist(sum))


# what is  a list by Chr, then subregion, then MAR/EXP for marginal lik
#what <- lapply(all.results,function(Y) { lapply(Y,do.bic,dif=3) } )

# ALTERNATIVE VERSION OF VS.LIK.MAX when only 1 level of list
vs.max.liks <- lapply(all.results,do.bic.max,dif=3) 



## load in support data for all SNPs so can match which names are on each chromosome
print(load("all.support.RData"))

#print(load("finalMetaTopHits15JAN.RData")) # new top list since mon 15/1/14, some removals curated manually 
print(load("finalMetaTopHitsi6FEB.RData"))
topsnplist <- paste(bonf.snps)
topsnplist[topsnplist=="rs3842727"] <- "rs689"
topsnplist[!topsnplist %in% qc.cloud.fail]

#topsnplist <- reader(sl,work.dir)
## Add correct names to the vs.max.liks list object
#topsnplist <- tab1snps   #top.snps   #readLines(cat.path(fn=sl,dir=work.dir))

tsl <- gsub("b","",ic.to.rs(rownames(spec)))
grp <- list()
#for (cc in 1:22) { grp[[cc]] <- topsnplist[topsnplist %in% all.support$dbSNP[all.support$Chr==cc]] }
for (cc in 1:22) { 
  next.cand <- tsl[Chr(tsl)==cc]
  ll <- length(vs.max.liks[[cc]])
  if(ll>0) {
    nn <- character(ll)
    for(dd in 1:ll) {
      nlist <- gsub("b","",names(vs.max.liks[[cc]][[dd]]))
      nn[dd] <- nlist[which(nlist %in% tsl)[1]] 
    }
    grp[[cc]] <- nn
  }
}


for (cc in 1:22) { names(vs.max.liks[[cc]]) <- grp[[cc]]  }
#vs.max.liks <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { X[!duplicated(names(X))] })}) #remove duplicates
vs.max.liks <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { names(X) <- rs.to.ic(names(X)); return(X) })}) 

vml <- do.call(c,vs.max.liks)
snp <- names(to.move.to.prime)[2]
grep(snp,names(vml))
vs.max.liks[[2]][[4]] <- to.move.to.prime[[1]]
vs.max.liks[[2]] <- vs.max.liks[[2]][c(1,4,2,3)]
names(vs.max.liks[[2]])[2] <- "rs35667974"

# generate a table of summary counts for snps passing bayes factor thresholds
bf3 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-3]) })})
bf5 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-5.3]) })})
bf10 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-10]) })})
bf20 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-20]) })})
bf30 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-30]) })})
bf100 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-100]) })})
bf3 <- unlist(bf3); bf5 <- unlist(bf5); bf10 <- unlist(bf10)
bf20 <- unlist(bf20); bf30 <- unlist(bf30); bf100 <- unlist(bf100)

cbind(bf3,bf5,bf10,bf20,bf30,bf100)



### Bayes factors should be relative to the top SNP from the meta analysis which isn't always
# the top SNP, so this loop fixes this, means sometimes log-BFs will be positive
backup <- vs.max.liks; vs.max.liks <- backup
if(F) {
  # hand editing to fix probs in table due to multiple p values of 0
  #names(vs.max.liks[[1]])[1] <- "rs2476601"
  names(vs.max.liks[[11]])[1] <- "rs689"
  vs.max.liks[[16]][[3]] <- NULL  # because list was the same for rs193778 as for rs12927355
}
## recalc BFs based on top SNP from meta, recalculate MAX comparison log BF values
for(cc in 1:22) {
  if(length(vs.max.liks[[cc]])<1) { next }
  for(dd in 1:length(vs.max.liks[[cc]])) {
      nm <- names(vs.max.liks[[cc]])[[dd]]
      mxmx <- vs.max.liks[[cc]][[dd]][nm] 
      vs.max.liks[[cc]][[dd]] <- vs.max.liks[[cc]][[dd]] + (vs.max.liks[[cc]][[dd]][1]-mxmx)
  }
}

# for indistinguishable list (previously) needed to manually insert rs689, with diff of 102.97
if(F) {
  vs.max.liks[[11]][[1]] <- c(102.97,vs.max.liks[[11]][[1]])
  names(vs.max.liks[[11]][[1]])[1] <- "rs689"
  vs.max.liks[[11]][[1]] <- vs.max.liks[[11]][[1]]-vs.max.liks[[11]][[1]][1]
}

# for indistinguishable list need to manually insert rs34536443, with diff of 0.735469 to rs74956615
if(T) {
  vs.max.liks[[19]][[1]] <- c(vs.max.liks[[19]][[1]][1],-0.735469,vs.max.liks[[19]][[1]][-1])
  names(vs.max.liks[[19]][[1]])[2] <- "rs34536443"
}

## save the final main list generated
save(vs.max.liks,file="IndistinguishableList_Jan24.RData")

## save the final conditional list generated
#save(vs.max.liks,file="IndistinguishableCondList_Jan24.RData")
