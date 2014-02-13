#setwd("~/Downloads/")
setwd("/chiswick/data/ncooper/iChipData/")
#print(load("all.resultsSat_Nov_30_13.23.02_2013.RData"))
print(load("all.resultsFri_Dec__6_19.47.18_2013.RData"))
#print(load("all.resultsTue_Dec__3_20.16.49_2013.RData"))
print(load("topSNPs.RData"))
print(load("finalMetaTopHits.RData"))

#print(load("all.resultsSat_Nov_30_13.23.02_2013.RData"))
top.snps <- gsub("imm_19_10324118","rs34536443",top.snps)
top.snps.dat[,"whichSNP"] <- gsub("imm_19_10324118","rs34536443",top.snps.dat[,"whichSNP"])
top.snps <- gsub("rs3842752","rs3842727",top.snps)
top.snps.dat[,"whichSNP"] <- gsub("rs3842752","rs3842727",top.snps.dat[,"whichSNP"])


tab1snps <- c("rs2476601","rs3024493","rs35667974","rs3087243","rs17630466","rs75793288","rs12212193","rs6907898","rs1574285","rs61839660","rs12416116","rs3842727","rs917911","rs1701704","rs653178","rs9585056","rs1350275","rs56994090","rs72729320","rs34593439","rs12927355","rs151234","rs3826110","rs36038753","rs16939895","rs1615504","rs34536443","rs516246","rs6043409","rs80054410","rs4820830","rs229533","rs10865035","rs113010081","rs11954020","rs12922409","rs6691977","rs2611215","rs722988")

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
                       "rs1574285","rs1615504","rs16939895","rs17630466","rs3024493","rs653178","rs80054410")

badz <- c("imm_19_10324118","rs1016431","rs1160544","rs12927355","rs151234","rs1788097","rs202535","rs229533","rs3024493","rs3826110","rs4759229","rs56994090","rs5763751","rs5763779","rs61839660","rs6691977","rs6907898","rs72729320","rs8074179","rs9585056")
goodz <- c("rs74956615","rs7790800","rs2309837","rs9746695","rs151233","rs17207042","rs2250261","rs229536","rs3024505","rs8054218","rs705704","rs8022656","rs1989870","rs4820830","rs12722496","rs7511678","rs2326451","rs55791667","rs757411","rs9557217")
ok <- c("rs3024493","rs61839660","rs56994090","rs12927355","rs151234","rs229533")
goodz <- goodz[!badz %in% ok]
badz <- badz[!badz %in% ok]

to.remove.and.recalc.t1 <- bonf.snps[bonf.snps %in% badz]

qc.cloud.fail <- qc.cloud.fail[!qc.cloud.fail %in% ok]

# these table 1's must fail
mustfail <- tab1snps[which((tab1snps %in% qc.cloud.fail) & (!tab1snps %in% ok))]
uu <- tab1snps[tab1snps %in% badz]
tab1snps[tab1snps %in% uu] <- goodz[match(tab1snps[tab1snps %in% uu],badz)]
top.snps[top.snps %in% badz]

# p <10^-5 but not passing bonferroni
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
to.chuck <- c(bands1[duplicated(b1)][-10],"EXT_2q24.2")  # after inspection, these are those in the bonf table to remove
duplicates <- bands1[duplicated(b1) | rev(duplicated(rev(b1)))] # show table of both members of dup pairs
superceeded <- bands2[b2 %in% b1]  # these have a passing hit passing bonferroni on same arm
non.bonfs.filt <- non.bonfs[-(match(superceeded,rownames(non.bonfs))),]
bonfs.filt <- bonfs[-(match(to.chuck,rownames(bonfs))),]

bonf.snps <- bonfs.filt[,1]
non.bonf.snps <- non.bonfs.filt[,1]

# these haven't had cloud QC
list.to.check <- paste(bonf.snps[(!bonf.snps %in% qc.cloud.fail) & (!bonf.snps %in% good.snps.checked)])
# these failed cloud qc and need next best for replacement
list.to.get.new <- paste(bonf.snps[bonf.snps %in% qc.cloud.fail])

# save curated table1 and p<10^-5 table, plus snp-lists to file (binary)
save(bonfs.filt,non.bonfs.filt,bonf.snps,non.bonf.snps,file="~/Downloads/finalMetaTopHits.RData")

# for table 1 snps, where are they in my CC analysis:
# these are top in mine: c(1,3,5,9,10,15,16,17,19,22,24,26,27,28,30,31,34,39)
# these are #2 in mine: c(2,7,21,27,30,33)
# these are #3 in mine: c(8,13,35,36,37,38)
# these are #4 in mine: 11
# these are #5 in mine: 4

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

##' @author Claudia Giambartolomei
logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

suck.bic <- function(X,dif=3) {
  bic <- X$BIC;
}

do.bic.marg <- function(X,dif=3) {
  bic <- X$BIC;
  #which.min <- (which(bic==min(bic)))
  denom <- logsum(-bic/2)
  #logsum((-bic/2)[-which.min])
  next.marg <- next.exp <- numeric(length(bic))
  for (cc in 1:length(bic)) {
    which.min <- which(bic==sort(bic)[cc])[1]
    next.min <- (-bic/2)[which.min]
    next.marg[cc] <- next.min - denom
    next.exp[cc] <- (exp(next.marg[cc])); #print(next.exp[cc])
  }
  names(next.marg) <- names(next.exp) <- names(bic)[order(bic)]
  #print(paste("SNP:",names(bic)[top.grp]))
  #lapply(X$GLM[top.grp],function(Z) { print(tail(Z,1)[4]) })
  return(list(MAR=next.marg,EXP=next.exp))
}

do.bic.max <- function(X,dif=3) {
  bic <- X$BIC;
  max.bic <- -min(bic/2)
  dif.bic <- (-bic/2)-max.bic
  return(rev(sort(dif.bic)))
}

clearly.suck <- function(X,thresh=100) {
  bic <- X$BIC;
  max.bic <- -min(bic/2)
  dif.bic <- (-bic/2)-max.bic
  failers <- (rev(sort(dif.bic)))
  failers <- failers[abs(failers)>=thresh]
  return(names(failers))
}

sum <- lapply(all.results,function(Y) { lapply(Y,suck.bic,dif=3) } )
marg.liks <- lapply(all.results,function(Y) { lapply(Y,do.bic.marg,dif=3) } )
#what <- lapply(all.results,function(Y) { lapply(Y,do.bic,dif=3) } )
vs.max.liks <- lapply(all.results,function(Y) { lapply(Y,do.bic.max,dif=3) } )
clear.not.tops <- lapply(all.results,function(Y) { lapply(Y,clearly.suck,thresh=20) } )
names(vs.max.liks) <- paste("chr",1:22)
#summary(unlist(sum))


# what is  a list by Chr, then subregion, then MAR/EXP for marginal lik
#what <- lapply(all.results,function(Y) { lapply(Y,do.bic,dif=3) } )


if(F) { 
  ## plot bayes against marginal likelihoods
  wh <- unlist(what,recursive=FALSE)
  wh <- unlist(wh,recursive=FALSE)
  wh1 <- wh[1+c(0:38)*2]
  wh2 <- wh[2+c(0:38)*2]
  wh1 <- unlist(wh1)
  wh2 <- unlist(wh2)
  plot(wh1[order(wh1)],wh2[order(wh1)],type="l",xlim=c(-11,0),ylim=c(0,0.01))
  abline(v=-1,col="green")
  abline(v=-2,col="blue")
  abline(v=-3,col="yellow")
  abline(v=-4,col="orange")
  abline(v=-5,col="red")
  abline(v=-6,col="purple")
  abline(v=-8,col="brown")
  abline(v=-10,col="black")
  abline(h=0.05,col="black")
  abline(h=0.001,col="grey")
  abline(h=0.00001,col="red")
}
# 3= .05, 0.01=4.6, .005=5.29, .001=6.9, .0005=7.59, 
# .0001 = 9.2, .00005=9.89, 0.00004543 = 10,  .00001=11.48



sl <- "nonbonf.txt"
sl <- "table1snplist.txt"
print(load("all.support.RData"))
work.dir <- getwd()
#topsnplist <- reader(sl,work.dir)
topsnplist <- readLines(cat.path(fn=sl,dir=work.dir))
grp <- list()
for (cc in 1:22) { grp[[cc]] <- topsnplist[topsnplist %in% all.support$dbSNP[all.support$Chr==cc]] }
for (cc in 1:22) { names(vs.max.liks[[cc]]) <- grp[[cc]]  }

bf3 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-3]) })})

bf3 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-3]) })})
bf5 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-5.3]) })})
bf10 <- lapply(vs.max.liks,function(Y) { lapply(Y,function(X) { length(X[X>-10]) })})
bf3 <- unlist(bf3); bf5 <- unlist(bf5); bf10 <- unlist(bf10)
cbind(bf3,bf5,bf10)

data.frame.to.SnpMatrix <- function(X){
  NN <- as.matrix(X)
  NN <- round(NN)
  SS <- as(NN,"SnpMatrix")
  return(SS)
}



if(!covs) {
  fm <- paste("Pheno ~ ",grp.labs[[grp]][hero])
} else {
  fm <- paste("Pheno ~ sex + region +",grp.labs[[grp]][hero])
}
snps.in.cond <- grp.labs[[grp]][hero]
keepgoing <- T
while(keepgoing) {
  cat("Testing",ncol(myDat)-length(snps.in.cond),"Snps against base model:",fm,"\n")
  result <- snp.rhs.tests(as.formula(fm),snp.data=myDatSnp,data=myDat)
  minp <- min(p.value(result),na.rm=T)
  if(minp<=bonf) { bst <- which(p.value(result)==minp) } else { bst <- NA }
  if(length(bst)>1) { warning("multiple best conditional snps:",paste(names(result)[bst],collapse=",")) }
  if(!is.na(bst)) { 
    best.snp <- names(result)[bst[1]] 
    snps.in.cond <- c(snps.in.cond,best.snp)
    fm <- paste(fm,best.snp,sep=" + ")
  } else { 
    best.snp <- paste(NA) 
    keepgoing <- FALSE
  }
}
      



#top snps not NA
top.snps1 <- lapply(vs.max.liks,function(Y) { 
  lapply(Y,function(X) { u <- names(X)[!is.na(X)][1] })})
# top snps not in cloud fail list
top.snps <- lapply(vs.max.liks,function(Y) { 
  lapply(Y,function(X) { u <- names(X)[1] ; cc<-1; while(u %in% qc.cloud.fail) { cc <- cc+1 ; u <- names(X)[cc] }; return(u) })})

top.snps <- lapply(fail.bonf.liks,function(Y) { 
  lapply(Y,function(X) { u <- names(X)[1]; cc<-1; while(u %in% qc.cloud.fail) { cc <- cc+1 ; u <- names(X)[cc] }; return(u) })})

  
## remove QC failing SNPs from list and recalculate MAX comparison log BF values
for(cc in 1:22) {
  if(length(vs.max.liks[[cc]])<1) { next }
  for(dd in 1:length(vs.max.liks[[cc]])) {
    poz <- 1
    while(names(vs.max.liks[[cc]][[dd]])[poz] %in% qc.cloud.fail) {
      nm <- names(vs.max.liks[[cc]][[dd]])[poz]
      vs.max.liks[[cc]][[dd]][poz] <- NA
      vs.max.liks[[cc]][[dd]] <- vs.max.liks[[cc]][[dd]] - vs.max.liks[[cc]][[dd]][poz+1]
      cat("removed",nm,"from grp",dd,"chr",cc,"\n")
      poz <- poz + 1
    }
  }
}
  
print(load("passAndFailBonfLiks.RData"))
#fail.bonf.liks <- vs.max.liks
#pass.bonf.liks <- vs.max.liks
#save(fail.bonf.liks,pass.bonf.liks,file="passAndFailBonfLiks.RData")
save(vs.max.liks,top.snps,file="finalEquivs.RData")
