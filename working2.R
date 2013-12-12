source("~/github/iChip/iFunctions.R")

bonf <- 3.23*(10^-7)
first.deg <- "no" # "yes" "no" "none"
covs <- FALSE # whether to use covariates


plot.one.chr <- function(result,chr=1,new.plot=TRUE,text=TRUE,x.off=0,line.col="blue",pt.col="red",text.pos=1,label.offset=3, text.cex=1) {
  wone <- result
  tab <- sapply(wone[[chr]],"[",13)  # or no number
  #sorter <- order(tab[1,])
  sorter <- order(tab)
  wone[[chr]] <- wone[[chr]][sorter]
  n.in.chr <- length(wone[[chr]])
  xo <- x.off
  l1 <- numeric(n.in.chr)
  if(new.plot) {
    plot(c(0,(n.in.chr+1)),c(0,2.5),xlab="marker(s)",ylab="Odds Ratio",col="white",main=paste("Chr",chr),bty="l",xaxt="n")
    abline(h=1,col="grey",lty="dashed",lwd=1.25)
  }
  for (jj in (1:n.in.chr)) {
    l1[jj] <- y <- tail(wone[[chr]][[jj]][,1],1)
    yl <- tail(wone[[chr]][[jj]][,2],1)
    yh <- tail(wone[[chr]][[jj]][,3],1)
    symbols(jj+xo, y, circles= 0.05,add=T,fg="black",inches=F,bg=pt.col)
    arrows(jj+xo,y,jj+xo,yl,angle=90,length=0.05)
    arrows(jj+xo,y,jj+xo,yh,angle=90,length=0.05)
    if(text) {
      text(jj,yl,labels=names(wone[[chr]])[jj],pos=text.pos,cex=text.cex,offset=label.offset)
    }
  }
  lines((1:n.in.chr)+xo,l1,col=line.col)
  return()
}


highlights <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { cat("more than one top hit: ",paste(wh,collapse=","),"\n") }
  next.row <- c(wh,top,length(X),length(which(X<bonf)))
  return(next.row) 
}


conditional <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { cat("more than one top hit: ",paste(wh,collapse=","),"\n") }
  next.row <- c(wh,top,length(X),length(which(X<bonf)))
  return(next.row) 
}


setwd("/chiswick/data/ncooper/iChipData")
library(reader)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")

id.switcher <- reader("~/Documents/necessaryfilesICHIP/alt.id.lookup.txt")
first.degree <- reader("/ipswich/data/Immunochip/support/casecontrol/subject-lookup-2013-11-08.tab")

sib.only <- (c("no","yes")[(1+as.numeric(first.degree$t1d_sib=="yes" & first.degree$t1d_mother=="no" & first.degree$t1d_father=="no"))])
sib.incl <- first.degree$t1d_sib

fsi <- match(first.degree$X.subjectid,id.switcher$subjectid)
new.id <- id.switcher$sampleid[fsi]
if(any(is.na(new.id))) { cat(length(which(is.na(new.id))),"/",length(new.id)," IDs were not found for first degree status\n")}
FD <- first.degree$t1d_first_degree[!is.na(new.id)]
FDSO <- sib.only[!is.na(new.id)]
FDSI <- sib.incl[!is.na(new.id)]
names(FD) <- names(FDSO) <- names(FDSI) <- narm(new.id)
cat("summary of first degree status:\n")
print(table(FD))
cat("\n")

if(F) {
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
  poz <- as.numeric(table1a[,3])
}
  
  # iChip regions
  if(F) {
    print(load("/chiswick/data/ncooper/iChipData/ichip.regions.RData"))
  }
  topsnplist <- reader(fn="/chiswick/data/ncooper/iChipData/topsnplist.txt")
  gg <- grep("rs689",topsnplist)
  if(length(gg)>0) { topsnplist <- topsnplist[-gg] } # need to remove this SNP as not actually on iChip, just Taqman
  
  ofn <- "/chiswick/data/ncooper/iChipData/temp.ichip-data.RData"

#  cyto <- get.cyto(); cyto[["gene"]] <- rownames(cyto)
  
  chr.dat <- cat.path(fn="temp.ichip-data",suf=paste(1:22),ext="RData")
  chr.dat <- as.list(chr.dat)
  if(!exists("ms")) { ms <- list.rowsummary(chr.dat) }
  chrz <- 1:22
  cr.filt <- ms$Call.rate>=0.953
  hz.filt <- ms$Heterozygosity>=0.19 & ms$Heterozygosity<=0.235
  sample.filt <- cr.filt & hz.filt
  excl.ids <- rownames(get.SnpMatrix.in.file(chr.dat[[22]]))[!sample.filt]
  all.results <- vector("list",22)
  for(next.chr in chrz) {
    Header(paste("Chromosome",next.chr))
    chr <- next.chr
    print(load(cat.path(fn=ofn,suf=next.chr,ext="RData")))
    #annotated.snp.support, t1d.data, t1d.support, control.data, control.support
    snp.support <- clean.snp.support(annotated.snp.support)
    Ph <- rep(c(0,1),times=c(nrow(control.data),nrow(t1d.data)))
    myData <- rbind(t1d.data,control.data)
    #myData <- fix.rownames(myData)
    if(covs) {
      nms <- c(rownames(control.support),rownames(t1d.support) )
      bdz <- substr(nms,1,1)
      nms[bdz=="X"] <- substr(nms,2,100000)[bdz=="X"]
      region13 <- c(control.support$b58cregion,t1d.support$b58cregion) 
      the.sex <- c(control.support$sex,t1d.support$sex)
      cov.dat <- data.frame(sex=the.sex,region=region13)
      rownames(cov.dat) <- nms
     ## cov.dat <- fix.rownames(cov.dat)
      cat(" using covariates for region and sex:\n")
      prv(cov.dat)
    }
    snpid.list <- topsnplist[topsnplist %in% snp.support$dbSNP]
    snpic.list <- snp.support$SNP[match(snpid.list,snp.support$dbSNP)]
    if(any(is.na(snpic.list))) { warning("missing dbSNP entries in list:",snpid.list[is.na(snpic.list)])}
    snpid.list <- snpid.list[!is.na(snpic.list)]
    snpic.list <- narm(snpic.list)
    FD <- FDSI; #FDSO
    t1y <- rownames(t1d.data)[rownames(t1d.data) %in% names(FD[FD=="yes"])] # cases with first degree relative
    t1n <- rownames(t1d.data)[rownames(t1d.data) %in% names(FD[FD=="no"])] # cases without first degree relative
    snp.filt <- narm(match(snpic.list,colnames(myData)))
    if(first.deg!="none") {
      if(first.deg=="yes") {
        smp.filt <- which((!rownames(myData) %in% excl.ids) & (!rownames(myData) %in% t1n))
      } else {
        smp.filt <- which((!rownames(myData) %in% excl.ids) & (!rownames(myData) %in% t1y))
      }
    } else {
      smp.filt <- which(!rownames(myData) %in% excl.ids)
    }
    myDataFilt <- myData[smp.filt,snp.filt]
    myDat <- SnpMatrix.to.data.frame(myDataFilt)
    cat("Analysis dataset:\n")
    print(Dim(myDat))
    result <- vector("list",length(snpid.list))
    
    ## Collect phenotype information and (optional) covariates ##
    Pheno <- rep(NA,nrow(myDataFilt)) # missing (default)
    Pheno[rownames(myDataFilt) %in% rownames(t1d.data)] <- 1 # CASE if in the T1d dataset row (id) names
    Pheno[rownames(myDataFilt) %in% rownames(control.data)] <- 0 # CONTROL if in the Controls dataset row (id) names
    if(any(is.na(Pheno))) { 
      stop("some samples were neither in the T1D nor Control dataset") 
    } else {
      cat("successfully matched",length(which(Pheno==1)),"cases and",length(which(Pheno==0)),"controls\n\n")
      if(first.deg!="none") { cat("Cases with first degree relative:",length(t1y),"without:",length(t1n),"\n") }
    }
    # single.snp.tests(phenotype, stratum, data = sys.parent(), snp.data,
    #                 rules=NULL, subset, snp.subset, uncertain = FALSE, score=FALSE)
    if(covs) { 
      which.covs <- match(rownames(myDataFilt),rownames(cov.dat))
      if(any(is.na(which.covs))) { 
        warning("missing covariates, ignoring covariates") ; covs <- F 
      } else {
        sex <- cov.dat$sex[which.covs]; region <- as.factor(cov.dat$region[which.covs])
      }
    } 
    
    # run GLM analysis for each SNP in the list on the current chromosome
    for(ss in 1:length(snpic.list)) {

      if(!covs) {
        fm <- paste("Pheno ~ ",snpic.list[ss])
      } else {
        fm <- paste("Pheno ~ sex + factor(region) +",snpic.list[ss])
      }
      nxt <- glm(as.formula(fm), family = binomial(logit),data=myDat)
      result[[ss]] <- mysumfun(nxt,p.digits=250,ci=T)[[1]]
    }
    names(result) <- paste(snpic.list,snpid.list,sep="/")
    print(result)
    names(result) <- snpid.list
    all.results[[next.chr]] <- result
  }

 #cov.results <- all.results

 #plain.results =  case/control analysis, univariate
 #cov.results =  case/control analysis with region and sex
 #yes.results =  case/control analysis, univariate, only affect sibs
 #cov.yes.results =  case/control analysis with region and sex, only affect sibs
 #no.results =  case/control analysis, univariate, no affect sibs
 #cov.no.results =  case/control analysis with region and sex, no affect sibs
#save(plain.results,cov.results,yes.results,cov.yes.results,no.results,cov.no.results,file="/chiswick/data/ncooper/iChipData/results.RData")
#save(yes.results,no.results,file="/chiswick/data/ncooper/iChipData/results2.RData")

do.plots <- TRUE
if(do.plots) {  
  load("/chiswick/data/ncooper/iChipData/results2.RData")
  
  pdf("/chiswick/data/ncooper/iChipData/myChrResults2.pdf",width=10,height=7)
  for (chr in 1:22) {
    plot.one.chr(cov.yes.results,chr,text.pos=3,text=F,line.col="red",x.off=-.1)
    plot.one.chr(cov.no.results,chr,new.plot=FALSE,pt.col="darkgreen",line.col="green",x.off=.1)
    legend("topleft",legend=c("T1D: affected family member","T1D: no affected family member"),col=c("red","green"),lwd=2,bty="n")
  }
  dev.off()
  
  cov <- FALSE  #!TRUE
  if(cov) {
    YESYES <- cov.yes.results; NONO <- cov.no.results; indx <- 13
  } else {
    YESYES <- yes.results; NONO <- no.results; indx <- 1
  }
  all.yes <- lapply(YESYES,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,1)) }) } )
  all.no <- lapply(NONO,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,1)) }) } )
  all.yes.seL <- lapply(YESYES,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,2)) }) } )
  all.no.seL <- lapply(NONO,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,2)) }) } )
  all.yes.seH <- lapply(YESYES,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,3)) }) } )
  all.no.seH <- lapply(NONO,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,3)) }) } )
  

  yy <- unlist(all.yes) ; yyl <- unlist(all.yes.seL); yyh <- unlist(all.yes.seH)
  nn <- unlist(all.no); nnl <- unlist(all.no.seL); nnh <- unlist(all.no.seH)
  
  ord <- order(nn)
  YY <- yy[ord]; NN <- nn[ord]; YYL <- yyl[ord]; YYH <- yyh[ord]; NNL <- nnl[ord]; NNH <- nnh[ord]; 
  

if(F) {
  plot(x=(1:length(yy)),y=yy[ord],col="red",type="l",main="71 SNPs across iChip",ylab="odds ratio",xaxt="n",xlab="ordered SNPs")
  abline(h=1,col="grey",lty="dashed",lwd=1.25)
  lines(x=(1:length(yy)),y=nn[ord],col="green")
  legend("topleft",legend=c("T1D: affected family member","T1D: no affected family member"),col=c("red","green"),lwd=2,bty="n")
  
  #pdf("/chiswick/data/ncooper/iChipData/affectedVSnon.pdf")
  plot(YY[YY<1.4]~NN[YY<1.4],xlab="OR: without affected family members",
       ylab="OR: with affected family members",main="Odds Ratios (OR) > 1",bty="l")
  abline(lm(YY[YY<1.4]~0+NN[YY<1.4]),col="red")
  abline(a=-.003,b=1,lty="dashed",col="grey")
  legend("topleft",legend=c("null line","linear fit"),col=c("grey","red"),lty=c("dashed","solid"),bty="n")
  
  plot(YY[YY<=1],NN[YY<=1],xlab="OR: with affected family members",
       ylab="OR: without affected family members",main="Odds Ratios (OR) < 1",bty="l")
  abline(lm(YY[YY<=1]~NN[YY<=1]),col="red")
  abline(a=0,b=1,lty="dashed",col="grey")
  legend("topleft",legend=c("null line","linear fit"),col=c("grey","red"),lty=c("dashed","solid"),bty="n")
  
  plot(YY~NN,xlab="OR: without affected family members",
       ylab="OR: with affected family members",main="All Odds Ratios (OR)",bty="l")
  abline(lm(YY~NN),col="red")
  abline(a=0,b=1,lty="dashed",col="grey")
  legend("topleft",legend=c("null line","linear fit"),col=c("grey","red"),lty=c("dashed","solid"),bty="n")
}
  ####### CALULATE VERSION WHERE ORs>1 with Error Bars ########
  pdf("newnick.pdf")
  YYY <- YY; NNN <- NN; YYYL <- YYL; NNNL <- NNL ; YYYH <- YYH; NNNH <- NNH
  YYY[YY<1] <- 1/YYY[YY<1]; YYYH[YYH<1] <- 1/YYYH[YYH<1] ; YYYL[YYL<1] <- 1/YYYL[YYL<1]
  NNN[NN<1] <- 1/NNN[NN<1]; NNNH[NNH<1] <- 1/NNNH[NNH<1] ; NNNL[NNL<1] <- 1/NNNL[NNL<1]
  plot(YYY~NNN,xlab="OR: without affected family members",
       ylab="OR: with affected family members",main="All Odds Ratios (OR)",bty="l",xlim=c(1,2.25),ylim=c(1,2.6))
  arrows(NNN,YYYL,NNN,YYYH,code = 3,angle=90, length=0.05,col="grey")
  arrows(NNNL,YYY,NNNH,YYY,code = 3,angle=90, length=0.05,col="lightgrey")
  points(YYY~NNN)
  
  plot(YYY[NNN<1.4]~NNN[NNN<1.4],xlab="OR: without affected family members",
       ylab="OR: with affected family members",main="Odds Ratios (OR) > 1",bty="l",xlim=c(1,1.5),ylim=c(1,1.5))
  abline(lm(YYY[NNN<1.4]~0+NNN[NNN<1.4]),col="red")
  abline(a=.016,b=1,lty="dashed",col="grey")
  legend("topleft",legend=c("null line","linear fit"),col=c("grey","red"),lty=c("dashed","solid"),bty="n")
  
  
  
  YYY <- YYY-1; NNN <- NNN-1
  lll <- lm(YYY~0 + NNN); 
  abline(a=(1-lll$coefficients),b=lll$coefficients,col="red")
  cat("Fit line SE:",round(summary(lll)$coefficients[2],5),"\n")
#  abline(lll,col="orange")
  abline(a=0,b=1,lty="dashed",col="blue")
  regr <- paste("b",round(lll$coefficients,3),sep="=")
  legend("topleft",legend=c("null line",paste("linear fit (",regr,")",sep=""),"95% C.I"),col=c("blue","red","grey"),lty=c("dashed",rep("solid",2)),bty="n")
  
  ###################
  
  dev.off()
  wilcox.test(YYY, NNN, paired=TRUE) 
  
  wilcox.test(YY, NN, paired=TRUE) 
  wilcox.test(YY[YY<=1], NN[YY<=1], paired=TRUE) 
  wilcox.test(YY[YY>1], NN[YY>1], paired=TRUE) 
  
}
  #table1:
  # 1.Chr 2.causal  3.best  4.r2  5.all.M>m  6.MAF  7.OR  8.P  9.prev  10.All  11.OR  12.P  13.Ref
  # 1: band - function lookup from location (chr=3, pos =4)
  # 2: gene names - function lookup from location (chr=3, pos =4)
  # 3: rsid - snpid in col2, some of which will need rsid lookup
  # 4: r^2: could look in dataset at col2, versus c#9
  # 5: col 5 > col 6 ?
  # 6: col 7
  # 7: col 10
  # 8: col 12
  # 9:  prev snp?
  # 10: prev snp alleles ?
  # 11: prev snp OR ?
  # 12: prev snp P ?
  # 13: prev snp ref ?
