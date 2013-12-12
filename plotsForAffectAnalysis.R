source("~/github/iChip/iFunctions.R")
work.dir <- "/chiswick/data/ncooper/iChipData/"
setwd(work.dir)
library(reader)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")

only.sibs <- "any"  #"only"  "any" "no"
name.for.pdf <- "affectedVsNot"  # name for main plots file
name.for.saving.data <- "results" # name for binary results file
chr.res.fn <- "myChrResults" # name for chromosome results file
do.plots <- TRUE   # produce plots of ORs for affected vs non-affected
plot.cov <- FALSE  # whether to plot the covariate version
                   # i.e, (cov.no.results/cov.yes.results instead of no.results, yes.results)
regen.data <- FALSE # regenerate the data from scratch or load from RData file
x.cutoff <- 1.4

## build file names
suffix <- c("","Cov")[1+as.numeric(plot.cov)] # covariation [or not] suffix
if(only.sibs!="no") { suffix <- paste(suffix, only.sibs,sep="_") }
fnm <- cat.path(work.dir,name.for.saving.data,suf=suffix,ext="RData")

if(regen.data) {
  ## GENERATE DATA FOR PLOTS FROM SCRATCH ##
  # must run separately for yes/no (affected/unaffected)
  covs <- plot.cov # whether to use covariates
  first.deg <- "no" 
  source('~/github/iChip/affectedAnalysis.R')
  if(covs)  { cov.no.results <- all.results } else { no.results <- all.results }
  first.deg <- "yes" 
  source('~/github/iChip/affectedAnalysis.R')
  if(covs)  { cov.yes.results <- all.results } else { yes.results <- all.results }
  if(covs)  {
    save(cov.yes.results,cov.no.results,file=fnm)
  } else {
    save(yes.results,no.results,file=fnm)
  }
} else {
  ## ALTERNATIVELY, LOAD SAVED DATA FOR PLOTS ##
  load(fnm)
}



if(do.plots) {  
  
  if(plot.cov) {
    YESYES <- cov.yes.results; NONO <- cov.no.results; indx <- 13
  } else {
    YESYES <- yes.results; NONO <- no.results; indx <- 1
  }
  
  out.file <- cat.path(work.dir,chr.res.fn,suf=suffix,ext="pdf")
  pdf(out.file,width=10,height=7)
  for (chr in 1:22) {
    plot.one.chr(YESYES,chr,text.pos=3,text=F,line.col="red",x.off=-.1)
    plot.one.chr(NONO,chr,new.plot=FALSE,pt.col="darkgreen",line.col="green",x.off=.1)
    legend("topleft",legend=c("T1D: affected family member","T1D: no affected family member"),col=c("red","green"),lwd=2,bty="n")
  }
  dev.off()
  cat("wrote plot to:",out.file,"\n")
  
  # extract OR, and lower/upper CIs from list, to vectors
  all.yes <- lapply(YESYES,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,1)) }) } )
  all.no <- lapply(NONO,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,1)) }) } )
  all.yes.seL <- lapply(YESYES,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,2)) }) } )
  all.no.seL <- lapply(NONO,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,2)) }) } )
  all.yes.seH <- lapply(YESYES,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,3)) }) } )
  all.no.seH <- lapply(NONO,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,3)) }) } )
  yy <- unlist(all.yes) ; yyl <- unlist(all.yes.seL); yyh <- unlist(all.yes.seH)
  nn <- unlist(all.no); nnl <- unlist(all.no.seL); nnh <- unlist(all.no.seH)
  
  ## order the all of the vectors by the unaffected odds ratios ##
  ord <- order(nn)
  YY <- yy[ord]; NN <- nn[ord]; YYL <- yyl[ord]; YYH <- yyh[ord]; NNL <- nnl[ord]; NNH <- nnh[ord]; 
  
  if(F) {
    ## what are these plots?
    plot(x=(1:length(yy)),y=yy[ord],col="red",type="l",main="71 SNPs across iChip",ylab="odds ratio",xaxt="n",xlab="ordered SNPs")
    abline(h=1,col="grey",lty="dashed",lwd=1.25)
    lines(x=(1:length(yy)),y=nn[ord],col="green")
    legend("topleft",legend=c("T1D: affected family member","T1D: no affected family member"),col=c("red","green"),lwd=2,bty="n")
    
    plot(YY~NN,xlab="OR: without affected family members",
         ylab="OR: with affected family members",main="All Odds Ratios (OR)",bty="l")
    abline(lm(YY~NN),col="red")
    abline(a=0,b=1,lty="dashed",col="grey")
    legend("topleft",legend=c("null line","linear fit"),col=c("grey","red"),lty=c("dashed","solid"),bty="n")
  }
  
  
  #### CALULATE Data WHERE ORs>1 with Error Bars #####
  
  YYY <- YY; NNN <- NN; YYYL <- YYL; NNNL <- NNL ; YYYH <- YYH; NNNH <- NNH
  YYY[YY<1] <- 1/YYY[YY<1]; YYYH[YYH<1] <- 1/YYYH[YYH<1] ; YYYL[YYL<1] <- 1/YYYL[YYL<1]
  NNN[NN<1] <- 1/NNN[NN<1]; NNNH[NNH<1] <- 1/NNNH[NNH<1] ; NNNL[NNL<1] <- 1/NNNL[NNL<1]
  
  
  ## The plot file ##
  out.file <- cat.path(work.dir,name.for.pdf,suf=suffix,ext="pdf")
  pdf(out.file)
  ###################
  
  ## for all odds ratios ##
  plot(YYY~NNN,xlab="OR: without affected family members",
       ylab="OR: with affected family members",main="All Odds Ratios (OR)",bty="l",xlim=c(1,2.25),ylim=c(1,2.6))
  arrows(NNN,YYYL,NNN,YYYH,code = 3,angle=90, length=0.05,col="grey")
  arrows(NNNL,YYY,NNNH,YYY,code = 3,angle=90, length=0.05,col="lightgrey")
  points(YYY~NNN)
  
  yyy <- YYY-1; nnn <- NNN-1
  lll <- lm(yyy~0 + nnn); 
  abline(a=(1-lll$coefficients),b=lll$coefficients,col="red")
  cat("Fit line SE:",round(summary(lll)$coefficients[2],5),"\n")
  #  abline(lll,col="orange")
  abline(a=0,b=1,lty="dashed",col="blue")
  regr <- paste("b",round(lll$coefficients,3),sep="=")
  legend("topleft",legend=c("null line",paste("linear fit (",regr,")",sep=""),"95% C.I"),col=c("blue","red","grey"),lty=c("dashed",rep("solid",2)),bty="n")
  
  
  ## exclude larger odds ratios ##
  plot(YYY[NNN<x.cutoff]~NNN[NNN<x.cutoff],xlab="OR: without affected family members",
       ylab="OR: with affected family members",main=paste("Odds Ratios (OR) <",x.cutoff),bty="l",xlim=c(1,x.cutoff),ylim=c(1,x.cutoff))

  lll <- lm(yyy[NNN<x.cutoff]~0 + nnn[NNN<x.cutoff]); 
  regr <- paste("b",round(lll$coefficients,3),sep="=")
  abline(a=(1-lll$coefficients),b=lll$coefficients,col="red")
  abline(a=0,b=1,lty="dashed",col="grey")
  legend("topleft",legend=c("null line",paste("linear fit (",regr,")",sep="")),col=c("grey","red"),lty=c("dashed","solid"),bty="n")
  
  ###################
  
  dev.off()
  cat("wrote plots to:",out.file,"\n")
  print(wilcox.test(YYY, NNN, paired=TRUE) )
  print(wilcox.test(YY, NN, paired=TRUE) )
}

