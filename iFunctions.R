source("~/github/iChip/firstscriptFunctions.R")

mysumfun <- function(glmr,o.digits=3,p.digits=6,lab=TRUE,ci=FALSE)
{
  co <- summary(glmr)$coefficients
  predz <- rownames(co)[-1]
  label <- paste(summary(glmr)$call)[2]
  p <- co[2:nrow(co),4]; #print(p)
  o.r <- exp(co[2:nrow(co),1]); #print(o.r)
  p <- round(p,p.digits)
  
  # outlist <- list(round(o.r,o.digits),p)
  # names(outlist) <- c("OR","p-value")
  if(ci) {
    co1 <- co[2:nrow(co),1]-(1.96*co[2:nrow(co),2])
    co2 <- co[2:nrow(co),1]+(1.96*co[2:nrow(co),2])
    if(sum(co1)<=sum(co2)) { co.l <- co1; co.h <- co2 } else {  co.h <- co1; co.l <- co2 }
    co.l <- exp(co.l); co.h <- exp(co.h)
    out.fr <- cbind(round(o.r,o.digits),round(co.l,o.digits),round(co.h,o.digits),p)
    colnames(out.fr) <- c("OR","OR-low","OR-hi","p-value")
  } else {
    out.fr <- cbind(round(o.r,o.digits),p)
    colnames(out.fr) <- c("OR","p-value")
  }
  if(lab) { out.fr <- list(out.fr); names(out.fr) <- label }
  return(out.fr)
}

or.conv <- function(X) { 
  X[X<1] <- 1/(X[X<1])
  return(X)
}

make.pheno <- function(X,cases,controls) {
  Pheno <- rep(NA,nrow(X)) # missing (default)
  Pheno[rownames(X) %in% cases] <- 1 # CASE if in the T1d dataset row (id) names
  Pheno[rownames(X) %in% controls] <- 0 # CONTROL if in the Controls dataset row (id) names
  return(Pheno)
}


remove.X <- function(str) {
  bdz <- substr(str,1,1)
  str[bdz=="X"] <- substr(str,2,100000)[bdz=="X"]
  return(str)
}


add.x <- function(str) {
  bdz <- substr(str,1,1)
  numy <-(paste(bdz) %in% paste(c(0:9)))
  str[numy] <- paste("X",str[numy],sep="") 
  return(str)
}


impute.missing <- function (X, bp = 1:ncol(X), strata = NULL, numeric = FALSE, verbose=FALSE, ...) {
  N <- as(X, "numeric")
  if(any(Dim(X)!=Dim(N))) { stop("as numeric lost desired dimensionality") }
  if (!is.null(strata)) {
    strata <- as.factor(strata)
    if (length(levels(strata)) > 10) 
      stop("too many levels in strata\n")
    for (i in levels(strata)) {
      cat("\nstrata", i, "\n")
      wh <- which(strata == i)
      N[wh, ] <- impute.missing(X[wh, , drop = FALSE], 
                                bp, numeric = TRUE, ...)
    }
  }
  else {
    csumm <- col.summary(X)
    use <- csumm[, "Certain.calls"] == 1
    X2 <- X[, use]
    bp <- bp[use]
    imp <- (csumm[, "Call.rate"] < 1 & !is.na(csumm[, "Call.rate"]))[use]
    cat(sum(imp,na.rm=T), "to impute\n")
    mx <- max(which(imp),na.rm=T)
    for (i in which(imp)) {
      loop.tracker(i,mx)
      supres <- capture.output(rule <- snp.imputation(X2[, -i, drop = FALSE], X2[, 
                                                        i, drop = FALSE], bp[-i], bp[i]))
      
      if (is.null(rule@.Data[[1]])) 
        next
      imp <- impute.snps(rules = rule, snps = X2[, rule@.Data[[1]]$snps, drop = FALSE], ...) 
      wh.na <- which(is.na(N[, i]))
      N[wh.na, colnames(X2)[i]] <- imp[wh.na]
    }
    cat("\n")
  }
  if (numeric) {
    return(as.data.frame(N))
  }
  else {
    print(Dim(N))
    print(Dim(X@snps))
    return(new("aSnpMatrix", .Data = new("SnpMatrix", data = (round(N) + 1),
                                           , nrow = nrow(N), ncol = ncol(N), dimnames = dimnames(N)), 
               snps = X@snps, samples = X@samples, phenotype = X@phenotype, 
               alleles = X@alleles))
  }
}


data.frame.to.SnpMatrix <- function(X){
  NN <- as.matrix(X)
  NN <- round(NN)
  SS <- as(NN,"SnpMatrix")
  return(SS)
}


clean.snp.support <- function(X) {
  X[["dbSNP"]][is.na(X$dbSNP)] <- X$SNP[is.na(X$dbSNP)]
  X[["SNP"]] <- clean.snp.ids(X$SNP)
  X[["dbSNP"]] <- clean.snp.ids(X$dbSNP)
  rownames(X) <- clean.snp.ids(rownames(X))
  # specific rsIDs to add
  X$dbSNP[narm(match("seq_VH_2643",X$SNP))] <- "rs75793288"
  X$dbSNP[narm(match("imm_21_42709079",X$SNP))] <- "rs80054410"
  X$dbSNP[narm(match("imm_3_46432416",X$SNP))] <- "rs113010081"
  X$dbSNP[narm(match("imm_3_46488935",X$SNP))] <- "rs115288943"
  X$dbSNP[narm(match("imm_19_10288721",X$SNP))] <- "rs74956615"
  return(X)
}

clean.snp.ids <- function(snpid.list) {
  snpid.list <- gsub("-",".",snpid.list,fixed=T)
  snpid.list <- gsub("1kg_","X1kg_",snpid.list,fixed=T)
  snpid.list <- add.x(snpid.list)
  snpid.list <- gsub(":","_",snpid.list,fixed=T)
  snpid.list <- gsub(",","_",snpid.list,fixed=T)
  snpid.list <- gsub(".","_",snpid.list,fixed=T)
  for (dd in 1:4) { snpid.list <- gsub("XX1kg_","X1kg_",snpid.list,fixed=T) }  # repeat!
  for (nn in 0:9) {
    for (dd in 1:4) { snpid.list <- gsub(paste("XX",nn,"_",sep=""),paste("X",nn,"_",sep=""),snpid.list,fixed=T) }  # repeat!
  }
  return(snpid.list)
}

# fix.rownames <- function(data) {
#   nms <- rownames(data)
#   nms[substr(nms,1,2)=="NA"] <- nms[substr(nms,2,100000)=="NA"] # there are NAxxxx names that should be Axxxx
#   nms[nms=="R58201702_C02"] <- "58201702_C02" # this name has a preceeding R in one place, not in another
#   rownames(data) <- nms.
#   return(data)
# }

SnpMatrix.to.data.frame <- function(SnpMat) {
  cov.data <- as.data.frame(SnpMat)
  for(jj in 1:ncol(cov.data)) { 
    nuxt <- as.numeric(cov.data[,jj])-1
    nuxt[nuxt<0] <- NA
    cov.data[,jj] <- nuxt
    # assign(colnames(cov.data)[jj], nuxt)
  }
  return(cov.data)
}


## this is very specific ##
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



### not sure if these still useful
highlights <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { cat("more than one top hit: ",paste(wh,collapse=","),"\n") }
  next.row <- c(wh[1],top,length(X),length(which(X<bonf)))
  return(next.row) 
}


conditional <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { cat("more than one top hit: ",paste(wh,collapse=","),"\n") }
  next.row <- c(wh,top,length(X),length(which(X<bonf)))
  return(next.row) 
}


## Given a region, add window cM either side and return.
## Everything is HapMap v2, release 22, build 36.

## Options to plot the rates and window exist:
##   do.plot=TRUE -> produce plot
##   add.plot=TRUE -> add lines to existing plot
##   do.lines=TRUE -> use lines to indicate positions of window
## NB, add.plot and do.lines are ignored, unless do.plot=TRUE.

## function from David.  Used because when running this on the queue,
## sometimes need multiple connections to open the damn files.
multitry <- function(expr, times=5, silent=FALSE, message=""){
  warn <- options()$warn
  options(warn=-1)
  for (i in 1:times) {
    res <- try(expr, silent=TRUE)
    if (inherits(res, "try-error")){
      if (i==times) {
        options(show.error.messages = TRUE)
        stop(geterrmessage(), " (", times, " failed attempts)", message)
      }
      next
    }
    break
  }
  options(warn=warn)
  if (!silent && (i>1)) {
    warning(geterrmessage(), " (", times-1, " failed attempts)")
  }
  res
}



recwindow <- function(chr,st,en=st,window=0.1, # cM either side
                      do.plot=FALSE, # if wanted to plot
                      add.plot=FALSE,do.lines=TRUE,...) {
  rates <- read.table(gzfile(sprintf("/dunwich/scratch/chrisw/HapMap/rates_rel22/genetic_map_chr%s_b36.txt.gz",chr)),
                      header=TRUE)
  cm.st <- rates[which.min(abs(rates$position-st)),3]
  cm.en <- rates[which.min(abs(rates$position-en)),3]

  mx <- max(window,1)
  kk <- rates[which.min(abs(rates[,3]-(cm.st-window))) : which.min(abs(rates[,3]-(cm.en+window))),]
  cat("n hapmap snps in window =",nrow(kk),"\n")
  from <- min(kk[,1])
  to <- max(kk[,1])

  if(do.plot) {
    kk <- rates[abs(rates[,3]-cm.st)<mx | abs(rates[,3]-cm.en)<mx,]
    if(add.plot) {
      lines(kk[,1:2])
    } else {
      plot(kk[,1:2],type="l",main=paste("Recombination rates on chr",chr),
           xlab="chromosome position (bp)",ylab="rec rate (cM/Mb)",...)
    }
    if(window>0 & do.lines) {
      abline(v=c(from,to),col="red")
      abline(v=c(st,en),col="blue")
      legend("topleft",lty=c(1,1),col=c("red","blue"),legend=c("window","target"))
    }
  }
  
  cat("window size is\nleft: ",(st-from)/1000,"kb\tright: ",(to-en)/1000,"kb\ttotal: ",(to-from)/1000,"kb\n",sep="")
  return(c(from,to))
}



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

