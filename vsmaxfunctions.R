
if(F) {  
  nms.for.prime <- ic.to.rs(rownames(spec)[spec$MULTI==0])
  nms.for.cond <- ic.to.rs(rownames(spec)[spec$MULTI==1])


vs.max.liks[[2]][[4]] <- to.move.to.prime[[1]]
vs.max.liks[[2]] <- vs.max.liks[[2]][c(1,4,2,3)]
vs.max.liks$chr10$rs61839660 <- to.move.to.prime[[2]]
vs.max.liks$chr11$rs689 <- to.move.to.prime[[3]]
vs.max.liks[[16]] <- vs.max.liks[[16]][c(1,2,4)]
vs.max.liks$chr16[[2]] <- to.move.to.prime[[4]]-to.move.to.prime[[4]][3]
vs.max.liks$chr17[[1]] <- to.move.to.prime[[5]]
vs.max.liks$chr18[[1]] <- to.move.to.prime[[6]]
vs.max.liks$chr19[[1]] <- to.move.to.prime[[7]]
}


## just for me



add.names.to.big.vml <- function(vml,tsl=NULL,ret.nms=FALSE,...) {
  if(!detect.vml.big(vml)) { warning("function only for big vml object"); return(vml) }
  if(is.null(tsl)) {
    spec <- read.spec.file(...); tsl <- gsub("b","",ic.to.rs(rownames(spec)))
  } 
  grp <- list()
  for (cc in 1:length(vml)) {
    next.cand <- tsl[Chr(tsl)==cc]
    ll <- length(vml[[cc]])
    #prv(ll)
    if(ll>0) {
      nn <- character(ll)
      for(dd in 1:ll) {
        nlist <- gsub("b","",names(vml[[cc]][[dd]]))
        ww <- which(nlist %in% tsl)[1]
        if(is.na(ww)) { 
          warn.txt <- "could not find any of the tsl (top snp list) SNPs in current row"
          warning(warn.txt)
          cat(warn.txt)
          cat(" Chr ",cc," Grp ",dd," so have inserted top hit instead: ",nlist[1],"\n")
          ww <- 1
        }
        nn[dd] <- nlist[ww]  
        #prv(nn,nlist,ww)
      }
      grp[[cc]] <- nn
      names(vml[[cc]]) <- grp[[cc]] 
    }
  }
  if(ret.nms) { return(grp) } else { return(vml) }
}


results.to.bf <- function(X,add.names=TRUE,...) {
  # ... can be name of spec file and/or topsnplist (tsl)
  vml <- apply.vs(X,do.bic.max,dif=3) 
  ## add names
  if(detect.vml.big(X)) { 
    if(add.names) { 
      names(vml) <- paste("chr",1:22,sep="") 
      vml <- add.names.to.big.vml(vml,ret.nms=FALSE,...) 
    } 
  } else {
    condlist <- get.condlist(...)
    condon <- sapply(condlist,function(x) paste(x,collapse="_"))
    if(length(condon)==length(vml)) {  
      if(add.names) { names(vml) <- condon }
    } else { 
      warning("length of names in conditional file didn't match length of vml object")
    }
  }
  return(vml)
}


  

read.spec.file <- function(st.fn="/chiswick/data/ncooper/iChipData/spectable1bf.csv") {
  if(!file.exists(st.fn)) { stop("spec file did not exist at: ",st.fn) }
  spec <- reader(st.fn)
  needcols <- c("ROW","MULTI")
  if(!all(needcols %in% colnames(spec))) {
    stop("invalid spec file, must contain columns ", paste(needcols,collapse=","))
  }
  return(spec)
}

get.pure.conds <- function(vml,...) {
  ct <- read.spec.file(...)
  to.keep.as.cond <- vml[which(ct$COND!=0)]
  names(to.keep.as.cond) <- ct$TABLE1[which(ct$COND!=0)]
  return(to.keep.as.cond)
}

get.index.conds <- function(vml,...) {
  ct <- read.spec.file(...)
  to.move.to.prime <- vml[which(ct$PRIME==1)]
  names(to.move.to.prime) <- ct$TABLE1[which(ct$PRIME==1)]
  return(to.move.to.prime)
}



read.cond.tests <- function(ct.fn="conditionalTests.csv") {
  if(!file.exists(ct.fn)) { stop("conditional tests file did not exist at: ",ct.fn) }
  ct <- reader(ct.fn,stringsAsFactors=F)
  needcols <- c("SNP1","SNP2","CHR","GRP","PRIME","COND","TABLE1")
  if(!all(needcols %in% colnames(ct))) {
    stop("invalid conditional tests file, must contain columns ", paste(needcols,collapse=","))
  }
  ct[,1] <- ic.to.rs(ct[,1]); ct[,2] <- ic.to.rs(ct[,2])
  return(ct)
}



get.condlist <- function(...) {
  ct <- read.cond.tests(...)
  condlist <- apply(ct,1,function(X) { paste(c(narm(X[1:2]))) })
  return(condlist)
}




### for all ##


apply.vs <- function(X,FUN,...) {
  if(detect.vml.big(X)) {
    out <- lapply(X,function(Y) { lapply(Y,FUN,...) } )
  } else {
    out <- lapply(X,FUN,...)
  } 
  return(out)
}


detect.vml.big <- function(X) {
  if(is.list(X)) {
   return(any(sapply(lapply(X,is),"[",1)=="list"))
  } else {
   stop("invalid vs.max.liks object")
  }
}



set.list.names.to.zero <- function(X) {
  if(detect.vml.big(X)) {
    for(cc in 1:22) {
      if(length(X[[cc]])<1) { next }
      for(dd in 1:length(X[[cc]])) {
          nm <- names(X[[cc]])[[dd]]
          mxmx <- X[[cc]][[dd]][nm]
          X[[cc]][[dd]] <- X[[cc]][[dd]] - mxmx
      }
    }
  } else {
    for(dd in 1:length(X)) {
          nm <- names(X)[dd]
          mxmx <- X[[dd]][nm]
          X[[dd]] <- X[[dd]] - mxmx
    }
  }
  return(X)
}


remove.dups.vml <- function(X) {
  FUN <- function(X) { X[!duplicated(names(X))] }
  X <- apply.vs(X,FUN)
  return(X)
}


rs.to.ic.vml <- function(X) {
  FUN <- function(X) { 
    names(X) <- rs.to.ic(names(X)) 
    if(any(is.na(names(X)))) { warning("NAs produced") } 
    return(X) 
  }
  X <- apply.vs(X,FUN)
  return(X)
}


ic.to.rs.vml <- function(X,dups=TRUE,bs=FALSE) {
  FUN <- function(X) { 
    names(X) <- ic.to.rs(names(X)) 
    if(any(is.na(names(X)))) { warning("NAs produced") } 
    return(X) 
  }
  FUNb <- function(X) { names(X) <- gsub("b","",names(X)); return(X) }
  X <- apply.vs(X,FUN)
  if(bs) { X <- apply.vs(X,FUNb) }
  if(dups) { X <- remove.dups.vml(X) }
  return(X)
}



summarise.bf <- function(X,thresh=c(-3,-5.3,-10,-20,-30,-100)) {
	# make a summary table of bayes factors in SNP lists
	nr <- length(unlist(apply.vs(X,function(Y,thr) { length(Y[Y>thr]) },thr=thresh[1])))	
	rw <- vector("list",length(thresh)) #matrix(nrow=nr,ncol=length(thresh))

	for (cc in 1:length(thresh)) {
	  rw[[cc]] <- unlist(apply.vs(X,function(Y,thr) { length(Y[Y>thr]) },thr=thresh[cc]))
	}
        out <- do.call("cbind",rw)
        colnames(out) <- paste("BF<",abs(thresh),sep=" ")
	return(out)
}




flatten.vml <- function(X,max=F) {
 if(!is.list(X)) { return(X) }
 if(detect.vml.big(X)) {
   Y <- do.call(c,X); 
   if(max) { return(flatten.vml(Y)) } else { return(Y) } 
 } else { 
   names(X) <- NULL; return(unlist(X)) 
 }
}




add.entry.vml <- function(vml,new,where=NULL) {
  if(!is.null(where)) { if(length(where)==1) { where <- rep(where,times=length(new)) } }
  if(detect.vml.big(vml)) {
    # big one
    chrz <- as.numeric(Chr(names(new)))
    if(length(chrz)==0) { stop("length 0 for list of snp ids") }
    if(any(!is.numeric(chrz))) { stop("invalid new snp ids") }
    for(cc in 1:length(chrz)) {
     cat("placing new SNPs into chromosome",chrz[cc],"\n")
     next.rep <- lilbit(vml=vml[[chrz[cc]]],new=new[cc],where={if(is.null(where)) { NULL } else { where[cc] }})
     if(is.list(next.rep)) { vml[[chrz[cc]]] <- next.rep } else { stop("Error, a non-list returned by subfunction lilbit") }
    }
    return(vml)
  } else {
    # little one
    out <- lilbit(vml=vml,new=new,where=where)
    if(is.list(out)) { vml <- out } else { stop("Error, a non-list was returned by subfunction lilbit") }
    return(vml)
  }
}


  ## in CHR or little one
  lilbit <- function(vml,new,where=NULL) { 
    if(is.null(names(new))) { stop("new must be named") }
    if(!is.null(where)) { 
      if(is.numeric(where)) { 
        if(!all(where %in% 1:length(vml))) { stop("bad pos #") } 
      } else {
        full.match <- match(where,names(vml))
        vb.match <- match(where,paste(names(vml),"b",sep="")) 
        nb.match <- match(paste(where,"b",sep=""),names(vml))
        new.where <- full.match
        new.where[is.na(new.where)] <- vb.match[is.na(new.where)]
        new.where[is.na(new.where)] <- nb.match[is.na(new.where)]
        if(any(is.na(new.where))) { prv(new.where,nb.match,vb.match); stop("found non-matching names",where[is.na(new.where)]) }
        where <- new.where
      }
    } else {
      ## autodetect 
      chrz <- Chr(names(sapply(vml,"[[",1)))
      ranges <- lapply(names(vml),function(x) { pp <- range(Pos(x)) })
      chrz.new <- Chr(names(new)); pos.new <- Pos(names(new))
      if(!all(chrz.new %in% chrz)) { stop("outside chromosome range") }
      for(dd in 1:length(chrz.new)) { 
        pozez <- match(chrz.new[dd],chrz)
        if(length(pozez)>1) {
          rg <- ranges[pozez] 
          lf <- sapply(rg,"[",1)
          rt <- sapply(rg,"[",2)
          sub.pos <- which(lf <= pos.new[dd] & rt >= pos.new[dd])
          where[dd] <- pozez[sub.pos]
        } else {
          where[dd] <- pozez
        }
      }
    }    
    cat("placing ",names(new)," at list entry ",where,"\n")
    for(cc in 1:length(new)) {
      vml[[where[cc]]] <- insert.one(vml[[where[cc]]],new[cc])
    }
    return(vml)
  }

insert.one <- function(llist,new) {
  if(length(new)>1) { stop("new must be length 1") }
  nm <- names(llist)
  vl <- as.numeric(llist)
  if(is.character(names(new))) { nnm <- names(new) } else { nnm <- paste("V",1:length(new),sep="") }
  if(any(is.na(nnm))) { nnm[is.na(nnm)] <- paste("Unkwn",1:length( nnm[is.na(nnm)]),sep="") }
  nvl <- as.numeric(new)
  if(all(vl<=nvl)) {
      #first
      onm <- c(nnm,nm)
      ovl <- c(nvl,vl)
  } else {
      if(all(vl>nvl)) {
        #last
        onm <- c(nm,nnm)
        ovl <- c(vl,nvl)
      } else {
        #normal
        onm <- c(nm[vl>nvl],nnm,nm[vl<=nvl])
        ovl <- c(vl[vl>nvl],nvl,vl[vl<=nvl])
      }
  }
  names(ovl) <- onm
  return(ovl)
}


