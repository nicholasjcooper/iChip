


## just for me


# get the list of table1 snps
get.tsl <- function(...,chr=NULL,keep.b=TRUE,ind=TRUE,cond=TRUE) {
  spec <- read.spec.file(...)
  ct <- read.spec.file(...)
  indxs <- ct$TABLE1[which(ct$PRIME==1)]
  conds <- ct$TABLE1[which(ct$COND!=0)]
  tsl <- ic.to.rs(rownames(spec))
  if(!is.null(chr)) { 
    chrz <- Chr(tsl)
    ret <- chrz %in% chr
    tsl <- tsl[ret]
  }
  if(!keep.b) { tsl <- gsub("b","",tsl) }
  if(!ind) { tsl <- tsl[!tsl %in% indxs] }
  if(!cond) { tsl <- tsl[!tsl %in% conds] }
  return(tsl)
}

# add snp names for the lists in a vs.max.liks object [2 level]
add.names.to.big.vml <- function(vml,tsl=NULL,ret.nms=FALSE,keep.b=TRUE,...) {
  if(!detect.vml.big(vml)) { warning("function only for big vml object"); return(vml) }
  if(is.null(tsl)) {
    spec <- read.spec.file(...); 
    if(keep.b) { 
      tsl <- ic.to.rs(rownames(spec))
    } else { tsl <- gsub("b","",ic.to.rs(rownames(spec))) }
  } 
  grp <- list()
  for (cc in 1:length(vml)) {
    next.cand <- tsl[Chr(tsl)==cc]
    ll <- length(vml[[cc]])
    #prv(ll)
    if(ll>0) {
      nn <- character(ll)
      for(dd in 1:ll) {
        nlist <- names(vml[[cc]][[dd]])
        if(!keep.b) { nlist <- gsub("b","",nlist) }
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
      if(any(duplicated(nn))) {
        warning("default algorithm returned duplicates for Chr ",cc," names, using alternative")
        indx.snps.in.chr <- ids.by.pos(get.tsl(chr=cc)) # assume these were in genome order
        dupz <- which(pduplicated(nn))
        not.dupz <- which(!pduplicated(nn))
        isic <- indx.snps.in.chr[!indx.snps.in.chr %in% nn[not.dupz]]
        NN <- nn[dupz]
        if(length(NN)==length(isic)) {
          flag <- TRUE
          for(zz in dupz) {
            #cat("looking for",comma(rs.to.ic(isic[zz])),"in",rs.to.ic(names(vml[[cc]][[zz]])),"\n")
            if(!rs.to.ic(isic[zz]) %in% rs.to.ic(names(vml[[cc]][[zz]]))) { flag <- FALSE }
          }
          if(flag) { 
            nn[dupz] <- isic ; cat("replaced",comma(NN),"with",comma(isic),"\n")
          } else { warning("didn't find expected snps in lists, reverting")}
        } else { print(NN); print(isic); warning("didn't find expected number of index snps in chromosome, reverting") }
      }
      grp[[cc]] <- nn
      names(vml[[cc]]) <- grp[[cc]] 
    }
  }
  if(ret.nms) { return(grp) } else { return(vml) }
}


# convert an 'all.results' object from indistinguishable...Analysis.R files into a vs.max.liks BF summary object
results.to.bf <- function(X,add.names=TRUE,...) {
  # ... can be name of spec file and/or topsnplist (tsl)
  ## add names
  if(detect.ar.big(X)) { 
    vml <- lapply(X,function(Y) { lapply(Y,do.bic.max,dif=3) } )
    if(add.names) { 
      names(vml) <- paste("chr",1:22,sep="") 
      vml <- add.names.to.big.vml(vml,ret.nms=FALSE,...) 
    } 
  } else {
    vml <- lapply(X,do.bic.max,dif=3)
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


  
# read a 'spec' file, which basicaly just lists whether each snp in table 1 is
# conditional, has a conditional or hit, or neither
read.spec.file <- function(st.fn="/chiswick/data/ncooper/iChipData/spectable1bf.csv") {
  if(!file.exists(st.fn)) { stop("spec file did not exist at: ",st.fn) }
  spec <- reader(st.fn)
  needcols <- c("ROW","MULTI")
  if(!all(needcols %in% colnames(spec))) {
    stop("invalid spec file, must contain columns ", paste(needcols,collapse=","))
  }
  return(spec)
}

# get the list of conditional snps
get.pure.conds <- function(vml,...) {
  ct <- read.cond.tests(...)
  to.keep.as.cond <- vml[which(ct$COND!=0)]
  names(to.keep.as.cond) <- ct$TABLE1[which(ct$COND!=0)]
  return(to.keep.as.cond)
}

# get the list of index snps with conditional hits
get.index.conds <- function(vml,...) {
  ct <- read.cond.tests(...)
  to.move.to.prime <- vml[which(ct$PRIME==1)]
  names(to.move.to.prime) <- ct$TABLE1[which(ct$PRIME==1)]
  return(to.move.to.prime)
}


# read the file entered to run the indistinguishableCondAnalysis.R script on the right SNPs
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


# wrapper for read.cond.tests to return the list of conditional snp sets
get.condlist <- function(...) {
  ct <- read.cond.tests(...)
  condlist <- apply(ct,1,function(X) { paste(c(narm(X[1:2]))) })
  return(condlist)
}




### for all ##

# apply function for a vs.max.liks object
apply.vs <- function(X,FUN,...) {
  if(detect.vml.big(X)) {
    out <- lapply(X,function(Y) { lapply(Y,FUN,...) } )
  } else {
    out <- lapply(X,FUN,...)
  } 
  return(out)
}

# determine whether it is a 2 level (T) or 1 level (F) vs.max.liks object
detect.vml.big <- function(X) {
  if(is.list(X)) {
   return(any(sapply(lapply(X,is),"[",1)=="list"))
  } else {
   stop("invalid vs.max.liks object")
  }
}


# determine whether it is a 3 level (T) or 2 level (F) all.results object
detect.ar.big <- function(X) {
  if(is.list(X)) {
    return(!all(sapply(X,length)==2))
  } else {
    stop("invalid all.results object")
  }
}

# BF lists are with respect to the index snp, so the BF for this snp should be zero,
#  if not, this script will perform the necessary transformation to make it so, it
# assumes the index snp in each case is defined by the sublist name
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


# remove duplicate names from vs.max.liks sublists
remove.dups.vml <- function(X) {
  FUN <- function(X) { X[!duplicated(names(X))] }
  X <- apply.vs(X,FUN)
  return(X)
}

# convert the ids in a vs.max.liks object to ichip ids
rs.to.ic.vml <- function(X) {
  FUN <- function(X) { 
    names(X) <- rs.to.ic(names(X)) 
    if(any(is.na(names(X)))) { warning("NAs produced") } 
    return(X) 
  }
  X <- apply.vs(X,FUN)
  return(X)
}

# convert the ids in a vs.max.liks object to rs ids
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


## produce a summary table of which BFs in a vs.max.liks object pass different thresholds
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



# flatten a vs.max.liks object either completely, or by one level
flatten.vml <- function(X,max=F) {
 if(!is.list(X)) { return(X) }
 if(detect.vml.big(X)) {
   Y <- do.call(c,X); 
   if(max) { return(flatten.vml(Y)) } else { return(Y) } 
 } else { 
   names(X) <- NULL; return(unlist(X)) 
 }
}



# add a set of 'new' SNPs and BFs into an existing vs.max.liks object 'vml'
# can specify which list name (where) or let the function detect the best
# insertion point based on a chr, pos lookup from the snp id 
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

# workhorse function used to insert new BF/Snp entries into a vs.max.liks object
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


# function to insert a single BF/snp into a vs.max.liks object
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








make.final.lists <- function(ind.fn,icnd.fn,rs.ids=FALSE,keep.b=TRUE,out.i="",out.c=out.i,save=FALSE) {
  ind.cond <- reader(icnd.fn)
  ind.list <- reader(ind.fn)
  vml.c <- results.to.bf(ind.cond)
  vml.i <- results.to.bf(ind.list)
  vml.i <- set.list.names.to.zero(vml.i)
  vml.i <- remove.dups.vml(vml.i)
  #vml.i <- rs.to.ic.vml(vml.i)
  vml.c <- remove.dups.vml(vml.c)
  #vml.c <- rs.to.ic.vml(vml.c)
  to.move.to.indx <- get.index.conds(vml.c)
  vml.c <- get.pure.conds(vml.c)
  ntmti <- names(to.move.to.indx)
  chr.list <- lapply(vml.i,names)
  for (cc in 1:length(ntmti)) {
    nxt.chr <- Chr(ntmti[cc])
    cat("Ready to transfer",ntmti[cc],"into vml.i\n")
    sub.list.ids <- chr.list[[nxt.chr]]
    cat("Chr",nxt.chr," has groups",comma(sub.list.ids),"\n")
    n.to.repl <- which(sub.list.ids==ntmti[cc])
    cat("so now will replace index ",n.to.repl,"with",names(to.move.to.indx)[cc],"\n")
    vml.i[[nxt.chr]][[n.to.repl]] <- to.move.to.indx[[cc]]
  }
  vml.i <- rs.to.ic.vml(vml.i)
  vml.c <- rs.to.ic.vml(vml.c)
  if(rs.ids) {
    vml.i <- ic.to.rs.vml(vml.i,T,!keep.b)
    vml.c <- ic.to.rs.vml(vml.c,T,!keep.b) 
  }
  if(save) {
    if(out.i=="") { out.i <- simple.date(time=F) }
    if(out.c=="") { out.c <- simple.date(time=F) }
    out.i <- cat.path(getwd(),"IndistinguishableList_",suf=out.i,ext="RData")
    out.c <- cat.path(getwd(),"IndistinguishableCondList_",suf=out.c,ext="RData")
    vs.max.liks <- vml.i;  save(vs.max.liks,file=out.i); cat("Saved index BFs to",out.i,"\n")
    vs.max.liks <- vml.c;  save(vs.max.liks,file=out.c); cat("Saved conditional BFs to",out.c,"\n")
  }
  out <- list(index=vml.i,conditional=vml.c)
  return(out)
}

