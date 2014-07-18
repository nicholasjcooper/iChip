

# from iFunctions.R, but only used by me for the ichip analysis, not actually in package or of general use.

## SPECIFIC TO ME ##


#' @examples
#' hwe.fn <- "snpdataout.hwe"
#' HH <- plink.allele.coding(hwe.fn)
#' head(HH)
plink.allele.coding <- function(hwe.fn) {
  hh <- read.table(hwe.fn,header=T)
  FF <- strsplit(paste(hh[[6]]),split="/",fixed=T)
  hh[["F1"]] <- as.numeric(sapply(FF,"[",1))
  hh[["F2"]] <- as.numeric(sapply(FF,"[",2))
  hh[["F3"]] <- as.numeric(sapply(FF,"[",3))
  hh[["Min"]] <- with(hh,((2*F1)+F2)/(2*(F1+F2+F3)))
  hh[["Maj"]] <- with(hh,((2*F3)+F2)/(2*(F1+F2+F3)))
  hh[["flip"]] <- apply(cbind(hh$A1,hh$A2),1,function(x) { all(x!=sort(x)) })
  hh[["RAF"]] <- hh$Min
  hh$RAF[!hh$flip] <- hh$Maj[!hh$flip]
  return(hh)
}


## imputation subscript ## 
# Nick
# used by conditional analysis and indistinguishable analyses scripts, does the imputation,
# saving imputed results as we go so don't need to keep recalculating the same SNPs
# will automatically proceed in the most efficient way possible
# for my ichip analysis, prv.file="allImputed.RData"
#' @param myData SnpMatrix object with missing data that you wish to impute
#' @param imp.file name to give the file where the result will be saved
#' @param smp.filt a list of sample (indexes or labels) to include
#' @param snp.filt a list of SNPs (indexes or labels) to include
#' @param prv.file character, if applicable the name of a previous imputation file, which
#' means you don't need to impute the same SNPs again if you've imputed them before
#' @param by integer, a parameter that is passed to impute.missing() to determine
#' what sized subsets to impute by (smaller subsets makes imputation faster)
#' @param ret.object logical, if TRUE then return the a SnpMatrix with missing values
#' filled in, if false, return the file name this object was written to (same as imp.file)
#' @return See ret.object. If ret.object=TRUE, then returns the same SnpMatrix with missing
#' values replaced.
#' @examples
#' test.mat <- rSnpMatrix(nsnp=5,nsamp=10)
#' ichip.imputation(test.mat, "myImp.RData")
#' ichip.imputation(test.mat, imp.file="myImp.RData",prv.file="myImp.RData")
#' test.mat2 <- rSnpMatrix(nsnp=500,nsamp=10)
#' ichip.imputation(test.mat2, imp.file="myImp.RData",prv.file="myImp.RData")
#' test.mat3 <- cbind(test.mat2[,1:200],rSnpMatrix(nsnp=500,nsamp=10))
#' ichip.imputation(test.mat3, imp.file="myImp.RData",prv.file="myImp.RData")
ichip.imputation <- function(myData, imp.file, smp.filt=NULL, snp.filt=NULL) {
  if(is.null(smp.filt)) { smp.filt <- 1:nrow(myData) }
  if(is.null(snp.filt)) { snp.filt <- 1:ncol(myData) }
  sample.names <- rownames(myData)[smp.filt]
  #print("1");print(length(smp.filt)); print(Dim(myData))
  if(!file.exists(imp.file)) {
    myDataFilt <- myData[smp.filt,snp.filt]
    myDat <- impute.missing(myDataFilt,numeric=T)
    myDat <- randomize.missing(myDat)
    myDatSnp <- data.frame.to.SnpMatrix(myDat)
    save(myDat,myDatSnp,cov.dat,file=imp.file)
    cat("wrote imputed to:",imp.file,"\n")
  } else { 
    cat("loaded imputed data from:",imp.file,"\n")
    print(load(imp.file))
    ### note this row/col check is not 100% foolproof!
    #print("2");print(length(smp.filt)); print(Dim(myDat))
    if(nrow(myDat)!=length(smp.filt)) { 
      if(!all(sample.names %in% rownames(myDat))) {
        use.exist <- FALSE
        stop("missing samples in loaded file, suggest deleting existing imputation and re-run") 
      } else {
        cat("re-arranging samples in file to match smp.filt\n")
        indxz <- match(sample.names,rownames(myDat))
        if(length(which(is.na(indxz)))<1) {  myDat <- myDat[indxz,] } else { stop("rearranging failed") }
        if(nrow(myDat)!=length(smp.filt)) { stop("something went wrong with sample matching") }
        if(any(rownames(myDat)!=sample.names)) { stop("selection went wrong with sample matching") }
      }
    }
    if(ncol(myDat)!=length(snp.filt)) { 
      ## if this has changed a bit, try to ressurrect without recalculating the whole thing
      cat("mismatching number of snps in loaded file\n") 
      if(!exists("bigDat")) { print(load("allImputed.RData")) }
      targs <- colnames(myData)[snp.filt]
      gotem <- narm(match(targs,colnames(bigDat)))
      aintgotem <- targs[!targs %in% colnames(bigDat)]
      if(length(aintgotem)>0) {
        use.exist <- length(gotem)>0
        if(use.exist) {
          cat("combining",length(gotem),"previously imputed with",length(aintgotem),"from scratch\n")
          myDat.part1 <- impute.missing(myData[smp.filt,match(aintgotem,colnames(myData))],numeric=T)
          if(!exists("bigDat")) { print(load("allImputed.RData")) }
          myDat.part2 <- bigDat[,gotem]
          if(any(rownames(myDat.part1)!=rownames(myDat.part2))) { 
            if(!all(rownames(myDat.part1) %in% rownames(myDat.part2))) {
              warning("samples missing from existing dataset") 
              cat("samples were missing from existing dataset, ")
              use.exist <- FALSE
            } else {
              indzx <- match(rownames(myDat.part1),rownames(myDat.part2))
              if(length(which(is.na(indzx)))<1) {  myDat.part2 <- myDat.part2[indzx,] } else { stop("rearranging samps failed") }
            }
          }
        }
        if(use.exist) {
          # ie, if still true [becomes false if samples were missing]
          myDat.cbind <- cbind(myDat.part1,myDat.part2)
          indz <- match(targs,colnames(myDat.cbind))
          if(any(is.na(indz))) { stop("combined data still missing target SNPs") }
          myDat <- myDat.cbind[,indz]
        } else { cat("no existing valid imputation, ") }
        if(!use.exist) {
          ## do all from scratch
          cat("imputing from scratch\n")
          myDataFilt <- myData[smp.filt,snp.filt]
          myDat <- impute.missing(myDataFilt,numeric=T)
        }
      } else { 
        ## have all the snps, but need to prune some
        cat("trimming loaded data to subset needed\n")
        indz <- match(targs,colnames(bigDat))
        if(any(is.na(indz))) { stop("not sure why loaded data is missing target SNPs") }
        if(length(indz)==0) { stop("no target snps in list") }
        myDat <- bigDat[,indz]
      }
      #prv(myDat)
      myDatSnp <- data.frame.to.SnpMatrix(myDat)
      myDat <- impute.missing(myDatSnp,numeric=T)
      myDat <- randomize.missing(myDat)
      myDatSnp <- data.frame.to.SnpMatrix(myDat)
      save(myDat,myDatSnp,cov.dat,file=imp.file)
      cat("wrote imputed to:",imp.file,"\n")
    } else { 
      myDat <- randomize.missing(myDat)
      myDatSnp <- data.frame.to.SnpMatrix(myDat)
      save(myDat,myDatSnp,cov.dat,file=imp.file)
      cat("loaded file is a good match to requested SNP-set\n")
    }
  }
  return(imp.file)
}



# make sure support file naming and symbol convention match internal conventions, see clean.snp.ids
#internal
clean.snp.support <- function(X) {
  if(is(X)[1]!="ChipInfo") { stop("this function only allows X to be a ChipInfo object") }
  mcols(X)[["dbSNP"]][is.na(X$dbSNP)] <- rownames(X)[is.na(mcols(X)[["dbSNP"]])]
  rownames(X) <- clean.snp.ids(rownames(X))
  mcols(X)[["dbSNP"]] <- clean.snp.ids(mcols(X)[["dbSNP"]])
  rownames(X) <- clean.snp.ids(rownames(X))
  return(X)
}


## FIX UVA XCHR name screw up ##
#internal
rmv.uva.dup.rows <- function(table1) {
  print(Dim(table1))
  bad.tab <- table1[which(is.na(id.to.rs(rownames(table1)))),]
  nr <- substr(rownames(bad.tab),1,nchar(rownames(bad.tab))-1)
  table1 <- table1[which(!rownames(table1) %in% paste(nr,"1",sep="")),] # remove duplicate rows
  print(Dim(table1))
  return(table1)
}


#internal
# this finds the full set of rows in a dataframe that share a duplicated column value (coln)
# i.e, differs to duplicated() which only returns 2..n duplicates and ignores the first members of the sets
# appends row numbers to allow to easily keep track of where these were in the original dataframe
dup.tracker <- function(X,col="dbSNP") {  
  keeps <- (X[(which(duplicated(X[[col]],fromLast=T))),])
  kills <- (X[(which(duplicated(X[[col]],fromLast=F))),])
  dupz <- rbind(keeps,kills)
  wkp <- which(duplicated(X[[col]],fromLast=F))
  wkl <- which(duplicated(X[[col]],fromLast=T))
  rn <- c(wkp,wkl)
  dupz[["rn"]] <- rn
  return(dupz[order(dupz[[col]]),])
}


#internal
## read in downloaded annotation from immunobase and return implied support file as dataframe
immunobase.support <- function(src="/chiswick/data/ncooper/iChipData/iChipAnnot.txt") {
  ici <- reader(src)
  ici[["rsid"]] <- rownames(ici)  
  ici[["A1"]] <- substr(ici[,"Alleles"],1,1)
  ici[["A2"]] <- substr(ici[,"Alleles"],3,3)
  info <- ici[,c("rsid","Immunochip.ID","Chromosome","Start","A1","A2")]
  colnames(info) <- c("rsID","iChipID","Chr","Pos","A1","A2")
  if(!anyDuplicated(info[,2])) { rownames(info) <- info[,2] } else { warning("ichip ids were not unique") }
  return(info)
}




#internal OLD
# remove leading X from variable names (e.g, if original name started with a number and changed by make.names)
remove.X <- function(str,char="X") {
  bdz <- substr(str,1,1)
  str[bdz==char] <- substr(str,2,100000)[bdz==char]
  return(str)
}

#internal OLD
# add X to the start of any string which has a digit as the first character
add.x <- function(str) {
  bdz <- substr(str,1,1)
  numy <-(paste(bdz) %in% paste(c(0:9)))
  str[numy] <- paste("X",str[numy],sep="") 
  return(str)
}


# trivial and specific ... internal for some function somewhere?
# flip odds ratios to always be positive (e.g, to change whether with respect to minor/major alleles)
or.pos <- function(X) { 
  x <- X; sel <- !is.na(X)
  X <- X[sel]
  X[X<1] <- 1/(X[X<1])
  x[sel] <- X
  return(x)
}


# flip odds ratios 
# true/false vector of whether the odds ratio should be greater than 1
# GWAS convention is that ORs>1 when cases have more of the minor allele than controls
or.flip <- function(X,greater.than.1=rep(T,length(X))) {
  gt1 <- greater.than.1 # save space
  if(!is.numeric(X)) { stop("X should be numeric (ie, odds ratios)") }
  if(!is.logical(gt1) | (length(gt1)!=length(X) & length(gt1)!=1) ) { 
    stop("'greater.than.1' should be logical (ie, TRUE/FALSE for each member of X)") }
  xx <- X; sel <- !is.na(gt1) & !is.na(X)
  X <- X[sel]; gtone <- gt1[sel]
  gt1s <- X>1
  lt1s <- X<1
  X[gt1s & !gtone] <- 1/X[gt1s & !gtone]
  X[lt1s & gtone] <- 1/X[lt1s & gtone]
  xx[sel] <- X
  return(xx)
}



##INTERNAL
# create a phenotype vector for a dataframe with rownames that are subject ids, where
# cases and controls are text vectors of which IDs are that category; ctrls coded 0, cases 1
# or can also enter X as a vector of all ids, e.g, X = rownames(someDataFrame)
make.pheno <- function(X,cases,controls) {
  if(length(Dim(X))==1) {
    Pheno <- rep(NA,length(X)) # missing (default)
    Pheno[X %in% cases] <- 1 # CASE if in the T1d dataset row (id) names
    Pheno[X %in% controls] <- 0 # CONTROL if in the Controls dataset row (id) names
  } else {
    Pheno <- rep(NA,nrow(X)) # missing (default)
    Pheno[rownames(X) %in% cases] <- 1 # CASE if in the T1d dataset row (id) names
    Pheno[rownames(X) %in% controls] <- 0 # CONTROL if in the Controls dataset row (id) names
  }  
  return(Pheno)
}


# to extract info from 22 separate chromosomes
# pheno <- make.pheno(rawdata,rownames(t1d.data),rownames(control.data)) # assume all are same length
# for (cc in 1:22) {
#   
#   load(paste("temp.ichip-data",cc,".RData",sep=""))
#   rawdata <- rbind(control.data,t1d.data)
#   tt <- fix.OR.directions(tt, OR.col=c(4,6,8), snp.data=rawdata, pheno=pheno, partial=T,verbose=F)
#   loop.tracker(cc,22)
# }
### specific to ichip paper analysis
fix.OR.directions <- function(results, OR.col=1, snp.data=NULL, pheno=NULL, partial=FALSE,
                              case.list=NULL, control.list=NULL, effect.labels=TRUE, verbose=TRUE, alt.return=FALSE) {
  nameslist <- colnames(snp.data)
  if(length(OR.col)>0) {
    if(is.character(OR.col)) {
      if(any(!OR.col %in% colnames(results))) { stop("OR.col must be column names or numbers of 'results'") } else { do.or.col <- T }
    }
    if(is.numeric(OR.col)) {
      if(any(OR.col > ncol(results))) { stop("OR.col was greater than the number of columns in 'results'") } else { do.or.col <- T }
    } 
  } else { effect.labels <- TRUE }
  rr_in_ss <- rownames(results) %in% nameslist
  if(!all(rr_in_ss)) { 
    if(!partial) { warning(length(which(!rr_in_ss))," rownames from results not found in colnames of snp.data") };
    incomplete <- T  
  } else { incomplete <- F }
  if(!any(rr_in_ss)) { stop("no rownames from results were found in colnames of snp.data") }
  if(is.null(pheno) & !is.character(case.list) & is.character(control.list)) {  
    pheno <- make.pheno(snp.data,case.list,control.list)
  }
  if(length(pheno)!=nrow(snp.data)) { 
    stop("pheno (or case.list+control.list) must be ",
         "the same length as the number of rows in snp.data") }
  indx <- match(rownames(results)[rr_in_ss],nameslist)
  majminlist <- majmin(snp.data[,indx])
  reflist <- caseway(snp.data[,indx],pheno)
  if(alt.return) { return(list(majminlist=majminlist,reflist=reflist,
                               nameslist=nameslist,indx=indx,rr_in_ss=rr_in_ss)) }
  # get logical vectors
  l1 <- (majminlist=="minor" & reflist=="CasesRef-")
  l2 <- (majminlist=="minor" & reflist=="CasesRef+")
  l3 <- (majminlist=="major" & reflist=="CasesRef-")
  l4 <- (majminlist=="major" & reflist=="CasesRef+")
  OR.dir <- rep(as.logical(NA),length(indx)); 
  caseshave <- rep("nodiff",length(indx))
  # assign non-default values based on logical vectors above
  OR.dir[l1] <- FALSE; caseshave[l1] <- "MoreMajor"
  OR.dir[l2] <- TRUE;  caseshave[l2] <- "MoreMinor"
  OR.dir[l3] <- TRUE;  caseshave[l3] <- "MoreMinor"
  OR.dir[l4] <- FALSE; caseshave[l4] <- "MoreMajor"
  if(do.or.col) {
    for(cc in 1:length(OR.col)) { 
      results[rr_in_ss,OR.col[cc]] <- or.flip(results[rr_in_ss,OR.col[cc]],OR.dir)
      if(verbose) { cat("inverted any inconsistent odds-ratios in column",OR.col[cc],"\n") }
    }
  }
  if(effect.labels) {
    # if(incomplete) { results[["CasesHave"]][!rr_in_ss] <- NA }
    results[["CasesHave"]][rr_in_ss] <- caseshave
  }
  return(results)
}


# to extract info from 22 separate chromosomes
#  pheno <- make.pheno(rawdata,rownames(t1d.data),rownames(control.data)) # assume all are same length
#  for (cc in 1:22) {
#    load(paste("temp.ichip-data",cc,".RData",sep=""))
#    rawdata <- rbind(control.data,t1d.data)
#    tt2 <- add.allele.to.result(tt2, OR.col=8, snp.data=rawdata, pheno=pheno, partial=T)
#    loop.tracker(cc,22)
#  }
### specific to ichip paper analysis
add.allele.to.result <- function(results, OR.col=1, snp.data=NULL, pheno=NULL, partial=FALSE,
                                 case.list=NULL, control.list=NULL) {
  lll <- fix.OR.directions(results, OR.col=OR.col[1], snp.data=snp.data, pheno=pheno, partial=partial,verbose=FALSE,
                           case.list=case.list, control.list=control.list, effect.labels=FALSE,alt.return=TRUE)
  list.to.env(lll) # gets: majminlist,reflist,nameslist,indx,rr_in_ss
  al.names <- c("allele.A","allele.B")
  if(all(al.names %in% colnames(results))) {
    tt <- results
  } else {
    tt <- cbind(results,AB(rownames(results)))
    colnames(tt)[(ncol(tt)-c(1,0))] <- al.names
  }
  ## T>A   #  cases h ave more T, less A
  ## T is ref because ref is allele.B
  ## if cases ref+ then  allele.B >  allele.A
  ## if cases ref- then  allele.A >  allele.B
  ## else allele.A ~ allele.B
  
  # get logical vectors
  l1 <- (majminlist=="minor")
  l2 <- (majminlist=="major")
  relship <- paste((tt$allele.A[rr_in_ss]),"~",(tt$allele.B[rr_in_ss]))
  relship[l1] <- paste((tt$allele.A[rr_in_ss]),">",(tt$allele.B[rr_in_ss]))[l1]
  relship[l2] <- paste((tt$allele.B[rr_in_ss]),">",(tt$allele.A[rr_in_ss]))[l2]
  if(!"effect" %in% colnames(tt)) { tt[["effect"]] <- NA }
  #prv(reflist,l1,l2,relship,indx,rr_in_ss)
  tt[rr_in_ss,"effect"] <- relship
  return(tt)
}


# specific to me
# calculate and print the overall lambdas reflecting inflation
calculate.overall.lambdas <- function(surround=FALSE) {
  print(load("lambdaTableVars.RData")) #"c2"      "p.m.out" "p.m.reg"
  # c2 is a table of results including armitage trend test 'Y2' column, calc derived using lambda.R
  work.dir <- "/chiswick/data/ncooper/iChipData/"
  print(load("/chiswick/data/ncooper/iChipData/finalMetaTopHitsFEB17.RData"))
  # ^file contains objects from getMetaTable1.R
  # "bonfs.filt"     "non.bonfs.filt" "bonf.snps"      "non.bonf.snps" 
  print(load("/chiswick/data/ncooper/iChipData/compiledTableAllResultsPassingQC.RData"))
  # "tt"
  print(load("pvaluesforqq.RData"))
  # "p.meta" "p.cc"
  near.region <- grep("EXT",tt$gene)
  out.region <- grep("OTHER",tt$gene)
  in.region <- which(!1:nrow(tt) %in% c(near.region,out.region))
  potential.novel.region <- which(1:nrow(tt) %in% c(near.region,out.region))
  p.m.reg <- names(p.meta[names(p.meta) %in% rs.to.id(tt$names[in.region])])
  p.m.out <- names(p.meta[names(p.meta) %in% rs.to.id(tt$names[out.region])])
  p.m.nov <- names(p.meta[names(p.meta) %in% rs.to.id(tt$names[potential.novel.region])])
  
  if(surround) {
    t1dsnps <- calibrate.cond.bonf(topsnplist,cm.window=0.2,bp.ext=0,build=37,qclist="snpsExcluded.txt",ret.snps=T)
  }
  t1region <- tt$gene %in% rownames(bonfs.filt)
  t1region.novel <- (tt$gene %in% rownames(bonfs.filt)) & (rs.to.id(tt$names) %in% rs.to.id(p.m.nov))
  t1regsnps <- rs.to.id(tt$names[t1region])
  t1regsnps.novel <- rs.to.id(tt$names[t1region.novel])
  p.m.out.not.t1 <- p.m.out[!p.m.out %in% t1regsnps.novel]
  cat("there were",length(t1regsnps.novel),"SNPs in novel regions\n")
  if(surround) { cat("raw lambda for the cM regions surrounding t1dhits:",median(c2[t1dsnps,"Y2"],na.rm=T)/.456,"\n") }
  cat("raw lambda for the regions containing t1dhits:",median(c2[t1regsnps,"Y2"],na.rm=T)/.456,"\n")
  lam0 <- median(c2[,"Y2"],na.rm=T)/.456
  cat("raw lambda for all SNPs passing QC:",lam0,"\n")  
  lam1000.0 <- lambda_nm(lam0,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T))
  cat("lambda_1000 for all SNPs passing QC:",lam1000.0,"\n")
  if(surround) { lam1 <- median(c2[-which(rownames(c2) %in% t1dsnps),"Y2"],na.rm=T)/.456 }
  lam2 <- median(c2[-which(rownames(c2) %in% t1regsnps),"Y2"],na.rm=T)/.456
  if(surround) { cat("raw lambda excluding the 0.2 cM regions surrounding t1dhits:",lam1,"\n") }
  cat("raw lambda excluding the dense regions containing t1dhits or other table-1 novel regions:",lam1,"\n")
  if(surround) { lam1000.1 <- lambda_nm(lam1,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T)) }
  lam1000.2 <- lambda_nm(lam2,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T))
  if(surround) { cat("lambda_1000 excluding the cM regions surrounding t1dhits:",lam1000.1,"\n") }
  cat("lambda_1000 excluding the dense regions containing t1dhits or other table-1 novel regions:",lam1000.2,"\n")
  lam3 <- median(c2[p.m.out,"Y2"],na.rm=T)/.456
  lam4 <- median(c2[p.m.reg,"Y2"],na.rm=T)/.456
  lam5 <- median(c2[p.m.out.not.t1,"Y2"],na.rm=T)/.456
  cat("raw lambda for outside ichip dense regions:",lam3,"\n")
  cat("raw lambda for inside ichip dense regions:",lam4,"\n")
  cat("raw lambda for outside ichip dense regions and excluding table1 novel regions:",lam5,"\n")
  lam1000.3 <- lambda_nm(lam3,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T))
  lam1000.4 <- lambda_nm(lam4,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T))
  lam1000.5 <- lambda_nm(lam5,1000,1000,median(c2$R,na.rm=T),median(c2$S,na.rm=T))
  cat("lambda_1000 for outside ichip dense regions:",lam1000.3,"\n")
  cat("lambda_1000 for inside ichip dense regions:",lam1000.4,"\n")
  cat("lambda_1000 for outside ichip dense regions and excluding table1 novel regions:",lam1000.5,"\n")
  
  #L = 1.641179
  #L1000 = 1.074309
  do.the.meta.ones <- function(p.meta,suf="") {
    pdf(cat.path("","chiqq",ext="pdf",suf=suf))
    #prv(p.meta)
    #p.meta <- narm(p.meta)
    txt <- paste("Meta Analysis (ChiSq) [",toheader(paste(gsub("."," ",suf,fixed=T),"regions")),"]",sep="")
    LL <- length(which(!is.na(qchisq(p.meta,1))))
    xx <- qchisq(1-((1:LL)/LL),1); yy <- qchisq(1-sort(p.meta),1)
    cond <- (is.na(xx) | is.na(yy) | !is.finite(yy) | !is.finite(xx))
    xx <- xx[!cond]; yy <- yy[!cond]
    three4 <- function(x) { mean(c(min(x,na.rm=T), rep(max(x,na.rm=T),3))) }
    plot(xx,yy,type="l",main=txt,xlab="expected",ylab="observed",xlim=c(0,20))
    text(three4(xx),three4(yy),paste("slope =",round(coefficients(lm(yy~xx))[2],3)))
    lines(x=qchisq(1-((1:LL)/LL),1),y=qchisq(1-((1:LL)/LL),1),col="red",lty="dotted")
    dev.off()
  }
  
  # QQ for dense versus non dense regions #
  p.meta.reg <- p.meta[names(p.meta) %in% rs.to.id(tt$names[in.region])]
  #p.meta.out <- p.meta[names(p.meta) %in% rs.to.id(tt$names[out.region])]
  p.meta.out <- p.meta[names(p.meta) %in% rs.to.id(p.m.out.not.t1)]
  do.the.meta.ones(p.meta.reg,"inside.dense")
  do.the.meta.ones(p.meta.out,"outside.dense")
  # QQ for t1d versus non dense regions #
  p.meta.reg <- p.meta[names(p.meta) %in% rs.to.id(t1regsnps)]
  p.meta.out <- p.meta[!names(p.meta) %in% rs.to.id(t1regsnps)]
  do.the.meta.ones(p.meta.reg,"inside.t1d")
  do.the.meta.ones(p.meta.out,"outside.t1d")
  
  return(c(lam1000.1,lam1000.2,lam1000.3,lam1000.4,lam1000.5))
}


#send cw: bonferonnis, numbers of snps


# function specific to getMetatable.R script
#  gets a list of equivalent SNPs to the current snp
get.equivs <- function(id,iden.list,na.fail=T) {
  id <- id.to.rs(id) 
  do.one <- function(id,iden.list,na.fail) {
    if(id %in% id.to.rs(unlist(iden.list))) {
      lnum <- which(sapply(iden.list,function(X) { id %in% id.to.rs(X) }))
      idens <- id.to.rs(iden.list[[lnum]])
      return(idens[-which(idens %in% id)])
    } else {
      if(na.fail) { return(NA) } else { return(id) }
    }
  }
  return(sapply(id,do.one,iden.list=iden.list,na.fail=na.fail))
}


# specific
### GLOBAL QC STATS ###
samp.summ <- function(ms,CR=0.953,HZlo=0.19,HZhi=0.235,by.pheno=FALSE) {
  if(by.pheno) {
    if(exists("control.data") & exists("t1d.data")) {
      ph <- make.pheno(ms,rownames(t1d.data),rownames(control.data))
      cat("Sample QC for CASES (phenotype=1)\n")
      sample.filt1 <- samp.summ(ms=ms[ph==1,],CR=CR,HZlo=HZlo,HZhi=HZhi,by.pheno=FALSE)
      cat("Sample QC for CONTROLS (phenotype=0)\n")
      sample.filt0 <- samp.summ(ms=ms[ph==0,],CR=CR,HZlo=HZlo,HZhi=HZhi,by.pheno=FALSE)
      return(c(sample.filt0,sample.filt1))
    }
  }
  nsamp <- nrow(ms)
  cr.filt <- ms$Call.rate>=CR
  hz.filt <- ms$Heterozygosity>=HZlo & ms$Heterozygosity<=HZhi
  sample.filt <- cr.filt & hz.filt
  cat(out.of(length(which(!cr.filt)),nsamp)," samples fail on call rate < ",CR,"\n",sep="")
  cat(out.of(length(which(!hz.filt)),nsamp)," samples fail on ",HZlo,"<Hz<",HZhi,"\n",sep="")
  cat(out.of(length(which(!sample.filt)),nsamp)," samples fail on call rate or HWE\n",sep="")
  return(sample.filt)
}

#specific
snp.summ <- function(MAF=0.005,CR=.99,HWE=3.8905,qc.file="snpqc.RData") {
  if(file.exists(qc.file)) { SNPQC <- reader(qc.file) } else { stop("couldn't find",qc.file) }
  mini.snp.qc <- function(SNPQC,MAF=0.005,CR=.99,HWE=3.8905) {
    maf <- SNPQC$MAF>MAF
    clr <- SNPQC$Call.rate>CR
    hwe <- abs(SNPQC$z.HWE)<HWE
    mono <- SNPQC$MAF==0 | SNPQC$MAF<0.0005
    maf.not.mono <- SNPQC$MAF>=0.0005 & SNPQC$MAF<MAF
    nsnp <- nrow(SNPQC)
    excl.rules <- rownames(SNPQC)[which(!maf | !clr | !hwe)]
    cat(out.of(length(which(!maf)),nsnp)," snps fail on MAF < ",MAF," ... \n",sep="")
    cat("... of which ",length(which(mono))," were monomorphic and ",length(which(maf.not.mono))," were just rare\n",sep="")
    cat(out.of(length(which(!clr & maf)),nsnp)," snps fail on call rate < ",CR,"\n",sep="")
    cat(out.of(length(which(!hwe & clr & maf)),nsnp)," snps fail on HWE < ",HWE," [p=",round(2*(1-pnorm(HWE)),8),"]\n",sep="")
    cat(out.of(length(excl.rules),nsnp)," snps fail on MAF, Callrate or HWE \n",sep="")
    return(excl.rules)
  }
  nsnp <- nrow(SNPQC)
  Header("DIL")
  mono <- SNPQC$MAF==0 | SNPQC$MAF<0.0005
  #excl.rules <- mini.snp.qc(SNPQC[!mono,],MAF=MAF,CR=CR,HWE=HWE)
  excl.rules <- mini.snp.qc(SNPQC,MAF=MAF,CR=CR,HWE=HWE)
  num.excl.rules <- length(excl.rules)
  if(!file.exists("snpsExcluded.txt")) { warning("couldn't find exclusions file") } else {
    which.cut <- readLines(cat.path(getwd(),"snpsExcluded.txt"))
    total.excl <- length((unique(which.cut)))
    cat(out.of(total.excl-num.excl.rules,nsnp),"snps fail for other reasons\n")
    cat(out.of(total.excl,nsnp),"snps failed altogether\n")
    doc.path <- "ImChip_T1D_20130913"
    docs <- cat.path(doc.path,list.files(doc.path))
    if(file.exists(docs[2])) {  
      Header("UVA")
      excl.uva <- rownames(reader(docs[2])) 
      cat(out.of(length(excl.uva),nsnp),"snps failed according to UVA list\n")
      cat(out.of(length(excl.uva %in% which.cut),length(excl.uva),oo=F),"of these also failed according to DIL list\n")
      uva.lookup <- match(rs.to.id(excl.uva),rs.to.id(rownames(SNPQC)))
      print(length(which(is.na(uva.lookup))))
      excl.rules.uva <- mini.snp.qc(SNPQC[uva.lookup,],MAF=.005,CR=.95,HWE=4.89)
      num.excl.rules.uva <- length(excl.rules.uva)
      cat(out.of(length(excl.uva)-num.excl.rules.uva,length(excl.uva)),"UVA failing snps fail for other reasons\n")
    }
  }
  return(excl.rules)
}
####################



# what does this do?
find.overlapping.regions  <- function(ranged) {
  if(is(ranged)[1]==("GRanges")) { ranged <- as(ranged,"RangedData") }
  if(!is(ranged)[1]==("RangedData")) { stop("Need RangedData or GRanges object to proceed")} else {
    ranged <- toGenomeOrder(ranged,strict=T) }  
  rr <- NULL
  chrz <- narm(rownames(chrInfo(ranged)))
  LenC <- length(chrz)
  if(LenC>1) { 
    ret <- vector("list",LenC)
    for (ccc in 1:LenC) { 
      cat("Chr",chr2(ranged[ccc])[1],":\n")
      ret[[ccc]] <- find.overlapping.regions(ranged[ccc]) 
      if(!is.null(ret[[ccc]])) { 
        print(paste("chr",ccc))
        uu <- ranged[ccc][sort(unique(as.numeric(ret[[ccc]]))),]
        print(uu); if(ccc==1) { rr <- uu } else { rr <- rbind(rr,uu) }
      }
    }
    names(ret) <- paste(chrz)[1:LenC]
    isnull <- sapply(ret,function(X) { all(is.na(X)) })
    return(rr) #ret[!isnull])
  }
  ov <- findOverlaps(ranged)
  qh <- queryHits((ov))
  sh <- subjectHits((ov))
  ovs <- which(qh!=sh)
  if(length(ovs)>0) {
    overlaps <- cbind(qh,sh)[ovs,]
    overlaps <- t(apply(overlaps,1,sort))
    #s return(overlaps)
    otxt <- apply(overlaps,1,function(x) { paste(x,collapse=",") })
    overlaps <- overlaps[!duplicated(otxt),]
    otxt <- otxt[!duplicated(otxt)]
    otxt <- paste("overlap found for regions: ",otxt,"\n")
    cat(otxt,"\n",sep="")
    if(length(Dim(overlaps))==1) { dim(overlaps) <- c(1,2) }
    if(!is.null(dim(overlaps))) {
      colnames(overlaps) <- c("region.a","region.b")
    }
    return(overlaps)
  } else {
    cat("no overlaps\n")
    return(NA)
  }
}

### this is used by getMetaTable1.R to get the top snp in each region-group ###
highlights <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { cat("more than one top hit: ",paste(wh,collapse=","),"\n")  }
  # if(length(wh)==2) {  wh <- wh[-1] } ## wh[-1] is hack to adjust for 1p13.2, although don't know why i bothered as it doesn't carry through
  next.row <- c(wh[1],top,length(X),length(which(X<bonf)))
  return(next.row) 
}


multihit.check <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { return(T) } else { return(F) }
}

multihit.return <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  return(wh)
}


## this function gets marginal likelihoods from the BIC, similar to how do.bic.max works on all.results


## my GWAS with SNPstats will not have ORs/betas with the right directions to be consistent 
# with UVA's TDT family analysis. In order to ensure the OR/betas are the right way around,
# flip each if necessary to have log(sign) consistent with the UVA case-control analysis
# e.g. 
#    new.table.fn <- "nick.meta.table2.RData"
#    table.nick <- convert.OR.directions(new.table.fn)
#    save(table.nick,file="nick.meta.table3.RData")
convert.OR.directions <- function(new.table.fn) {
  nick.table <- reader(new.table.fn)
  print(Dim(nick.table))
  rownames(nick.table) <- clean.snp.ids(rownames(nick.table))
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  excl <- reader(docs[2])
  uva.table <- reader(docs[1])
  rownames(uva.table) <- clean.snp.ids(rownames(uva.table))
  uva.table <- uva.table[rownames(nick.table),]
  wwu <- sign(log(uva.table[,"OR_CC"])) 
  names(wwu) <- rownames(uva.table)
  wwn <- sign(log(nick.table[,"OR_CC"])) 
  names(wwn) <- rownames(nick.table)
  wwun <- wwu[match(names(wwn),names(wwu))]
  length(which(is.na(wwun)))
  wwun[is.na(wwun)] <- wwn[is.na(wwun)] 
  print(head(cbind(wwn,wwun)))
  print(length(which(wwn==wwun)))
  print(length(which(wwn!=wwun)))
  print(Dim(nick.table))
  nick.table[["OR_CC_Raw"]] <- nick.table[,"OR_CC"]
  nick.table[,"OR_CC"][wwn!=wwun] <- 1/(nick.table[,"OR_CC"][wwn!=wwun])
  cat("removing columns",paste(colnames(nick.table)[15:17],collapse=","),"\n")
  cat("replacing with:",paste(colnames(meta.me(nick.table)),collapse=","),"\n")
  table.nick <- cbind(nick.table[,-15:-17],(meta.me(nick.table)))
  colnames(table.nick)[17:ncol(table.nick)] <- c("OR_Meta","b_Meta","SE_Meta","Z_Meta","P_Meta")
  return(table.nick)
}

do.bic.marg <- function(X,dif=3) {
  if(!"BIC" %in% names(X)) { return(NA) }
  bic <- X$BIC;
  if(is.list(bic)) { bic <- unlist(bic) }
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


## this is the function that extracts 'BIC' from the all.results list and converts to relative Bayes factors
# against the lowest BIC in the list
do.bic.max <- function(X,dif=3) {
  if(!"BIC" %in% names(X)) { return(NA) }
  bic <- X$BIC;
  if(is.list(bic)) { bic <- unlist(bic) }
  max.bic <- -min(bic/2)
  dif.bic <- (-bic/2)-max.bic
  return(rev(sort(dif.bic)))
}


# extract the BIC information from 'all.results' list which is returned by the indistinguishableAnalysis.R file
suck.bic <- function(X,dif=3) {
  bic <- X$BIC;
}

# not used in the end?
clearly.suck <- function(X,thresh=100) {
  if(!"BIC" %in% names(X)) { return(NA) }
  bic <- X$BIC;
  if(is.list(bic)) { bic <- unlist(bic) }
  max.bic <- -min(bic/2)
  dif.bic <- (-bic/2)-max.bic
  failers <- (rev(sort(dif.bic)))
  failers <- failers[abs(failers)>=thresh]
  return(names(failers))
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




conditional <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { cat("more than one top hit: ",paste(wh,collapse=","),"\n") }
  next.row <- c(wh,top,length(X),length(which(X<bonf)))
  return(next.row) 
}




data.frame.to.ranged2 <- function(dat,ids=NULL,start="start",end="end",width=NULL,
                                  chr="chr",exclude=NULL,build="hg18") 
{
  ## convert any data frame with chr,start,end, or pos data into a RangedData object
  # not case sensitive
  ## abandon longer names as they clash with function names
  st <- paste(start); en <- paste(end); ch <- paste(chr); wd <- paste(width)
  must.use.package(c("genoset","IRanges"),T)
  if(is.matrix(dat)) { dat <- as.data.frame(dat) }
  if(!is.data.frame(dat)) { stop("Error: not a dataframe")}
  key.nms <- c(ids,st,en,ch,wd)
  tries <- 0
  while(!all(key.nms %in% colnames(dat))) { 
    colnames(dat) <- tolower(colnames(dat)); key.nms <- tolower(key.nms)
    st <- tolower(st); en <- tolower(en); ch <- tolower(ch); wd <- tolower(wd)
    if(tries>1) {
      if((tolower(st)=="pos" | tolower(en)=="pos") & !(tolower(st)=="pos" & tolower(en)=="pos")) {
        st <- en <- "pos"
      } else {
        if(tolower(st)=="start" & tolower(en)=="end") { st <- en <- "pos" }
      }
    }
    key.nms <- c(ids,st,en,ch,wd)
    tries <- tries+1
    if(tries > 3) { warning("columns not found") ; break }
  }
  if(!is.null(ids)) { 
    if(anyDuplicated(dat[[ids]])==0) { 
      id <- dat[[ids]] 
    } else { 
      key.nms <- key.nms[-match(ids,key.nms)] # allow non-unique ids as regular
      ids <- NULL
      warning("id must be unique to form rownames, will insert as a separate column") 
    }
  }
  if(is.null(ids)) { 
    if(!is.null(rownames(dat)) & all(rownames(dat)!=paste(1:nrow(dat)))) { 
      id <- rownames(dat)
    } else { 
      id <- paste(1:nrow(dat)) 
    }
  }
  if(length(ch)>0) { ch1 <- gsub("X","chrX",gsub("chr","",dat[[ch]],ignore.case=T)) } else { ch1 <- NULL }
  if(length(st)>0) { st1 <- as.numeric(dat[[st]]) } else { st1 <- NULL }
  if(length(en)>0) { en1 <- as.numeric(dat[[en]]) } else { en1 <- NULL }
  if(length(wd)>0) { en1 <- st1+as.numeric(dat[[wd]]) } # { en1 <- st1+dat[[wd]] }
  #print(length(st1)); print(head(st1))
  outData <- RangedData(ranges=IRanges(start=st1,end=en1,names=id),space=ch1,universe=build[1])
  outData <- toGenomeOrder(outData,strict=T)
  # note when adding data subsequently that 'RangedData' sorts by genome order, so need
  # to resort any new data before adding.
  if(is.null(rownames(outData))) { rownames(outData) <- paste(1:nrow(outData)) }
  reorder <- match(rownames(outData),id)
  more.cols <- colnames(dat)[!colnames(dat) %in% key.nms]
  more.cols <- more.cols[!more.cols %in% exclude]
  if(length(more.cols)>0) {
    for (cc in 1:length(more.cols)) {
      outData[[more.cols[cc]]] <- dat[[more.cols[cc]]][reorder]
    }
  }
  return(outData)
}


## usage 
# export.all.or.pv(TR,qc.excluded.snps,qc.cloud.fail,fn="forolly.RData")
# function takes the object 'TR' tidies and makes a clean dataframe of results ready for export
export.all.or.pv <- function(TR,qc.excluded.snps,qc.cloud.fail,fn="forolly.RData") {
  jj <- order(TR[["meta.p.value"]])
  tt <- as.data.frame(TR)[jj,]
  kk2 <- which(!(rs.to.id(tt[,5]) %in% unique(rs.to.id(c(qc.excluded.snps,qc.cloud.fail)))))
  tt <- tt[kk2,]
  tt <- tt[,-3:-4]
  colnames(tt)[1:2] <- c("Chr","Pos")
  aabb <- AB(tt$names)
  forolly <- cbind(tt,aabb)
  colnames(forolly)[11:12] <- c("allele.A","allele.B")
  colnames(forolly)[10] <- "band"
  colnames(forolly)[3] <- "rsid"
  rownames(forolly) <- rs.to.id(forolly[,3])
  
  save(forolly,file=fn)
  return(forolly)
}


# Create a table comparing the meta and case-control p-values for the old version of the table (uva)
# versus the new version (DIL, includes extra ~2,500 CBR samples) 
# e.g, 
#  > compare.to.previous.results(bonfs.filt,first=T)
#  > compare.to.previous.results(bonfs.filt,first=F,true.reps=c(6,20,23,39))
compare.to.previous.results <- function(bonfs.filt,first=T,true.reps=NULL,
                                        prv.fn="finalMetaTopHits7FEB.RData") {
  new.tab <- bonfs.filt
  print(load(prv.fn))
  old.tab <- bonfs.filt
  all.reg <- unique(c(rownames(new.tab),rownames(old.tab)))
  nr <- length(all.reg); nc <- 5; com.tab <- as.data.frame(matrix(nrow=nr,ncol=nc))
  rownames(com.tab) <- all.reg
  colnames(com.tab) <- c("oldSnp","old_p","newSnp","new_p","equivalent.to.NEW")
  sel1 <- which(rownames(com.tab) %in% rownames(old.tab))
  sel2 <- which(rownames(com.tab) %in% rownames(new.tab))
  com.tab[sel1,1:2] <- old.tab[match(rownames(com.tab)[sel1],rownames(old.tab)),1:2]
  com.tab[sel2,3:4] <- new.tab[match(rownames(com.tab)[sel2],rownames(new.tab)),1:2]
  com.tab[sel2,5] <- sapply(get.equivs(com.tab[sel2,3],iden.list),paste,collapse=",")
  ## by inspection first time, checked whether any old snps in the equivalent list (shouldn't be):
  # if they are, put the row numbers into the vector 'true.reps'
  if(first) { 
    iii <- (which("NA"!=(com.tab[,5]) & (com.tab[,1]!=com.tab[,3])))
    print(com.tab[iii,]) ; stop() 
  }
  # then once vector is checked, first ==FALSE
  if(length(true.reps)>0) {  com.tab[true.reps,3] <- com.tab[true.reps,1] }
  com.tab[sel2,5] <- sapply(get.equivs(com.tab[sel2,3],iden.list),paste,collapse=",")
  
  com.tab[,4] <- paste(substr(com.tab[,4],1,6),substr(com.tab[,4],nchar(com.tab[,4])-3,nchar(com.tab[,4])),sep="")
  chrzz <- apply(cbind(Chr(com.tab[,1]),Chr(com.tab[,3])),1,mean,na.rm=T)
  com.tab <- com.tab[order(chrzz),]
  
  com.tab[["equiv.SNPs"]] <- com.tab[,5]
  com.tab[,5] <- Pos(com.tab[,3],build=37)
  
  colnames(com.tab)[5] <- "New.snp.pos"
  
  return(com.tab)
}





# returns list of snps from old table, not in new
# prints summary of new vs old pvalues and odds ratios for those snps
# e.g, examine.no.longer.t1.snps(table1snpsfinaljan30,bonf.snps,table1a,identicals,T,T)
examine.no.longer.t1.snps <- function(table1snpsfinaljan30,bonf.snps,table1a,identicals,do.OR=FALSE,do.p=TRUE) {
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  excl <- reader(docs[2])
  uva.table <- reader(docs[1])
  why <- table1snpsfinaljan30[!id.to.rs(table1snpsfinaljan30) %in% id.to.rs(bonf.snps)]
  ct.fn <- "conditionalTests.csv"
  ct <- reader(ct.fn,stringsAsFactors=F)
  conditionals <- ct$TABLE1[ct$COND!=0]
  if(length(why)>0 & exists("conditionals")) {  why <- why[!id.to.rs(why) %in% id.to.rs(conditionals)] }
  if(length(why)>0) {   why <- why[!id.to.rs(why) %in% id.to.rs(identicals)] }
  if(length(why)==0) { return("none found a mystery")}  
  
  if(do.p) {
    # compare pvalues new vs old for SNPs in the old table 1, no longer in the new table 1
    ww <- cbind(table1a[rs.to.id(why),"P_Meta"],
                uva.table[rs.to.id(why),"P_Meta"],Chr(why),Pos(why))[order(Chr(why)),]
    xx <- cbind(table1a[rs.to.id(why),"P_CC"],
                uva.table[rs.to.id(why),"P_CC"])[order(Chr(why)),]
    zz<-cbind(xx,ww)
    colnames(zz) <- c("New_CC","Old_CC","New_Meta","Old_Meta","Chr","Pos")
    zz <- zz[,c(5,6,1,2,3,4)]
    rownames(zz) <- why[order(Chr(why))]
    print(zz,digits=4)
  }
  if(do.OR) {
    # compare odds ratios new vs old for SNPs in the old table 1, no longer in the new table 1
    ww2 <- cbind(table1a[rs.to.id(why),"OR_Fam"],
                 uva.table[rs.to.id(why),"OR_Fam"],Chr(why),Pos(why))[order(Chr(why)),]
    xx2 <- cbind(table1a[rs.to.id(why),"OR_CC"],
                 uva.table[rs.to.id(why),"OR_CC"])[order(Chr(why)),]
    zz2 <-cbind(xx2,ww2)
    colnames(zz2) <- c("New_CC","Old_CC","New_Fam","Old_Fam","Chr","Pos")
    zz2 <- zz2[,c(5,6,1,2,3,4)]
    rownames(zz2) <- why[order(Chr(why))]
    print(zz2,digits=4)  
  }
  return(why) # returns list of snps from old table, not in new
}


# T/F see whether snps in a list have been checked for signal clouds previous
unchecked <- function(snps) {
  source("~/github/iChip/hardCodedSnpLists.R",local=TRUE,echo=FALSE)
  list.to <- ((!snps %in% qc.cloud.fail) & (!snps %in% ok) & (!snps %in% good.snps.checked))
  return(list.to)
}

# turn the 'condit.res' object returned by conditionalAnalysis.R into a 
# clean table ready to paste into table 1
condit.to.res <- function(condit.res) {
  newt1rown <- unlist(lapply(condit.res,names))
  newt1p <- unlist(lapply(condit.res,function(X) { lapply(X,"[",1) }))
  newt1simp <- unlist(unlist(lapply(condit.res,function(X) { lapply(X,"[","glm") }),recursive=F),recursive=F)
  OR_SE <- lapply(newt1simp,function(X) { Y <- tail(X,1); return(Y[c(1,5)]) })
  out <- cbind(newt1rown,newt1p,sapply(OR_SE,"[",1),sapply(OR_SE,"[",2))
  colnames(out) <- c("rsid","aov.P","CC_OR","CC_SE.beta")
  out <- cbind(out,Chr(out[,"rsid"]),Pos(out[,"rsid"],37))
  colnames(out)[5:6] <- c("Chr","Pos")
  out <- as.data.frame(shift.rownames(out))
  rn <- rownames(out)
  out <- as.data.frame(lapply(out,as,"numeric"))
  rownames(out) <- rn
  out[["allele.A"]] <- out[["allele.B"]] <- rep("A",nrow(out))
  out[,6:7] <- AB(rownames(out))
  return(out)
}




##### SPECIFIC TO ME ######




# calibrate the bonferroni threshold for conditional analyses
calibrate.cond.bonf <- function(snplist,cm.window=0.1,bp.ext=0,build=37,qclist="snpsExcluded.txt",ret.snps=FALSE) {
  all.snps.tested <- NULL
  qc.excluded.snps <- reader(qclist)
  qc.excluded.snps <- qc.excluded.snps[!rs.to.id(qc.excluded.snps) %in% rs.to.id(snplist)]
  
  for (cc in 1:22) {
    cat("chr",cc,"\n")
    snpid.list <- snplist[Chr(snplist) %in% cc]
    if(length(snpid.list)<1) { next }
    grp.labs <- get.nearby.snp.lists(snpid.list,cM=cm.window,bp.ext=bp.ext,build=37,excl.snps=qc.excluded.snps,name.by.bands=FALSE)
    all.snps.tested <- c(all.snps.tested,unlist(grp.labs))  
  }
  
  ast <- all.snps.tested[!rs.to.id(all.snps.tested) %in% rs.to.id(c(qc.excluded.snps,qc.cloud.fail))]
  bcf <- length(ast)
  bcfu <- length(unique(ast))
  cat("implied bonferroni (count all tests) threshold is:",.05/bcf,"\n")
  cat("implied bonferroni (count unique snps only) threshold is:",.05/bcfu,"\n")
  if(ret.snps) {
    return(unique(ast))
  } else {
    return(list(all.tests=bcf,unique.snps=bcfu))
  }
}




if(F) {
  all.support <- reader("/chiswick/data/ncooper/iChipData/ImmunoChip_ChipInfo.RData")
  
  LL <- length(which(!is.na(qchisq(p.cc,1)))) ; plot(qchisq(1-((1:LL)/LL),1),qchisq(1-sort(p.cc),1),xlim=c(0,6),type="l",main="Case-control Analysis (ChiSq)",xlab="expected",ylab="observed"); lines(x=qchisq(1-((1:LL)/LL),1),y=qchisq(1-((1:LL)/LL),1),col="red",lty="dotted")
  ww <- (which(!is.na(-log10(p.cc)))); LL <- length(ww)
  plot(-log10(((1:LL)/LL)),-log10(sort(p.cc)),xlim=c(0,6),type="l",main="Case-control Analysis -Log10(p)",xlab="expected",ylab="observed"); lines(x=-log10(((1:LL)/LL)),y=-log10(((1:LL)/LL)),col="red",lty="dotted")
  
  do.the.meta.ones <- function(p.meta) {
    LL <- length(which(!is.na(qchisq(p.meta,1)))) ; plot(qchisq(1-((1:LL)/LL),1),qchisq(1-sort(p.meta),1),xlim=c(0,6),type="l",main="Meta Analysis (ChiSq)",xlab="expected",ylab="observed"); lines(x=qchisq(1-((1:LL)/LL),1),y=qchisq(1-((1:LL)/LL),1),col="red",lty="dotted")
    ww <- (which(!is.na(-log10(p.meta)))); LL <- length(ww)
    plot(-log10(((1:LL)/LL)),-log10(sort(p.meta)),xlim=c(0,6),type="l",main="Meta Analysis -Log10(p)",xlab="expected",ylab="observed"); lines(x=-log10(((1:LL)/LL)),y=-log10(((1:LL)/LL)),col="red",lty="dotted")
  }
  
  print(load("/chiswick/data/ncooper/iChipData/compiledTableAllResultsPassingQC.RData"))
  near.region <- grep("EXT",tt$gene)
  out.region <- grep("OTHER",tt$gene)
  in.region <- which(!1:nrow(tt) %in% c(near.region,out.region))
  
  print(load("pvaluesforqq.RData"))
  p.m.reg <- p.meta[names(p.meta) %in% rs.to.id(tt$names[in.region])]
  p.m.out <- p.meta[names(p.meta) %in% rs.to.id(tt$names[out.region])]
  do.the.meta.ones(p.m.reg)
  do.the.meta.ones(p.m.out)
}





# wut?
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

