## imputation subscript ## 
# used by conditional analysis and indistinguishable analyses scripts, does the imputation,
# saving imputed results as we go so don't need to keep recalculating the same SNPs

imp.file <- cat.path(work.dir,pref="Imputed_chr",fn=next.chr,suf=paste("grp",grp,sep=""),ext="RData")
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
  if(nrow(myDat)!=length(smp.filt)) { stop("mismatching number of samples in loaded file") }
  if(ncol(myDat)!=length(snp.filt)) { 
    ## if this has changed a bit, try to ressurrect without recalculating the whole thing
    cat("mismatching number of snps in loaded file\n") 
    if(!exists("bigDat")) { print(load("allImputed.RData")) }
    targs <- colnames(myData)[snp.filt]
    gotem <- narm(match(targs,colnames(bigDat)))
    aintgotem <- targs[!targs %in% colnames(bigDat)]
    if(length(aintgotem)>0) {
      if(length(gotem)>0) {
        cat("combining",length(gotem),"previously imputed with",length(aintgotem),"from scratch\n")
        myDat.part1 <- impute.missing(myData[smp.filt,match(aintgotem,colnames(myData))],numeric=T)
        if(!exists("bigDat")) { print(load("allImputed.RData")) }
        myDat.part2 <- bigDat[,gotem]
        if(any(rownames(myDat.part1)!=rownames(myDat.part2))) { stop("samples did not match up") }
        myDat.cbind <- cbind(myDat.part1,myDat.part2)
        indz <- match(targs,colnames(myDat.cbind))
        if(any(is.na(indz))) { stop("combined data still missing target SNPs") }
        myDat <- myDat.cbind[,indz]
      } else {
        ## do all from scratch
        cat("none were previously imputed, calculating from scratch\n")
        myDataFilt <- myData[smp.filt,snp.filt]
        myDat <- impute.missing(myDataFilt,numeric=T)
      }
    } else { 
      ## have all the snps, but need to prune some
      cat("trimming loaded data to subset needed\n")
      indz <- match(targs,colnames(bigDat))
      if(any(is.na(indz))) { stop("not sure why loaded data is missing target SNPs") }
      myDat <- bigDat[,indz]
    }
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

