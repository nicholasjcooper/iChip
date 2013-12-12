#R --no-save --slave < ~/github/iChip/conditionalAnalysis.R > condanalFri29.log

source("~/github/iChip/iFunctions.R")
work.dir <- "/chiswick/data/ncooper/iChipData/"
setwd(work.dir)
library(reader)
library(annotSnpStats)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")

### output: 'all.results'; a list containing the OR,OR_L,OR_H,P-value for each snp in 'topsnplist', by chromosome

bonf <- .05/23236      #3.23*(10^-7) # bonferroni threshold #
covs <- TRUE

if(!exists("covs")) { covs <- c(TRUE,FALSE)[2] }

load.uva.table <- TRUE  # results of UVA analyses for each SNP
load.ichip.regions <- TRUE # ichip regions as RangedData

### naughty SNPs ###
qc.cloud.fail <- c("imm_19_10324118","rs1016431","rs1160544","rs12743005","rs12927355","rs151234","rs1788097","rs202535","rs229533","rs2793108","rs3024493","rs34934650","rs3826110","rs4656327","rs4759229","rs4935063","rs56994090","rs5763751","rs5763779","rs61955086","rs666619","rs6691977","rs6710297","rs6907898","rs72729320","rs8074179","rs9585056","rs34536443","rs36038753") #,"rs61839660" last 2 new
ok <- c("rs3024493","rs61839660","rs56994090","rs12927355","rs151234","rs229533")
qc.cloud.fail <- qc.cloud.fail[!qc.cloud.fail %in% ok]

#OLD and outdated list
#table1passtopsnps <- c("rs2476601","rs3024505","rs7511678","rs2111485","rs3087243","rs4850921","rs115288943","rs34046593","rs75793288","rs2611215","rs11950831","rs72928038","rs2326451","rs7867224","rs61839660","rs12416116","rs722988","rs6357","rs917911","rs877636","rs3184504","rs9517719","rs1350276","rs1054000","rs61573680","rs34593439","rs741172","rs151233","rs8056814","rs36038753","rs1893217","rs71178827","rs74956615","rs485186","rs2250261","rs9981624","rs1548389","rs229536")

print(load("finalMetaTopHits09DEC.RData")) # new top list since mon 9/12/13, curated manually using clouds

## load the results of the UVA meta, c/c and family analyses ##
if(load.uva.table) {
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  excl <- reader(docs[2])
}

### get list of snps to test
sl <- "table1snplist.txt"
sl <- "nonbonf.txt"
#sl <- "topsnplist.txt"
#sl <- "allthesnps.txt"
#topsnplist <- reader(sl,work.dir)
#topsnplist <- table1passtopsnps
grpzo <- names(bonf.snps)
topsnplist <- paste(bonf.snps)
opsnplist <- topsnplist[!topsnplist %in% c("rs4378452","rs12924112","rs11066320","rs1989870")] # remove overlapping regions Snps
gg <- grep("rs689",topsnplist)
if(length(gg)>0) { topsnplist <- topsnplist[-gg] } # need to remove this SNP as not actually on iChip, just Taqman

print(load(cat.path(work.dir,"all.support.RData")))   ## load object: all.support [snp support for whole chip]
#all.support <- clean.snp.support(all.support) # this clean should have recently been done
topsnplistic <- all.support$SNP[match(topsnplist,all.support$dbSNP)]

## load the results of the UVA meta, c/c and family analyses ##

## load the ichip dense mapping regions  
if(load.ichip.regions) {
  # iChip regions
  print(load("/chiswick/data/ncooper/iChipData/ichip.regions.RData"))
  if(!exists("cyto"))  { cyto <- get.cyto(); cyto[["gene"]] <- rownames(cyto) }
}


  ######################
  #### MAIN PROGRAM ####
  ######################

  ### GET LIST OF FILE NAMES FOR RDATA FILES FOR EACH CHROMOSOME ###
  ofn <- "/chiswick/data/ncooper/iChipData/temp.ichip-data.RData"
  chr.dat <- cat.path(fn="temp.ichip-data",suf=paste(1:22),ext="RData")
  chr.dat <- as.list(chr.dat)

  ## do sample QC on each chrosome and recombine ##
  if(!exists("ms")) { ms <- list.rowsummary(chr.dat) } # if re-running the script, save repeating the QC
  chrz <- 1:22
  ## create filters for call rate and Heterozygosity ##
  cr.filt <- ms$Call.rate>=0.953
  hz.filt <- ms$Heterozygosity>=0.19 & ms$Heterozygosity<=0.235
  sample.filt <- cr.filt & hz.filt  # logical filter
  excl.ids <- rownames(get.SnpMatrix.in.file(chr.dat[[22]]))[!sample.filt]  # ID exclusion list
  excl.snps <- rownames(excl) # from UVA, above
  #not.top.snps = non significant ids from #clear.not.tops <- lapply(all.results,function(Y) { lapply(Y,clearly.suck,thresh=20) } )
  excl.snps <- unique(c(excl.snps,qc.cloud.fail)) #  ,not.top.snps))
  sample.excl1 <- reader("sample-exclusions-2011-08-09.tab")[,2]
  sample.excl2 <- reader("unacknowledged-dups-2011-08-11.tab")[,5]
  excl.ids <- unique(c(excl.ids,sample.excl1,sample.excl2))
  ###

  all.results <- vector("list",22)
  file.out <- cat.path(work.dir,"all.results",suf=gsub(":",".",gsub(" ","_",date())),ext="RData")

  # run separately for each chromosome #
  for(next.chr in chrz) {
    Header(paste("Chromosome",next.chr))
    chr <- next.chr
    print(load(cat.path(fn=ofn,suf=next.chr,ext="RData")))
    #loads: annotated.snp.support, t1d.data, t1d.support, control.data, control.support

    # gather data and support for next chromosome
    rownames(control.support) <- remove.X(rownames(control.support))
    rownames(t1d.support) <- remove.X(rownames(t1d.support))
    snp.support <- clean.snp.support(annotated.snp.support) # fix little issues with SNP vs dbSNP names
    snp.support <- clean.snp.support(snp.support)
    #ioio <- substr(snp.support$SNP,1,2)
    #print("any with bad prefix after clean?")
    #print(snp.support[grep("2_",ioio),])
    myData <- rbind(control.data,t1d.data) # combine the cases and controls
    colnames(myData) <- clean.snp.ids(colnames(myData))
    snp.qc <- col.summary(myData[sample.filt,])
    maf <- snp.qc$MAF>0.005
    clr <- snp.qc$Call.rate>.99
    hwe <- abs(snp.qc$z.HWE)<3.75
    qc.excl.snps <- rownames(snp.qc)[which(!maf | !clr | !hwe)]
    qc.excl.snps <- qc.excl.snps[!qc.excl.snps %in% c("rs3842727","imm_11_2141424","rs12212193","imm_6_91053490")]
    excl.snps <- unique(c(excl.snps,qc.excl.snps))
    ## extract covariates as 'cov.dat' ##
    if(covs) {
      nms <- c(rownames(control.support),rownames(t1d.support))
      region13 <- c(control.support$b58cregion,t1d.support$b58cregion) 
      the.sex <- c(control.support$sex,t1d.support$sex)
      the.pheno <- make.pheno(myData,rownames(t1d.support),rownames(control.support))
      cov.dat <- data.frame(sex=the.sex,region=region13,pheno=the.pheno)
      rownames(cov.dat) <- nms
     ## cov.dat <- fix.rownames(cov.dat)
      cat(" using covariates for region and sex:\n")
      prv(cov.dat)
    }
    #AmyData <- new("aSnpMatrix", .Data = myData, 
     #              snps = snp.support, samples = cov.dat, phenotype = "pheno", 
      #             alleles = c("allele.A","allele.B"))
    
    ## get the list of SNP/dbSNP ids that are in the current chromosome
    snpid.list <- topsnplist[topsnplist %in% snp.support$dbSNP]
    snpic.list <- snp.support$SNP[match(snpid.list,snp.support$dbSNP)]
    if(any(is.na(snpic.list))) { warning("missing dbSNP entries in list:",snpid.list[is.na(snpic.list)])}
    snpid.list <- snpid.list[!is.na(snpic.list)]
    snpic.list <- narm(snpic.list)
    if(length(snpic.list)<1) { next } # skip this chromosome if no SNPs were found to analyse
    
    ### READY TO ADD RECWINDOW FROM 'WORKING.R' ##
    #recwindow(chr,st,en=st,window=0.1,bp=0)
    if(T) {
      which.snps <- match(snpid.list,snp.support$dbSNP)
      if(any(is.na(which.snps))) { stop(paste("NAs in dbSNP match:",paste(snpid.list[is.na(which.snps)],collapse=","))) }
      #snps.next <- snpid.list
      #snps.next.SNP <- snpic.list
      snps.locs <- snp.support$Pos[which.snps]
      if(all(snp.support$SNP[which.snps]==snpic.list)) { cat("it's true!\n") }
      warning("dup SNPs:",snp.support$SNP[which.snps][duplicated(snp.support$SNP[which.snps])],"\n")
      snp.rd <- RangedData(ranges=IRanges(start=snps.locs,
                                          end=snps.locs,names=snp.support$SNP[which.snps]),
                           space=rep(next.chr,length(snps.locs)))
      snp.rd <- annot.cnv(snp.rd,gs=cyto); colnames(snp.rd) <- "band"
      bands <- snp.rd$band
      nxt.window <- lapply(snps.locs, function(X,...) { recwindow(st=X,...) },chr=next.chr,window=1,bp=50000)
      print(nxt.window)
      st.window <- lapply(nxt.window, "[",1)
      en.window <- lapply(nxt.window, "[",2)
      n.snps <- vector("list",length(st.window))
      for(cc in 1:length(st.window)) {
        n.snps[[cc]] <- which(snp.support$Chr==next.chr &
                              snp.support$Pos>=st.window[cc] & 
                              snp.support$Pos<=en.window[cc] &
                            (!snp.support$SNP %in% excl.snps) &
                            (!snp.support$dbSNP %in% excl.snps) 
        )
      }
      grp.labs <- lapply(n.snps,function(X) { snp.support$SNP[X] })
      grp.snps <- lapply(n.snps,function(X) { snp.support$dbSNP[X] })
      for (cc in 1:length(grp.labs)) { 
        if(any(is.na(grp.snps))) { stop(paste("NAs in grp.labs match:",paste(grp.snps[is.na(grp.snps)],collapse=","))) }
        grp.snps[[cc]][is.na(grp.snps[[cc]])] <- grp.labs[[cc]][is.na(grp.snps[[cc]])]
        grp.snps[[cc]][duplicated(grp.snps[[cc]])] <- grp.labs[[cc]][duplicated(grp.snps[[cc]])]
      }
      if(length(unique(bands))!=length(bands)) { warning("these bands are not unique ==> ",paste(bands[duplicated(bands)],collapse=",")) }
      grpz <- 1:length(bands)
      names(grp.snps) <- names(grp.labs) <- paste(grpz,bands,sep=":")
    } 
    
    n.grps <- length(grp.labs) ; chr.results <- vector("list",n.grps)
    if(length(snpic.list)!=n.grps) { 
      stop("why is length(snpic.list)=",length(snpic.list),"but length(grp.labs)=",n.grps) }
    #### MAIN LOOP OF TOP-SNP GROUPINGS ####
    for (grp in 1:n.grps) {
      snp.filt <- (match(grp.labs[[grp]],colnames(myData)))
      if(any(is.na(snp.filt))) { 
        warning("These were missing:",paste(grp.labs[[grp]][is.na(snp.filt)],collapse=","),"\n") ; 
        X <- grp.labs[[grp]][is.na(snp.filt)]
        Y <- clean.snp.ids(grp.labs[[grp]][is.na(snp.filt)]); print(Y)
        snp.filt <- narm(snp.filt) 
      } 
      #1_113034078,1_113044737,1_113056157
      smp.filt <- which(!rownames(myData) %in% excl.ids)
      
      ## select desired SNPs and samples, imputed and make a dataframe and SnpMatrix version
      source("~/github/iChip/imputationBit.R")

      ## summarise dimension of working dataset ##
      cat("Analysis dataset:\n")
      print(Dim(myDat))
      
      ## Collect phenotype information  ##
      Pheno <- make.pheno(myDat,rownames(t1d.support),rownames(control.support))

      if(any(is.na(Pheno))) { 
        stop("some samples were neither in the T1D nor Control dataset") 
      } else {
        cat("successfully matched",length(which(Pheno==1)),"cases and",length(which(Pheno==0)),"controls\n\n")
      }
      
      ## check covariates, and if ok/present, put into sex and region variables ##
      if(covs) { 
        which.covs <- match(rownames(myDat),rownames(cov.dat))
        if(any(is.na(which.covs))) { 
          warning("missing covariates, ignoring covariates") ; covs <- F 
        } else {
          sex <- cov.dat$sex[which.covs]; region <- as.factor(cov.dat$region[which.covs])
        }
      } 
      
      ## MAIN ANALYSIS ##
      # run GLM analysis for each SNP in the list on the current chromosome
      n.in.grp <- length(grp.labs[[grp]])
      #n.in.grp <- min(200,n.in.grp)
      #result <- list() 
      #bic <- numeric()
      ## decide between model with and without covariates ##
      cat("running conditional analysis on",n.in.grp,"SNPs in surrounding 1 centimorgan + 50Kb region of marker:",snpid.list[grp],"\n")
      hero <- which(grp.labs[[grp]] %in% snpic.list[grp]) #topsnplistic)
           
      if(!covs) {
        fm <- paste("Pheno ~ ",grp.labs[[grp]][hero])
      } else {
        fm <- paste("Pheno ~ sex + region +",grp.labs[[grp]][hero])
      }
      snps.in.cond <- grp.labs[[grp]][hero]
      keepgoing <- T; ccc <- 0
      result.glm <- bic <- ano <- list()
      basemod <- glm(as.formula(fm), family = binomial(logit),data=myDat)
      bic[[1]] <- BIC(basemod)
      while(keepgoing) {
        ccc <- ccc+1
        cat("Testing",ncol(myDat)-length(snps.in.cond),"Snps against base model:",fm,"\n")
        result <- snp.rhs.tests(as.formula(fm),snp.data=myDatSnp,data=myDat)
        legal.ps <- p.value(result)[!names(result) %in% snps.in.cond]
        minp <- min(legal.ps,na.rm=T)
        bst <- which(p.value(result)==minp)
        #if(minp<=bonf) { bst <- which(p.value(result)==minp) } else { bst <- NA }
        if(length(bst)>1) { warning("multiple best conditional snps:",paste(names(result)[bst],collapse=",")) }
        #if(!is.na(bst)) { 
          best.snp <- names(result)[bst[1]] 
          fm <- paste(fm,best.snp,sep=" + ")
        #} else { 
        #  best.snp <- paste(NA) 
        #  keepgoing <- FALSE
        #}
        nxt.mod <- glm(as.formula(fm), family = binomial(logit),data=myDat)
        nxt.anova <- anova(basemod,nxt.mod,test="Chisq")
        ano[[ccc]] <- nxt.anova
        panova <- tail(nxt.anova[["Pr(>Chi)"]],1)
        if(panova<=bonf) {
          snps.in.cond <- c(snps.in.cond,best.snp)
          result.glm[[ccc]] <- mysumfun(nxt.mod,p.digits=250,ci=T)[[1]]
          bic[[ccc+1]] <- BIC(nxt.mod)
          basemod <- nxt.mod
        } else {
          best.snp <- paste(NA)
          keepgoing <- FALSE
        }
      }
      #nxt <- glm(as.formula(fm), family = binomial(logit),data=myDat)
      #result.glm <- mysumfun(nxt,p.digits=250,ci=T)[[1]]
      #bic <- BIC(nxt)
      chr.results[[grp]] <- list(GLM=result.glm,BIC=bic,SNPs=snps.in.cond,AOV=ano)  # add result, just with dbSNP name
    }
    all.results[[next.chr]] <- chr.results
    save(all.results,file=file.out)
    cat("updated object in: ",file.out,"\n")
  }


  ##


##################################################################################
if(F) {
  if(covs) { indx <- 13 } else { indx <- 1 }
  
  ORs <- lapply(all.results,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,1)) }) } )
  Ps <- lapply(all.results,function(Y) { lapply(Y,function(X) { do.call("[",args=list(X,indx,4)) }) } )
  
  all.fn <- cat.path(work.dir,"totalGWAS.RData")
  save(ORs,Ps,all.results,file=all.fn)
  cat("wrote file",all.fn,"\n")
  
  indxr <- match(topsnplist,table1a$SNP_ID)
  
  
  ss <- cbind(unlist(ORs),unlist(Ps))
  rownames(ss) <- gsub(".OR","",rownames(ss))
  tt <- cbind(ss[topsnplist,],table1a[indxr,c("OR_CC","P_CC")])[,c(1,3,2,4)]
  tt[,1] <- or.conv(tt[,1])
  tt[,2] <- or.conv(tt[,2])
  colnames(tt) <- c("ncOR","uvaOR","ncPv","uvaPv")
  print(tt)
  
}



#######PLOTS

#print(load("~/Downloads/all.resultsSat_Nov_30_13.23.02_2013.RData"))
print(load("all.resultsSat_Nov_30_13.23.02_2013.RData"))

##' @author Claudia Giambartolomei
logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}

suck.bic <- function(X,dif=3) { 
  bic <- X$BIC; 
}

do.bic <- function(X,dif=3) { 
  bic <- X$BIC; 
  -min(bic)
  which.min <- (which(bic==min(bic)))
  logsum((-bic/2)[-which.min])
  which.min <- which(bic==sort(bic)[2])
  logsum((-bic/2)[-which.min])
  #print(X$GLM)
  min.n <- min(bic)
  top.grp <- which(bic<(min.n+dif))
  print(paste("SNP:",names(bic)[top.grp]))
  lapply(X$GLM[top.grp],function(Z) { print(tail(Z,1)[4]) })
}
sum <- lapply(all.results,function(Y) { lapply(Y,suck.bic,dif=3) } )
summary(unlist(sum))
bic

