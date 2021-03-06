## This file is used to get the conditional results list 
# iterates through each table 1 hit and tests whether any additional SNPs are able to significantly contribute
# to the model


#R --no-save --slave < ~/github/iChip/conditionalAnalysis.R > condanalThu16.log

source("~/github/iChip/iFunctions.R")
work.dir <- "/chiswick/data/ncooper/iChipData/"
setwd(work.dir)
library(reader)
library(annotSnpStats)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")

### output: 'all.results'; a list containing the OR,OR_L,OR_H,P-value for each snp in 'topsnplist', by chromosome

#bonf <- .05/23236      #3.23*(10^-7) # bonferroni threshold #
thresholds <- list(MAF=0.005,CR=.99,HWE=3.8905,SCR=0.953,HZlo=.19,HZhi=.235,
                                      bonf=3.68*(10^-7),bonfcond=.05/18718) #23236 ) 16782) 18718 #
redo.bonf <- FALSE # whether to recalculate the bonferroni threshold
current.qc.file <- "SNPQC_FEB14.RData"
cm.window <- 1; bp.ext <- 0  #50000  # 0
use.se.OR <- FALSE  # whether to report SE of the beta or the odds ratio (log/exp units)
covs <- TRUE   # use covariates of sex and region in analysis
store.top.snps <-  FALSE #T  #FALSE  store the topsnplist snps in a single matrix for later use
proper.run <- FALSE    # whether to save results to be used in further calcs (or just testing, etc)

if(!exists("covs")) { covs <- c(TRUE,FALSE)[2] }

load.uva.table <- TRUE  # results of UVA analyses for each SNP
load.ichip.regions <- TRUE # ichip regions as RangedData
use.imputation <- TRUE

#print(load("finalMetaTopHits09DEC.RData")) # new top list since mon 9/12/13, curated manually using clouds
# get resulting Snps from 'getMetaTable1.R'
#print(load("/chiswick/data/ncooper/iChipData/finalMetaTopHits6FEB.RData")) # new top list since mon 15/1/14, curated manually 
print(load("/chiswick/data/ncooper/iChipData/finalMetaTopHitsFEB17.RData"))


# get hard coded lists of SNPs passing, failing QC, appearing in table 1, etc #tab1snps
source('~/github/iChip/hardCodedSnpLists.R', echo=FALSE)


## load the results of the UVA meta, c/c and family analyses ##
if(load.uva.table) {
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  excl <- reader(docs[2])
}

### get list of snps to test
grpzo <- names(bonf.snps)
topsnplist <- paste(bonf.snps)
#topsnplist <- c("rs16939895","rs2111485","rs10795791","rs66718203","rs757411","rs12720356")

#just removed these!
##topsnplist <- topsnplist[!topsnplist %in% c("rs4378452","rs12924112","rs11066320","rs1989870")] # remove overlapping regions Snps
##topsnplist <- tab1snps ### ??? last use was line above, but maybe this is better?


#gg <- grep("rs689",topsnplist)
#if(length(gg)>0) { topsnplist <- topsnplist[-gg] } # need to remove this SNP as not actually on iChip, just Taqman
topsnplist[topsnplist %in% "rs3842727"] <- "rs689" # do now with rs689 instead
topsnplist <- topsnplist[!topsnplist %in% c(qc.cloud.fail,not.top.snps)]
topsnplist <- ids.by.pos(topsnplist) # order by genome order

if(redo.bonf) {
  thresholds$bonfcond <- calibrate.cond.bonf(topsnplist,cm.window=cm.window,bp.ext=bp.ext,build=37,qclist="snpsExcluded.txt")
}


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
  ## create filter for call rate and Heterozygosity ##
  sample.filt <- samp.summ(ms,CR=thresholds$SCR,HZlo=thresholds$HZlo,HZhi=thresholds$HZhi) # + print summary 
  # prints SNPQC overall summary, but actually for analysis this is done chr by chr
  ignore.for.now <- snp.summ(MAF=thresholds$MAF,CR=thresholds$CR,HWE=thresholds$HWE,qc.file=current.qc.file)
  excl.ids <- 
    suppressWarnings(rownames(get.SnpMatrix.in.file(chr.dat[[22]]))[!sample.filt])  # ID exclusion list
  excl.snps <- clean.snp.ids(rownames(excl)) # from UVA, above
  #not.top.snps = non significant ids from #clear.not.tops <- lapply(all.results,function(Y) { lapply(Y,clearly.suck,thresh=20) } )
  excl.snps <- unique(c(excl.snps,qc.cloud.fail)) #  ,not.top.snps))
  sample.excl1 <- reader("sample-exclusions-2011-08-09.tab")[,2]
  sample.excl2 <- reader("unacknowledged-dups-2011-08-11.tab")[,5]
  excl.ids <- unique(c(excl.ids,sample.excl1,sample.excl2,nonconsent))
  ###

  all.results <- vector("list",22)
  file.out <- cat.path(work.dir,"cond.results",suf=gsub(":",".",gsub(" ","_",date())),ext="RData")

  # run separately for each chromosome #
  for(next.chr in chrz[chrz>=16 & chrz<=16]) {
    #chrz[chrz %in% c(2,10,16,17,19)]) {
 # for(next.chr in c(19)) {
    Header(paste("Chromosome",next.chr))
    chr <- next.chr
    print(load(cat.path(fn=ofn,suf=next.chr,ext="RData")))
    #loads: annotated.snp.support, t1d.data, t1d.support, control.data, control.support

    # gather data and support for next chromosome
    rownames(control.support) <- remove.X(rownames(control.support))
    rownames(t1d.support) <- remove.X(rownames(t1d.support))
    ## moved from snp.support to all.support ##
    ##snp.support <- clean.snp.support(annotated.snp.support) # fix little issues with SNP vs dbSNP names
    ##snp.support <- clean.snp.support(snp.support)
    #ioio <- substr(snp.support$SNP,1,2)
    #print("any with bad prefix after clean?")
    #print(snp.support[grep("2_",ioio),])
    myData <- rbind(control.data,t1d.data) # combine the cases and controls
    colnames(myData) <- clean.snp.ids(colnames(myData))
    if(next.chr==11){
      rsmtch <- which( colnames(control.data) %in% "imm_11_2138800")
      myData[,rsmtch] <- rep(as.raw("00"),times=nrow(myData)) # "imm_11_2138800" = "rs689"
      rs689.snpmatrix <- get(paste(load("taqman_rs689.RData"))) #get.taqman.snp.from.text.file()
      indzz <- match(rownames(rs689.snpmatrix),rownames(myData))
      myData[narm(indzz),rsmtch] <- rs689.snpmatrix[which(!is.na(indzz)),1]
    }
    snp.qc <- col.summary(myData[sample.filt,])
    maf <- snp.qc$MAF>thresholds$MAF
    clr <- snp.qc$Call.rate>thresholds$CR
    hwe <- abs(snp.qc$z.HWE)<thresholds$HWE
    qc.excl.snps <- rownames(snp.qc)[which(!maf | !clr | !hwe)]
    excl.snps <- unique(c(excl.snps,qc.excl.snps))
    excl.snps <- excl.snps[!excl.snps %in% unique(c(unexclude,topsnplist))]
    
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
      #prv(cov.dat)
    }
    
    ## get the list of SNP/dbSNP ids that are in the current chromosome
    snpid.list <- topsnplist[topsnplist %in% all.support$dbSNP[all.support$Chr==next.chr]]
    snpic.list <- rs.to.ic(snpid.list) #snp.support$SNP[match(snpid.list,snp.support$dbSNP)]
    if(any(is.na(snpic.list))) { warning("missing dbSNP entries in list:",snpid.list[is.na(snpic.list)])}
    snpid.list <- ic.to.rs(snpid.list[!is.na(snpic.list)])
    snpic.list <- narm(snpic.list)
    if(length(snpic.list)<1) { next } # skip this chromosome if no SNPs were found to analyse
    
    ### get snps surrounding top snps by 0.1cM ##
    grp.labs <- get.nearby.snp.lists(snpid.list,cM=cm.window,bp.ext=bp.ext,build=37,excl.snps=excl.snps)
    
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
      if(use.imputation) {
        imp.file <- cat.path(work.dir,pref="Imputed_chr",fn=next.chr,suf=paste("grp",grp,sep=""),ext="RData")
        print(load(ichip.imputation(myData, imp.file, smp.filt, snp.filt)))
      } else {
        myDatSnp <- myData[smp.filt,snp.filt]
        myDat <- SnpMatrix.to.data.frame(myDatSnp)
      }
      # NB: old way was to source("~/github/iChip/imputationBit.R")

      ## summarise dimension of working dataset ##
      cat("Analysis dataset:\n")
      print(Dim(myDat))
      
      ## store top snps in new data object ##
      if(store.top.snps) {
        tsl <- topsnplist #table1snpsfinaljan30
        if(next.chr==1) {
          newfr <- matrix(numeric(),nrow=nrow(myDat),ncol=length(tsl))
          colnames(newfr) <- tsl #table1snpsfinaljan30
          rownames(newfr) <- rownames(myDat)
        }
        to.grab <- rs.to.ic(tsl[Chr(tsl) %in% next.chr])
        #if(rs.to.ic(c("rs12927355","rs193778")) %in% to.grab) { print("THey're HERE!!!!") }
        found <- colnames(myDat)[colnames(myDat) %in% to.grab]
        #if(rs.to.ic(c("rs12927355","rs193778")) %in% found) { print("THey're FOUND!!!!") } else { print("BOO!")}
        if(length(found)>0) {
          wrn <- which(rownames(newfr) %in% rownames(myDat))
          if(length(wrn)!=nrow(newfr)) { warning("mismatch in rownames between newfr and myDat")}
          toget <- rownames(newfr)[wrn]
          if(length(wrn)>0) {
            for(dd in 1:length(found)) {
              newfr[rownames(newfr)[wrn],ic.to.rs(found[dd])] <- myDat[rownames(newfr)[wrn],found[dd]]
            }
            cat("grabbed snps",paste(ic.to.rs(found),collapse=","),"from myDat\n")
          } else {
            warning("no rownames to match")
          }
          #suppressWarnings(prv(newfr))
        }
      }
      
      ## Collect phenotype information  ##
      Pheno <- make.pheno(myDat,rownames(t1d.support),rownames(control.support))

      ## do test of difference in genotype distribution between sanger and uva controls for each SNP ##
      sanger.not.uva <- substr(rownames(myDat[Pheno==0,]),1,1)=="7"
      ppp <- numeric(ncol(myDat)); names(ppp) <- colnames(myDat)
      for (cc in 1:ncol(myDat)) {
        ppp[cc] <- suppressWarnings(chisq.test(table(sanger.not.uva[Pheno==0],round(myDat[Pheno==0,cc])))$p.value)
      }
      baduns <- which(ppp<thresholds$bonfcond) #(0.05/length(ppp)))
      toppers <- (names(ppp)[baduns] %in% rs.to.ic(topsnplist))
      if(any(toppers)) { warning("top snp ",names(ppp)[baduns][toppers]," failed sanger/uva qc") ; baduns <- baduns[-which(toppers)] }
      grp.labs[[grp]] <- grp.labs[[grp]][!grp.labs[[grp]] %in% rs.to.ic(names(ppp)[baduns])]
      cat(length(baduns),"/",length(ppp)," SNPs failed on UVA vs Sanger genotype frequencies\n",sep="")
      if(length(baduns)>0) { cat(paste(names(ppp)[baduns],collapse=","),"\n") }
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

      ## decide between model with and without covariates ##
      cat("running conditional analysis on",n.in.grp,"SNPs in surrounding",cm.window,"centimorgan region of marker:",snpid.list[grp],"\n")
      hero <- which(grp.labs[[grp]] %in% snpic.list[grp]) #topsnplistic)
      if(length(hero)==0) { cat("The following top snps were not found: ",
                                paste(snpic.list[grp],collapse=",","\n check QC lists\n",sep="")); next }
      if(!covs) {
        fm <- paste("Pheno ~ ",grp.labs[[grp]][hero])
      } else {
        fm <- paste("Pheno ~ sex + region +",grp.labs[[grp]][hero])
      }
      snps.in.cond <- grp.labs[[grp]][hero]
      keepgoing <- T; ccc <- 0
      result.glm <- ano <- RES <- list(); bic <- numeric()
      basemod <- glm(as.formula(fm), family = binomial(logit),data=myDat)
      bic[1] <- BIC(basemod)
      while(keepgoing) {
        ccc <- ccc+1
        cat("Testing",ncol(myDat)-length(snps.in.cond),"Snps against base model:",fm,"\n")
        RES[[ccc]] <- result <- snp.rhs.tests(as.formula(fm),snp.data=myDatSnp,data=myDat)
        legal.ps <- p.value(result)[!names(result) %in% c(snps.in.cond)]
        minp <- min(legal.ps,na.rm=T)
        bst <- which(p.value(result)==minp)
        #if(minp<=bonf) { bst <- which(p.value(result)==minp) } else { bst <- NA }
        if(length(bst)>1) { warning("multiple best conditional snps:",paste(names(result)[bst],collapse=",")) }
        best.snp <- names(result)[bst[1]] 
        fm <- paste(fm,best.snp,sep=" + ")      
        nxt.mod <- glm(as.formula(fm), family = binomial(logit),data=myDat)
        nxt.anova <- anova(basemod,nxt.mod,test="Chisq")
        ano[[ccc]] <- nxt.anova
        panova <- tail(nxt.anova[["Pr(>Chi)"]],1)
        if(panova<=thresholds$bonfcond) {
          cw <- caseway(myDat[,best.snp],Pheno)
          mm <- majmin(myDat[,best.snp])
          if((cw=="CasesRef+" & mm=="major") | (cw=="CasesRef-" & mm=="minor") ) { cat("T1d more Major\n") }
	        if((cw=="CasesRef-" & mm=="major") | (cw=="CasesRef+" & mm=="minor") ) { cat("T1d more Minor\n") }
          print(ic.to.rs(best.snp))
          print(cw);print(mm)
          snps.in.cond <- c(snps.in.cond,best.snp)
          result.glm[[ccc]] <- mysumfun(nxt.mod,p.digits=250,ci=T)[[1]]
          bic[ccc+1] <- BIC(nxt.mod)
          basemod <- nxt.mod
        } else {
          cw <- mm <- NA
          best.snp <- paste(NA)
          keepgoing <- FALSE
        }
      }
      chr.results[[grp]] <- list(GLM=result.glm,BIC=bic,SNPs=snps.in.cond,AOV=ano,SNPRHS=RES,QC=ppp,CW=cw,MM=mm)  # add result, just with dbSNP name
    }
    all.results[[next.chr]] <- chr.results
    save(all.results,file=file.out)
    cat("updated object in: ",file.out,"\n")
  }


if(store.top.snps) {
  if(proper.run) {
   save(newfr,file="topsnpmatrix.RData")
   cat("Saved top snps as new SnpMatrix to: topsnpmatrix.RData\n")
  }
  
  ## calculate major / minor alleles and direction of odds ratio for top snp list
  print(load("topsnpmatrix.RData"))
  #topsnpmatrix <- newfr
  Pheno <- make.pheno(newfr,rownames(t1d.support),rownames(control.support))
  mm <- majmin(newfr[Pheno==0,]); cc <- caseway(newfr,Pheno)
  print(cbind(names(cc),paste(cc),paste(mm)))
  
}

# extract lists of conditional snps (as rsids)
qclist <- unlist(lapply(all.results,function(X) { lapply(X,function(Y) { Y$QC })}))
ano.results <- lapply(all.results,function(Y) { lapply(Y,"[",4) })
p.values <- lapply(ano.results,function(Y) { lapply(Y,function(X) { X[[1]][[1]][2,5] })})
cond.results <- lapply(all.results,function(Y) { lapply(Y,"[",3) })
cond <- lapply(cond.results,function(X) { lapply(X,function(Y) { if(length(Y[[1]])>1) { Y[[1]] } else { NULL } }) })
sapply(cond.results,function(X) { 1:length(X) })
sapply(cond.results,function(X) { which(is.null(X)) })
out <- vector("list",length(cond))
for (cc in 1:length(cond)){
  XX <- sapply(cond[[cc]],is.null)
  if(length(XX)>0) {
    out[[cc]] <- rep((which(!XX)),times=sapply(cond[[cc]],length)[which(!XX)]-1)
  } else { out[[cc]] <- NULL }
}
grpz <- unlist(out)
cond <- unlist(cond,recursive=F)
cond <- cond[sapply(cond,length)>0]
cond <- lapply(cond,ic.to.rs)
names(cond) <- paste("chr",sapply(lapply(cond,Chr),"[",1),sep="")

# generate the list used to configure 'indistinguishableCondAnalysis.R'
n.conds <- sum(sapply(cond,length))-length(cond)
condlist <- vector("list",n.conds)
cnt <- 1
for (cc in 1:length(cond)) {
  for (dd in 1:(length(cond[[cc]])-1)) {
    condlist[[cnt]] <- cond[[cc]][1:dd]
    cnt <- cnt + 1
  }
}

# extract p values, etc
glm.results <- lapply(all.results,function(Y) { lapply(Y,"[",1) })
condit.res <- list(); cnt <- 0; nmz <- NULL

# print conditional results where there is a significant n>1 SNP
for (cc in 1:length(glm.results)) {
  if(length(glm.results[[cc]])>0) {
  for (dd in 1:length(glm.results[[cc]])) {
    if(length(glm.results[[cc]][[dd]])>0) {
      for (ee in 1:length(glm.results[[cc]][[dd]])) {
        if(length(glm.results[[cc]][[dd]][[ee]])>0) {
          for(ff in 1:length(glm.results[[cc]][[dd]][[ee]])) {
            if(ee==1 & ff==1) {  cat("Chromosome",cc,"\n")  ; cnt <- cnt + 1; nmz <- c(nmz,cc) ;  condit.res[[cnt]] <- list() }
            rowz <- glm.results[[cc]][[dd]][[ee]][[ff]][-12:-1,]
            ps <- ano.results[[cc]][[dd]][[1]][[ff]][2,5]
            print(ps)
            se.OR <- rowMeans(abs(cbind((((rowz[,2])-(rowz[,1]))/1.96), ((((rowz[,3])-(rowz[,1]))/1.96)))))
            se.beta <- rowMeans(abs(cbind(((log(rowz[,2])-log(rowz[,1]))/1.96), (((log(rowz[,3])-log(rowz[,1]))/1.96)))))
            if(use.se.OR) { se <- se.OR } else { se <- se.beta }
            rowz <- cbind(rowz,se)
            RSs <- ic.to.rs(tail(rownames(rowz),1))
            rownames(rowz) <- ic.to.rs(rownames(rowz))
            condit.res[[cnt]][[paste(RSs)]] <- list()
            condit.res[[cnt]][[paste(RSs)]][["aov"]] <- ps
            condit.res[[cnt]][[paste(RSs)]][["glm"]] <- rowz
            print(RSs)
            print(rowz)
          }
        }
      }
    }
  } }
}

names(condit.res) <- paste("chr",nmz,sep="")

cat("Info to add to table 1 for conditional rows:\n")
if(length(condit.res)>1) {
  print(sum.tab <- condit.to.res(condit.res))
} else { print(condit.res) }

if(proper.run) {
 save(sum.tab,condit.res, grpz,cond,condlist,glm.results,ano.results,p.values,qclist,file="conditionalSummaryObjects.RData")
}
## A PROBLEM= OFTEN CC BEST != BEST - how did i deal with this before???


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
#print(load("all.resultsSat_Nov_30_13.23.02_2013.RData"))

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
#sum <- lapply(all.results,function(Y) { lapply(Y,suck.bic,dif=3) } )
#summary(unlist(sum))
#bic

