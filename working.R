bonf <- 3.23*(10^-7)


highlights <- function(X) {
  top <- min(X,na.rm=T)
  wh <- which(X==top)
  if(length(wh)>1) { cat("more than one top hit: ",paste((X)[wh],collapse=","),"\n") }
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


setwd("/chiswick/data/ncooper/iChipData")
library(reader)
require(genoset); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")


if(T) {
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
  
  
  # iChip regions
  if(F) {
    print(load("/chiswick/data/ncooper/iChipData/dense.ic.regions.b36.RData"))
    ## gives  regions.gr
    cyto <- get.cyto(); cyto[["gene"]] <- rownames(cyto)
    cy.ranges <- toGenomeOrder(as(regions.gr,"RangedData"))
    cy.rangesX <- annot.cnv(cy.ranges,gs=cyto)  ### YES!!!! ###
    ichip.regions <- cy.rangesX
    save(ichip.regions,file="/chiswick/data/ncooper/iChipData/ichip.regions.RData")
  } else {
    print(load("/chiswick/data/ncooper/iChipData/ichip.regions.RData"))
  }
  
  
  table.ranges <- RangedData(ranges=IRanges(start=poz,end=poz,names=table1a[,1]),
                             space=table1a[,2],OR=table1a[,9],p.value=table1a[,11],fam.OR=table1a[,12],
                             fam.p.value=table1a[,14],meta.OR=table1a[,15],meta.p.value=table1a[,17])
  
  cyto <- get.cyto(); cyto[["gene"]] <- rownames(cyto)
  
  gs <- ichip.regions
  table.ranges1 <- annot.cnv(table.ranges,gs=gs)
  gs[["gene"]] <- paste("EXT",gs[["gene"]],sep="_")
  start(gs) <- pmax(1,(start(gs)-50000))
  end(gs) <- end(gs)+50000
  table.ranges2 <- annot.cnv(table.ranges,gs=gs)
  gs2 <- cyto
  gs2[["gene"]] <- paste("OTHER",gs2[["gene"]],sep="_")
  table.ranges3 <- annot.cnv(table.ranges,gs=gs2)
  table.ranges1[["gene"]][which(table.ranges1$gene=="")] <- table.ranges2[["gene"]][which(table.ranges1$gene=="")]
  table.ranges1[["gene"]][which(table.ranges1$gene=="")] <- table.ranges3[["gene"]][which(table.ranges1$gene=="")]
  
  table.ranges <- table.ranges1
  MHC <- RangedData(ranges=IRanges(start=25000000,end=35000000),space=6)
  remv <- queryHits(findOverlaps(table.ranges[6],MHC))
  #new6 <- table.ranges[6][-remv,]
  st6 <- (chrIndices(table.ranges)[6,1]-1)
  TR <- table.ranges[-(remv+st6),]
  if(length(queryHits(findOverlaps(TR[6],MHC)))) { cat("MHC not removed\n") } else { cat("MHC successfully removed\n") }
  
  
  # chr(cy.ranges) <- gsub("chr","",chr(cy.ranges))
  # cy.chr <- gsub("chr","",chr(cy.ranges))
  # annot.cnv(table.ranges,gs=cy.ranges)
  # table.rangesX <- annot.cnv(table.ranges,gs=cy.ranges)
  
  topsnplist <- reader(fn="/chiswick/data/ncooper/iChipData/topsnplist.txt")
  #table.ranges <- annot.cnv(table.ranges,gs=gs)
  jj <- order(TR[["meta.p.value"]])
  tt <- as.data.frame(TR)[jj,]
  print(load("all.support.RData"))
  kk2 <- all.support$SNP[match(topsnplist,all.support$dbSNP)]
  tt <- tt[kk2,]
  kk <- which(as.numeric(tt[["meta.p.value"]])<(10^-5))
  tt <- tt[kk,]
  prv.large(tt,rows=100)
}



if(T) {
  ## extract infuva-controls.txt	3	txt	1	2	5	6	12	1nplist.txt")
  #table.ranges <- annot.cnv(table.ranges,gs=gs)
#do for each 'row' (regional summary)
  out.list <- tapply(tt$meta.p.value,tt$gene,highlights) # main stats
  out.snps <- tapply(tt$names,tt$gene,"[",1) #  top snp (1st because sorted)
  grp.snps <- tapply(tt$names,tt$gene,c) # snp list
  # convert list to data.frame and format
  out.frame <- cbind(sapply(out.list,"[",1),sapply(out.list,"[",2),sapply(out.list,"[",3),sapply(out.list,"[",4))
  colnames(out.frame) <- c("whichSNP","best.p.value","region.hits","hits<bonferroni")
  out.frame <- cbind(out.frame,out.frame[,3]-out.frame[,4])
  colnames(out.frame)[5] <- c("hits>=bonferroni")
  out.frame[,1] <- as.character(out.frame[,1]) 
  out.frame[,1] <- out.snps
  ###
  top.snps <- out.frame[as.numeric(out.frame[,"best.p.value"])<(10^-5),"whichSNP"]
  names(top.snps)[1] <- "ChrX"
  top.snps.dat <- out.frame[as.numeric(out.frame[,"best.p.value"])<(10^-5),]
  save(top.snps, top.snps.dat, file="topSNPs.RData")

  # add genes for each region
  genes <- get.gene.annot()
  bandz <- paste(chr(genes),genes[["band"]],sep="")
  nmz <- rownames(out.frame)
  nmz <- gsub("OTHER_","",nmz)
  nmz <- gsub("EXT_","",nmz)
  genz <- lapply(nmz,function(nmz) { genes[["gene"]][which(bandz %in% nmz)] })
  ## actually this is too many genes to look at!
  
  
  ## GATHER THE REGIONS ##
  all.nms <- names(grp.snps)
  three.list <- all.chr <- vector("list",3) # hold the normal, extended and other lists
  all.pref1 <- substr(all.nms,1,2)
  three.list[[1]] <- all.nms[!((all.pref1 %in% "OT") | (all.pref1 %in% "EX"))]
  three.list[[2]] <- gsub("EXT_","",all.nms[all.pref1 %in% "EX"])
  three.list[[3]] <- gsub("OTHER_","",all.nms[all.pref1 %in% "OT"])
  for (cc in 1:length(three.list)) {
    all.chr[[cc]] <- substr(three.list[[cc]],1,2)
    selec <- which(substr(three.list[[cc]],2,2) %in% c("p","q"))
    all.chr[[cc]][selec] <- substr(three.list[[cc]],1,1)[selec]
  }
  names(three.list) <- c("","EXT","OTHER")
  ## NB must choose rule to account for OTHER_ and EXT_
  
  ## NB: original code source was: "/home/chrisw/local/R/scripts/hapmap-rates.R"
  #recwindow <- function(chr,st,en=st,window=0.1,bp=0,do.plot=FALSE,
  #                      add.plot=FALSE,do.lines=TRUE,...)
  
}


stop()

### FOLLOWUP ANALYSIS FOR SNPs ###
topsnplist <- reader(fn="/chiswick/data/ncooper/iChipData/topsnplist.txt")
topsnplist[topsnplist=="rs6679677"] <- "rs2476601"

######## do each chr ##########
#chrz <- 12  # add all chrs here
followup <- TRUE
if(followup) {
  chrz <- 1:22
} else {
  chrz <- narm(as.numeric(unique(unlist(all.chr))))
}

for(next.chr in chrz) {
  # load data for one chromosome at a time
  Header(paste("Chromosome",next.chr))
  ofn <- "/chiswick/data/ncooper/iChipData/temp.ichip-data.RData"
  chr <- next.chr; #st <- 1; en <- get.chr.lens()[next.chr]
  ###system(paste("~/github/iChip/load-ichip.R chr=",chr," file=",out.file=ofn,sep=""))
  #system("~/github/iChip/load-ichip.R --chr",chr,"--start",st,"--end",en,"--file",out.file=ofn)
  
  print(load(cat.path(fn=ofn,suf=next.chr,ext="RData")))
  #annotated.snp.support, t1d.data, t1d.support, control.data, control.support
  
  Ph <- rep(c(1,0),times=c(nrow(t1d.data),nrow(control.data)))
  myData <- rbind(t1d.data,control.data)
  covs <- FALSE
  
  if(covs) {
    region13 <- c(t1d.support$b58cregion,control.support$b58cregion) 
    the.sex <- c(t1d.support$sex + control.support$sex)
  }
  
  ### NOW CYCLE THROUGH SNP GROUPS ON THIS CHR ####
#  grpz <- 4 # work out which groups to analyse
  if(followup) {
    which.snps <- which(annotated.snp.support$dbSNP %in% topsnplist)
    snps.next <- annotated.snp.support$dbSNP[which.snps]
    snps.next.SNP <- annotated.snp.support$SNP[which.snps]
    snps.locs <- annotated.snp.support$Pos[which.snps]
    snp.rd <- RangedData(ranges=IRanges(start=snps.locs,
                                        end=snps.locs,names=annotated.snp.support$SNP[which.snps]),
                                        space=rep(next.chr,length(snps.locs)))
    snp.rd <- annot.cnv(snp.rd,gs=cyto); colnames(snp.rd) <- "band"
    bands <- snp.rd$band
    nxt.window <- lapply(snps.locs, function(X,...) { recwindow(st=X,...) },chr=next.chr,window=1,bp=50000)
    st.window <- lapply(nxt.window, "[",1)
    en.window <- lapply(nxt.window, "[",2)
    n.snps <- vector("list",length(st.window))
    for(cc in 1:length(st.window)) {
      n.snps[[cc]] <- which(annotated.snp.support$Chr==next.chr &
                            annotated.snp.support$Pos>=st.window[cc] & 
                            annotated.snp.support$Pos<=en.window[cc] &
                            (!annotated.snp.support$SNP %in% rownames(excl)) 
                            )
    }
    grp.labs <- lapply(n.snps,function(X) { annotated.snp.support$SNP[X] })
    grp.snps <- lapply(n.snps,function(X) { annotated.snp.support$dbSNP[X] })
    for (cc in 1:length(grp.labs)) { 
      grp.snps[[cc]][is.na(grp.snps[[cc]])] <- grp.labs[[cc]][is.na(grp.snps[[cc]])]
      grp.snps[[cc]][duplicated(grp.snps[[cc]])] <- grp.labs[[cc]][duplicated(grp.snps[[cc]])]
    }
    names(grp.snps) <- names(grp.labs) <- bands    
    grpz <- 1:length(bands)
  } else {
    mainhit <- three.list[[1]][which(all.chr[[1]]==paste(next.chr))]
    nearhit <- paste("EXT",three.list[[2]][which(all.chr[[2]]==paste(next.chr))],sep="_")
    outhit <- paste("OTHER",three.list[[3]][which(all.chr[[3]]==paste(next.chr))],sep="_")
    grpz <- which(names(grp.snps) %in% mainhit)
  }
  
  if(length(grpz)==0) { next }
  for(grp in grpz) {
    the.sig <- NULL
    cat("Testing group:",grp,", around SNP",snps.next.SNP[grp],"[",snps.next[grp],"]"," band",bands[grp],"\n")
    if(!followup) {  
      snpid.list <- annotated.snp.support$SNP[match(grp.snps[[grp]],annotated.snp.support$dbSNP)] 
    } else {
      snpid.list <- grp.labs[[grp]]
    }
    snpid.list <- gsub("-",".",snpid.list,fixed=T)
    snpid.list <- gsub("1kg_","X1kg_",snpid.list,fixed=T)
    snpid.list <- gsub(",","_",snpid.list,fixed=T)
    grp.snps[[grp]] <- gsub("-",".",grp.snps[[grp]],fixed=T)
    grp.snps[[grp]] <- gsub("1kg_","X1kg_",grp.snps[[grp]],fixed=T)
    grp.snps[[grp]] <- gsub(",","_",grp.snps[[grp]],fixed=T)
    #gsub("1kg_","X1kg_",snpid.list,fixed=T)
    were.missing <- which(!snpid.list %in% colnames(myData))
   # prv(snpid.list)
    
    if(length(were.missing)>0){  
      cat("missing SNPs from grp list:\n"); print(grp.snps[[grp]][were.missing]) ; cat("\n") 
      snpid.list <- snpid.list[-were.missing]
      rename <- grp.snps[[grp]][-were.missing]
    } else { 
      rename <- grp.snps[[grp]]
    }
    
    # conditional analysis p value?
    #sum(sapply(grp.snps,length))
    #[1] 2319
    #> .05/2319
    #[1] 2.156102e-05
   # prv(snpid.list)
    if(length(snpid.list)<1) { warning("no snps left to test in this set"); next }
    myDataFilt <- myData[,snpid.list]
    colnames(myDataFilt) <- rename
    if(followup) {
      top.snp <- as.numeric(myDataFilt[,snps.next[grp]])-1
    } else {
      top.snp <- as.numeric(myDataFilt[,1])-1
    }
    #snp.rhs.tests(as.formula("Ph ~ 1"),snp.data=myDataFilt,allow.missing=.1,tests=rename)
    #snp.rhs.tests(as.formula("Ph ~ top.snp"),snp.data=myDataFilt[,-1],allow.missing=.1,tests=rename[-1])
    
    ### now loop through until no further covariates are found ##
    #if(is.null(the.sig)) { the.sig <- top.snp } # else it will be all the conditionally significant so far
    found.conds <- TRUE
    if(followup) {
      excl.cols <- which(colnames(myDataFilt)==snps.next[grp])[1]
    } else {
      excl.cols <- 1
    }
    first <- T
    sig.reg.snps <- snps.next[grp]   #"top.snp"   # snpid.list[1]
    while(found.conds) {
##      ## hack to convert selection to a matrix ##
##      top.snp <- as.data.frame(myDataFilt[,excl.cols])
##      for(jj in 1:ncol(top.snp)) { 
##        top.snp[,jj] <- as.numeric(top.snp[,jj])-1 
##      }
##      top.snp <- as.matrix(top.snp)

      ##
#      prv(top.snp)
      terms <- paste(sig.reg.snps,collapse=" + ")
      #print(colnames(myDataFilt)); print(terms)
      if(covs) {
        cond.res <- snp.rhs.tests(as.formula(paste("Ph ~ top.snp + strata(region13) + the.sex")),snp.data=myDataFilt[,],allow.missing=.1,tests=rename[-excl.cols])
      } else {
        cat("testing",terms,"[",excl.cols,"]\n")
        top.snp[top.snp<0] <- NA
       # prv.large(top.snp)
##      print(summary(apply(top.snp,1,function(X) { length(which(is.na(X))) })))
        mod.txt <- paste("Ph ~ top.snp")

        for (dd in 1:length(sig.reg.snps)) {
          if(!exists(sig.reg.snps[dd])) {
            nuxt <- as.numeric(myDataFilt[,sig.reg.snps[dd]])-1
            nuxt[nuxt<0] <- NA
            assign(sig.reg.snps[dd],nuxt)
          }
        }
        mod.txt <- paste("Pheno ~",terms)
      ##  cond.res <- snp.rhs.tests(as.formula(mod.txt),snp.data=myDataFilt[,-excl.cols]) #,allow.missing=.1,uncertain=TRUE,tests=rename[-excl.cols])
        cov.data1 <- myDataFilt[,-excl.cols]; 
        cond.res <- list()
        cov.data <- as.data.frame(cov.data1)
        for(jj in 1:ncol(cov.data)) { 
          nuxt <- as.numeric(cov.data[,jj])-1
          nuxt[nuxt<0] <- NA
          cov.data[,jj] <- nuxt
         # assign(colnames(cov.data)[jj], nuxt)
        }
        if(first) {
	  row.mis <- (apply(cov.data,1,function(X) { length(which(is.na(X)))/length(X) }))
        }
        col.mis <- apply(cov.data,2,function(X) { length(which(is.na(X)))/length(X) })
        cov.data <- cov.data[row.mis<.03,col.mis<.03]
        if(first){
          Pheno <- Ph[row.mis<.03]
          for(jj in 1:length(sig.reg.snps)) { 
            assign(terms[jj],get(sig.reg.snps[jj])[row.mis<.03])
          }
          for(jj in 1:ncol(cov.data)) {
            #print(colnames(cov.data)[jj])
            assign(colnames(cov.data)[jj], cov.data[,jj])
          }
        }
        #cov.data <- as.matrix(cov.data)
        cat("testing",ncol(cov.data),"snps\n")
      ## cond.res <- snp.rhs.tests(as.formula(paste(mod.txt)), snp.data=cov.data1)
        #print(cond.res)
#        print(head(p.value(cond.res)[order(p.value(cond.res))],20) )
        for (dd in 1:ncol(cov.data)) {
           nxt <- glm(as.formula(paste(mod.txt,colnames(cov.data)[dd],sep=" + ")), family = binomial(logit),data=cov.data)
           cond.res[[dd]] <- mysumfun(nxt,p.digits=250)[[1]][,2]
        }
        #print(cond.res)
      }
      
      if(is.list(cond.res)) {   p.cond <- do.call("rbind",args=cond.res)  } else { p.cond <- p.value(cond.res) }
      if(!is.null(dim(p.cond))) { p.cond <- p.cond[,ncol(p.cond)] }
      cond.sigs <- which(p.cond<=bonf)
      if(length(cond.sigs)>0) {
        naz <- is.na(p.cond)
        the.p <- min(p.cond[!naz],na.rm=T)
        the.min <- which(p.cond[!naz]==the.p)
        if(length(the.min)>1) {
  #        cat("min length",length(the.min),"\n"); 
  #        the.chi <- max(chi.squared(cond.res)[!naz],na.rm=T)
  #        the.min <- which(chi.squared(cond.res)[!naz]==the.chi)
          if(length(the.min)>1) { cat("min chi length",length(the.min),"\n"); the.min <- the.min[1]  }
        }
        the.sig.lab <- rename[-excl.cols][the.min]
        cat("new snp:",the.sig.lab," with conditional p value:",the.p,"\n")
        the.sig <- myDataFilt[,the.sig.lab]
        excl.cols <- c(excl.cols,which(colnames(myDataFilt)==the.sig.lab))
      } else {
        found.conds <- FALSE
      }
      sig.reg.snps <- colnames(myDataFilt)[excl.cols]
      glm.result <- glm(as.formula(mod.txt), family = binomial(logit))
      ii <- mysumfun(glm.result,p.digits=250,o.digits=3)
      print(ii)
      first <- FALSE
      if(length(excl.cols)>5) { found.conds <- FALSE }
    }

    
    print(sig.reg.snps) # append this to the list of top/conditional snps
  }
  ## add chr result
}



### OTHER STUFF

if(F) {
  
  tab <- all.summary[order(all.summary[,"CHR"]),]
  
  
  colnames(table1a)
  cn <- c("CHR","MAF_UKcontrol","OR_CC","P_CC","OR_Fam","P_Fam","OR_Meta","P_Meta")
  snp.list <- outframe$whichSNP
  more.data <- table1a[match(snp.list,table1a[,1]),cn]
  colnames(more.data)[2] <- "MAF"
  all.summary <- cbind(out.frame,more.data)
  print(all.summary[order(all.summary[,"CHR"]),])
  
  do.one.chr <- function(chr.tab) {
    grp <- tapply(chr.tab$gene,chr.tab$gene,function(X) X)
    maxs <- tapply(chr.tab$p.value,chr.tab$gene,min,na.rm=T)
    #names(maxs) <- grp
  }
  
  #chr.results.list <- lapply(tt,do.one.chr)
  results.list <- do.one.chr(tt)
  prv(results.list)
  gen.list <- names(results.list)
  gen.list.spl <- strsplit(gen.list,";",fixed=T)
  prv(gs)
  bands <- lapply(gen.list.spl,
                  function(X) { 
                    ind <- match(X,gs$gene)
                    OO <- paste(unique(paste(chr(gs)[ind],gs$band[ind],sep="")),collapse=";")
                    return(OO)  })
  res2 <- data.frame(genes=substr(gen.list,1,20),p.value=as.numeric(results.list),band=unlist(bands))
  res2 <- res2[order(res2[,2]),]
  prv.large(res2,rows=50,cols=3)
  unique(res2$band)
  
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
