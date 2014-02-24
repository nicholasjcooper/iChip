
# to import the ichip files in the first place:
ofn <- "/chiswick/data/ncooper/iChipData/temp.ichip-data.RData"
for(chr in 1:22) {
  system(paste("~/github/iChip/load-ichip.R chr=",chr," file=",out.file=ofn,sep=""))
}
system(paste("~/github/iChip/load-ichip.R chr=X file=",out.file=ofn,sep=""))



ambigs <- (which(all.support$dbSNP=="AMBIG"))
ambo <- (all.support[ambigs,])
amb <- character(nrow(ambo))
for (cc in 1:nrow(ambo)) { amb[cc] <- (paste(ambo[cc,c(1,2,3)],collapse="_")) }
all.support[ambigs,"dbSNP"] <- amb
save(all.support,file="all.support.RData")

#table(all.support$dbSNP[(which(duplicated(all.support$dbSNP)))])
# this entry has 3 listings! remove one with same pos
all.support <- all.support[-which(rownames(all.support)=="rs3754055"),]



length(which(dupz$SNP==dupz$dbSNP))

as1 <- all.support[which(!all.support$Pos %in% start(sample.info)),]
si1 <- sample.info[which(!start(sample.info) %in% all.support$Pos),]
# do within dups
ii <- numeric(nrow(as1))
for (cc in 1:nrow(as1)) { uu <- (as1$Pos[cc] - start(si1)); ii[cc] <- which(abs(uu)==min(abs(uu)))[1]  }

# do for entire array
ii <- numeric(nrow(all.support))
for (cc in 1:nrow(all.support)) { uu <- (all.support$Pos[cc] - start(sample.info)); ii[cc] <- which(abs(uu)==min(abs(uu)))[1]  ; loop.tracker(cc,192000)}

# these two have no nearby in dup list c("rs2844750","rs9481155") # these two aren't duplicated
closest <- cbind(as1[,1:3],ii,start(si1)[ii])
closest[["diff"]] <- closest[,5]-closest[,3]

closest <- cbind(all.support[,1:3],ii,start(sample.info)[ii])

print(load("/chiswick/data/ncooper/immunochipRunTest/ANNOTATION/snpinfo.RData"))
colnames(support) <- c("rsID","A1","A2","Pos37")
support[["rsID"]] <- ic.to.rs(rownames(support))

## load the results of the UVA meta, c/c and family analyses ##
if(load.uva.table) {
  doc.path <- "ImChip_T1D_20130913"
  docs <- cat.path(doc.path,list.files(doc.path))
  excl <- reader(docs[2])
  table1 <- reader(docs[1])
}

support[["Pos37"]] <- table1$Position[match(rownames(support),rownames(table1))]
ibs <- immunobase.support()
indx <- match(rownames(support),rownames(ibs))
support[["Pos37"]][!is.na(indx)] <- ibs$Pos[indx][!is.na(indx)]
support[["rsID"]][!is.na(indx)] <- ibs$rsID[indx][!is.na(indx)]
 
aybe <- AB(rownames(support))
support[["A1"]] <- aybe[,1]
support[["A2"]] <- aybe[,2]
save(support,file="enhancedInfo.RData")
print(load("all.support.RData"))

support[["rsID"]][is.na(support[["rsID"]])] <- rownames(support)[is.na(support[["rsID"]])]
support.df <- as.data.frame(support[,1:4])
colnames(support.df)[1:5] <- c("Chr","Pos","p","w","icID")
support.df <- support.df[,c(5,6,1,2,7,8,9)]
support.df <- shift.rownames(support.df)

dupz2 <- dup.tracker(support.df,"rsID")
unq <- unique(dupz2$rsID)
for(cc in 1:length(unq)){
  this <- dupz2[dupz2$rsID==unq[cc],]
  if(nrow(this)==2) {
    if(rownames(this)[1]==this$rsID[1]) {
      if(rownames(this)[2]!=this$rsID[2]) {
        this$rsID[1] <- paste(this$rsID[1],"b",sep="")
      } else {
        stop("unforseen duplicate pattern in",this$rsID[1])
      }
    } else {
      if(rownames(this)[2]==this$rsID[2]) {
        this$rsID[2] <- paste(this$rsID[2],"b",sep="")
      } else {
        warning("strange pattern, adding 'b' to second"); print(this)
        this$rsID[2] <- paste(this$rsID[2],"b",sep="")
      }
    }
  } else {
    stop("not expecting sets of ",nrow(this))
  }
  dupz2[dupz2$rsID==unq[cc],] <- this
}
support.df$rsID[dupz2$rn] <- dupz2$rsID
dothem <- is.na(support$Pos37) & (!chr(support) %in% c("MT","XY","YX"))
b37 <- (conv.36.37(support[dothem,]))
support[["Pos37"]][match(rownames(b37),rownames(support))] <- b37[,"Pos"]
rownames(support.df) <- clean.snp.ids(rownames(support.df))

support.df[match(clean.snp.ids(rownames(b37)),rownames(support.df)),"Pos37"] <- b37[,"Pos"]

as.extras <- all.support[match(clean.snp.ids(rownames(support.df)),clean.snp.ids(rownames(all.support))),c("type","gene.annotation","cSNP","consequence")]

all.support2 <- cbind(support.df,as.extras)
all.support <- all.support2
save(all.support,file="all.support.RData")

