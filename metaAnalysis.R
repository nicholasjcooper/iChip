meta.me <- function(X) {
 OR.CC <- X[,"OR_CC"]
 beta.CC  <- log(X[,"OR_CC"])
 se.CC <- X[,"SE_CC"]
 OR.family <- X[,"OR_Fam"]
 beta.family  <- log(X[,"OR_Fam"])
 se.family <- X[,"SE_Fam"]
 z.CC <- beta.CC/se.CC
 z.family <- beta.family/se.family

 inv.CC <- 1 / (se.CC^2)
 inv.family <- 1 / (se.family^2)
 var.meta <- 1 / (inv.CC+inv.family)
 weight.CC <- inv.CC * var.meta
 weight.family <- inv.family * var.meta
 
 famN <- 3819*2  #3509*2   #  3819*2   #  10796
 ccN <- 9416+6670
 WeightFam = sqrt(famN)/(sqrt(famN)+sqrt(ccN))
 #WeightFam = wf
 WeightCC <- 1-WeightFam
 
 beta.meta <- round((weight.CC * beta.CC) + (weight.family * beta.family),digit=3)
 z.metaW1 <- round((weight.CC * z.CC) + (weight.family * z.family),digit=6)
 z.metaW2 <- round((WeightCC * z.CC) + (WeightFam * z.family),digit=6)
 se.meta <- round(sqrt(var.meta), digit=3)
 z.meta <- beta.meta/se.meta
 OR.meta <- exp(beta.meta)
 p.meta <- 2*pnorm(-abs(z.meta))
 p.metaW1 <- 2*pnorm(-abs(z.metaW1))
 p.metaW2 <- 2*pnorm(-abs(z.metaW2))
 out <- (cbind(OR.meta,beta.meta,se.meta,z.meta,p.meta)) #,z.metaW1,p.metaW1,z.metaW2,p.metaW2))
 rownames(out) <- rownames(X)
 return(out)
}

if(F) {
for(cc in 35) {
print(cc/100)
mm <- meta.me(uva.table,cc/100)
iii <- mm[,"p.metaW2"]
kkk <- mm[,"p.meta"]
jjj <- uva.table[,"P_Meta"]
print(Cor(-l10(iii[iii<1.1]),-l10(jjj[iii<1.1])))
}
}

#table1a[,"OR_Meta"] <- OR.meta
#table1a[,"SE_Meta"] <- se.meta
#table1a[,"P_Meta"] <- p.meta

#cor(abs(log(table1[rownames(all.results),"OR_Meta"])),abs(log(table1a[rownames(all.results),"OR_Meta"])),use="pairwise.complete")
#[1] 0.9819474

