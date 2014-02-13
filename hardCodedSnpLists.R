## These are manually entered lists of SNPs passing failing manual QC, and reaching table 1,
# that are manually hard-coded, and used in various scripts 

#nonconsent <- "R58201702_C02"

#rs74956615 <- rs34536443b; ==> rs12720356, ?

table1snpsfinaljan30 <- c("rs34536443b","rs2476601","rs3024493","rs6691977","rs10865035","rs35667974","rs2111485","rs3087243","rs113010081","rs17630466","rs75793288","rs2611215","rs11954020","rs12212193","rs1538171","rs1574285","rs61839660","rs10795791","rs41295121","rs12416116","rs689","rs72853903","rs917911","rs1701704","rs653178","rs9585056","rs1350275","rs56994090","rs72727394","rs34593439","rs151234","rs12927355","rs66718203","rs8056814","rs8067378","rs757411","rs16939895","rs1615504","rs74956615","rs12720356","rs516246","rs6043409","rs80054410","rs4820830","rs229533")
#rs193778b  not included as in a seperate region to merged


# failures hand checked by me
qc.cloud.fail <- c("rs34536443","rs1016431","rs1160544","rs12743005","rs12927355",
                   "rs151234","rs1788097","rs202535","rs229533","rs2793108","rs3024493",
                   "rs34934650","rs3826110","rs4656327","rs4759229","rs4935063","rs56994090",
                   "rs5763751","rs5763779","rs61955086","rs666619","rs6691977","rs6710297",
                   "rs6907898","rs72729320","rs8074179","rs9585056","rs36038753",
                   "rs4142967", "rs62097857", "rs62323883","rs722988") #,"rs61839660" last 3 new (last one not so bad)

## this is a list of snps initially failing manual QC, but later decided were ok
ok <- c("rs3024493","rs61839660","rs56994090","rs12927355","rs151234",
         "rs229533","rs1788097","rs9585056","rs6691977")
# most of these clouds not so bad, or rs9585056 for instance, clouds not great, but replicated in previous
# studies including crohns and UC

qc.cloud.fail <- qc.cloud.fail[!qc.cloud.fail %in% ok]

## SNPs hand checked that pass QC
good.snps.checked <- c("rs193778b","rs6679677","rs2476601","rs11552449","rs2269240","rs12568515","rs7511678","rs3024505","rs2111485",
                       "rs1534422","rs2309837","rs984971","rs3087243","rs62625034","rs75793288","rs6827756","rs6840119","rs4505848",
                       "rs2611215","rs1445898","rs71624119","rs12212193","rs2326451","rs1010474","rs7753008","rs10272724","rs7790800",
                       "rs10758594","rs61839660","rs12416116","rs12722563","rs689","rs3842727","rs7122407","rs2651830","rs11215766",
                       "rs705704","rs4378452","rs917911","rs3176775","seq-t1d-12-56449572-T-C","rs3184504","rs11066320","rs9557217","rs6573857","rs174213",
                       "rs4900384","rs8022656","rs55791667","rs34593439","rs151233","rs9746695","rs12922409","rs12924112","rs4238595",
                       "rs6498184","rs8054218","rs2290438","rs757411","rs2435200","rs1893217","rs17207042","rs74956615","rs12978105",
                       "rs10408526","rs516246","rs6043409","rs2250261","rs9981624","rs1989870","rs4820830","rs229536",
                       "rs1701704","rs10865035","rs113010081","rs11954020","rs12927355","rs1350275","rs151234",
                       "rs1574285","rs1615504","rs16939895","rs17630466","rs3024493","rs653178","rs80054410",
                       "rs35667974","rs8067378","rs11069349","rs72727394","rs8056814","rs1538171","rs722988","rs56994090",
                       "rs229533","rs7960225","rs4767000","rs9585056","rs193778","rs679574")  #rs722988 singleton to check
#rs6691977 slightly dodgy
# a list of bad, then good, where the bad failed QC, then the good was the next hit down on the list that passed QC
badz <- c("imm_19_10324118","rs1016431","rs1160544","rs12927355","rs151234","rs1788097","rs202535","rs229533","rs3024493","rs3826110","rs4759229","rs56994090","rs5763751","rs5763779","rs61839660","rs6691977","rs6907898","rs72729320","rs8074179","rs9585056")
goodz <- c("rs74956615","rs7790800","rs2309837","rs9746695","rs151233","rs17207042","rs2250261","rs229536","rs3024505","rs8054218","rs705704","rs8022656","rs1989870","rs4820830","rs12722496","rs7511678","rs2326451","rs55791667","rs757411","rs9557217")

## ok is a list of snps initially failing manual QC, but later decided were ok

goodz <- goodz[!badz %in% ok]
badz <- badz[!badz %in% ok]


## these were SNPs ollie noticed were in the analysis list and the QC exclusion list, in the 
# end this was because of the underscores and match() function failure
weird <- c("ccc_2_100113638_G_A","ccc_2_100118573_G_C","ccc_2_100126398_C_A",
           "ccc_2_100127085_C_T","ccc_2_100174734_C_T","ccc_2_100175954_C_T",
           "ccc_2_100176881_A_G","ccc_2_100177022_A_G","ccc_22_35871557_C_A",
           "imm_12_54649933","rs1690500","rs181581",
           "rs2304237","rs35397078_rs8078409","rs62212954_rs66733675",
           "rs72687017","rs72687027","rs73966411",
           "seq_NOVEL_11243","seq_NOVEL_11254","seq_NOVEL_11256",
           "seq_NOVEL_3186","seq_NOVEL_7020","seq_NOVEL_99")



# this is a list of SNPs that may be excluded in QC, not to exclude as they are known locii
unexclude <- c("rs1788097","rs9585056","rs13119723","rs9268645","imm_18_65694668","X1kg_13_98879767",
                               "rs3842727","imm_11_2141424","rs12212193","imm_6_91053490",
                                "imm_11_2138800","rs689")   #,"rs34536443")

dif.mis <- c("rs9585056", "X1kg_13_98879767") # have differential missingness but kept as t1d locii

# t1d locii, but appalling signal clouds so must exclude
appalling <- c("rs6916742","rs601338","imm_14_100375798","rs941576","imm_20_1614253","rs202536")

## after inspection of table, these are duplicate regions not removed automatically
not.top.snps <- c("rs11600952","rs7960225","rs4378452","rs4767000")

