###NAMESPACE ADDITIONS###
# Depends: R (>= 2.10), grDevices, graphics, stats, utils, genoset, reader, NCmisc, snpStats
# Imports: 
# Suggests:
# importFrom(proftools, readProfileData, flatProfile)
# importFrom(tools, toHTML)
# importFrom(BiocInstaller, biocVersion)
# import(grDevices, graphics, stats, utils)
###END NAMESPACE###



## options ##
options(chip.info="~/github/iChip/ImmunoChip_ChipInfo_New.RData") # if you can access the file, you won't need to change this path
options(ucsc="hg19") # depends on which analysis, need to set something though
options(save.annot.in.current=1)  # 1 = TRUE, stores annotation in current folder to speed up subsequent lookups


if(Sys.info()[["user"]]=="ncooper")
{
  #source("~/github/plumbCNV/generalCNVFunctions.R", echo=FALSE)
  if(getwd()!= "/home/ncooper"){
    source("/home/ncooper/github/iChip/ChipInfoClass.R", echo=FALSE)
  } else {
    warning("You need to source the file 'ChipInfoClass.R to get some of these functions to work properly")
  }
  internal.analyses <- FALSE
  if(internal.analyses) {
    source("/home/ncooper/github/iChip/firstscriptFunctions.R", echo=FALSE) # only needed for internal analyses
    source('/home/ncooper/github/iChip/specificIFunctions.R', echo=FALSE) # only needed for internal analyses
  }
}

if(getwd()!= "/home/ncooper"){
  must.use.package("snpStats",T)
  must.use.package("reader")
  must.use.package("genoset",T)
}
#


##file includes the generally useful functions: simple.date, out.of, randomize.missing
 # simple.date, out.of should now be in NCmisc

