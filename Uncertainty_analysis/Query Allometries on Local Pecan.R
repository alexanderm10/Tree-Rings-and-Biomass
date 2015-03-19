# Script Querying allometries
setwd("~/Desktop/pecan/modules/allometry/R")
outdir <- "~/Desktop/PalEON CR/Tree Rings/Tree-Rings-and-Biomass/Uncertainty_analysis/AllomFiles"

source("AllomAve.R")
source("query.allom.data.R")
source("allom.BayesFit.R")
source("read.allom.data.R")



# PIPO -- the new species
pfts = list(test = data.frame(spcd=122,acronym="PIPO"))

# Example with just 1 eq
#pfts = list(test = data.frame(spcd=91,acronym="PIAB"))

# Example with multiple Pecan eq
pfts = list(test = data.frame(spcd=202,acronym="PSME"))

AllomAve(pfts,2,outdir="~/Desktop",parm="../data/Table3_GTR-NE-319.v2_RossAdendum.csv",ngibbs=500)