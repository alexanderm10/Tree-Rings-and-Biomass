####################################################
# Script to Query Allometry from local Pecan R scripts (not Pecan VM)
# Christy Rollinson, crollinson@gmail.com
# --------------------------------------------------
# Notes: 
# CR has added a feature in Pecan to eliminate equations or change weighting 
#   of equations based on diameter distribution of interest
# As of 24 March 2015, this version is not on the Pecan Mainline, but can be pulled 
#   for CR's github repository: https://github.com/crollinson/pecan
# See Pecan_Size_Testing folder for updated scripts to pull allom data using this new method
####################################################

# Script Querying allometries
setwd("~/Desktop/pecan/modules/allometry/R")

outdir <- "~/Dropbox/PalEON CR/Tree Rings/Tree-Rings-and-Biomass/Uncertainty_analysis/AllomFiles" # CR Office
#outdir <- "~/Desktop/PalEON CR/Tree Rings/Tree-Rings-and-Biomass/Uncertainty_analysis/AllomFiles" # CR Laptop

source("AllomAve.R")
source("query.allom.data.R")
source("allom.BayesFit.R")
source("read.allom.data.R")


# # 
# # PIPO -- the new species
# # pipo = list(test = data.frame(spcd=122,acronym="PIPO"))

# # Example with just 1 eq
# #pfts = list(test = data.frame(spcd=91,acronym="PIAB"))

# # Example with multiple Pecan eq
# #pfts = list(test = data.frame(spcd=202,acronym="PSME"))

# # AllomAve(pipo,2,outdir=file.path(outdir, "TysonOnly"),parm="../data/Table3_GTR-NE-319.v2_RossAdendum.csv",ngibbs=500)


# Example with just 1 eq
picea = list(Picea.Bad.Eq = data.frame(spcd=90,acronym="PICEA"))

AllomAve(picea,2,outdir=file.path(outdir, "Picea"),parm="../data/Table3_GTR-NE-319.v2_RossAdendum.csv",ngibbs=500)

AllomAve(picea,2,outdir=file.path(outdir, "Picea"),parm="../data/Table3_GTR-NE-319.v2_RossAdendum.csv",ngibbs=500, dmin=2, dmax=250)


# # # Species level equation "if available" for Ross and Marcy Valles Data'
# valles.sp = list(PIPO = data.frame(spcd=122, acronym="PIPO"),
              # PIEN = data.frame(spcd=93, acronym="PIEN"),
              # PSME = data.frame(spcd=202, acronym="PSME"),
              # ABCO = data.frame(spcd=15, acronym = "ABCO"))
# AllomAve(valles.sp,2,outdir=outdir,parm="../data/Table3_GTR-NE-319.v2_RossAdendum.csv",ngibbs=1000)

# valles.sp2 = list(PIPO = data.frame(spcd=122, acronym="PIPO"))
# AllomAve(valles.sp2,2,outdir= outdir,parm="../data/Table3_GTR-NE-319.v2_RossAdendum.csv",ngibbs=1000)

# # Genus level equations for Ross's valles data'
# valles.genus = list(pinus.sp = data.frame(spcd=100, acronym="PINUS"),
                 # picea.sp = data.frame(spcd=90, acronym="PICEA"),
                 # abies.sp = data.frame(spcd=10, acronym="ABIES"))
# AllomAve(valles.genus,2,outdir= outdir,parm="../data/Table3_GTR-NE-319.v2_RossAdendum.csv",ngibbs=1000)

# valles.genus2 = list(picea.sp = data.frame(spcd=90, acronym="PICEA"))
# AllomAve(valles.genus2,2,outdir= outdir,parm="../data/Table3_GTR-NE-319.v2_RossAdendum.csv",ngibbs=1000)

# #making up our own PFT's, just need a species--query simon's database adn cross list based on christy's code
# #pfts = list(PFT.TITLE = data.frame(spcd=PFT.TITLE$spcd,acronym=PFT.TITLE$acronym))

