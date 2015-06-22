# Transforming DBH reconstructions from Tree-ring Widths into biomass (kg/m2)

library(dplR)
library(ggplot2)
se <- function(x){
  sd(x, na.rm=TRUE) / sqrt((length(!is.na(x))))}


# Load in diameter reconstrcutions generated from step 3.

g.filled.diam <- read.csv("processed_data/GapFilling_DBHrecon_ALL.csv", header=T, row.names=1)
summary(g.filled.diam)

# read in tree data
tree.data <- read.csv("processed_data/TreeData.csv", header=T)
summary(tree.data)

#select sites for which to run the biomass reconstruction
# NOTE: the North Carolina coast doesn't have a density estimate, so we're going to ignore it for now
trees.use <- tree.data[tree.data$TreeID %in% names(g.filled.diam),] # If you want to do this later, it'll be a special case
# trees.use <- tree.data
summary(trees.use)

plot.data <- read.csv("raw input files/DOE_plus_Valles.csv")
plot.data$Year.Sample <- as.numeric(substr(plot.data$date.sample,7,10))
summary(plot.data)


##########################################################################
# Allometric Equations
##########################################################################


#Convert to biomass with the allometric equation
#using the PECAN generated bayesian equations
library(car)
load("processed_data/allometries_list.Rdata")

# Find out what species we don't have allometries for that need to be renames
unique(trees.use$Species)[!(unique(trees.use$Species) %in% names(allometries))]

trees.use$spp.allom <- recode(trees.use$Species, " 'PIEN'='picea.sp'; 'FRAX'='FRAM'; 'ASTR'='e.hard'; 'PRSE'='e.hard'; 'ULRU'='e.hard'")
summary(trees.use)
plots <- unique(trees.use$PlotID) # You had the right idea, but it was throwing errors because you were trying to evaluate plots you haven't gotten to yet


# making the generalized allometry function to run
# log(AGB) = mu0 + mu1*log(DBH) --equaton form of PECAN allometrics
allom.eq <- function(mu0, mu1, DBH) { exp(mu0 + mu1 * log(DBH) )}

# A temporary allometry file that will get filled
allom.temp <- g.filled.diam
allom.temp[,] <- NA

# dbh=0 causes problems, so we're going to make those NA
g.filled.diam[g.filled.diam==0] <- 1e-6
min(g.filled.diam, na.rm=T)
summary(g.filled.diam)
dim(g.filled.diam)

# This is where the PLOT-LEVEL allometry will be stored
bm.array <- array(NA, dim=c(nrow(g.filled.diam), length(unique(trees.use$PlotID)), nrow(allometries[[1]])))
row.names(bm.array) <- row.names(g.filled.diam)
dimnames(bm.array)[[2]] <- unique(trees.use$PlotID)
summary(bm.array[,,1])
#--------------------------------------------------
# INSERT i LOOP HERE to go through each iteration of randomness from MCMC
# This is one big loop that goes through each layer of the 500 iterations
#--------------------------------------------------
for(i in 1:nrow(allometries[[1]])){
  
# Species loop for calculating tree biomass
for(j in unique(trees.use$spp.allom)){
  cols <- which(names(g.filled.diam) %in% trees.use[trees.use$spp.allom==j, "TreeID"])

  # Pulling coefficients from the randomly pulled estimates from Pecan; 
  # Need to use Bg0 because the heierarchical means were being weird
  mu0=allometries[[j]][i,"Bg0"]
  mu1=allometries[[j]][i,"Bg1"]
  allom.temp[,cols] <- allom.eq(mu0=mu0, mu1 = mu1, DBH = g.filled.diam[,cols])
}
# summing to the plot level

allom.temp[is.na(allom.temp)] <- 0

# biomass loop for summing trees to plots in kg/m2
for(p in 1:length(plots)){
  cols <- which(names(allom.temp) %in% trees.use[trees.use$PlotID==plots[p], "TreeID"])
  if(substr(plots[p],1,1)=="V"){
    bm.array[,p,i] <- rowMeans(allom.temp[,cols])*plot.data[plot.data$PlotID==paste(plots[p]), "Density.Total..stems.ha."]/10000 #mean tree * trees/ha (do for Valles only bc sum of trees != plot density; different sampling method than Neil)
  } else {
    temp <- allom.temp[,cols]
    for(n in names(temp)){ # Convert biomass/tree to biomass/ha
      temp[,n] <- temp[,n] * tree.data[tree.data$TreeID==n,"Density..stems.ha."]/10000
      }
    bm.array[,p,i] <- rowSums(temp) #sum biomass/ha
    }
}
}
#--------------------------------------------------

#bm.array[,,1]
summary(bm.array[,,1])


plots <- paste(unique(tree.data[tree.data$TreeID %in% names(g.filled.diam),"PlotID"]))
years <- as.numeric(row.names(g.filled.diam))
for(p in 1:length(plots)){ # will go by plot
	yr.sample <- plot.data[plot.data$PlotID==plots[p], "Year.Sample"]
	if(max(years) > yr.sample){ # if tree wasn't cored in the last year of sampling, make it's outer years NA
	  rows.na <- which(years > yr.sample)
	  bm.array[rows.na,p,] <- NA
	}
}
summary(bm.array[,,1])

save(bm.array, file="processed_data/Plot_Biomass_Array.Rdata")


### OUTSIDE of all LOOPs (iteration + species + plots)
# You should now have a 3-dimensional array with plots as columns, years as rows, and iterations as layers
biomass.list <- list(mean.biom=data.frame(apply(bm.array[,,], c(1,2), mean, na.rm=T)), 
				     ci.lower =data.frame(apply(bm.array[,,], c(1,2), quantile, 0.025, na.rm=T)),
				     ci.upper =data.frame(apply(bm.array[,,], c(1,2), quantile, 0.975, na.rm=T)))

summary(biomass.list$mean.biom)

biomass <- stack(biomass.list$mean.biom)[,c(2,1)]
names(biomass) <- c("PlotID", "Mean.Biom")
biomass$Year <- as.numeric(paste(row.names(bm.array))) 
biomass$Plot <- as.factor(substr(biomass$PlotID, 3,3))
biomass$Site <- as.factor(substr(biomass$PlotID, 1,2))
biomass$ci.lower <- stack(biomass.list$ci.lower)[,1]
biomass$ci.upper <- stack(biomass.list$ci.upper)[,1]
summary(biomass)

write.csv(biomass, "processed_data/Plot_Biomass_Summary.csv")



pdf("figures/Biomass_AllPlots.pdf")
ggplot(biomass[,]) + facet_wrap(~Site, scales="fixed") +
  geom_ribbon(aes(x=Year, ymin=ci.lower, ymax=ci.upper, fill=Plot), alpha=0.5) +
  geom_line(aes(x=Year, y=Mean.Biom, color=Plot)) +
  theme_bw()
dev.off()




