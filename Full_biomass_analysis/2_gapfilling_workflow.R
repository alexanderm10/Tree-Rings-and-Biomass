##################################################################################
## Basic Components necessary for data management for doing dendro modelling
##################################################################################

# clear memory
rm(list=ls())

# importing libraries
library(dplR)
library(lattice)

# Getting Libraries
library(reshape)
library(car)
library(mgcv)
library(nlme)
library(lmeSplines)
#library(lme4)
library(splines)
library(MASS)
library(MuMIn)
library(ggplot2)
library(grid)
se <- function(x){
	sd(x, na.rm=TRUE) / sqrt((length(!is.na(x))))}

q.blank <- theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=14, face="bold"), axis.text.y=element_text(color="black", size=12, face="bold"), axis.title.x=element_text(face="bold", size=14),  axis.title.y=element_text(face="bold", size=14))


#################################################################################################
# STEP 1: Gap-filling measured trees
# STEP 1b: Pith-correction in measured trees (useful for stand dynamics; won't do unless you ask for it)
# STEP 2: Gap-filling missing trees
#
# BIG SELLING POINT OF THIS APPROACH: we can quantify different levels of uncertainty & variability
# Caveats: fitting the initial GAMM is very time-intensive (it may take hours with your full data) because current form fits a spline to each tree in a mixed model framework

#################################################################################################
# Previous workflow
# 1) Read in RWL, QA/QC
# 3) Aggregate to tree (factoring in whether cores were dated) level & decide if an entire tree is dated or not -- WRITE THIS AS A FILE!
# 4) Stack RWL & Merge with metadata (becomes "ring.data" or named equivalent)

# Ring.data format: stack all of the core BAI, so that data frame with a SIGNLE BAI column, and then all of the factors in other columns
ring.data <- read.csv("TreeRWL_AllSites_stacked.csv")
ring.data$Tree <- as.factor(ring.data$Tree) 
summary(ring.data)

# Tree Data
tree.data <- read.csv("TreeData.csv")
summary(tree.data)

# Site Data (for year cored) 
Site.data <- read.csv("raw input files/DOE_plus_valles.csv", na.strings="")
Site.data$Year.sample <- as.numeric(substr(Site.data$date.sample,7,10))
summary(Site.data)

# merging in the year sampled into the tree data & calculating age
tree.data <- merge(tree.data, Site.data[,c("PlotID", "Year.sample")], all.x=T, all.y=F)
tree.data$PlotID
tree.data$Age <- tree.data$Year.sample - tree.data$Pith
summary(tree.data)

# We're going to run 2 sets of fillin models:
# 1) Model based on only DATED trees (m1d)
# 2) Model based on both DATED and UNDATED trees (m1u)

trees.dated <- ring.data[ring.data$Dated=="Y","TreeID"]

# using the gamm allows us to fit a smoothing spline to each tree, which allows us to basically gapfill taking into account unique tree temporal trends
#	current spline parameter: shrinkage version of cubic spline with 3 knots (stiff CR spline)
#	when we fit a generalized version for missing trees, we'll have to decide what to fit it to instead of TreeID; I think probably species|plot
# m1 <- gamm(RW ~ s(Year, bs="cs", k=3) + species + DBH..cm., random=list(Site=~1, PlotID=~1, TreeID=~1), data=ring.data, na.action=na.omit)

# ----------------------------------------------------------------
# IDEAL MODEL FORM (it won't work for many reasons)
#	-- won't predict outside range of years observed on each core
#	-- end up with singularity issues
#	-- would take FOREVER to fit even if it did work
# m1d <- gamm(RW ~ s(Year, bs="cs", k=3, by=TreeID) + species*DBH..cm.*canopy.class, random=list(Site=~1, PlotID=~1, TreeID=~1), data=trees.dated.full, na.action=na.omit)
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# The spline doesn't fit outside the range of observed values, so we need to give it a "null" guess
# As a very very rough guess right now, filling missing with the measurement from the oldest ring
# 1) Create rough size-age relationships to give a narrow window of rings to fill (i.e. don't fill a 10 cm oak back to 1900 if our rings stop in 1980)
# 2) fill the modeling window with non-0 values... perhaps mean growth from last decade or past observed trend?

# ---------------------------------------
# 1) Create rough size-age relationships to give a narrow window of rings to fill (i.e. don't fill a 10 cm oak back to 1900 if our rings stop in 1980)
# ---------------------------------------
# What we need: core summary data (DBH..cm., estimated pith date)

# Ignoring all the Sites we don't have and doing some exploratory graphing
tree.data2 <- tree.data[tree.data$PlotID %in% unique(ring.data$PlotID),]
summary(tree.data2)

# Need to remove species for which we have no pith estimates for the time being
Species.pith <- unique(tree.data2[!is.na(tree.data2$Pith), "Species"])
tree.data3 <- tree.data2[tree.data2$Species %in%  Species.pith,]
summary(tree.data3)

qplot(DBH..cm., Age, color=Species, data=tree.data3) + facet_wrap(~Species) +
	stat_smooth(method="lm", alpha=0.5, size=1.5) +
	theme_bw()


# Making a very basic linear model looking at Site-specific species-DBH..cm.-age relationships
dbh.age <- lm(Age ~ Species*DBH..cm.*Site-1, data=tree.data3)
summary(dbh.age)
summary(dbh.age)$r.squared # Note, this very basic model works pretty well!

# Using the prediction interval to get us a higher upper bound
age.pi <- predict(dbh.age, newdata=tree.data3, interval="predict")
summary(age.pi)
dim(age.pi); dim(tree.data3) # Making sure we didn't lose any rows along the way

tree.data3 <- cbind(tree.data3, age.pi)
summary(tree.data3)

# Setting the filling window to the upper p.i. limit
tree.data3$fill.year <- ifelse(is.na(tree.data3$Pith), tree.data3$Year.sample-tree.data3$upr, tree.data3$Pith)
summary(tree.data3)


# Merging this back into a data frame that contains info for all the trees we're modeling right now
tree.data.model <- merge(tree.data3, tree.data2, all.x=T, all.y=T)
summary(tree.data.model)

#--------------------- 
# QUESTION: what to do about the species with no pith estimates?
#	For now, I'm just going to do a species-naive fit within each plot
#--------------------- 
dbh.age.plot <- lm(Age ~ DBH..cm.*PlotID-1, data=tree.data.model)
summary(dbh.age.plot) 

age.pi.plot <- predict(dbh.age.plot, newdata= tree.data.model, interval="predict")
summary(age.pi.plot)
dim(age.pi.plot); dim(tree.data.model)

summary(tree.data.model)
# Filling missing fill.year with one caluclated form age.pi.plot fill 
tree.data.model[is.na(tree.data.model$fill.year),"fill.year"] <- tree.data.model[is.na(tree.data.model$fill.year),"Year.sample"] - age.pi.plot[which(is.na(tree.data.model$fill.year)),3]
summary(tree.data.model)
#--------------------- 
# ---------------------------------------

# ---------------------------------------
# ORIGINAL USAGE: Provide dummy variable
# 2) fill the modeling window with non-0 values... perhaps mean growth from last decade or past observed trend?
#		QUESTION: What dummy value do we want to give it?
#			For now, I'm going to do the mean of the last decade
# ---------------------------------------
# summary(ring.data)

# ring.data$RW0 <- ring.data$RW # Making a column of dummy-filled ring widths
# for(i in unique(ring.data$TreeID)){
	# #------------------------------
	# #year.fill = the oldest year to fill based on above step (prediction interval, size-species-Site relationships)
	# #year.min = the oldest year measured
	# #rw.fill = the dummy ring width to feed the model; right now this is for the last measured decade 
	# #			-- (yr.min + 10) means the 10 years following year.min (
	# #			-- e.g. if min year = 1950, we'll fill with the mean value from 1950:1960
	# #------------------------------
	# yr.fill <- tree.data.model[tree.data.model$TreeID==i,"fill.year"]
	# yr.min <- min(ring.data[ring.data$TreeID==i & !is.na(ring.data$RW), "Year"])
	# rw.fill <- mean(ring.data[ring.data$TreeID==i & ring.data$Year>=yr.min & ring.data$Year<=(yr.min+10), "RW"],na.rm=T)
	# #------------------------------

	# #------------------------------
	# # The actual insertion of the dummy fil value into the fill range
	# #------------------------------
	# ring.data[ring.data$TreeID==i & is.na(ring.data$RW) & ring.data$Year>=yr.fill, "RW0"] <- rw.fill
	# #------------------------------
# }
# summary(ring.data)
# # ---------------------------------------

# ---------------------------------------
# New Usage: delete NAs for years way outside what we think we should be fitting
# 2b) I still don't think we want to be fi
# ---------------------------------------
summary(tree.data.model)
summary(ring.data)
dim(ring.data)

for(i in unique(ring.data$TreeID)){
	#------------------------------
	#year.fill = the oldest year to fill based on above step (prediction interval, size-species-Site relationships)
	#------------------------------
	yr.fill <- tree.data.model[tree.data.model$TreeID==i,"fill.year"]
	#------------------------------

	#------------------------------
	# The actual insertion of the dummy fil value into the fill range
	#------------------------------
	ring.data <- ring.data[!(ring.data$TreeID==i & ring.data$Year<yr.fill),]
	#------------------------------
}
summary(ring.data)
dim(ring.data)
# ---------------------------------------
# ----------------------------------------------------------------

# ################################################################
# ################################################################
# RUNNING THE GAMM!!
# ################################################################
# ################################################################
# A generalized additive mixed model (gamm) allows us to fit splines in a mixed model context
# we can let these splines vary by tree which essentially detrends the core
# here's we're using our dummy-filled ring widths as a response so that the spline will fit over the whole time period of interest
# ################################################################

# ----------------------------------------------------------------
# Gapfilling trees for which we have at least some measurements
# ----------------------------------------------------------------
source("0_gapfill_gamm_function.R") # This is where the generalized function that runs the gamm & saves the diagnostic plots is


# ---------------------------------------
# Loop through by Site using measured trees only
# ---------------------------------------
site.codes <- unique(substr(ring.data$PlotID,1,2))

# Making dead trees have a special canopy class
ring.data$RW.modeled <- NA # making a placeholder vector that we're going to fill in
for(s in site.codes){
	rows.site <- which(substr(ring.data$PlotID,1,2)==s) # figuring out which rows belong to a given site

	data.use <- ring.data[rows.site,] # subset the data for each plot
	data.use <- droplevels(data.use) # git rid of levels for factors like Species that no longer exist

	out.path <- paste0("GapFill_Metadata/", s)
	if(substr(s,1,1)=="V" | substr(s,1,1)=="N"){
		gamm.out <- gapfill.gamm(data= data.use, DBH="DBH..cm.", Species.Use="Species", Canopy.Class="Canopy.Class", canopy=F, out.prefix=out.path)
	} else {
		gamm.out <- gapfill.gamm(data= data.use, DBH="DBH..cm.", Species.Use="Species", Canopy.Class="Canopy.Class", canopy=T, out.prefix=out.path)

	}
	ring.data[rows.site, "RW.modeled"] <- gamm.out$data$RW.modeled
}



# ----------------------------------------------------------------
# Gapfilling trees that we have no measurements for (punky, dead...)
# ----------------------------------------------------------------
# Going to gapfill live trees only (won't try and figure out when dead trees welcome)
# We can't fit a TreeID spline for trees we don't have any measurements for, so we need something more generalizeable
#	The best options are probably species or PlotID; Ross votes PlotID because of variation in the Valles
# ----------------------------------------------------------------
# creating a data frame of blank ring widths to fill with the gamm
trees.missing <- tree.data.model[!(tree.data.model$TreeID %in% unique(ring.data$TreeID)), c("Species", "Site", "PlotID", "TreeID", "DBH..cm.", "Canopy.Class")]
summary(trees.missing) 
dim(trees.missing) # number pf rows = number of missing trees

# Creating a vector of years and dummy ring widths we're going to want to fill
rw.dummy <- data.frame(Year=min(ring.data$Year):max(ring.data$Year), RW=NA)
summary(rw.dummy)

# making a data frame with blank rings widths for the full range of years for missing trees
trees.missing <- merge(trees.missing, rw.dummy, all.x=T, all.y=T)
summary(trees.missing)

# ---------------------------------------
# delete NAs for years way outside what we think we should be fitting (this was done earlier)
# ---------------------------------------
summary(tree.data.model)
summary(trees.missing)
dim(trees.missing)

for(i in unique(trees.missing$TreeID)){
	#------------------------------
	#year.fill = the oldest year to fill based on above step (prediction interval, size-species-Site relationships)
	#------------------------------
	yr.fill <- tree.data.model[tree.data.model$TreeID==i,"fill.year"]
	#------------------------------

	#------------------------------
	# The actual insertion of the dummy fil value into the fill range
	#------------------------------
	trees.missing <- trees.missing[!(trees.missing$TreeID==i & trees.missing$Year<yr.fill),]
	#------------------------------
}
summary(trees.missing)

# ----------------------------------------------------------------
# Creating a gamm for missing trees
#	Because we don't have any ring measurements for missing trees, we can't use a TreeID spline (because we can't predict what we can't fit)
# 	After talking with Ross, we decided to fit the spline by plot, which should get us the general plot dynamics and the hierarchical effects should help fill in the species; The fixed DBH..cm. effect will also help adjust growth based on size; the alterative would be to fit by species or something like that if we think the signal is more regional rather than driven by gap dynamics
# ----------------------------------------------------------------
# We had a couple species for which we have no measured trees, so nothign to train the model on, so we have to call it something else

# putting missing & measured trees together
trees.missing$Measured <- as.factor("NO")
ring.data$Measured     <- as.factor("YES")
data.all <- merge(ring.data, trees.missing, all.x=T, all.y=T)


# Do the missing tree gapfilling based off of dated trees only or species that have no dated trees
species.keep <- unique(data.all$Species)[!(unique(data.all$Species) %in% unique(data.all[!data.all$Dated=="N","Species"]))] # exceptions to the dated only rule
fill.missing <- data.all[!data.all$Dated=="N" | data.all$Species %in% species.keep,] # This will subset completely missing trees and the dataed
fill.missing$Species.Model <- recode(fill.missing$Species, "'POTR'='POGR'") # there are no measurements for these species, so we have to pretend they're something else
summary(fill.missing)

out.path <- paste0("GapFill_Metadata/", "MissingTrees")

# Right now we're going to run this with Canopy off because the Valles trees still don't have a Canopy Class
gamm.missing <- gapfill.gamm(data=fill.missing, DBH="DBH..cm.", Species.Use="Species.Model", Canopy.Class="Canopy.Class", smooth.by="PlotID", canopy=T, out.prefix=out.path)

data.all[data.all$Measured=="NO", "RW.modeled"] <- gamm.missing$data[gamm.missing$data$Measured=="NO", "RW.modeled"]



##################################################################################
# Create a gapfilled data set
# Now we have a lot of different piece meal data sets that we need to put back together to have estimated or measured ring widths for everything
##################################################################################
data.all$RW.gapfilled <- ifelse(is.na(data.all$RW), data.all$RW.modeled, data.all$RW)
summary(data.all)

write.csv(data.all, "RingData_All_Gapfilled.csv", row.names=F)
# ----------------------------------------------------------------



