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
ring.data$tree <- as.factor(ring.data$tree) 
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
# NOTE: for this to work with canopy class, we'll need to figure out what to do about dead trees 
# NOTE: Right now this is set up for each Site separately.  If you want to borrow strength from other Sites to help gap fill certain species, run them together
# ################################################################

# ----------------------------------------------------------------
# Gapfilling trees for which we have at least some measurements
# ----------------------------------------------------------------
# ---------------------------------------
# Morgan Monroe Forest
# ---------------------------------------
ring.data.mmf <- ring.data[substr(ring.data$PlotID,1,2)=="MM",]
ring.data.mmf <- droplevels(ring.data.mmf) 
# Note: the log link with the gaussian family to prevent fitting negative values requires that you can't actually have 0 growth so we're going to make it really really small instead
ring.data.mmf$RW <- ifelse(ring.data.mmf$RW==0, 1e-6, ring.data.mmf$RW)
summary(ring.data.mmf)


gamm.mmf <- gamm(log(RW) ~ s(Year, bs="cs", k=3, by=TreeID) + DBH..cm., random=list(Species=~1, Site=~1, PlotID=~1), data= ring.data.mmf, na.action=na.omit)

# Saving the GAMM From above so we can load it without having to refit it
save(gamm.mmf, file="GapFilling_gamm_mmf_02.2015.rData")

# # This isn't a great way of doing it, but it'll let you see the splines for each of the cores
# par(mfrow=c(4,5), mar=c(2,2,0,0)+0.1)
# plot(gamm.mmf$gam)

# Making Predicted ring widths
ring.data.mmf$RW.modeled <- exp(predict(gamm.mmf, ring.data.mmf))
summary(ring.data.mmf)
write.csv(ring.data.valles, "GapFilling_MMF_Measured_Filled.csv", row.names=F)

# Graphing of the modeled rings we will use to fill the data (note this isn't truncating ones that are past where we think pith actually is)
pdf("gamm_gapfill_mmf.pdf", height=7.5, width=10)
ggplot() + facet_wrap(~Species) +
	geom_path(aes(x=Year, y=RW, color=Species), data=ring.data.mmf, size=0.5) +
	geom_point(aes(x=Year, y=RW.modeled), ring.data.mmf[is.na(ring.data.mmf$RW),], size=0.5) +
	scale_y_continuous(limits=c(0,1.25)) +
	theme_bw()
dev.off()
# ---------------------------------------


# ---------------------------------------
# Valles Caldera (Upper & Lower modeled together)
# ---------------------------------------
ring.data.valles <- ring.data[substr(ring.data$PlotID,1,2)=="VL" | substr(ring.data$PlotID,1,2)=="VU",]
ring.data.valles <- droplevels(ring.data.valles)
ring.data.valles$RW <- ifelse(ring.data.valles$RW==0, 1e-6, ring.data.valles$RW)
summary(ring.data.valles)

gamm.valles <- gamm(log(RW) ~ s(Year, bs="cs", k=3, by=TreeID) + DBH..cm., random=list(Species=~1, Site=~1, PlotID=~1), data= ring.data.valles, na.action=na.omit)

# Saving the GAMM From above so we can load it without having to refit it
save(gamm.valles, file="GapFilling_gamm_valles_02.2015.rData")

# # This isn't a great way of doing it, but it'll let you see the splines for each of the cores
# par(mfrow=c(4,5), mar=c(2,2,0,0)+0.1)
# plot(gamm.valles$gam)

ring.data.valles$RW.modeled <- exp(predict(gamm.valles, ring.data.valles))
summary(ring.data.valles)
write.csv(ring.data.valles, "GapFilling_Valles_Measured_Filled.csv", row.names=F)

# VERY rough graphing of the modeled rings we will use to fill the data
pdf("gamm_gapfill_valles.pdf", height=7.5, width=10)
ggplot() + facet_grid(Species~Site) +
	geom_path(aes(x=Year, y=RW, color=Species), data=ring.data.valles, size=0.25) +
	geom_point(aes(x=Year, y=RW.modeled), ring.data.valles[is.na(ring.data.valles$RW),], size=0.5) +
	scale_y_continuous(limits=c(0,1.25)) +
	theme_bw()
dev.off()
# ---------------------------------------
# ----------------------------------------------------------------


# ################################################################
# ################################################################

# ----------------------------------------------------------------
# Gapfilling trees that we have no measurements for (punky, dead...)
# ----------------------------------------------------------------
# Going to gapfill live trees only (won't try and figure out when dead trees welcome)
# We can't fit a TreeID spline for trees we don't have any measurements for, so we need something more generalizeable
#	The best options are probably species or PlotID; Ross votes PlotID because of variation in the Valles
# ----------------------------------------------------------------
# Ring.data format: stack all of the core BAI, so that data frame with a SIGNLE BAI column, and then all of the factors in other columns
ring.data <- read.csv("TreeRWL_AllSites_stacked.csv")
ring.data$tree <- as.factor(ring.data$tree) 
summary(ring.data)

# Tree Data
tree.data <- read.csv("TreeData.csv")
summary(tree.data)

# Site Data (for year cored) 
Site.data <- read.csv("input files/DOE_plus_valles.csv", na.strings="")
Site.data$Year.sample <- as.numeric(substr(Site.data$date.sample,7,10))
summary(Site.data)

# merging in the year sampled into the tree data & calculating age
tree.data <- merge(tree.data, Site.data[,c("PlotID", "Year.sample")], all.x=T, all.y=F)
tree.data$Age <- tree.data$Year.sample - tree.data$Pith
summary(tree.data)

# ---------------------------------------
# ---------------------------------------
# NOTE: Re-run the size-age section above so that we can create a data frame of NAs to fill for the missing trees
# ---------------------------------------
# ---------------------------------------

# creating a data frame of blank ring widths to fill with the gamm
summary(tree.data.model)
trees.missing <- tree.data.model[!(tree.data.model$TreeID %in% unique(ring.data$TreeID)), c("Species", "Site", "PlotID", "TreeID", "DBH..cm.")]
summary(trees.missing)# 29 missing trees: 4 from MM, 22 VUF, 3 VLF
dim(trees.missing)

# Creating a vector of years and dummy ring widths we're going to want to fill
rw.dummy <- data.frame(Year=min(ring.data$Year):max(ring.data$Year), RW=NA)
summary(rw.dummy)

trees.missing <- merge(trees.missing, rw.dummy, all.x=T, all.y=T)
summary(trees.missing)

# ---------------------------------------
# delete NAs for years way outside what we think we should be fitting
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
# 	After talking with Ross, we decided to fit the spline by plot, which should get us the general plot dynamics and the hierarchical effects should help fill in the species; The fixed DBH..cm. effect will also help adjust growth based on size
# ----------------------------------------------------------------
# We're just fitting a single model for all missing trees, but it'll be okay since we're doing the spline by PlotID and have the random effects structure; this is what we could do for everything if the spline fitting by tree didn't take so long

# the log RW doesn't like 0s, so lets just make 0 really small
ring.data$RW <- ifelse(ring.data$RW==0, 1e-6, ring.data$RW)
summary(ring.data)

gamm.missing <- gamm(log(RW) ~ s(Year, bs="cs", k=3, by=PlotID) + DBH..cm., random=list(Species=~1, Site=~1, PlotID=~1), data= ring.data, na.action=na.omit)

par(mfrow=c(3,3))
plot(gamm.missing$gam)

trees.missing$RW.modeled <- exp(predict(gamm.missing, trees.missing))
summary(trees.missing)
write.csv(trees.missing, "GapFilling_Missing_Filled.csv", row.names=F)

# ----------------------------------------------------------------

##################################################################################
# Create a gapfilled data set
##################################################################################
# ---------------------------------------
# Data Formatting
# ---------------------------------------
# putting the Sites back into 1 file
ring.data.valles <- read.csv("GapFilling_Valles_Measured_Filled.csv")
ring.data.mmf <- read.csv("GapFilling_MMF_Measured_Filled.csv")
ring.data.missing <- read.csv("GapFilling_Missing_Filled.csv")

dim(ring.data.valles); dim(ring.data.mmf)
ring.data <- rbind(ring.data.valles, ring.data.mmf)
ring.data <- merge(ring.data, ring.data.missing, all.x=T, all.y=T)

# putting the modeled data where there are no ring widths (will cause negative DBH..cm.)
ring.data$RW.gapfilled <- ifelse(is.na(ring.data$RW), ring.data$RW.modeled, ring.data$RW)
summary(ring.data)

# making a data frame with trees as columns and years as ros
ring.data$Year <- as.factor(ring.data$Year)
trees.gapfilled <- recast(ring.data[,c("Year", "TreeID", "RW.gapfilled")], Year ~ TreeID)
summary(trees.gapfilled)

row.names(trees.gapfilled) <- trees.gapfilled$Year
trees.gapfilled <- trees.gapfilled[,2:ncol(trees.gapfilled)]
trees.gapfilled[(nrow(trees.gapfilled)-10):nrow(trees.gapfilled), 1:10]
# ---------------------------------------


##################################################################################
# DBH Reconstruction
##################################################################################

# ---------------------------------------
# DBH..cm. Reconstruction
# ---------------------------------------
# Tree Data
tree.data <- read.csv("TreeData.csv")
summary(tree.data)

# Site Data (for year cored) 
Site.data <- read.csv("input files/DOE_plus_valles.csv", na.strings="")
Site.data$Year.sample <- as.numeric(substr(Site.data$date.sample,7,10))
summary(Site.data)

# merging in the year sampled into the tree data & calculating age
tree.data <- merge(tree.data, Site.data[,c("PlotID", "Year.sample")], all.x=T, all.y=F)
tree.data$Age <- tree.data$Year.sample - tree.data$Pith
summary(tree.data)

core.data <- read.csv("Core_data_DOE_summer_2014.csv", na.strings=c("", "NA", "#VALUE!", "*"), header=T)
#adding a column include which plot at the Site the trees belong to
names(core.data)
core.data$plot <- substr(core.data$plot.id, 3, 3)
core.data$plot <- as.factor(core.data$plot)
summary(core.data)

# Ordering the data
trees.gapfilled <- trees.gapfilled[order(row.names(trees.gapfilled), decreasing=T),order(names(trees.gapfilled))]
trees.gapfilled[1:10, 1:10]
trees.gapfilled[1:10, (ncol(trees.gapfilled)-10):ncol(trees.gapfilled)]

dbh.recon <- trees.gapfilled
trees.check <- vector() # trees with negative DBH..cm.
for(j in names(dbh.recon)){

	# Step 1: Replace filled years beyond the year in which a tree was sampled with NA (both trees.gapfilled & DBH..cm. recon); 
	# 	Gapfilled: filling the years where I changed 0 to 1e-6 back to 0
	#	DBH..cm.recon: fill year of sample with DBH..cm. when sampled
	trees.gapfilled[as.numeric(row.names(trees.gapfilled))>tree.data[tree.data$TreeID==j, "Year.sample"],j] <- NA
	trees.gapfilled[,j] <- ifelse(trees.gapfilled[,j]==1e-6, 0, trees.gapfilled[,j])

	dbh.recon[as.numeric(row.names(dbh.recon))>tree.data[tree.data$TreeID==j, "Year.sample"],j] <- NA
	dbh.recon[as.numeric(row.names(dbh.recon))==tree.data[tree.data$TreeID==j, "Year.sample"],j] <- tree.data[tree.data$TreeID==j, "DBH..cm."]
	
	# Doing the DBH..cm. reconstruction	
	for(i in 2:(length(dbh.recon[,j]))){
		dbh.recon[i,j] <- ifelse(!is.na(trees.gapfilled[i-1,j]), dbh.recon[i-1,j] - trees.gapfilled[i-1,j]*2, dbh.recon[i,j]) # subtracting the previous year's growth from DBH..cm. to get that year's DBH..cm.
	}
	
	# Get rid of DBH..cm. past the guestimated pith dates -- both DBH..cm. recon & gapfilled
	if(!is.na(tree.data[tree.data$TreeID==j, "Pith"])){
		DBH..cm..recon[as.numeric(row.names(DBH..cm..recon))<tree.data[tree.data$TreeID==j, "Pith"],j] <- NA
		trees.gapfilled[as.numeric(row.names(trees.gapfilled))<tree.data[tree.data$TreeID==j, "Pith"],j] <- NA
	} else 
	if(is.na(tree.data[tree.data$TreeID==j, "Pith"])){ # Get rid of negative modeled DBH..cm.
		dbh.recon[,j] <- ifelse(dbh.recon[,j]<0, NA, dbh.recon[,j]) 
		trees.gapfilled[,j] <- ifelse(dbh.recon[,j]<0, NA, trees.gapfilled[,j]) 		
	}
	# also getting rid of these rings in the stacked ring data too
	years.na <- row.names(trees.gapfilled)[which(is.na(trees.gapfilled[,j]))]
	ring.data[ring.data$TreeID==j & ring.data$Year %in% years.na,"RW.gapfilled"] <- NA

	if(min(DBH..cm..recon[,j], na.rm=T)<0) trees.check <- c(trees.check, j)
}
DBH..cm..recon[1:20, 1:10]
DBH..cm..recon[1:20, (ncol(DBH..cm..recon)-20):ncol(DBH..cm..recon)]
min(DBH..cm..recon, na.rm=T)
trees.check
summary(DBH..cm..recon[,trees.check])

# for trees with negative DBH..cm., working from the inside out
	# If a tree has a negative DBH..cm., we're just going to add from the inside out (at time )
for(j in trees.check){
	if(min(DBH..cm..recon[,j], na.rm=T)<0){
		DBH..cm..recon[,j] <- trees.gapfilled[,j]
		for(i in (nrow(DBH..cm..recon)-1):1){
			DBH..cm..recon[i,j] <- sum(DBH..cm..recon[i+1, j], trees.gapfilled[i,j]*2, na.rm=T)
			}
	}
}
min(DBH..cm..recon, na.rm=T)
DBH..cm..recon[1:10,trees.check]
summary(DBH..cm..recon[, trees.check])
#trees.gapfilled[, trees.check]
tree.data[tree.data$TreeID %in% trees.check,]
# ---------------------------------------
write.csv(dbh.recon, "GapFilling_DBHrecon_ALL.csv", row.names=T)
write.csv(trees.gapfilled, "GapFilling_RingWidths_ALL.csv", row.names=T)

summary(ring.data)
write.csv(ring.data, "TreeRWL_AllSites_stacked_gapfilled_ALL.csv", row.names=F)

##################################################################################
# see next script for reconstructing basal area of trees with no samples 
##################################################################################
