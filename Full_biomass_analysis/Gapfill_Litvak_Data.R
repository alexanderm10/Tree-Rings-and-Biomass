# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# Workflow for modeling generalized tree growth for the Valles Caldera
# Note: This is developed so that we can use our tree rings to gap-fill the trees
#       in data from Marcy Litvak.  This means we have to remove at least PlotID &
#       either generalize species or drop that as well. 
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# Loading in libraries and some other useful things
# --------------------------------------------------------------------------------
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
# -----------------------------------------------


# --------------------------------------------------------------------------------
# Loading in & comparing our tree ring data & Marcy's 
#    data so we know how to the structure the gamm
# --------------------------------------------------------------------------------
# -----------------------------------------------
# Tree Ring Data
# Format: stack all of the core measurements (BAI or RW) to create a single data frame with a single column of measurements and all predictors in other columns
# -----------------------------------------------
ring.data <- read.csv("TreeRWL_AllSites_stacked.csv")
ring.data$Tree <- as.factor(ring.data$Tree) 
ring.data <- ring.data[,!names(ring.data)=="Density..stems.ha."]
summary(ring.data)

# Tree Data
tree.data <- read.csv("TreeData.csv")
tree.data <- tree.data[,!names(tree.data)=="Density..stems.ha."]
summary(tree.data)

# Site Data (for year cored) 
site.data <- read.csv("raw input files/DOE_plus_valles.csv", na.strings="")
site.data$Year.sample <- as.numeric(substr(site.data$date.sample,7,10))
summary(site.data)

# merging in the year sampled into the tree data & calculating age
# also adding the correct stem density here
tree.data <- merge(tree.data, site.data[,c("PlotID", "Year.sample")], all.x=T, all.y=F)
tree.data$PlotID
tree.data$Age <- tree.data$Year.sample - tree.data$Pith
summary(tree.data)

# -----------------------------------------------

# -----------------------------------------------
# subsetting only data for the Valles & Formatting
# -----------------------------------------------
ring.data <- ring.data[substr(ring.data$PlotID,1,1)=="V",]
summary(ring.data)
tree.data <- tree.data[substr(tree.data$PlotID,1,1)=="V",]

# merge in the modern stem density here
tree.data <- merge(tree.data, site.data[,c("PlotID", "Density.Total..stems.ha.")], all.x=T, all.y=F)
ring.data <- merge(ring.data, site.data[,c("PlotID", "Density.Total..stems.ha.")], all.x=T, all.y=F)
tree.data$Density.Total..stems.ha. <- as.numeric(paste(tree.data$Density.Total..stems.ha.))
ring.data$Density.Total..stems.ha. <- as.numeric(paste(ring.data$Density.Total..stems.ha.))
summary(tree.data)
summary(ring.data)

# Recoding the sites to match Marcy's
tree.data$Site <- recode(tree.data$Site, "'Valles Caldera Upper'='PPINE'; 'Valles Caldera Lower'='MCON'")
ring.data$Site <- recode(ring.data$Site, "'Valles Caldera Upper'='PPINE'; 'Valles Caldera Lower'='MCON'")
summary(ring.data)
summary(tree.data)

# Just for reference: find the trees that were dated
trees.dated <- ring.data[ring.data$Dated=="Y","TreeID"]
# -----------------------------------------------

# -----------------------------------------------
# Marcy Litvak's data
#   ***ADD WHATEVER WE NEED TO FOR FORMATTING HERE***
#	***This might include recoding species in both Marcy's & our data
#   ***Make sure site codes line up, etc.
# -----------------------------------------------
marcy.ppine <- read.csv("raw input files/marcy_ppine_2013.csv")
summary(marcy.ppine)
marcy.mcon <- read.csv("raw input files/marcy_mcon_2012.csv")
summary(marcy.mcon)

names(marcy.ppine)
names(marcy.mcon)

# Merging the two into 1 data frame and doing a bit of formatting
marcy <- rbind(marcy.ppine[,names(marcy.ppine) %in%  names(marcy.mcon)], marcy.mcon)
marcy$PlotID <- as.factor(paste(marcy$Site, marcy$Plot_Name, sep="."))
names(marcy)[16] <- "DBH..cm." 
marcy$Tree_Tag_Number <- as.factor(marcy$Tree_Tag_Number)
summary(marcy)


# Finding Marcy's plot Densities
#### ------------
#### ROSS -- Double check that my assumptions about Marcy's sampling design are correect
#### ------------
marcy.density <- aggregate(marcy[,c("Plot_Radius")], by=list(marcy[,"PlotID"]), FUN=length)
names(marcy.density) <- c("PlotID", "n.stems")
marcy.density$Density.Total..stems.ha. <- marcy.density$n.stems/(pi*10^2)*10000 # stems/ha
summary(marcy.density)

# Merging the plot densities back into marcy's data
marcy <- merge(marcy, marcy.density)
summary(marcy)
# -----------------------------------------------

# -----------------------------------------------
# Recode Species
#### ------------
#### ROSS -- I'll leave it to you to decide what to do with these species
#### ------------
# -----------------------------------------------
unique(ring.data$Species) # PIPO, PIEN, PSME
unique(tree.data$Species) # PIPO, PIEN, PSME, POTR
unique(marcy$Species)     # PIPO, PIEN,     , POTR, ABCO


# Just for now I created 3 categories: PIPO, PIEN, and OTHER
ring.data$Species.Model <- recode(ring.data$Species, "'PIPO'='PIPO'; 'PIEN'='PIEN'; 'PSME'='OTHER'; 'POTR'='OTHER'; 'ABCO'='OTHER'")
tree.data$Species.Model <- recode(tree.data$Species, "'PIPO'='PIPO'; 'PIEN'='PIEN'; 'PSME'='OTHER'; 'POTR'='OTHER'; 'ABCO'='OTHER'")
marcy$Species.Model <- recode(marcy$Species, "'PIPO'='PIPO'; 'PIEN'='PIEN'; 'PSME'='OTHER'; 'POTR'='OTHER'; 'ABCO'='OTHER'")
# -----------------------------------------------

# -----------------------------------------------
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# Pre-GAMM formatting
# --------------------------------------------------------------------------------
# -----------------------------------------------
# Create rough size-age relationships to give a narrow window of rings to fill (i.e. don't fill a 10 cm oak back to 1900 if our rings stop in 1980)
# -----------------------------------------------
# Ignoring all the Sites we don't have and doing some exploratory graphing
tree.data2 <- tree.data[tree.data$PlotID %in% unique(ring.data$PlotID),]
summary(tree.data2)

# Need to remove species for which we have no pith estimates for the time being
Species.pith <- unique(tree.data2[!is.na(tree.data2$Pith), "Species"])
tree.data3 <- tree.data2[tree.data2$Species %in%  Species.pith,]
summary(tree.data3)

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


summary(tree.data.model)
summary(ring.data)
dim(ring.data)

# Do the actual trimming of years way outside the time frame we're interested in
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

# -----------------------------------------------

# --------------------------------------------------------------------------------



# --------------------------------------------------------------------------------
# Generate the gamm for gap-filling purposes
# Thought: We could add a factor for initial plot density to get some of the plot effects
#          that were being captured by the random PlotID effect earlier
#
# Note: We need to think about at what level we want to fit the spline
#       -- For the tree ring data, we moved to fitting by plot, but that won't work here
#       -- My reccomendation would be to fit by species or just by site (i.e. PIPO & MCON
#          have different patterns, but within the site there's a common history/dynamic)
#
# Note: In the gamm below, DBH and Density are the modern DBH and density because we can't fit
# --------------------------------------------------------------------------------
# The actual gamm; note: I added a fixed modern density effect to try and account for one of the major differences among plots
gamm.valles.general <- gamm(log(RW) ~ s(Year, bs="cs", k=3, by=Site) + DBH..cm. + Density.Total..stems.ha., random=list(Species.Model=~1, Site=~1), data= ring.data, na.action=na.omit)

# This plots the overall shape of growth curves by site
#     Take a look. This actually makes a good argument for the sites sharing a common history
plot(gamm.valles.general$gam)

# Saving the GAMM From above so we can load it without having to refit it
save(gamm.valles.general, file="GapFilling_gamm_valles_generalized_2015.04.02.rData")

# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Apply the gamm to Mary's data
# --------------------------------------------------------------------------------
# -----------------------------------------------
# Create a blank data set from Marcy's data
# -----------------------------------------------
summary(marcy)
marcy.fill <- marcy[,c("Tree_Tag_Number", "Species.Model", "Site", "PlotID", "DBH..cm.", "Density.Total..stems.ha.")]
summary(marcy.fill)


# Creating a vector of years and dummy ring widths we're going to want to fill
rw.dummy <- data.frame(Year=min(ring.data$Year):max(ring.data$Year), RW=NA)
summary(rw.dummy)

# merge the dummy frame with Marcy's trees
marcy.fill <- merge(marcy.fill, rw.dummy, all.x=T, all.y=T)
summary(marcy.fill)
# -----------------------------------------------

# -----------------------------------------------
# Note: I'm skipping the guessing on age a priori step; see christy_gapfilling.R for details
# -----------------------------------------------

# -----------------------------------------------
# Apply the gamm to Marcy's data
# -----------------------------------------------
marcy.fill$RW <- exp(predict(gamm.valles.general, marcy.fill))
summary(marcy.fill)

write.csv(marcy.fill, "GapFilling_Litvak_Raw.csv", row.names=F)
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# DBH reconstruction on Marcy's data (remove negative DBH)
# --------------------------------------------------------------------------------
summary(marcy)

# Get rid of trees that had no DBH & making year a factor for recasting
marcy.fill <- marcy.fill[!is.na(marcy.fill$RW),]
marcy.fill$Year <- as.factor(marcy.fill$Year)
summary(marcy.fill)

# Some trees are in the dataframe twice for some reason.  
# Going to aggregate to take the mean ring width just to be safe
marcy.fill2 <- aggregate(marcy.fill[,"RW"], by=list(marcy.fill$Tree_Tag_Number, marcy.fill$Year), FUN=mean)
names(marcy.fill2) <- c("Tree_Tag_Number", "Year", "RW")
summary(marcy.fill2)

# Re-casting the data so that trees are columns and years are rows
#   Then deleting the Year column & sorting the data so that most recent is tat the top
marcy.fill3 <- recast(marcy.fill2[,c("Year", "Tree_Tag_Number", "RW")], Year ~ Tree_Tag_Number, FUN=mean)
row.names(marcy.fill3) <- marcy.fill3$Year
marcy.fill3 <- marcy.fill3[,2:ncol(marcy.fill3)]
marcy.fill3 <- marcy.fill3[order(row.names(marcy.fill3), decreasing=T),]
summary(marcy.fill3[,1:20])
marcy.fill3[1:10,1:10]


marcy.dbh <- marcy.fill3
for(j in names(marcy.dbh)){
	# put NA in any rows older than when the tree was last measured
	marcy.fill3[as.numeric(row.names(marcy.fill3))>max(marcy[marcy$Tree_Tag_Number==j, "Year"]),j] <- NA
	marcy.dbh[as.numeric(row.names(marcy.dbh))>max(marcy[marcy$Tree_Tag_Number==j, "Year"]),j] <- NA

	# Put DBH in the year dbh last measured in the dield
	marcy.dbh[as.numeric(row.names(marcy.dbh))==max(marcy[marcy$Tree_Tag_Number==j, "Year"]),j] <- max(marcy[marcy$Tree_Tag_Number==j, "DBH..cm."])
	
	# Doing the actual DBH reconstruction
	for(i in 2:(length(marcy.dbh[,j]))){
		marcy.dbh[i,j] <- ifelse(!is.na(marcy.fill3[i-1,j]), marcy.dbh[i-1,j] - marcy.fill3[i-1,j]*2, marcy.dbh[i,j]) # subtracting the previous year's growth from DBH..cm. to get that year's DBH..cm.
	}
	
	# We don't have any estimated pith dates, so we're just going to stop when DBH goes negative
	marcy.dbh[,j] <- ifelse(marcy.dbh[,j]<0, NA, marcy.dbh[,j]) 
	marcy.fill3[,j] <- ifelse(marcy.dbh[,j]<0, NA, marcy.fill3[,j]) 
}

write.csv(marcy.dbh, "Gapfilling_DBHrecon_Litvak.csv", row.names=T)
write.csv(marcy.fill3, "Gapfilling_RingWidths_Modeled_Final_Litvak.csv", row.names=T)
