library(ggplot2)

# exploratory analysis of Marcy Litvak's Valles Caldera Data
ross.trees <- read.csv("raw input files/tree_metadata_DOE_plus_valles.csv", na.strings=c("", "NA", "#VALUE!", "*"),header=T)

# subsetting only the valles trees from the overall trees list
summary(ross.trees)
ross.valles <- ross.trees[substr(ross.trees$PlotID,1,1)=="V",]
summary(ross.valles)

# reducing the number of columns in the data frame
ross.valles<- ross.valles[,c("TreeID", "PlotID", "Species","DBH..cm.", "Density..stems.ha.")]
names(ross.valles)<- c("tree.id", "plot.id", "species", "dbh", "density")
ross.valles$density <- ross.valles$density/10000
summary(ross.valles)

# settign site as a factor
ross.valles$site<- as.factor(substr(ross.valles$plot.id, 1, 2))
summary(ross.valles)

ross.pipo <- ross.valles[substr(ross.valles$plot.id,1,2)=="VL",]
ross.mcon <- ross.valles[substr(ross.valles$plot.id,1,2)=="VU",]
summary(ross.pipo)
summary(ross.mcon)

# read in Marcy's tree tree data
marcy.ppine <- read.csv("raw input files/marcy_ppine_2013.csv")
summary(marcy.ppine)
marcy.mcon <- read.csv("raw input files/marcy_mcon_2012.csv")
summary(marcy.mcon)

names(marcy.ppine)
names(marcy.mcon)

# Merging the two into 1 data frame and doing a bit of formatting
marcy <- rbind(marcy.ppine[,names(marcy.ppine) %in%  names(marcy.mcon)], marcy.mcon)
# marcy$PlotID <- as.factor(paste(substr(marcy$Site, 1,3), marcy$Transect, sep="."))
marcy$PlotID <- as.factor(paste(substr(marcy$Site, 1,3), marcy$Plot_Name, sep="."))
names(marcy)[16] <- "DBH..cm." 
marcy$Tree_Tag_Number <- as.factor(paste0("X", marcy$Tree_Tag_Number))
#marcy <- marcy[!marcy$Tree_Tag_Number=="X476",]
summary(marcy)

#### ------------
marcy.density <- aggregate(marcy[,c("Plot_Radius")], by=list(marcy[,"PlotID"]), FUN=length)
names(marcy.density) <- c("PlotID", "n.stems")
marcy.density$Density.Total..stems.m2. <- marcy.density$n.stems/(pi*10^2)
summary(marcy.density)

# Merging the plot densities back into marcy's data
marcy <- merge(marcy, marcy.density)
summary(marcy)

marcy.valles <- marcy[, c("Tree_Tag_Number", "PlotID", "Species", "DBH..cm.", "Density.Total..stems.m2.", "Site")]
names(marcy.valles) <- c("tree.id", "plot.id", "species", "dbh", "density", "site")
summary(marcy.valles)
marcy.valles$year <- ifelse(marcy.valles$site == "PPINE", 2013, 2012)

save(marcy.valles, file="marcy_valles_metadata.Rdata")
dim(marcy.valles)



# Creating one object with Ross and Marcy data
all.valles<-merge(marcy.valles, ross.valles, all.x=T, all.y=T)
summary(all.valles)




all.valles$site<- ifelse(substr(all.valles$site, 1, 2)=="VL", "ppine.ross", 
                         ifelse(substr(all.valles$site,1,2)=="VU", "mcon.ross", 
                                ifelse(all.valles$site=="PPINE", "ppine.marcy", "mcon.marcy")))
all.valles$site<-as.factor(all.valles$site)
summary(all.valles)

dbh.bins1 <- seq(0, max(all.valles$dbh, na.rm=T), 2)


qplot(x=dbh, data=all.valles, geom="histogram", breaks=dbh.bins1, fill=species)+
  facet_grid(site ~ .)+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_continuous(name="DBH")+scale_y_continuous(name="Stems") + ggtitle("Size Distribution") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))

qplot(x=dbh, data=all.valles, geom="histogram", breaks=dbh.bins1, fill=species, weight=density)+
  facet_grid(site ~ .)+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_continuous(name="DBH")+scale_y_continuous(name="Stems/m2") + ggtitle("Size Distribution weighted by Density") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))

qplot(x=dbh, data=all.valles, geom="histogram", breaks=dbh.bins1, fill=species, weight=density)+
  facet_grid(site ~ .)+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_continuous(name="DBH")+scale_y_continuous(name="Stems/m2") + ggtitle("Size Distribution weighted by Density") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))

# -----------------------
# Ignoring species for a minute
unique(all.valles$density)
valles.plot <- aggregate(all.valles$density, by=list(all.valles$site, all.valles$plot.id), FUN=mean, na.rm=T)
names(valles.plot) <- c("site", "plot.id", "density")
unique(valles.plot$density)

valles.plot

valles.site <- aggregate(valles.plot$density, by=list(valles.plot$site), FUN=mean)
names(valles.site) <- c("site", "density")
valles.site

ggplot()+
  geom_histogram(data=valles.site, aes(x=site, weight=density, fill=site)) +
  geom_histogram(data=valles.plot[valles.plot$plot.id=="VLB",], aes(x=site, weight=density))+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_discrete(name="Site")+scale_y_continuous(name="Stems/m2") + ggtitle("Size Distribution weighted by Density; VLB overlay in black") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))
  
# -----------------------

all.valles$dbh.2 <- round(all.valles$dbh)
summary(all.valles)


valles.plot2 <- aggregate(all.valles$density, by=list(all.valles$site, all.valles$plot.id, all.valles$species, all.valles$dbh.2), FUN=mean, na.rm=T)
names(valles.plot2) <- c("site", "plot.id", "species", "dbh", "density")
summary(valles.plot2)
unique(valles.plot2$density)

library(car)
valles.plot2$site2 <- as.factor(ifelse(substr(valles.plot2$site,1,4)=="mcon", "mcon", "ppine"))

unique(valles.plot2$plot.id) 

valles.plot2$plot.id2 <- as.factor(ifelse(substr(valles.plot2$plot.id,1,3)=="PPI", paste("Marcy", substr(valles.plot2$plot.id,5,5), sep = "."),
                                ifelse(substr(valles.plot2$plot.id,1,3)=="MCO", paste("Marcy",substr(valles.plot2$plot.id,5,5), sep="."),
                                       ifelse(substr(valles.plot2$plot.id,3,3)=="A","Ross.A", "Ross.B"))))


summary(valles.plot2) 

unique(valles.plot2$density)

qplot(x=dbh, data=valles.plot2, geom="histogram", breaks=dbh.bins1, fill=species, weight=density)+
  facet_grid(plot.id2 ~ site2)+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_continuous(name="DBH")+scale_y_continuous(name="Stems/m2") + ggtitle("Size Distribution weighted by Density") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))


valles.site2 <- aggregate(valles.plot2$density, by=list(valles.plot2$site, valles.plot2$species, valles.plot2$dbh), FUN=mean)
names(valles.site2) <- c("site", "species", "dbh", "density")
unique(valles.site2$density)
summary(valles.site2)

qplot(x=dbh, data=valles.site2, geom="histogram", breaks=dbh.bins1, fill=species, weight=density)+
  facet_grid(site ~ .)+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_continuous(name="DBH")+scale_y_continuous(name="Stems/m2") + ggtitle("Size Distribution weighted by Density") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))
