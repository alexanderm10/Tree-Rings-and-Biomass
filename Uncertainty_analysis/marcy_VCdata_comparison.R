library(ggplot2)

# exploratory analysis of Marcy Litvak's Valles Caldera Data

marcy.pipo <- read.csv("raw input files/marcy_ppine_2013.csv", header=T)
marcy.mcon <- read.csv("raw input files/marcy_mcon_2012.csv", header=T)
ross.trees <- read.csv("raw input files/tree_metadata_DOE_plus_valles.csv", na.strings=c("", "NA", "#VALUE!", "*"),header=T)

# subsetting only the valles trees from the overall trees list
summary(ross.trees)
ross.valles <- ross.trees[substr(ross.trees$PlotID,1,1)=="V",]
summary(ross.valles)

# reducing the number of columns in the data frame
ross.valles<- ross.valles[,c("TreeID", "PlotID", "Species","DBH..cm.", "Density..stems.ha.")]
names(ross.valles)<- c("tree.id", "plot.id", "species", "dbh", "density")
summary(ross.valles)

# settign site as a factor
ross.valles$site<- as.factor(substr(ross.valles$plot.id, 1, 2))
summary(ross.valles)

ross.pipo <- ross.valles[substr(ross.valles$plot.id,1,2)=="VL",]
ross.mcon <- ross.valles[substr(ross.valles$plot.id,1,2)=="VU",]
summary(ross.pipo)
summary(ross.mcon)

# Calculating density per plot for Marcy's data; she used 10m plots
summary(marcy.pipo)

marcy.pipo$density <- 1/(pi*marcy.pipo$Plot_Radius^2)/1e-4
marcy.pipo$plot <- as.factor(paste(marcy.pipo$Site, marcy.pipo$Plot_Name, sep=""))

summary(marcy.pipo)


marcy.mcon$density <- 1/(pi*marcy.mcon$Plot_Radius^2)/1e-4
marcy.mcon$plot <- as.factor(paste(marcy.mcon$Site, marcy.mcon$Plot_Name, sep=""))
summary(marcy.mcon)
#now density is in trees per hectare to match what christy calculated in the DOE spreadsheets.

marcy.valles<- merge(marcy.pipo, marcy.mcon, all.x=T, all.y=T)
summary(marcy.valles)

marcy.valles <- marcy.valles[,c("Site", "Tree_Tag_Number","Species", "DBH_1", "density", "plot")]
summary(marcy.valles)
names(marcy.valles)<- c("site", "tree.id", "species", "dbh", "density", "plot.id")
summary(marcy.valles)

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
  scale_x_continuous(name="DBH")+scale_y_continuous(name="Stems/ha") + ggtitle("Size Distribution weighted by Density") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))

qplot(x=dbh, data=all.valles, geom="histogram", breaks=dbh.bins1, fill=species, weight=density)+
  facet_grid(site ~ .)+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_continuous(name="DBH")+scale_y_continuous(name="Stems/ha") + ggtitle("Size Distribution weighted by Density") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))

# -----------------------
# Ignoring species for a minute
valles.plot <- aggregate(all.valles$density, by=list(all.valles$site, all.valles$plot.id), FUN=sum, na.rm=T)
names(valles.plot) <- c("site", "plot.id", "density")
valles.plot

valles.site <- aggregate(valles.plot$density, by=list(valles.plot$site), FUN=mean)
names(valles.site) <- c("site", "density")
valles.site

ggplot()+
  geom_histogram(data=valles.site, aes(x=site, weight=density, fill=site)) +
  geom_histogram(data=valles.plot[valles.plot$plot.id=="VLB",], aes(x=site, weight=density))+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_discrete(name="Site")+scale_y_continuous(name="Stems/ha") + ggtitle("Size Distribution weighted by Density; VLB overlay in black") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))
  
# -----------------------

all.valles$dbh.2 <- round(all.valles$dbh)
summary(all.valles)


valles.plot2 <- aggregate(all.valles$density, by=list(all.valles$site, all.valles$plot.id, all.valles$species, all.valles$dbh.2), FUN=sum, na.rm=T)
names(valles.plot2) <- c("site", "plot.id", "species", "dbh", "density")
summary(valles.plot2)

library(car)
valles.plot2$site2 <- as.factor(ifelse(substr(valles.plot2$site,1,4)=="mcon", "mcon", "ppine"))

unique(valles.plot2$plot.id)
valles.plot2$plot.id2 <- as.factor(ifelse(substr(valles.plot2$plot.id,1,5)=="PPINE", substr(valles.plot2$plot.id,6,6),
                                ifelse(substr(valles.plot2$plot.id,1,4)=="MCON", substr(valles.plot2$plot.id,5,5),
                                       ifelse(substr(valles.plot2$plot.id,3,3)=="A",5, 6))))
summary(valles.plot2)

qplot(x=dbh, data=valles.plot2, geom="histogram", breaks=dbh.bins1, fill=species, weight=sum(density)/density)+
  facet_grid(plot.id2 ~ site2)+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_continuous(name="DBH")+scale_y_continuous(name="Stems/ha") + ggtitle("Size Distribution weighted by Density") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))


valles.site2 <- aggregate(valles.plot2$density, by=list(valles.plot2$site, valles.plot2$species, valles.plot2$dbh), FUN=mean)
names(valles.site2) <- c("site", "species", "dbh", "density")
summary(valles.site2)

qplot(x=dbh, data=valles.site2, geom="histogram", breaks=dbh.bins1, fill=species, weight=density)+
  facet_grid(site ~ .)+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_continuous(name="DBH")+scale_y_continuous(name="Stems/ha") + ggtitle("Size Distribution weighted by Density") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))




qplot(x=dbh, data=valles.site2, geom="histogram", breaks=dbh.bins1, fill=species, weight=sum(density)/density)+
  facet_grid(site ~ .)+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_continuous(name="DBH")+scale_y_continuous(name="Stems/ha") + ggtitle("Size Distribution weighted by Density") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))

test <- t.test(marcy.pipo$DBH_1, ross.pipo$dbh, na.rm=T)
summary(test)
test


