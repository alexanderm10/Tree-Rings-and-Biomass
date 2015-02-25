library(ggplot2)

# exploratory analysis of Marcy Litvak's Valles Caldera Data

marcy.pipo <- read.csv("raw input files/marcy_ppine_2013.csv", header=T)
marcy.mcon <- read.csv("raw input files/marcy_mcon_2012.csv", header=T)
ross.trees <- read.csv("raw input files/tree_metadata_DOE_plus_valles.csv", na.strings=c("", "NA", "#VALUE!", "*"),header=T)

summary(ross.trees)
ross.valles <- ross.trees[substr(ross.trees$PlotID,1,1)=="V",]


# ross.pipo <- ross.trees[c("TreeID", "Species", "DBH..cm."),]
# names(ross.pipo) <- c("tree.id", "species", "dbh")

summary(ross.valles)

ross.valles<- ross.valles[,c("TreeID", "PlotID", "Species","DBH..cm.", "Density..stems.ha.")]
names(ross.valles)<- c("tree.id", "plot.id", "species", "dbh", "density")
summary(ross.valles)

ross.valles$site<- as.factor(substr(ross.valles$plot.id, 1, 2))
summary(ross.valles)

ross.pipo <- ross.valles[substr(ross.valles$plot.id,1,2)=="VL",]
ross.mcon <- ross.valles[substr(ross.valles$plot.id,1,2)=="VU",]
summary(ross.pipo)
summary(ross.mcon)

marcy.valles<- merge(marcy.pipo, marcy.mcon, all.x=T, all.y=T)
summary(marcy.valles)

marcy.valles <- marcy.valles[,c("Site", "Tree_Tag_Number","Species", "DBH_1")]
summary(marcy.valles)
names(marcy.valles)<- c("site", "tree.id", "species", "dbh")
summary(marcy.valles)

all.valles<-merge(marcy.valles, ross.valles, all.x=T, all.y=T)
summary(all.valles)




all.valles$site<- ifelse(substr(all.valles$site, 1, 2)=="VL", "ross.ppine", 
                         ifelse(substr(all.valles$site,1,2)=="VU", "ross.mcon", 
                                ifelse(all.valles$site=="PPINE", "marcy.ppine", "marcy.mcon")))
all.valles$site<-as.factor(all.valles$site)
summary(all.valles)

dbh.bins1 <- seq(0, max(all.valles$dbh, na.rm=T), 2)



qplot(x=dbh, data=all.valles, geom="histogram", breaks=dbh.bins1, fill=species)+
  facet_grid(site ~ .)+
  theme(axis.line=element_line(color="black", size=0.5), panel.grid.major=element_blank(), panel.grid.minor= element_blank(), panel.border= element_blank(), panel.background= element_blank(), axis.text.x=element_text(angle=0, color="black", size=12), axis.text.y=element_text(color="black", size=12)) + 
  scale_x_continuous(name="DBH") + ggtitle("Size Distribution") #+ scale_fill_manual(values=as.vector(spp.col.tree$Color))


test <- t.test(marcy.pipo$DBH_1, ross.pipo$dbh, na.rm=T)
summary(test)
test
