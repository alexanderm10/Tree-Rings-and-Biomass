library(ggplot2)
##############################################################################
# wanting to make a composite figure showing the contribution of uncertainties
##############################################################################

# loading in the various datasets that will be needed

# allometric uncertainty of BM at the site level
load("valles_allometry_uncertainty.Rdata")
allom.uncert$LB.dev <- allom.uncert$Mean - allom.uncert$LB
allom.uncert$UB.dev <-  allom.uncert$UB - allom.uncert$Mean
allom.uncert <- allom.uncert[order(allom.uncert$Year),]
allom.uncert <- allom.uncert[order(allom.uncert$SiteID),]
summary(allom.uncert)

allom.uncert$range <- allom.uncert$UB - allom.uncert$LB

summary(allom.uncert[allom.uncert$SiteID=="VLF",])
summary(allom.uncert[allom.uncert$SiteID=="VUF",])

# density BM--uses mean allometric eqtn. and accounts for differences in density with just ROSS plots
load("valles_density_uncertainty.Rdata")
dens.uncert$LB.dev <- dens.uncert$Mean - dens.uncert$LB
dens.uncert$UB.dev <-  dens.uncert$UB - dens.uncert$Mean
dens.uncert <- dens.uncert[order(dens.uncert$Year),]
dens.uncert <- dens.uncert[order(dens.uncert$SiteID),]
summary(dens.uncert)

dens.uncert$range <- dens.uncert$UB - dens.uncert$LB

summary(dens.uncert[dens.uncert$SiteID=="VLF",])
summary(dens.uncert[dens.uncert$SiteID=="VUF",])

# mortality Uncertainty of BM at the site level
load("valles_mortality_uncertainty.Rdata")
dummy.year <- data.frame(Year=dens.uncert$Year, SiteID=dens.uncert$SiteID)
mort.uncert <- merge(mort.uncert, dummy.year, all.x=T, all.y=T)
mort.uncert <- mort.uncert[order(mort.uncert$Year),]
mort.uncert <- mort.uncert[order(mort.uncert$SiteID),]

mort.uncert$LB.dev <- dens.uncert$Mean - mort.uncert$LB
mort.uncert$UB.dev <-  mort.uncert$UB - dens.uncert$Mean
summary(mort.uncert)

mort.uncert$range <- mort.uncert$UB - mort.uncert$LB

summary(mort.uncert[mort.uncert$SiteID=="VLF",])
summary(mort.uncert[mort.uncert$SiteID=="VUF",])

# uncertainty in the increment
load("dated_v_all_valles.inc.stack.Rdata")
valles.inc.stack <- valles.inc.stack[valles.inc.stack$group=="all",]
valles.inc.stack <- valles.inc.stack[order(valles.inc.stack$Year),]
valles.inc.stack <- valles.inc.stack[order(valles.inc.stack$SiteID),]

valles.inc.stack$LB.dev <- valles.inc.stack$Mean.inc - valles.inc.stack$inc.LB
valles.inc.stack$UB.dev <- valles.inc.stack$inc.UB - valles.inc.stack$Mean.inc
summary(valles.inc.stack)

valles.inc.stack$range <- valles.inc.stack$inc.UB - valles.inc.stack$inc.LB

summary(valles.inc.stack[valles.inc.stack$SiteID=="VLF",])
summary(valles.inc.stack[valles.inc.stack$SiteID=="VUF",])







#combine the different areas into one figure
# 
# Order of operations
# 1)variability in the tree rings
# 2)allometric uncertainty
# 3)density uncertainty
# 4)mortality adustment


# Creating a dataframe that adds the uncertainties together
# will use the mean from the allometric uncertainty as our root

# doing a test first with just VLF--it worked--- expanding to both sites at Valles Caldera.
bm.final <- data.frame(Year=unique(allom.uncert$Year), SiteID=allom.uncert$SiteID, 
                        Base=allom.uncert[,"Mean"])
bm.final$LB.inc <- bm.final$Base - valles.inc.stack[,"LB.dev"]
bm.final$UB.inc <- bm.final$Base + valles.inc.stack[,"UB.dev"]
bm.final$LB.inc <- ifelse(bm.final$LB.inc<0, 0, bm.final$LB.inc)

bm.final$LB.allom <- bm.final$LB.inc - allom.uncert[, "LB.dev"]
bm.final$UB.allom <- bm.final$UB.inc + allom.uncert[, "UB.dev"]
bm.final$LB.allom <- ifelse(bm.final$LB.allom<0, 0, bm.final$LB.allom)

bm.final$LB.dens <- bm.final$LB.allom - dens.uncert[, "LB.dev"]
bm.final$UB.dens <- bm.final$UB.allom + dens.uncert[, "UB.dev"]
bm.final$LB.dens <- ifelse(bm.final$LB.dens<0, 0, bm.final$LB.dens)

bm.final$LB.mort <- bm.final$LB.dens - mort.uncert[, "LB.dev"]
bm.final$UB.mort <- bm.final$UB.dens + mort.uncert[, "UB.dev"]
bm.final$LB.mort <- ifelse(bm.final$LB.dens<0, 0, bm.final$LB.dens)

summary(bm.final)

summary(bm.final[bm.final$SiteID=="VLF",])   
summary(bm.final[bm.final$SiteID=="VUF",])   

# adding the uncertainties together in quadarture

vlf.final <- bm.final[bm.final$SiteID=="VLF",]
vuf.final <- bm.final[bm.final$SiteID=="VUF",]

vlf.final$range <- vlf.final$UB.mort - vlf.final$LB.mort
vuf.final$range <- vuf.final$UB.mort - vuf.final$LB.mort

summary(vlf.final)
summary(vuf.final)




ggplot(bm.final[bm.final$Year >= 1925 & bm.final$Year <=2011,]) + facet_grid(SiteID ~ .) +
  geom_line(aes(x=Year, y=Base), size=1.5, color="black") +

  #1) Increment Uncertainty
  geom_ribbon(aes(x=Year, ymin=LB.inc, ymax=UB.inc), alpha=0.5, fill="green3") +

  #2) Allometric Uncertainty -- separate for upper & lower to make things clearer
  geom_ribbon(aes(x=Year, ymin=LB.allom, ymax=LB.inc), alpha=0.5, fill="blue") +
  geom_ribbon(aes(x=Year, ymin=UB.allom, ymax=UB.inc), alpha=0.5, fill="blue") +
  
  #3) Density Uncertainty -- separate for upper & lower to make things clearer
  geom_ribbon(aes(x=Year, ymin=LB.dens, ymax=LB.allom), alpha=0.5, fill="red") +
  geom_ribbon(aes(x=Year, ymin=UB.dens, ymax=UB.allom), alpha=0.5, fill="red") +
  
  #4) Mortality Uncertainty -- separate for upper & lower to make things clearer
  geom_ribbon(aes(x=Year, ymin=LB.mort, ymax=LB.dens), alpha=0.5, fill="orange2") +
  geom_ribbon(aes(x=Year, ymin=UB.mort, ymax=UB.dens), alpha=0.5, fill="orange2") +
  
  # Reiterate mean line for clarity
  geom_line(aes(x=Year, y=Base), size=1.5, color="black") +
  labs(title= "Stacked Uncertainties", x="Year", y="Aboveground Biomass (kg/m2)")

  



########################################################
# Uncertainty percentages
########################################################



# Looking at relative magnitudes of uncertainties

# Allometry
vlf.allom.rel <- allom.uncert[allom.uncert$SiteID=="VLF", "range"]/allom.uncert[allom.uncert$SiteID=="VLF", "Mean"]
vuf.allom.rel <- allom.uncert[allom.uncert$SiteID=="VUF", "range"]/allom.uncert[allom.uncert$SiteID=="VUF", "Mean"]
summary(vlf.allom.rel)
summary(vuf.allom.rel)

mean(vlf.allom.rel, na.rm=T); sd(vlf.allom.rel, na.rm=T)
mean(vuf.allom.rel, na.rm=T); sd(vuf.allom.rel, na.rm=T)


# Density
vlf.dens.rel <- dens.uncert[dens.uncert$SiteID=="VLF", "range"]/allom.uncert[allom.uncert$SiteID=="VLF", "Mean"]
vuf.dens.rel <- dens.uncert[dens.uncert$SiteID=="VUF", "range"]/allom.uncert[allom.uncert$SiteID=="VUF", "Mean"]
summary(vlf.dens.rel)
summary(vuf.dens.rel)

mean(vlf.dens.rel, na.rm=T); sd(vlf.dens.rel, na.rm=T)
mean(vuf.dens.rel, na.rm=T); sd(vuf.dens.rel, na.rm=T)


# Mortality
vlf.mort.rel <- mort.uncert[mort.uncert$SiteID=="VLF", "range"]/allom.uncert[allom.uncert$SiteID=="VLF", "Mean"]
vuf.mort.rel <- mort.uncert[mort.uncert$SiteID=="VUF", "range"]/allom.uncert[allom.uncert$SiteID=="VUF", "Mean"]
summary(vlf.mort.rel)
summary(vuf.mort.rel)

mean(vlf.mort.rel, na.rm=T); sd(vlf.mort.rel, na.rm=T)
mean(vuf.mort.rel, na.rm=T); sd(vuf.mort.rel, na.rm=T)


# TR increment
vlf.inc.rel <- valles.inc.stack[valles.inc.stack$SiteID=="VLF", "range"]/allom.uncert[allom.uncert$SiteID=="VLF", "Mean"]
vuf.inc.rel <- valles.inc.stack[valles.inc.stack$SiteID=="VUF", "range"]/allom.uncert[allom.uncert$SiteID=="VUF", "Mean"]
summary(vlf.inc.rel)
summary(vuf.inc.rel)

mean(vlf.mort.rel, na.rm=T); sd(vlf.mort.rel, na.rm=T)
mean(vuf.mort.rel, na.rm=T); sd(vuf.mort.rel, na.rm=T)

# Overall uncertainty

vlf.final.rel  <- vlf.final$range/vlf.final$Base
vuf.final.rel <- vuf.final$range/vuf.final$Base

mean(vlf.final.rel, na.rm=T); sd(vlf.final.rel, na.rm=T)
mean(vuf.final.rel, na.rm=T); sd(vuf.final.rel, na.rm=T)
  
  