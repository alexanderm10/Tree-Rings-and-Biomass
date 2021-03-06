# Mortality data gathered from van Mantgem 2009
# used interior numbers start year 1979 and end year 2011
# used to change densities for the valles 
load("site_density.Rdata")
summary(site.density)
ross.density <- read.csv("raw input files/ross_density.csv", header=T)

ross.density.site <- ross.density[substr(ross.density$PlotID, 1, 1)=="V",]


start.m <- rnorm(5000, mean=0.4843, sd= 0.2823)
rate.m <- rnorm(5000, mean= 0.024, sd= 0.027)


start.yr.vm = 1979 # First year for which the increasing mortality rate is valid
start.yr = 1905
end.yr = 2011

mort.rate <- data.frame(array(dim=c(5000, length(start.yr:end.yr))))
names(mort.rate) <- start.yr:end.yr

for(j in 1:5000){
  for(i in 1:ncol(mort.rate)){
    if(as.numeric(names(mort.rate)[i])<=start.yr.vm){
      mort.rate[j,i] = sample(start.m, size=1, replace=T)
    } else {
    mort.rate[j,i] <- mort.rate[j,i-1] + mort.rate[j,i-1]*sample(rate.m, size=1, replace=T)
    }
  }
}
summary(mort.rate)
ci <- apply(mort.rate, 2, FUN=quantile, c(0.025, 0.975))
mort.mean <- apply(mort.rate, 2, FUN=mean)

plot(mort.mean, type="l", ylim=range(ci))
lines(ci[1,], lty="dashed")
lines(ci[2,], lty="dashed")

#ci<- as.data.frame(ci)

mort.rate <- mort.mean[sort(names(mort.mean), decreasing=T)]
ci.mort.rate <- ci[,sort(colnames(ci), decreasing=T)]


vla.density <- vector(length=length(end.yr:start.yr))
vla.density[] <- 0.590278
names(vla.density) <- end.yr:start.yr

for(i in 2:length(vla.density)){
  vla.density[i] <- vla.density[i-1] + vla.density[i-1]*mort.mean[i]/100
}
vla.density

vla.density.ci <-  data.frame(array(dim=c(nrow(ci), ncol(ci)))) 
rownames(vla.density.ci) <- row.names(ci)
row.names(vla.density.ci) <- names(ci.mort.rate)
vla.density.ci <-  ci.mort.rate 
vla.density.ci[] <- 0.590278

for(j in 1:nrow(vla.density.ci)){
  for(i in 2:ncol(vla.density.ci)){
    vla.density.ci[j,i] <-vla.density.ci[j,i-1] + vla.density.ci[j,i-1]*ci.mort.rate[j,i]/100
  }
}
plot(vla.density~names(vla.density), type="l", ylim=range(vla.density.ci))
lines(vla.density.ci[1,]~names(vla.density), lty="dashed")
lines(vla.density.ci[2,]~names(vla.density), lty="dashed")

vlb.density <- vector(length=length(end.yr:start.yr))
vlb.density[] <- 0.089744
names(vlb.density) <- end.yr:start.yr

for(i in 2:length(vlb.density)){
  vlb.density[i] <- vlb.density[i-1] + vlb.density[i-1]*mort.mean[i]/100
}
vlb.density

vlb.density.ci <-  data.frame(array(dim=c(nrow(ci), ncol(ci)))) 
rownames(vlb.density.ci) <- row.names(ci)
row.names(vlb.density.ci) <- names(ci.mort.rate)
vlb.density.ci <-  ci.mort.rate 
vlb.density.ci[] <- 0.089744

for(j in 1:nrow(vlb.density.ci)){
  for(i in 2:ncol(vlb.density.ci)){
    vlb.density.ci[j,i] <-vlb.density.ci[j,i-1] + vlb.density.ci[j,i-1]*ci.mort.rate[j,i]/100
  }
}
plot(vlb.density~names(vlb.density), type="l", ylim=range(vlb.density.ci))
lines(vlb.density.ci[1,]~names(vlb.density), lty="dashed")
lines(vlb.density.ci[2,]~names(vlb.density), lty="dashed")


vua.density <- vector(length=length(end.yr:start.yr))
vua.density[] <- 0.111111
names(vua.density) <- end.yr:start.yr

for(i in 2:length(vua.density)){
  vua.density[i] <- vua.density[i-1] + vua.density[i-1]*mort.mean[i]/100
}
vua.density

vua.density.ci <-  data.frame(array(dim=c(nrow(ci), ncol(ci)))) 
rownames(vua.density.ci) <- row.names(ci)
row.names(vua.density.ci) <- names(ci.mort.rate)
vua.density.ci <-  ci.mort.rate  
vua.density.ci[] <- 0.111111

for(j in 1:nrow(vua.density.ci)){
  for(i in 2:ncol(vua.density.ci)){
    vua.density.ci[j,i] <-vua.density.ci[j,i-1] + vua.density.ci[j,i-1]*ci.mort.rate[j,i]/100
  }
}
plot(vua.density~names(vua.density), type="l", ylim=range(vua.density.ci))
lines(vua.density.ci[1,]~names(vua.density), lty="dashed")
lines(vua.density.ci[2,]~names(vua.density), lty="dashed")

vub.density <- vector(length=length(end.yr:start.yr))
vub.density[] <- 0.194444
names(vub.density) <- end.yr:start.yr

for(i in 2:length(vub.density)){
  vub.density[i] <- vub.density[i-1] + vub.density[i-1]*mort.mean[i]/100
}
vub.density

vub.density.ci <-  data.frame(array(dim=c(nrow(ci), ncol(ci)))) 
rownames(vub.density.ci) <- row.names(ci)
row.names(vub.density.ci) <- names(ci.mort.rate)
vub.density.ci <-  ci.mort.rate 
vub.density.ci[] <- 0.194444

for(j in 1:nrow(vub.density.ci)){
  for(i in 2:ncol(vub.density.ci)){
    vub.density.ci[j,i] <-vub.density.ci[j,i-1] + vub.density.ci[j,i-1]*ci.mort.rate[j,i]/100
  }
}
plot(vub.density~names(vub.density), type="l", ylim=range(vub.density.ci))
lines(vub.density.ci[1,]~names(vub.density), lty="dashed")
lines(vub.density.ci[2,]~ names(vub.density), lty="dashed")

# vuf.density <- vector(length=length(end.yr:start.yr))
# # vuf.density[] <- 0.1177709
# vuf.density[] <- 0.1527775
# names(vuf.density) <- end.yr:start.yr
# 
# for(i in 2:length(vuf.density)){
#   vuf.density[i] <- vuf.density[i-1] + vuf.density[i-1]*mort.mean[i]/100
# }
# vuf.density
# 
# vuf.density.ci <-  data.frame(array(dim=c(nrow(ci), ncol(ci)))) 
# rownames(vuf.density.ci) <- row.names(ci)
# row.names(vuf.density.ci) <- names(ci.mort.rate)
# vuf.density.ci <-  ci.mort.rate 
# # vuf.density.ci[] <- 0.1177709
# vuf.density.ci[] <- 0.1527775
# 
# for(j in 1:nrow(vuf.density.ci)){
#   for(i in 2:ncol(vuf.density.ci)){
#     vuf.density.ci[j,i] <-vuf.density.ci[j,i-1] + vuf.density.ci[j,i-1]*ci.mort.rate[j,i]/100
#   }
# }
# plot(vuf.density~names(vuf.density), type="l", ylim=range(vuf.density.ci))
# lines(vuf.density.ci[1,]~names(vuf.density), lty="dashed")
# lines(vuf.density.ci[2,]~ names(vuf.density), lty="dashed")
# 
# vlf.density <- vector(length=length(end.yr:start.yr))
# # vlf.density[] <- 0.1568394
# vlf.density[] <- 0.3400110
# names(vlf.density) <- end.yr:start.yr
# 
# for(i in 2:length(vlf.density)){
#   vlf.density[i] <- vlf.density[i-1] + vlf.density[i-1]*mort.mean[i]/100
# }
# vlf.density
# 
# vlf.density.ci <-  data.frame(array(dim=c(nrow(ci), ncol(ci)))) 
# rownames(vlf.density.ci) <- row.names(ci)
# row.names(vlf.density.ci) <- names(ci.mort.rate)
# vlf.density.ci <-  ci.mort.rate 
# # vlf.density.ci[] <- 0.1568394
# vlf.density.ci[] <- 0.3400110
# 
# for(j in 1:nrow(vlf.density.ci)){
#   for(i in 2:ncol(vlf.density.ci)){
#     vlf.density.ci[j,i] <-vlf.density.ci[j,i-1] + vlf.density.ci[j,i-1]*ci.mort.rate[j,i]/100
#   }
# }
# plot(vlf.density~names(vlf.density), type="l", ylim=range(vlf.density.ci))
# lines(vlf.density.ci[1,]~names(vlf.density), lty="dashed")
# lines(vlf.density.ci[2,]~ names(vlf.density), lty="dashed")

vla.density <- vla.density[order(names(vla.density), decreasing=T)]
vla.density.ci <- vla.density.ci[,sort(colnames(vla.density.ci), decreasing=T)]

vlb.density <- vlb.density[order(names(vlb.density), decreasing=T)]
vlb.density.ci <- vlb.density.ci[,sort(colnames(vlb.density.ci), decreasing=T)]

vua.density <- vua.density[order(names(vua.density), decreasing=T)]
vua.density.ci <- vua.density.ci[,sort(colnames(vua.density.ci), decreasing=T)]

vub.density <- vub.density[order(names(vub.density), decreasing=T)]
vub.density.ci <- vub.density.ci[,sort(colnames(vub.density.ci), decreasing=T)]

# vuf.density <- vuf.density[order(names(vuf.density), decreasing=T)]
# vuf.density.ci <- vuf.density.ci[,sort(colnames(vuf.density.ci), decreasing=T)]
# 
# vlf.density <- vlf.density[order(names(vlf.density), decreasing=T)]
# vlf.density.ci <- vlf.density.ci[,sort(colnames(vlf.density.ci), decreasing=T)]


mort.list <- array(list(vua.density, vub.density, vla.density, vlb.density))
mort.ci.list <- array(list(vua.density.ci, vub.density.ci, vla.density.ci, vlb.density.ci))

names(mort.list) <- names(mort.ci.list) <- c("VUA", "VUB", "VLA", "VLB")
summary(mort.list)

summary(mort.ci.list)

# mort.site <- array(list(vlf.density, vuf.density))
# mort.site.ci <-array(list(vlf.density.ci, vuf.density.ci))
# 
# names(mort.site) <- names(mort.site.ci) <- c("VLF", "VUF")
# summary(mort.site)
# summary(mort.site.ci)

#####################################################
# applying changing densities to the biomass estimate
#####################################################
library(dplR)
library(ggplot2)
se <- function(x){
  sd(x, na.rm=TRUE) / sqrt((length(!is.na(x))))}

g.filled.diam <- read.csv("GapFilling_DBHrecon_ALL.csv", header=T, row.names=1)
g.filled.diam <- g.filled.diam[,substr(names(g.filled.diam),1,1)=="V"]
summary(g.filled.diam)
head(g.filled.diam)

#need to truncate the diameter time series to be the same as the mortality curve time series
g.filled.diam <- g.filled.diam[row.names(g.filled.diam)>=1905 & row.names(g.filled.diam)<=2011,]


# read in tree data
tree.data <- read.csv("TreeData.csv", header=T)
summary(tree.data)

trees.use <- tree.data[substr(tree.data$PlotID, 1, 1)=="V",]
summary(trees.use)

plot.data <- read.csv("raw input files/DOE_plus_Valles.csv")
summary(plot.data)



##########################################################################
# Allometric Equations
##########################################################################


#Convert to biomass with the allometric equation
#using the PECAN generated bayesian equations
library(car)
load("allometries_list.Rdata")
# Getting rid of POTR for now for conceptual figure purposes
#trees.use <- trees.use[!(trees.use$Species=="POTR"),]
summary(trees.use)
unique(trees.use$Species)

trees.use$spp.allom <- recode(trees.use$Species, " 'PIEN'='picea.sp'; 'PIPO'='pipo'; 'PSME'='psme'; 'POTR' = 'potr'")
summary(trees.use)
plots <- unique(trees.use$PlotID) # You had the right idea, but it was throwing errors because you were trying to evaluate plots you haven't gotten to yet


# will want to do general equations and pft level equations as well, but later
# log(AGB) = mu0 + mu1*log(DBH) --equaton form of PECAN allometrics

#allom.eq <- function(mu0, mu1, DBH) { mu0 * DBH^mu1}
allom.eq <- function(mu0, mu1, DBH) { exp(mu0 + mu1 * log(DBH) )}

allom.temp <- g.filled.diam
allom.temp[,] <- NA

# dbh=0 causes problems, so we're going to make those NA
g.filled.diam[g.filled.diam==0] <- 1e-6
min(g.filled.diam, na.rm=T)
summary(g.filled.diam)
dim(g.filled.diam)

#want to set up an array that has the layers 1) the time series of biomass 2) mean mortality numbers at the plot level 3) confidence intervals for the plots
bm.array.mean <- array(NA, dim=c(nrow(g.filled.diam), length(unique(trees.use$PlotID)), nrow(allometries[[1]])))
bm.array.high <- array(NA, dim=c(nrow(g.filled.diam), length(unique(trees.use$PlotID)), nrow(allometries[[1]])))
bm.array.low <- array(NA, dim=c(nrow(g.filled.diam), length(unique(trees.use$PlotID)), nrow(allometries[[1]])))



row.names(bm.array.mean) <- row.names(g.filled.diam)  #CRR Added
row.names(bm.array.high) <- row.names(g.filled.diam)  #CRR Added
row.names(bm.array.low) <- row.names(g.filled.diam)  #CRR Added

dimnames(bm.array.mean)[[2]] <- unique(trees.use$PlotID)
dimnames(bm.array.high)[[2]] <- unique(trees.use$PlotID)
dimnames(bm.array.low)[[2]] <- unique(trees.use$PlotID)

summary(bm.array.high[,,1])

dim(bm.array.low)
dim(bm.array.mean)
dim(bm.array.high)

#####################################################################
# Starting the Loop
#####################################################################

for(i in 1:nrow(allometries[[1]])){
  allom.temp <- g.filled.diam
  allom.temp[,] <- NA
  
  # Species loop for calculating tree biomass
  for(j in unique(trees.use$spp.allom)){
    cols <- which(names(g.filled.diam) %in% trees.use[trees.use$spp.allom==j, "TreeID"])
    # Note: we'll have to make this a bit fancier in the future for species with mu0==0
    #   allom.temp[,cols] <- allom.eq(mu0= -3.5185,
    #                          mu1 = 2.6909,
    #                         #DBH = seq(from=30, to=1, length=nrow(g.filled.diam)))
    #                          DBH = g.filled.diam[,cols])
    # test <- allom.eq(mu0=ifelse(!(allometries[[j]][i,"mu0"]==0 & allometries[[j]][i,"mu1"]==0),allometries[[j]][i,"mu0"], allometries[[j]][i,"Bg0"]),
    #                               mu1 =ifelse(!(allometries[[j]][i,"mu0"]==0 & allometries[[j]][i,"mu1"]==0),allometries[[j]][i,"mu1"], allometries[[j]][i,"Bg1"]),
    #                               DBH = g.filled.diam[,cols])
    # mu0 = ifelse(!(allometries[[j]][i,"mu0"]==0 & allometries[[j]][i,"mu1"]==0),allometries[[j]][i,"mu0"], allometries[[j]][i,"Bg0"])
    # mu1 = ifelse(!(allometries[[j]][i,"mu0"]==0 & allometries[[j]][i,"mu1"]==0),allometries[[j]][i,"mu1"], allometries[[j]][i,"Bg1"])
    mu0=allometries[[j]][i,"Bg0"]
    mu1=allometries[[j]][i,"Bg1"]
    allom.temp[,cols] <- allom.eq(mu0=mu0, mu1 = mu1, DBH = g.filled.diam[,cols])
  }
  # summing to the plot level
  
  allom.temp[is.na(allom.temp)] <- 0
  
  # biomass loop for summing trees to plots
  # We're doing the unit conversions here; we had calculated density in stems/ha, but Christy wants to look at Biomass in kg/m2, so we're putting everything in kg/m2 here
  for(p in 1:length(plots)){
    cols <- which(names(allom.temp) %in% trees.use[trees.use$PlotID==plots[p], "TreeID"])
#     if(substr(plots[p],1,1)=="V"){
#       bm.array[,p,i] <- rowMeans(allom.temp[,cols])*plot.data[plot.data$PlotID==paste(plots[p]), "Density.Total..stems.ha."]/10000 #mean tree * trees/ha (do for Valles only bc sum of trees != plot density; different sampling method than Neil)
#     } else {

#         if(substr(plots[p],1,1)=="V"){
                bm.array.mean[,p,i] <- rowMeans(allom.temp[,cols])*mort.list[[p]] #mean tree * trees/ha (do for Valles only bc sum of trees != plot density; different sampling method than Neil)
              #  bm.array.mean[,p,i] <- rowMeans(allom.temp[,cols])*plot.data[plot.data$PlotID==paste(plots[p]), "Density.Total..stems.ha."]/10000
                
                bm.array.low[,p,i] <- rowMeans(allom.temp[,cols])*mort.ci.list[[p]][1,] #mean tree * trees/ha (do for Valles only bc sum of trees != plot density; different sampling method than Neil)
                bm.array.high[,p,i] <- rowMeans(allom.temp[,cols])*mort.ci.list[[p]][2,] #mean tree * trees/ha (do for Valles only bc sum of trees != plot density; different sampling method than Neil)
            }

#                 } else {

#       temp1 <- temp2 <- temp3 <- allom.temp[,cols]
#       for(n in names(allom.temp[,cols])){ # Convert biomass/tree to biomass/ha
#         temp1[,n] <- temp1[,n] * mort.list[[plots[p]]]
#         temp2[,n] <- temp2[,n] * mort.ci.list[[plots[p]]][2,]
#         temp3[,n] <- temp3[,n] * mort.ci.list[[plots[p]]][1,]
#       
#       bm.array.mean[,p,i] <- rowSums(temp1) #sum biomass/ha
#       bm.array.hi[,p,i] <- rowSums(temp2) #sum biomass/ha
#       bm.array.lo[,p,i] <- rowSums(temp3) #sum biomass/ha

}

#--------------------------------------------------
summary(bm.array.mean[,,1])

biom.mort.mean <- apply(bm.array.mean[,,], c(1,2), mean)
biom.mort.mean <- as.data.frame(biom.mort.mean)

biom.mort.UB <- apply(bm.array.high[,,], c(1,2), mean)
biom.mort.UB <- as.data.frame(biom.mort.UB)
names(biom.mort.UB) <- paste(plots, "UB", sep=".")

biom.mort.LB <- apply(bm.array.low[,,], c(1,2), mean)
biom.mort.LB <- as.data.frame(biom.mort.LB)
names(biom.mort.LB) <- paste(plots, "LB", sep=".")

summary(biom.mort.mean)
summary(biom.mort.UB)
summary(biom.mort.LB)

# put the three separate DF into one DF to have mean, upper bound, and lower bound
biom.mort.valles <- as.data.frame(c(biom.mort.mean, biom.mort.LB, biom.mort.UB))
row.names(biom.mort.valles) <- row.names(biom.mort.mean)

summary(biom.mort.valles)
write.csv(biom.mort.valles, "valles_bm_mort_corrected.csv")


# Stack biom.mort.valles 
biom.mort.valles.stack <- stack(biom.mort.valles[1:4])
names(biom.mort.valles.stack) <- c("Biom.Mean", "PlotID")
biom.mort.valles.stack$Year <- as.numeric(paste(row.names(biom.mort.valles)))
biom.mort.valles.stack$Plot <- as.factor(substr(biom.mort.valles.stack$PlotID, 3,3))
biom.mort.valles.stack$Site <- as.factor(substr(biom.mort.valles.stack$PlotID, 1,2))
summary(biom.mort.valles.stack)

biom.mort.valles.stack.lb <- stack(biom.mort.valles[5:8])
names(biom.mort.valles.stack.lb) <- c("Biom.LB", "PlotID")

biom.mort.valles.stack.ub <- stack(biom.mort.valles[9:12])
names(biom.mort.valles.stack.ub) <- c("Biom.UB", "PlotID")

biom.mort.valles.stack$Biom.LB <- biom.mort.valles.stack.lb[,1]
biom.mort.valles.stack$Biom.UB <- biom.mort.valles.stack.ub[,1]
summary(biom.mort.valles.stack)


save(biom.mort.valles.stack, file="valles_bm_mortality_corrected_stack.Rdata")
load("valles_bm_recon_stack.Rdata")
load("valles_bm_recon.Rdata")

summary(biom.valles.stack)

ggplot(data=biom.mort.valles.stack[biom.mort.valles.stack$Year<2012 & (biom.mort.valles.stack$Site=="VL"),])  + facet_grid(Plot ~ Site) +
  geom_line(data=biom.valles.stack[biom.valles.stack$Year<2012 & (biom.valles.stack$Site=="VL"),], aes(x=Year, y=Biom.Mean, color=PlotID), linetype="dashed", size=0.5) +
  # plotting total site basal area  
  geom_ribbon(aes(x=Year, ymin=Biom.LB, ymax=Biom.UB, fill=PlotID), alpha=0.5) +
  geom_line(aes(x=Year, y=Biom.Mean, color=PlotID)) +
  scale_y_continuous(name=expression(bold(paste("Biomass (kg m" ^ "-2 ", ")"))))+
  scale_x_continuous(name=expression(bold(paste("Year"))))+
  ggtitle("Valles Caldera Lower (PIPO)")

ggplot(data=biom.mort.valles.stack[biom.mort.valles.stack$Year<2012 & (biom.mort.valles.stack$Site=="VU"),])  + facet_grid(Plot ~ Site) +
  #ggplot(data=biom.mort.valles.stack[biom.mort.valles.stack$Year<2012 & (biom.mort.valles.stack$PlotID=="VUB"),])  + facet_grid(Plot ~ Site) +
  # plotting total site basal area  
  geom_ribbon(aes(x=Year, ymin=Biom.LB, ymax=Biom.UB, fill=PlotID), alpha=0.5) +
  geom_line(data=biom.valles.stack[biom.valles.stack$Year<2012 & (biom.valles.stack$Site=="VU"),], aes(x=Year, y=Biom.Mean, color=PlotID), linetype="dashed", size=0.5) +
  geom_line(aes(x=Year, y=Biom.Mean, color=PlotID)) +
  scale_y_continuous(name=expression(bold(paste("Biomass (kg m" ^ "-2 ", ")"))))+
  scale_x_continuous(name=expression(bold(paste("Year"))))+
  ggtitle("Valles Caldera Upper (MCON)")



# valles.cum.plot<- 
#pdf("Valles_Biomass_19March_Christy_ChojnackyOnly.pdf", height=8.5, width=11)
ggplot(data=biom.mort.valles.stack[biom.mort.valles.stack$Year<2012,])  + facet_grid(.~Site) +
  geom_line(data=biom.valles.stack[biom.valles.stack$Year<2012,], aes(x=Year, y=Biom.Mean, color=PlotID), linetype="dashed", size=0.5) +
  # plotting total site basal area  
  geom_ribbon(aes(x=Year, ymin=Biom.LB, ymax=Biom.UB, fill=PlotID), alpha=0.5) +
  geom_line(aes(x=Year, y=Biom.Mean, color=PlotID)) +
  scale_y_continuous(name=expression(bold(paste("Biomass (kg m" ^ "-2 ", ")"))))+
  scale_x_continuous(name=expression(bold(paste("Year"))))+
  ggtitle("Valles Caldera Motality")

#scale_x_continuous(name=expression(bold(paste("Mean Annual Temperature (" ^"o", "C)"))))
# Aggregating the individual plot densities to a site level density by taking the mean of the CI's and the mean of the means

ci.vu <- data.frame(Year=as.numeric(row.names(biom.mort.mean)), SiteID= "VUF",
                    Mean = rowMeans(biom.mort.mean[,c("VUA", "VUB")]),
                    LB=apply(biom.mort.LB[,c("VUA.LB", "VUB.LB")],1, mean), 
                    UB=apply(biom.mort.UB[,c("VUA.UB", "VUB.UB")],1, mean))

ci.vl <- data.frame(Year=as.numeric(row.names(biom.mort.mean)), SiteID= "VLF",
                    Mean = rowMeans(biom.mort.mean[,c("VLA", "VLB")]),
                    LB=apply(biom.mort.LB[,c("VLA.LB", "VLB.LB")],1, mean), 
                    UB=apply(biom.mort.UB[,c("VLA.UB", "VLB.UB")],1, mean))


mort.uncert <- rbind(ci.vl, ci.vu)
summary(mort.uncert)
load("valles_density_uncertainty.Rdata")

ggplot(mort.uncert[mort.uncert$Year <=2011,]) + facet_grid(SiteID ~ .) +
  geom_ribbon(aes(x=Year, ymin=LB, ymax=UB, fill=SiteID), alpha=0.5) +
  geom_line(aes(x=Year, y=Mean, color=SiteID))+
  geom_line(data=dens.uncert[dens.uncert$Year <= 2011,], aes(x=Year, y=Mean, color=SiteID), linetype="dashed")

save(mort.uncert, file="valles_mortality_uncertainty.Rdata")



summary(biom.mort.valles)
summary(biom.valles)

vua.diff<- biom.mort.valles[row.names(biom.mort.valles)>= 1979 & row.names(biom.mort.valles)<=2011, "VUA"] - biom.valles[row.names(biom.valles)>= 1979 & row.names(biom.valles)<=2011, "VUA"]
vua.diff<- as.data.frame(vua.diff); row.names(vua.diff)<- c(2011:1979)
vua.diff

vub.diff<- biom.mort.valles[row.names(biom.mort.valles)>= 1979 & row.names(biom.mort.valles)<=2011, "VUB"] - biom.valles[row.names(biom.valles)>= 1979 & row.names(biom.valles)<=2011, "VUB"]
vub.diff<- as.data.frame(vub.diff); row.names(vub.diff)<- c(2011:1979)
vub.diff

vla.diff<- biom.mort.valles[row.names(biom.mort.valles)>= 1979 & row.names(biom.mort.valles)<=2011, "VLA"] - biom.valles[row.names(biom.valles)>= 1979 & row.names(biom.valles)<=2011, "VLA"]
vla.diff<- as.data.frame(vla.diff); row.names(vla.diff)<- c(2011:1979)
vla.diff

vlb.diff<- biom.mort.valles[row.names(biom.mort.valles)>= 1979 & row.names(biom.mort.valles)<=2011, "VLB"] - biom.valles[row.names(biom.valles)>= 1979 & row.names(biom.valles)<=2011, "VLB"]
vlb.diff<- as.data.frame(vlb.diff); row.names(vlb.diff)<- c(2011:1979)
vlb.diff



