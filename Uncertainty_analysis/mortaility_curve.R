# Mortality data gathered from van Mantgem 2009
# used interior numbers start year 1979 and end year 2011
# used to change densities for the valles 

start.m <- rnorm(1000, mean=0.4843, sd= 0.2823)
rate.m <- rnorm(1000, mean= 0.024, sd= 0.027)


start.yr = 1979
end.yr = 2011

mort.rate <- data.frame(array(dim=c(1000, 33)))
names(mort.rate) <- start.yr:end.yr

for(j in 1:1000){
  mort.rate[j,1] = sample(start.m, size=1, replace=T)
  
  for(i in 2:length(mort.rate)){
    mort.rate[j,i] <- mort.rate[j,i-1] + mort.rate[j,i-1]*sample(rate.m, size=1, replace=T)
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
vla.density.ci <-  ci 
vla.density.ci[,] <- 0.590278

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
vlb.density.ci <-  ci 
vlb.density.ci[,] <- 0.089744

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
vua.density.ci <-  ci 
vua.density.ci[,] <- 0.111111

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
vub.density.ci <-  ci 
vub.density.ci[,] <- 0.194444

for(j in 1:nrow(vub.density.ci)){
  for(i in 2:ncol(vub.density.ci)){
    vub.density.ci[j,i] <-vub.density.ci[j,i-1] + vub.density.ci[j,i-1]*ci.mort.rate[j,i]/100
  }
}
plot(vub.density~names(vub.density), type="l", ylim=range(vub.density.ci))
lines(vub.density.ci[1,]~names(vub.density), lty="dashed")
lines(vub.density.ci[2,]~ names(vub.density), lty="dashed")

mort.list <- list(vla.density, vlb.density, vua.density, vub.density)
mort.ci.list <- list(vla.density.ci, vlb.density.ci, vua.density.ci, vub.density.ci)
names(mort.list) <- names(mort.ci.list) <- c("VLA", "VLB", "VUA", "VUB")
summary(mort.list)

summary(mort.ci.list[[1]])
names(mort.ci.list) <- names(mort.ci.list) <- c("VLA", "VLB", "VUA", "VUB")


#####################################################
# applying changing densities to the biomass estimate
#####################################################
library(dplR)
library(ggplot2)
se <- function(x){
  sd(x, na.rm=TRUE) / sqrt((length(!is.na(x))))}

g.filled.diam <- read.csv("gap_filled_dbh.recon.csv", header=T, row.names=1)
g.filled.diam <- g.filled.diam[,substr(names(g.filled.diam),1,1)=="V"]
summary(g.filled.diam)

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

# Getting rid of POTR for now for conceptual figure purposes
trees.use <- trees.use[!(trees.use$Species=="POTR"),]
summary(trees.use)
unique(trees.use$Species)

trees.use$spp.allom <- recode(trees.use$Species, " 'PIEN'='picea.sp'; 'PIPO'='pipo'; 'PSME'='psme'")
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
bm.array <- array(NA, dim=c(nrow(g.filled.diam), length(unique(trees.use$PlotID)), ncol(allometries[[1]])))
row.names(bm.array) <- row.names(g.filled.diam)  #CRR Added

summary(bm.array[,,1])
  
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
    mu0=allometries[[j]][i,mean("Bg0")]
    mu1=allometries[[j]][i,mean("Bg1")]
    allom.temp[,cols] <- allom.eq(mu0=mu0, mu1 = mu1, DBH = g.filled.diam[,cols])
  }
  # summing to the plot level
  
  allom.temp[is.na(allom.temp)] <- 0
  
  # biomass loop for summing trees to plots
  # We're doing the unit conversions here; we had calculated density in stems/ha, but Christy wants to look at Biomass in kg/m2, so we're putting everything in kg/m2 here
  for(p in 1:length(plots)){
    cols <- which(names(allom.temp) %in% trees.use[trees.use$PlotID==plots[p], "TreeID"])
    if(substr(plots[p],1,1)=="V"){
      bm.array[,p,i] <- rowMeans(allom.temp[,cols])*plot.data[plot.data$PlotID==paste(plots[p]), "Density.Total..stems.ha."]/10000 #mean tree * trees/ha (do for Valles only bc sum of trees != plot density; different sampling method than Neil)
    } else {
      temp1 <- temp2 <- temp3 <- allom.temp[,cols]
      for(t in names(temp)){ # Convert biomass/tree to biomass/ha
        temp1[,t] <- temp1[,t] * mort.list[[plots[p]]]
        temp2[,t] <- temp2[,t] * mort.list[[plots[p]]][2,]
        temp3[,t] <- temp3[,t] * mort.list[[plots[p]]][1,]
      }
      bm.array.mean[,p,i] <- rowSums(temp1) #sum biomass/ha
      bm.array.hi[,p,i] <- rowSums(temp2) #sum biomass/ha
      bm.array.lo[,p,i] <- rowSums(temp3) #sum biomass/ha
    }
  }

#--------------------------------------------------

