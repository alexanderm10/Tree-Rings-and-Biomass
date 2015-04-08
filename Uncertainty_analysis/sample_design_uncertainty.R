################################################################################
# assessing the differences between Marcy BM recon and Valles BM recon
################################################################################

load("marcy_bm_recon.Rdata") #loads in the biom.marcy
load("valles_bm_recon.Rdata") #loads in the biom.valles
summary(biom.marcy)

names(biom.marcy)

################################################################################
# making mean Site biomass for Marcy's PIPO data
################################################################################
names(biom.marcy)

marcy.pipo.bm <- data.frame(apply(biom.marcy[,c(5:8)], 1, mean))
marcy.pipo.ci <- apply(biom.marcy[,c(5:8)], 1, quantile, c(0.025, 0.975))
summary(marcy.pipo.ci)

marcy.pipo.lbound <- data.frame(marcy.pipo.ci[1,])
marcy.pipo.ubound<- data.frame(marcy.pipo.ci[2,])

pipo.marcy <- as.data.frame(c(marcy.pipo.bm, marcy.pipo.lbound, marcy.pipo.ubound))
summary(pipo.marcy)

names(pipo.marcy) <- c("bm.mean", "LB", "UB")
row.names(pipo.marcy) <- row.names(biom.marcy)
summary(pipo.marcy)

plot(pipo.marcy$bm.mean~ row.names(pipo.marcy), type="l", ylim=c(0, max(pipo.marcy$UB)))
  lines(pipo.marcy$LB~ row.names(pipo.marcy), lty="dashed")
  lines(pipo.marcy$UB~ row.names(pipo.marcy), lty="dashed")

################################################################################
# making mean Site biomass for Marcy's MCON data
################################################################################


marcy.mcon.bm <- data.frame(apply(biom.marcy[,c(1:4)], 1, mean))
marcy.mcon.ci <- apply(biom.marcy[,c(1:4)], 1, quantile, c(0.025, 0.975))
summary(marcy.mcon.ci)

marcy.mcon.lbound <- data.frame(marcy.mcon.ci[1,])
marcy.mcon.ubound<- data.frame(marcy.mcon.ci[2,])

mcon.marcy <- as.data.frame(c(marcy.mcon.bm, marcy.mcon.lbound, marcy.mcon.ubound))
summary(mcon.marcy)

names(mcon.marcy) <- c("bm.mean", "LB", "UB")
row.names(mcon.marcy) <- row.names(biom.marcy)
summary(mcon.marcy)


plot(mcon.marcy$bm.mean~ row.names(mcon.marcy), type="l", ylim=c(0, max(mcon.marcy$UB)))
lines(mcon.marcy$LB~ row.names(mcon.marcy), lty="dashed")
lines(mcon.marcy$UB~ row.names(mcon.marcy), lty="dashed")

################################################################################
# making mean Site biomass for Ross PIPO data
################################################################################
summary(biom.valles)
ross.pipo.bm <- data.frame(apply(biom.valles[,c(3:4)], 1, mean))
ross.pipo.ci <- apply(biom.valles[,c(3:4)], 1, quantile, c(0.025, 0.975))
summary(ross.pipo.ci)

ross.pipo.lbound <- data.frame(ross.pipo.ci[1,])
ross.pipo.ubound<- data.frame(ross.pipo.ci[2,])

pipo.ross <- as.data.frame(c(ross.pipo.bm, ross.pipo.lbound, ross.pipo.ubound))
summary(pipo.ross)

names(pipo.ross) <- c("bm.mean", "LB", "UB")
row.names(pipo.ross) <- row.names(biom.valles)
summary(pipo.ross)

plot(pipo.ross$bm.mean ~ row.names(pipo.ross), type="l", ylim=c(0, max(pipo.ross$UB)))
  lines(pipo.ross$LB~ row.names(pipo.ross), lty="dashed")
  lines(pipo.ross$UB~ row.names(pipo.ross), lty="dashed")


################################################################################
# making mean Site biomass for Ross MCON data
################################################################################

summary(biom.valles)
ross.mcon.bm <- data.frame(apply(biom.valles[,c(1:2)], 1, mean))
ross.mcon.ci <- apply(biom.valles[,c(1:2)], 1, quantile, c(0.025, 0.975))
summary(ross.mcon.ci)

ross.mcon.lbound <- data.frame(ross.mcon.ci[1,])
ross.mcon.ubound<- data.frame(ross.mcon.ci[2,])

mcon.ross <- as.data.frame(c(ross.mcon.bm, ross.mcon.lbound, ross.mcon.ubound))
summary(mcon.ross)

names(mcon.ross) <- c("bm.mean", "LB", "UB")
row.names(mcon.ross) <- row.names(biom.valles)
summary(mcon.ross)

plot(mcon.ross$bm.mean ~ row.names(mcon.ross), type="l", ylim=c(0, max(mcon.ross$UB)))
  lines(mcon.ross$LB~ row.names(mcon.ross), lty="dashed")
  lines(mcon.ross$UB~ row.names(mcon.ross), lty="dashed")


################################################################################
# Calculating a % difference between Ross Valles Data and Marcy Valles Data
################################################################################
plot(pipo.marcy$bm.mean ~ row.names(pipo.marcy), typ="l", ylim=range(pipo.diff))
  lines(pipo.ross$bm.mean ~ row.names(pipo.ross), type="l", col="red")

marcy.mean <- pipo.marcy[as.numeric(paste(row.names(pipo.marcy)))>=1920 & as.numeric(paste(row.names(pipo.marcy)))<=2011,"bm.mean"]
ross.mean <- pipo.ross[as.numeric(paste(row.names(pipo.ross)))>=1920 & as.numeric(paste(row.names(pipo.ross)))<=2011,"bm.mean"]

pipo.diff <- (ross.mean - marcy.mean)/((ross.mean))*100
names(pipo.diff)<- 2011:1920


plot(pipo.diff ~ names(pipo.diff), typ="l", ylim=range(pipo.diff))
  lines(pipo.ross$bm.mean ~ row.names(pipo.ross), type="l", col="red")

row.names(pipo.ross)
length(row.names(pipo.ross))
length(pipo.ross$bm.mean)

################################################################################
# graphing the different samplings at the two sites to show the differences and overlapping uncertainties
################################################################################

summary(pipo.ross)
summary(pipo.marcy)
summary(mcon.ross)
summary(mcon.marcy)

pipo.ross$Site <- as.factor("PIPO")
pipo.marcy$Site <- as.factor("PIPO")
mcon.ross$Site <- as.factor("MCON")
mcon.marcy$Site <- as.factor("MCON")
pipo.ross$Group <- as.factor("Ross")
pipo.marcy$Group <- as.factor("Marcy")
mcon.ross$Group <- as.factor("Ross")
mcon.marcy$Group <- as.factor("Marcy")
pipo.ross$Year <- as.numeric(paste(row.names(pipo.ross)))
pipo.marcy$Year <- as.numeric(paste(row.names(pipo.marcy)))
mcon.ross$Year <- as.numeric(paste(row.names(mcon.ross)))
mcon.marcy$Year <- as.numeric(paste(row.names(mcon.marcy)))

biomass.group <- rbind(pipo.ross, pipo.marcy, mcon.ross, mcon.marcy)
summary(biomass.group)


library(ggplot2)
ggplot(data=biomass.group[biomass.group$Year>=1920 & biomass.group$Year<=2011,]) + facet_grid(Site ~ .) +
  geom_ribbon(aes(x=Year, ymin=LB, ymax=UB, fill=Group), alpha=0.3)+
  geom_line(aes(x=Year, y=bm.mean, color=Group), size=1.5)