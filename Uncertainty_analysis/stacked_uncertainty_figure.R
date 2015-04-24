library(ggplot2)
##############################################################################
# wanting to make a composite figure showing the contribution of uncertainties
##############################################################################

# loading in the various datasets that will be needed

#allometric mean and cloud per site
load("valles_bm_recon_site_stack.Rdata")

ggplot(data=site.valles.stack[site.valles.stack$Year<2012 ,])  + facet_grid(SiteID~.) +
  # plotting total site basal area  
  geom_ribbon(aes(x=Year, ymin=site.LB, ymax=site.UB, fill=SiteID), alpha=0.5) +
  geom_line(aes(x=Year, y=site.Mean, color=SiteID)) +
  ggtitle("Valles Caldera mean Densities")



# density BM--uses mean allometric eqtn. and accounts for differences in density with just ROSS plots
load("allTrees_valles_sample_stack.Rdata")

ggplot(data=all.valles.sample.stack[all.valles.sample.stack$Year<2012,])  + facet_grid(SiteID ~.) +
  # plotting total site basal area  
  geom_ribbon(aes(x=Year, ymin=all.sample.LB, ymax=all.sample.UB, fill=SiteID), alpha=0.5) +
  geom_line(aes(x=Year, y=all.Mean.sample, color=SiteID))+
  labs(title= "All Trees BM sample", x="Year", y="Bm sample. (kg/m2)")


# Need to get an error for the difference between dated and undated trees (Christy Help!)



#combine the different areas into one figure