#Compare CTmax and Topt
library(ggplot2)

#Similar concept (plus fish data): https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12707

#Dataset notes
#ROHR et al.- Estimates Topt from Dell data but doesn't include CTmax, Can recaluculate or ask Jason
#BACTERIA- Knies et al.- Bacteria, has data and compiled other data, but data not published. Not geographic. Ask Joel?

#OTHER POTENTIAL DATA
#SEA URCHINS: https://www.int-res.com/articles/meps2018/589/m589p153.pdf
#MARINE DIATOM SELECTION: https://www.biorxiv.org/content/biorxiv/early/2017/07/24/167817.full.pdf
#FISH: https://core.ac.uk/download/pdf/51490125.pdf, https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2435.12618
#AUST LIZARDS: https://onlinelibrary.wiley.com/doi/full/10.1111/evo.13064
#SEA CUCUMBER: https://d-nb.info/1147545456/34#page=91
#INSECT DEVELOPMENT? https://onlinelibrary.wiley.com/doi/full/10.1111/eea.12693
#AMPHIBIANS TABLE 2: https://www.imperial.ac.uk/media/imperial-college/faculty-of-natural-sciences/department-of-life-sciences/public/postgraduate/masters/cmee/InvestigatingClimateChangeExtinctionRisksOfAmphibiansBySimulation.pdf
#REEF FISH: https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.12455
#SCELOPORUS POPS: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2435.12406
#FRESH WATER?: https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2427.2004.01317.x
#ABALONE: https://www.sciencedirect.com/science/article/pii/S0306456599000327
#LIZARDS: https://www.journals.uchicago.edu/doi/full/10.1086/660021

#---
#PLANKTON- Thomas
#Thomas data downloaded from his website, https://mridulkthomas.weebly.com/data--code.html
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/Phytoplankton_temperature_growth_rate_dataset_2016_01_29/")
plank= read.csv('traits_derived_2016_01_29.csv')
#Topt: plank$mu.g.opt.list
#CTmax: plank$tmax

#Plot Topt vs Ctmax
plot(plank$mu.g.opt.list, plank$tmax)

#TPC breadth
plot(abs(plank$isolation.latitude), plank$tmax-plank$tmin, log="y")

#Breadth of declining part of TPC: Distance from Topt to CTmax
plank$dbreadth= plank$tmax - plank$mu.g.opt.list
plot(abs(plank$isolation.latitude), plank$dbreadth, log="y")
#plot by CTmax, TMion, Topt
plot(plank$tmax, plank$dbreadth, log="y", ylab="CTmax-Topt",xlab="CTmax" )
plot(plank$tmin, plank$dbreadth, log="y", ylab="CTmax-Topt",xlab="CTmin" )
plot(plank$mu.g.opt.list, plank$dbreadth, log="y", ylab="CTmax-Topt",xlab="Topt" )

#vs TTB
plot(plank$tmax-plank$tmin, plank$dbreadth, log="y", ylab="CTmax-Topt")
#Calculate proportion
plot(abs(plank$isolation.latitude), plank$mu.g.opt.list/plank$tmax)

#Use max growth data to calculate as slope
#mu.g.opt.val.list = estimated maximum specific growth rate (per day) based on the thermal reaction norm model fit
slope= plank$mu.g.opt.val.list/plank$dbreadth
plot(abs(plank$isolation.latitude), slope, log="y")
#plot by CTmax
plot(plank$tmax, slope, log="y")

#---
#LIZARDS- Huey
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
liz= read.csv('Hueyetal2009.csv', na.strings ='-9999')
plot(liz$newTopt,liz$CTmax)

#Breadth of declining part of TPC: Distance from Topt to CTmax
liz$dbreadth= liz$CTmax - liz$newTopt
plot(liz$AbsLat, liz$dbreadth)
#by CTMax, CTmin, Topt
plot(liz$CTmax, liz$dbreadth, ylab="CTmax-Topt",xlab="CTmax" )
plot(liz$CTmin, liz$dbreadth, ylab="CTmax-Topt",xlab="CTmin" )
plot(liz$newTopt, liz$dbreadth, ylab="CTmax-Topt",xlab="Topt" )
#vs TTB
plot(liz$CTmax-liz$CTmin,liz$dbreadth)

#by family
liz$Family= as.factor(liz$Family)
ggplot(liz) + aes(x=CTmax, y = dbreadth, color=Family, 
      group=Family)+geom_point()+geom_smooth(method="lm", se=FALSE)

#---
#INSECTS- Deutsch et al.
#Load Deutsch Data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/DeutschData/")
ins= read.csv('Deutschetel.2008Insect.TPCdata.csv')

#Breadth of declining part of TPC: Distance from Topt to CTmax
ins$dbreadth= ins$CTmax - ins$Topt
plot(abs(ins$Lat), ins$dbreadth)
#plot by CTmax, CTMin
plot(ins$CTmax, ins$dbreadth)
plot(ins$Ctmin, ins$dbreadth)

#---
#LIZARD Tp
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
liz2= read.csv('SpeciesPerformanceData_SCT.csv')

plot(liz2$Tp,liz2$CTmax)

#Breadth of declining part of TPC: Distance from Topt to CTmax
liz2$dbreadth= liz2$CTmax - liz2$Tp
#by CTMax, CTmin, Topt
plot(liz2$CTmax, liz2$dbreadth, ylab="CTmax-Topt",xlab="CTmax" )
plot(liz2$CTmin, liz2$dbreadth, ylab="CTmax-Topt",xlab="CTmin" )
plot(liz2$Tp, liz2$dbreadth, ylab="CTmax-Topt",xlab="Tp" )

