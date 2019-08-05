#Compare CTmax and Topt
library(ggplot2)

#Similar concept (plus fish data): https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12707

#Dataset notes
#ROHR et al.- Estimates Topt from Dell data but doesn't include CTmax, Can recaluculate or ask Jason
#BACTERIA- Knies et al.- Bacteria, has data and compiled other data, but data not published. Not geographic. Ask Joel?

#OTHER POTENTIAL DATA
#SEA URCHINS: https://www.int-res.com/articles/meps2018/589/m589p153.pdf
#MARINE DIATOM SELECTION: https://www.biorxiv.org/content/biorxiv/early/2017/07/24/167817.full.pdf
#AUST LIZARDS: https://onlinelibrary.wiley.com/doi/full/10.1111/evo.13064
#SEA CUCUMBER: https://d-nb.info/1147545456/34#page=91
#INSECT DEVELOPMENT? https://onlinelibrary.wiley.com/doi/full/10.1111/eea.12693
#AMPHIBIANS TABLE 2: https://www.imperial.ac.uk/media/imperial-college/faculty-of-natural-sciences/department-of-life-sciences/public/postgraduate/masters/cmee/InvestigatingClimateChangeExtinctionRisksOfAmphibiansBySimulation.pdf
#REEF FISH: https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.12455
#SCELOPORUS POPS: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2435.12406
#FRESH WATER?: https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2427.2004.01317.x
#ABALONE: https://www.sciencedirect.com/science/article/pii/S0306456599000327
#LIZARDS: https://www.journals.uchicago.edu/doi/full/10.1086/660021
#ANTS (data available?): https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/15-1225

#---------------
#Assemble data in simplified form for easy plotting

#PLANKTON- Thomas
#Thomas data downloaded from his website, https://mridulkthomas.weebly.com/data--code.html
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/Phytoplankton_temperature_growth_rate_dataset_2016_01_29/")
plank= read.csv('traits_derived_2016_01_29.csv')
#Topt: plank$mu.g.opt.list
#CTmax: plank$tmax

#subset to good fits
plank= subset(plank, plank$minqual=="good" & plank$maxqual=="good" & plank$curvequal=="good" )

#subset
tpc= plank[,c("species","genus","family","tmin","tmax","mu.g.opt.list","habitat")]
tpc$taxa="plankton"
names(tpc)[4:6]<- c("CTmin", "CTmax", "Topt")

#---
#LIZARDS- Huey
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
liz= read.csv('Hueyetal2009.csv', na.strings ='-9999')
tpc2= liz[,c("Species","Genus","Family","CTmin","CTmax","newTopt")]
tpc2$habitat="terrestrial"
tpc2$taxa="lizards"
#bind
tpc= rbind(tpc, setNames(tpc2, names(tpc)))

#---
#INSECTS- Deutsch et al.
#Load Deutsch Data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/DeutschData/")
ins= read.csv('Deutschetel.2008Insect.TPCdata.csv')

tpc2= ins[,c("Species","genus","Order","Ctmin","CTmax","Topt")]
tpc2$habitat="terrestrial"
tpc2$taxa="insects"
#bind
tpc= rbind(tpc, setNames(tpc2, names(tpc)))

#---

#LIZARD Tp
#Clusella-trullas et al., https://www.journals.uchicago.edu/doi/abs/10.1086/660021
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
liz2= read.csv('SpeciesPerformanceData_SCT.csv')
liz2$family=NA

tpc2= liz2[,c("Species","genus","family","CTmin","CTmax","Tp")]
tpc2$habitat="terrestrial"
tpc2$taxa="lizards_Tp"
#bind
tpc= rbind(tpc, setNames(tpc2, names(tpc)))

#---
#FISH: https://core.ac.uk/download/pdf/51490125.pdf, https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2435.12618
#CHECK DATA
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
fish= read.csv('Fish_Payne_Smith_EL.csv')
fish$family=NA
fish$CTmin=NA

tpc2= fish[,c("Species","Genus","family","CTmin","CTmax...C.","Topt...C.")]
tpc2$habitat="aquatic"
tpc2$taxa="fish"
#bind
tpc= rbind(tpc, setNames(tpc2, names(tpc)))

#--------
#Write out
write.csv(tpc, "tpcs.csv")


