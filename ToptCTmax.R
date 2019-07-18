#Compare CTmax and Topt

#Dataset notes
#ROHR et al.- Estimates Topt from Dell data but doesn't include CTmax
#BACTERIA- Knies et al.- Bacteria, has data and compiled other data, but data not published. Not geographic. Ask Joel?

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

#Use max growth data to calculate as slope
#mu.g.opt.val.list = estimated maximum specific growth rate (per day) based on the thermal reaction norm model fit
slope= plank$mu.g.opt.val.list/plank$dbreadth
plot(abs(plank$isolation.latitude), slope, log="y")

#---
#LIZARDS- Huey
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
liz= read.csv('Hueyetal2009.csv', na.strings ='-9999')

#Breadth of declining part of TPC: Distance from Topt to CTmax
liz$dbreadth= liz$CTmax - liz$newTopt
plot(liz$AbsLat, liz$dbreadth)

#---
#INSECTS- Deutsch et al.
#Load Deutsch Data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/DeutschData/")
ins= read.csv('Deutschetel.2008Insect.TPCdata.csv')

#Breadth of declining part of TPC: Distance from Topt to CTmax
ins$dbreadth= ins$CTmax - ins$Topt
plot(abs(ins$Lat), ins$dbreadth)



