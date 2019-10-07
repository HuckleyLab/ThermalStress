#Compare CTmax and Topt
library(ggplot2)

#Similar concept (plus fish data): https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12707

#Dataset notes
#BACTERIA- Knies et al.- Bacteria, has data and compiled other data, but data not published. Not geographic. Ask Joelfor a,b,m of thermal reaction norm? https://www.journals.uchicago.edu/doi/full/10.1086/597224.
#ANTS (data available?): https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/15-1225, Ask Mike?

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

#AMPHIBIAN LARVAE? https://onlinelibrary.wiley.com/doi/full/10.1111/evo.12711
#PLANT POPULATIONS? https://academic.oup.com/icb/article/51/5/733/627422
#Algae: https://royalsocietypublishing.org/doi/full/10.1098/rspb.2018.1076


#---------------
#Assemble data in simplified form for easy plotting

#PLANKTON- Thomas
#Thomas data downloaded from his website, https://mridulkthomas.weebly.com/data--code.html, https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12387
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/Phytoplankton_temperature_growth_rate_dataset_2016_01_29/")
plank= read.csv('traits_derived_2016_01_29.csv')
#Topt: plank$mu.g.opt.list
#CTmax: plank$tmax

#subset to good fits
plank= subset(plank, plank$minqual=="good" & plank$maxqual=="good" & plank$curvequal=="good" )

#subset
tpc= plank[,c("species","genus","family","tmin","tmax","mu.g.opt.list","habitat", "isolation.latitude", "isolation.longitude")]
tpc$taxa="plankton"
names(tpc)[4:6]<- c("CTmin", "CTmax", "Topt")
names(tpc)[8:9]<- c("lat", "lon")

#analysis
#Use max growth data to calculate as slope
#mu.g.opt.val.list = estimated maximum specific growth rate (per day) based on the thermal reaction norm model fit
slope= plank$mu.g.opt.val.list/(plank$tmax -plank$mu.g.opt.list)
plot(plank$mu.g.opt.list, slope, log="y")
#plot by CTmax and optima
par(mfrow=c(1,2))
plot(plank$tmax, slope, log="y")

#---
#LIZARDS- Huey
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
liz= read.csv('Hueyetal2009.csv', na.strings ='-9999')
tpc2= liz[,c("Species","Genus","Family","CTmin","CTmax","newTopt")]
tpc2$habitat="terrestrial"
tpc2= cbind(tpc2, liz[,c("Lat","Long")])
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
tpc2= cbind(tpc2, ins[,c("Lat","Long")])
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
tpc2$lat= NA
tpc2$lon= NA
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
tpc2$lat= NA
tpc2$lon= NA
tpc2$taxa="fish"
#bind
tpc= rbind(tpc, setNames(tpc2, names(tpc)))

#---
#DROSOPHILA: https://royalsocietypublishing.org/doi/full/10.1098/rstb.2018.0548
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
fly= read.csv('DrosophilaTopt_Maclean.csv')
fly$family=NA
fly$Genus=NA

tpc2= fly[,c("Species","Genus","family","CTMin","CTMax","Fitness.Topt")] #also EgglayingTopt
tpc2$habitat="terrestrial"
tpc2= cbind(tpc2, fly[,c("latitude")])
tpc2$lon= NA
tpc2$taxa="flies"
#bind
tpc= rbind(tpc, setNames(tpc2, names(tpc)))

#---
#PHAGE, https://www.journals.uchicago.edu/doi/full/10.1086/597224

# a -height
# m- location of optimum
# b- width

phage.growth= function(temp, a,b,m,P=1) {a +(1/b)*P*(1/b*(temp-m))}

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/phage/")
phage= read.csv('TMV_output_Knies2009.csv')

ggplot(phage) + aes(x=1:60, y = phage.growth(1:60,a,b,m), color=Genotype, group=Genotype)+geom_line()

for(ind in 1:nrow(phage)){
if(ind==1) plot(1:60, phage.growth(1:60, phage[ind, "a"], phage[ind, "b"], phage[ind, "m"]), type="l")
  points(1:60, phage.growth(1:60, phage[ind, "a"], phage[ind, "b"], phage[ind, "m"]), type="l")
}
  
#---
#OTHERS GATHERED FROM LITERATURE
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/ToptAssembly/")
dat= read.csv('Topt_Plantsetc.csv')
#just add lizards
dat<- subset(dat, dat$Taxa=="Australian lizards")
dat$Taxa<- "lizards"

dat$family=NA
dat$habitat<- "terrestrial"
#dat$habitat[which(dat$Taxa=="Sea Urchins")]<-"marine"

tpc2= dat[,c("Species","Genus","family","CTmin","CTmax","Topt","habitat","Latitide","Longitude","Taxa")] #also EgglayingTopt

#bind
tpc= rbind(tpc, setNames(tpc2, names(tpc)))

#--------
#Write out
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
write.csv(tpc, "tpcs.csv")

#=====================
#DELL ET AL.: https://www.pnas.org/content/pnas/suppl/2011/05/19/1015178108.DCSupplemental/sapp.pdf
#ROHR et al.- Estimates Topt from Dell data but doesn't include CTmax, Can recaluculate or ask Jason
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
dell= read.csv('DelletalTPC.csv')

#use equation to estimate CTmin and CTmax
delleq= function(E,Ef,temp,Topt, c=1) {
  k=1.3806*10^(-23)
  c*exp(-E/(k*temp))/(1+exp(-1/(k*temp)*(Ef-(Ef/Topt+k*log(E/(Ef-E))*temp)) ))}

plot(1:50, delleq(E=dell[4,"Er"],Ef=dell[4,"Ef"],temp=1:50,Topt=dell[4,"Topt"]), type="l")

delleq(E=dell[3,"Er"],Ef=dell[3,"Ef"],temp=25,Topt=dell[3,"Topt"])

#-----------
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
dat.full<-read.csv('Delletal2013.csv')
#adjust names
dat.full$temperature= dat.full$AmbientTemp
dat.full$growth.rate=dat.full$TraitValueSI
dat.full$curve.id= dat.full$DataSeriesID
#drop records without temperatures
dat.full= dat.full[!is.na(dat.full$temperature),]

ids= unique(dat.full$curve.id)
#restrict to ids with at least 4 temperatures
tabs= table(dat.full$curve.id)
ids= as.numeric( names(tabs)[which(tabs>5)])

keep=NA
#check cases max is intermediate temperature
for(k in ids){
  dat=dat.full[which(dat.full$curve.id==k),]
  #order by temperature 
  dat= dat[order(dat$temperature),]
  
  do.keep=which.max(dat$TraitValueSI)>2 & which.max(dat$TraitValueSI)<(nrow(dat)-1) & length(unique(dat$temperature))>3
  
  if(do.keep) keep=c(keep, k)
}
#Drop NA and subset data
ids=keep[2:length(keep)]
dat.sub= dat.full[which(dat.full$curve.id %in% ids),]

#normalize max to 1
dat.sub$trait=NA
for(k in 1:length(ids)){
  inds=which(dat.sub$curve.id==ids[k])
  #normalize max to 1
  dat.sub$trait[inds]= dat.sub$TraitValueSI[inds]/(max(dat.sub$TraitValueSI[inds]))
 }

##PLOT
for(k in 1:100){ #length(ids)
  inds=which(dat.sub$curve.id==ids[k])
  
  dat.sub2= dat.sub[inds,]
  dat.sub2= dat.sub2[order(dat.sub$temperature),]
  if(k==1) plot(dat.sub$temperature[inds], dat.sub$trait[inds], type="l", xlim=c(0,60))
  points(dat.sub$temperature[inds], dat.sub$trait[inds], type="l")
  }

##plot curves
library(ggplot2)
ggplot(dat.sub) + aes(x=temperature, y = trait, color=curve.id, group=curve.id)+ylim(0,1)+xlim(-5,60)+geom_smooth(se=FALSE)

#write out
write.csv(dat.sub, "Delletal2013_forfitting.csv")



