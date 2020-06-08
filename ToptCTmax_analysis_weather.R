#Load packages
library(dplyr)
library(rerddap)
library(ncdf4)
library(tidyverse)
library(heatwaveR)
library(RNCEP)
library(lubridate) #date and time manipulation
library(tidyverse) #data manipulation and visualization
library(RColorBrewer) #color schemes
library(sf) #to import a spatial object and to work with geom_sf in ggplot2
library(raster)
library(zoo)
library(cowplot)
library(plyr)
library(patchwork)

#convert longituide
convert.lon= function(r0) ifelse(r0 > 180, -360 + r0, r0)

#----------------------
#LOAD DATA

#SEE PREVIOUS VERSIONS FOR OTHER DATA SOURCES

#Load combined NCEP data
#FROM https://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NCEP+Reanalysis&Variable=Surface+lifted+index&group=0&submit=Search
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/NCEP/")

ncep.nc= nc_open('air.sig995.2015.nc')
#ncep.nc= nc_open('lftx.sfc.2015.nc')
print(ncep.nc)

#extract metrics
ncep.temp=ncvar_get(ncep.nc,"air")
dim(ncep.temp) #dimensions lon, lat, time

#extract lon, lat, time
ncep.lons= ncvar_get(ncep.nc,"lon") #get info about long
dim(ncep.lons)
ncep.lats= ncvar_get(ncep.nc,"lat") #get info about latitude
dim(ncep.lats)
ncep.times= ncvar_get(ncep.nc,"time") #julian date, calendar day, since 0000-01-01
dim(ncep.times)
#change to dates
ncep.dates= as.POSIXct(ncep.times*3600,origin='1800-01-01 00:00') 
year= as.numeric(format(ncep.dates, "%Y"))
doy= as.numeric(format(ncep.dates, "%j"))
hours= as.numeric(format(ncep.dates, "%H"))

#-----------------
#load and combine data across years
years= 2000:2018

#make array to store data
ncep.temp.yrs= array(NA, dim=c(length(years),length(ncep.lons),length(ncep.lats),length(ncep.times)  ))

for(year.k in 1:length(years)){
  
  file.name= paste('air.sig995.', years[year.k], '.nc', sep='')
  ncep.nc= nc_open(file.name)
  
  #extract metrics
  ncep.temp=ncvar_get(ncep.nc,"air") -273.15
  dim(ncep.temp) 
  ncep.temp.yrs[year.k,,,]= ncep.temp[,,1:1460] 
}

#---
#to raster

ncep.lons.neg= convert.lon(ncep.lons)

ncep.r <- raster(t(ncep.temp[,,1]), xmn=min(ncep.lons), xmx=max(ncep.lons), ymn=min(ncep.lats), ymx=max(ncep.lats), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

#plot
ncep.r <- flip(ncep.r, direction='y')
ncep.r <- rotate(ncep.r)
plot(ncep.r)
#add country outline
map('world', fill = FALSE, col = "grey", add=TRUE)
#=====================================
#Functions
#Convert to thermodynamic scale
#E=0.757 #eV
#k= 1.38*10^-23 #J K^-1
#k=8.617* 10^-5 #eV K^-1
#function using temp in C
thermo.temp= function(t, E=0.757, k=8.617* 10^-5) exp(-E/(k*(t+273.15))) 

#ESTIMATE IN Seasonality_TTB.R
#TTB, based on GlobTherm
#terrestrial
TTB.terr= function(AbsLat) 29.15383 + 0.19897 * AbsLat
#marine
TTB.mar= function(AbsLat) 0.6945813 + 0.0020061 * AbsLat
#TTB= 26.25588 + 0.09722 * AbsLat

#Topt as percent of TTB, based on Huey data
ToptPer= function(AbsLat) 0.6945813 + 0.0020061 * AbsLat

#METABOLIC FUNTIONS
#duration of days above Topt
maxc.overTopt= function(x){ 
  # get rle object
  consec <- rle(x>0)
  # max consecutive values
  max(consec$lengths[consec$values==TRUE], na.rm=TRUE)
}
meanc.overTopt= function(x){ 
  # get rle object
  consec <- rle(x>0)
  # mean consecutive values
  mean(consec$lengths[consec$values==TRUE], na.rm=TRUE)
}

#count
count=function(x) length(na.omit(x))

#=====================================
#LOAD BIOLOGICAL DATA

#Topt data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
tol.h=read.csv("tpcs.csv")

#species with lat lon
tol.h= subset(tol.h, !is.na(tol.h$lat) & !is.na(tol.h$lon) )
#drop data without all metrics
tol.h= subset(tol.h, !is.na(tol.h$CTmin) & !is.na(tol.h$Topt) & !is.na(tol.h$CTmax) )

#add asym
tol.h$asym= (2*tol.h$Topt-tol.h$CTmax - tol.h$CTmin)/(tol.h$CTmax-tol.h$CTmin )
#===================================================
#Estimate thermal stress with increasing information

#fit assymetry vs Topt by taxa
models <- dlply(tol.h, "taxa", function(df) 
  lm(asym ~ Topt, data = df))

# Apply coef to each model and return a data frame
mods=ldply(models, coef)

# Print the summary of each model
#l_ply(models, summary, .print = TRUE)

#add Topt with no asymetry
tol.h$Topt.noasym= tol.h$CTmin + (tol.h$CTmax - tol.h$CTmin)/2

#add estimated Topt
mod.k= match(tol.h$taxa, mods$taxa)
asym= mods[mod.k,"(Intercept)"]+mods[mod.k,"Topt"]*tol.h[,"Topt"]
tol.h$Topt.aveasym= (asym*(tol.h$CTmax - tol.h$CTmin) +tol.h$CTmax +tol.h$CTmin)/2
 
#===================================================
#THERMAL STRESS ESTIMATES

ts= array(NA, dim= c(length(years), nrow(tol.h), 365,6) )

#calculate degree days above Topt
for(spec.k in 1:nrow(tol.h)){
  
  #find closest grid cell
  lon.ind= which.min(abs(ncep.lons.neg - tol.h[spec.k, "lon" ]))
  lat.ind= which.min(abs(ncep.lats - tol.h[spec.k, "lat" ]))
  
  #loop years
  for(year.k in 1:length(years)){
    
    #extract data
    #daily min max
    ncep.cell= cbind(year, doy, hours, ncep.temp.yrs[year.k,lon.ind,lat.ind,])
    colnames(ncep.cell)[4]="temp"
    ncep.cell= as.data.frame(ncep.cell)
    #make day factor
    ncep.cell$doy= as.factor(ncep.cell$doy)
    
    # #daily min and max
    # ncep.cmm= ncep.cell %>%
    #   group_by(doy) %>%
    #   summarise(min = min(temp), max= max(temp))
    # ncep.cmm= as.matrix(ncep.cmm)
    # tmax.k= as.numeric( ncep.cmm[,"max"] )
    # tmin.k= as.numeric( ncep.cmm[,"min"] )
    
    tmax.k= tapply(ncep.cell$temp, ncep.cell$doy, max)
    tmin.k= tapply(ncep.cell$temp, ncep.cell$doy, min)
    
    #daily safety margins
    ts[year.k,spec.k,,1]= tol.h[spec.k,'CTmax']-tmax.k
    #ts[year.k,spec.k,,2]= tol.h[spec.k,'Topt']-tmax.k 
    
    #metabolic scaled thermal stress
    inds= which(tmax.k > tol.h[spec.k,'Topt'])
    
    if(length(inds)>0){ 
      ts[year.k,spec.k,inds,2]= tmax.k[inds]- tol.h[spec.k,'Topt']  
      
      #performance detriment
      #Topt no asymetry
      ts[year.k,spec.k,inds,3]= 1- tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt.noasym'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax'])  
      
      #Topt average asymetry
      ts[year.k,spec.k,inds,4]= 1- tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt.aveasym'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax'])  
      
      #actual Topt
      ts[year.k,spec.k,inds,5]= 1- tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax'])  
      
      #number of days with >= 50% loss of performance
      perf= tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax']) 
      t50= rep(NA, length(inds))
      t50[perf<0.5]<-1
      ts[year.k,spec.k,inds,6]= t50
      
    } #end check >Topt
    
    # #thermodynamic scale
    # Topt.tt= thermo.temp(tol.h[spec.k,'topt'])
    # tmax.k.tt= thermo.temp(tmax.k)
    # 
    # inds= which(tmax.k.tt > Topt.tt)
    # 
    # if(length(inds>0)) ts[year.k,spec.k,inds,4]= tmax.k.tt[inds]- Topt.tt  
    
  } # end year loop
  
} # end loop species

#------------
#TSM
tsm<-  array(NA, dim= c(length(years),nrow(tol.h),6) )

#reorder ts
ts2 <- aperm(ts, c(2,1,3,4))

#loop years
for(year.k in 1:length(years)){
  
  #min TSM
  tsm[year.k,,1]=t(apply(ts2[,year.k,,1], MARGIN=1, FUN='min', na.rm=TRUE))
  
  #performace detriment
  #no asym
  tsm[year.k,,3]= apply(ts2[,year.k,,3], MARGIN=1, FUN=sum, na.rm=T)
  #ave asym
  tsm[year.k,,4]= apply(ts2[,year.k,,4], MARGIN=1, FUN=sum, na.rm=T)
  #actual asym
  tsm[year.k,,5]= apply(ts2[,year.k,,5], MARGIN=1, FUN=sum, na.rm=T)
  
  #Replace -Inf with NA
  tsm[year.k,which(is.infinite(tsm[,year.k,2])),2]=NA
  tsm[year.k,which(is.infinite(tsm[,year.k,2])),6]=NA
  
} #end year loop

#AGGREGATE ACROSS YEARS
tsm.yrs= apply(tsm, MARGIN=c(2,3), FUN=mean, na.rm=T)

#----
#number of days with 50% loss of performance
d1= apply(ts2[,,,6], MARGIN=c(1,2), FUN=sum, na.rm=TRUE)
#average across years
days_p50= rowMeans(d1)

#----------------------
#PLOT
#comparison plots
tol2= cbind(tol.h, days_p50, tsm.yrs)
colnames(tol2)[c(16,18:20)]=c('minTSM','Perf.noAsym','Perf.aveAsym','Perf')
#drop fish
tol2=tol2[-which(tol2$taxa=="fish"),] 

#drop unneeded columns
tol2s= tol2[,c("taxa","asym","days_p50","minTSM","Perf.noAsym","Perf.aveAsym","Perf")]
#change names
names(tol2s)[5:6]=c("without asymetry","fitted asymetry")

#-----
#compare to TSM
fig4a= ggplot(tol2, aes(x=log(Perf),y=minTSM, color=asym)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+scale_color_viridis(name="asymmetry")+ theme(legend.position = "bottom")+
  xlab("log annual performance detriment")+ylab("annual minimum of daily TSM (°C)")+
  ylim(-10,15)

#-----
#Compare performance estimates
#to long format
tol.l <- melt(tol2s, id=c("taxa","asym","days_p50","minTSM","Perf"))

fig4b= ggplot(tol.l, aes(x=log(Perf),y=log(value), color=variable)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+ theme(legend.position = "bottom", legend.title = element_blank())+
  ylab("Estimated log annual performance detriment")+xlab("log annual performance detriment")+
  scale_color_viridis(discrete=TRUE)+geom_abline(slope=1, intercept=0)

#-----
#proportion days with 50% performance loss

fig4c= ggplot(tol2, aes(x=minTSM,y=days_p50/365, color=asym)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+xlim(-10,10)+ylim(0,0.65) +
  scale_color_viridis(name="asymmetry")+ theme(legend.position = "bottom",legend.key.width = unit(2, "cm"))+
  ylab("proportion days with 50% performance loss")+xlab("annual minimum of daily TSM (°C)")
#+geom_smooth(method='loess',se=TRUE)

#----
#Plot
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/figures/")
pdf("Figs4_TSM.pdf", height = 8, width = 8)

fig4a +fig4b +plot_annotation(tag_levels = 'a') +plot_layout(nrow=2) 

dev.off()

#----
#LATITUDINAL PLOT
#latitudinal figure for plankton
tol.p= tol2[which(tol2$taxa=="plankton"), c(1:12,15:16,18:20) ]
#convert days_p50 to proportion
tol.p$days_p50= tol.p$days_p50/365
#convert performance
tol.p$Perf= log(tol.p$Perf)
tol.p$Perf.noAsym= log(tol.p$Perf.noAsym)
tol.p$Perf.aveAsym= log(tol.p$Perf.aveAsym)

#to long format
tol.pl<- tol.p %>%
  gather("metric", "value", c("minTSM","Perf","Perf.noAsym","Perf.aveAsym") ) #"days_p50"
#make labels
tol.pl$metric.lab<-NA
tol.pl$metric.lab[tol.pl$metric=="days_p50"]<- "proportion days with 50% performance loss"
tol.pl$metric.lab[tol.pl$metric=="minTSM"]<- "annual minimum of daily TSM"
tol.pl$metric.lab[tol.pl$metric=="Perf"]<- "observed" #"log annual performance detriment"
tol.pl$metric.lab[tol.pl$metric=="Perf.noAsym"]<- "without asymetry"
tol.pl$metric.lab[tol.pl$metric=="Perf.aveAsym"]<- "average asymetry"
tol.pl$metric.lab= factor(tol.pl$metric.lab, levels=c("annual minimum of daily TSM","observed","without asymetry","average asymetry"))

#TSM plot
fig5a= ggplot(tol.p, aes(x=abs(lat),y=minTSM, color=asym) ) +geom_point() +geom_smooth(method='loess',se=TRUE) +
  theme_bw()+scale_color_viridis(name="asymmetry")+ theme(legend.position = "bottom",legend.key.width = unit(2, "cm"))+
  xlab("absolute latitude (°)")+ylab("annual minimum of daily TSM (°C)")+
  ylim(-10,15)

#Perf plot
tol.p2= tol.pl[which(tol.pl$metric %in% c("Perf","Perf.noAsym")),]

fig5b= ggplot(tol.p2, aes(x=abs(lat),y=value, color=metric.lab) ) +geom_point()+geom_smooth(method='loess',se=TRUE) +
  theme_bw()+scale_color_viridis(name="", discrete=TRUE)+ theme(legend.position = "bottom",legend.key.width = unit(2, "cm"))+
  xlab("absolute latitude (°)")+ylab("log annual performance detriment")

pdf("Figs5_TSMlat.pdf", height = 8, width = 8)
fig5a +fig5b +plot_annotation(tag_levels = 'a') +plot_layout(nrow=2) 
dev.off()




#OLD VERSION***********************************
#THERMAL STRESS ESTIMATES

ts= array(NA, dim= c(length(years), nrow(tol.h), 365,6) )

#calculate degree days above Topt
for(spec.k in 1:nrow(tol.h)){
  
  #find closest grid cell
  lon.ind= which.min(abs(ncep.lons.neg - tol.h[spec.k, "lon" ]))
  lat.ind= which.min(abs(ncep.lats - tol.h[spec.k, "lat" ]))
  
  #loop years
  for(year.k in 1:length(years)){
    
    #extract data
    #daily min max
    ncep.cell= cbind(year, doy, hours, ncep.temp.yrs[year.k,lon.ind,lat.ind,])
    colnames(ncep.cell)[4]="temp"
    ncep.cell= as.data.frame(ncep.cell)
    #make day factor
    ncep.cell$doy= as.factor(ncep.cell$doy)
    
    #daily min and max
    ncep.cmm= ncep.cell %>%
      group_by(doy) %>%
      summarise(min = min(temp), max= max(temp))
    ncep.cmm= as.matrix(ncep.cmm)
    tmax.k= as.numeric( ncep.cmm[,"max"] )
    tmin.k= as.numeric( ncep.cmm[,"min"] )
    
    #daily safety margins
    ts[year.k,spec.k,,1]= tol.h[spec.k,'CTmax']-tmax.k
    ts[year.k,spec.k,,2]= tol.h[spec.k,'Topt']-tmax.k 
    
    #metabolic scaled thermal stress
    inds= which(tmax.k > tol.h[spec.k,'Topt'])
    
    if(length(inds)>0){ 
      ts[year.k,spec.k,inds,3]= tmax.k[inds]- tol.h[spec.k,'Topt']  
    
    #performance detriment
      ts[year.k,spec.k,inds,5]= 1- tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax'])  
    
    #number of days with >= 50% loss of performance
    perf= tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax']) 
    t50= rep(NA, length(inds))
    t50[perf<0.5]<-1
    ts[year.k,spec.k,inds,6]= t50
    
    } #end check >Topt
    
    #thermodynamic scale
    Topt.tt= thermo.temp(tol.h[spec.k,'topt'])
    tmax.k.tt= thermo.temp(tmax.k)
    
    inds= which(tmax.k.tt > Topt.tt)
    
    if(length(inds>0)) ts[year.k,spec.k,inds,4]= tmax.k.tt[inds]- Topt.tt  
    
  } # end year loop
  
} # end loop species

#------------
#TSM
tsm<-  array(NA, dim= c(length(years),nrow(tol.h),13,4) )

#reorder ts
ts <- aperm(ts, c(2,1,3,4))

#loop years
for(year.k in 1:length(years)){
  
  #TSM daily: min, 10%, 50%, 90%
  #CTmax-Tmax
  tsm[year.k,,1:4,1]=t(apply(ts[,year.k,,1], MARGIN=1, FUN='quantile', probs=c(0,0.1,0.5,0.9), na.rm=TRUE))
  #Topt-Tmax
  tsm[year.k,,1:4,2]= t(apply(ts[,year.k,,2], MARGIN=1, FUN='quantile', probs=c(0,0.1,0.5,0.9), na.rm=TRUE))
  
  #lowest 7 day average
  tsm[year.k,,5,1]= apply(rollmean(ts[,year.k,,1], 7, fill=NA), MARGIN=1, FUN='min', na.rm=TRUE)
  tsm[year.k,,5,2]= apply(rollmean(ts[,year.k,,2], 7, fill=NA), MARGIN=1, FUN='min', na.rm=TRUE)
  #lowest 14 day average
  tsm[year.k,,6,1]= apply(rollmean(ts[,year.k,,1], 14, fill=NA), MARGIN=1, FUN='min', na.rm=TRUE)
  tsm[year.k,,6,2]= apply(rollmean(ts[,year.k,,2], 14, fill=NA), MARGIN=1, FUN='min', na.rm=TRUE)
  
  #count of days <5
  count.sub5= function(x) sum(x<5)
  tsm[year.k,,7,1]= apply(rollmean(ts[,year.k,,1], 14, fill=NA), MARGIN=1, FUN='count.sub5')
  tsm[year.k,,7,2]= apply(rollmean(ts[,year.k,,2], 14, fill=NA), MARGIN=1, FUN='count.sub5')
  
  #Metabolic integration
  tsm[year.k,,8,3]= apply(ts[,year.k,,3], MARGIN=1, FUN='maxc.overTopt')
  tsm[year.k,,8,4]= apply(ts[,year.k,,4], MARGIN=1, FUN='maxc.overTopt')
  #Replace -Inf with NA
  tsm[year.k,which(is.infinite(tsm[,year.k,8,3])),8,3]=NA
  tsm[year.k,which(is.infinite(tsm[,year.k,8,4])),8,4]=NA
  
  #-----
  #aggregate
  #sum, count, mean
  tsm[year.k,,9,3]= apply(ts[,year.k,,3], MARGIN=1, FUN=sum, na.rm=T)
  tsm[year.k,,9,4]= apply(ts[,year.k,,4], MARGIN=1, FUN=sum, na.rm=T)
  
  tsm[year.k,,10,3]= apply(ts[,year.k,,3], MARGIN=1, FUN=count)
  tsm[year.k,,10,4]= apply(ts[,year.k,,4], MARGIN=1, FUN=count)
  
  tsm[year.k,,11,3]= apply(ts[,year.k,,3], MARGIN=1, FUN=mean, na.rm=T)
  tsm[year.k,,11,4]= apply(ts[,year.k,,4], MARGIN=1, FUN=mean, na.rm=T)
  
  #performance detriment
  tsm[year.k,,12,3]= apply(ts[,year.k,,5], MARGIN=1, FUN=sum, na.rm=T)
  
  #Replace -Inf with NA
  tsm[year.k,which(is.infinite(tsm[,year.k,9,3])),9,3]=NA
  tsm[year.k,which(is.infinite(tsm[,year.k,9,4])),9,4]=NA
  tsm[year.k,which(is.infinite(tsm[,year.k,10,3])),10,3]=NA
  tsm[year.k,which(is.infinite(tsm[,year.k,10,4])),10,4]=NA
  tsm[year.k,which(is.infinite(tsm[,year.k,11,3])),11,3]=NA
  tsm[year.k,which(is.infinite(tsm[,year.k,11,4])),11,4]=NA
  tsm[year.k,which(is.infinite(tsm[,year.k,12,3])),12,3]=NA
} #end year loop

#number of days with 50% loss of performance
d1= apply(ts[,,,6], MARGIN=c(1,2), FUN=sum, na.rm=TRUE)
#average across years
days_p50= rowMeans(d1)

dat_p50= cbind(tol.h,days_p50)

plot.p50= ggplot(dat_p50, aes(x=abs(lat),y=days_p50) ) +geom_point()+facet_wrap(~taxa) +geom_smooth(method='loess',se=TRUE)

#----------------------
#AGGREGATE ACROSS YEARS
tsm.yrs= apply(tsm, MARGIN=c(2,3,4), FUN=mean, na.rm=T)

#===============================
#PLOT
#comparison plots
tol2= cbind(tol.h, days_p50, tsm.yrs[,1:7,1],tsm.yrs[,8:12,3])
colnames(tol2)[14:20]=c('minTSM','TSM10p','TSM50p', 'TSM90p', 'low7d', 'low14d', 'count5' )
colnames(tol2)[21:25]=c('dTopt','sumI','countI','meanI','Perf')
#drop fish
tol2=tol2[-which(tol2$taxa=="fish"),] 

fig4a= ggplot(tol2, aes(x=minTSM,y=log(Perf), color=asym)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+scale_color_viridis()+xlim(-10,10)+ theme(legend.position = "none")+
  ylab("log annual performance detriment")+xlab("annual minimum of daily TSM (°C)")

fig4b= ggplot(tol2, aes(x=minTSM,y=days_p50/365, color=asym)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+xlim(-10,10)+ylim(0,0.65) +
  scale_color_viridis(name="asymmetry")+ theme(legend.position = "bottom",legend.key.width = unit(2, "cm"))+
  ylab("proportion days with 50% performance loss")+xlab("annual minimum of daily TSM (°C)")
#+geom_smooth(method='loess',se=TRUE)

#Plot
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/figures/")
pdf("Figs4_TSM.pdf", height = 8, width = 8)
plot_grid(fig4a,fig4b, labels = c("A", "B"), ncol = 1)
dev.off()

#----
#latitudinal figure for plankton
tol.p= tol2[which(tol2$taxa=="plankton"), c(1:12,13,14,25) ]
#convert days_p50 to proportion
tol.p$days_p50= tol.p$days_p50/365
#convert performance
tol.p$Perf= log(tol.p$Perf)

#to long format
tol.pl<- tol.p %>%
  gather("metric", "value", 13:15)
#make labels
tol.pl$metric.lab<-NA
tol.pl$metric.lab[tol.pl$metric=="days_p50"]<- "proportion days with 50% performance loss"
tol.pl$metric.lab[tol.pl$metric=="minTSM"]<- "annual minimum of daily TSM"
tol.pl$metric.lab[tol.pl$metric=="Perf"]<- "log annual performance detriment"
tol.pl$metric.lab= factor(tol.pl$metric.lab, levels=c("annual minimum of daily TSM","log annual performance detriment","proportion days with 50% performance loss"))

fig5= ggplot(tol.pl, aes(x=abs(lat),y=value, color=asym) ) +geom_point()+facet_wrap(~metric.lab,ncol=1,scales="free_y") +geom_smooth(method='loess',se=TRUE) +
  theme_bw()+scale_color_viridis(name="asymmetry")+ theme(legend.position = "bottom",legend.key.width = unit(2, "cm"))+
  xlab("absolute latitude (°)")

pdf("Figs5_TSMlat.pdf", height = 12, width = 8)
fig5
dev.off()




