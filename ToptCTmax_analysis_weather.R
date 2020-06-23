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

ts= array(NA, dim= c(length(years), nrow(tol.h), 365,13) )
#counts of Topt exceedences
ts.exceed= array(NA, dim= c(length(years), nrow(tol.h), 2) )

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
    
    #find Topt exceedences
    inds.s= which(tmax.k > tol.h[spec.k,'Topt.noasym'])
    inds= which(tmax.k > tol.h[spec.k,'Topt'])
    
      ts[year.k,spec.k,inds,2]= tmax.k[inds]- tol.h[spec.k,'Topt']  
      
      #find Topt exceedences
      inds.s= which(tmax.k > tol.h[spec.k,'Topt.noasym'])
      inds= which(tmax.k > tol.h[spec.k,'Topt'])
      
      #count exceed Topt
      ts.exceed[year.k,spec.k,1]= length(inds.s)
      ts.exceed[year.k,spec.k,2]= length(inds)
      
      #performance detriment
      #Topt no asymetry
      if(length(inds.s)>0) ts[year.k,spec.k,inds.s,3]= 1- tpc.plot(tmax.k[inds.s],tol.h[spec.k,'Topt.noasym'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax'])  
      
      #shift Topt
      Topt.shift= tol.h[spec.k,'Topt']-tol.h[spec.k,'Topt.noasym']
      if(length(inds)>0) ts[year.k,spec.k,inds,4]= 1- tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt'],tol.h[spec.k,'CTmin'], (tol.h[spec.k,'CTmax']+Topt.shift) )  
      
      #shift slope
      if(length(inds.s)>0) ts[year.k,spec.k,inds.s,5]= 1- tpc.plot(tmax.k[inds.s],tol.h[spec.k,'Topt.noasym'],tol.h[spec.k,'CTmin'], (tol.h[spec.k,'CTmax']-Topt.shift) )  
      
      #Topt average asymetry
      if(length(inds)>0) ts[year.k,spec.k,inds,6]= 1- tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt.aveasym'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax'])  
      
      #actual Topt
      if(length(inds)>0) ts[year.k,spec.k,inds,7]= 1- tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax'])  
      
      #number of days with >= 50% loss of performance
      #symetric
      if(length(inds.s)>0) {
      perf= tpc.plot(tmax.k[inds.s],tol.h[spec.k,'Topt.noasym'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax']) 
      t50= rep(NA, length(inds.s))
      t50[perf<0.5]<-1
      ts[year.k,spec.k,inds.s,8]= t50
     
      #shift Topt
      perf= tpc.plot(tmax.k[inds.s],tol.h[spec.k,'Topt'],tol.h[spec.k,'CTmin'], (tol.h[spec.k,'CTmax']+Topt.shift) )
      t50= rep(NA, length(inds.s))
      t50[perf<0.5]<-1
      ts[year.k,spec.k,inds.s,9]= t50}
      
      if(length(inds)>0) {
      #shift slope
      perf= tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt.noasym'],tol.h[spec.k,'CTmin'], (tol.h[spec.k,'CTmax']-Topt.shift) )  
      t50= rep(NA, length(inds))
      t50[perf<0.5]<-1
      ts[year.k,spec.k,inds,10]= t50
      
      #Topt ave sym
      perf= tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt.aveasym'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax']) 
      t50= rep(NA, length(inds))
      t50[perf<0.5]<-1
      ts[year.k,spec.k,inds,11]= t50
      
      #Actual Topt
      perf= tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax']) 
      t50= rep(NA, length(inds))
      t50[perf<0.5]<-1
      ts[year.k,spec.k,inds,12]= t50}
      
      #-----
      #performance
      if(length(inds)>0) ts[year.k,spec.k,inds,13]= tpc.plot(tmax.k[inds],tol.h[spec.k,'Topt'],tol.h[spec.k,'CTmin'], tol.h[spec.k,'CTmax'])  
      
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
tsm<-  array(NA, dim= c(length(years),nrow(tol.h),8) )

#reorder ts
ts2 <- aperm(ts, c(2,1,3,4))

#loop years
for(year.k in 1:length(years)){
  
  #min TSM
  tsm[year.k,,1]=t(apply(ts2[,year.k,,1], MARGIN=1, FUN='min', na.rm=TRUE))
  
  #performace detriment
  #no asym
  tsm[year.k,,3]= apply(ts2[,year.k,,3], MARGIN=1, FUN=sum, na.rm=T)
  #shift Topt
  tsm[year.k,,4]= apply(ts2[,year.k,,4], MARGIN=1, FUN=sum, na.rm=T)
  #shift slope
  tsm[year.k,,5]= apply(ts2[,year.k,,5], MARGIN=1, FUN=sum, na.rm=T)
  #ave asym
  tsm[year.k,,6]= apply(ts2[,year.k,,6], MARGIN=1, FUN=sum, na.rm=T)
  #actual asym
  tsm[year.k,,7]= apply(ts2[,year.k,,7], MARGIN=1, FUN=sum, na.rm=T)
  
  #perf
  tsm[year.k,,8]= apply(ts2[,year.k,,13], MARGIN=1, FUN=sum, na.rm=T)
  
  #Replace -Inf with NA
  tsm[year.k,which(is.infinite(tsm[,year.k,2])),2]=NA
  tsm[year.k,which(is.infinite(tsm[,year.k,2])),6]=NA
  
} #end year loop

#AGGREGATE ACROSS YEARS
tsm.yrs= apply(tsm, MARGIN=c(2,3), FUN=mean, na.rm=T)

#----
#number of days with 50% loss of performance
d1= apply(ts2[,,,8:12], MARGIN=c(1,2,4), FUN=sum, na.rm=TRUE)
#average across years
days_p50= apply(d1, MARGIN=c(1,3), FUN=mean, na.rm=TRUE)

#----------------------
#PLOT
#comparison plots
tol2= cbind(tol.h, tsm.yrs)
colnames(tol2)[c(15,17:21 )]=c('minTSM','Perf.noAsym','Perf.dTopt','Perf.dSlope','Perf.aveAsym','Perf')
#drop fish
tol2=tol2[-which(tol2$taxa=="fish"),] 

#drop unneeded columns
tol2s= tol2[,c("taxa","asym","minTSM","Perf.noAsym",'Perf.dTopt','Perf.dSlope',"Perf.aveAsym","Perf")]
#change names
#names(tol2s)[5:6]=c("without asymetry","fitted asymetry")

#-----
#compare to TSM

#max performance by taxa
max.perf= aggregate(tol2s$Perf, list(tol2s$taxa), FUN="max")
colnames(max.perf)=c("taxa","Perf")

#Scale
match1= match(tol2$taxa, max.perf$taxa)

tol2$PerfS= tol2$Perf/max.perf$Perf[match1]
tol2$PerfS[tol2$PerfS>0.5]=0.5

fig4a= ggplot(tol2, aes(x=PerfS,y=minTSM, color=asym)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+scale_color_viridis(name="asymmetry")+ theme(legend.position = "bottom")+
  xlab("annual proportional performance detriment")+ylab("annual minimum of daily TSM (°C)")+
  ylim(-10,15)

#-----
#Compare performance estimates
#to long format
tol.l <- melt(tol2s, id=c("taxa","asym","minTSM","Perf"))

#scale to max
match1= match(tol.l$taxa, max.perf$taxa)
tol.l$Perf= tol.l$Perf/max.perf$Perf[match1]
tol.l$value= tol.l$value/max.perf$Perf[match1]

#set max to 2.5
tol.l$Perf[tol.l$Perf>2.5]=2.5
tol.l$value[tol.l$value>2.5]=2.5

#make labels
tol.l$metric.lab<-NA
tol.l$metric.lab[tol.l$variable=="Perf.noAsym"]<- "no asymmetry"
tol.l$metric.lab[tol.l$variable=="Perf.dTopt"]<- "omit slope"
tol.l$metric.lab[tol.l$variable=="Perf.dSlope"]<- "omit Topt shift"
tol.l$metric.lab[tol.l$variable=="Perf.aveAsym"]<- "taxa asymmetry"
#drop symmetric
tol.l= subset(tol.l, tol.l$metric.lab %in% c("omit Topt shift","omit slope","taxa asymmetry"))
tol.l$metric.lab= factor(tol.l$metric.lab, levels=c("omit Topt shift","omit slope","taxa asymmetry")) #drop symmetric

fig4b= ggplot(tol.l, aes(x=(Perf),y=(value), color=metric.lab)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+ theme(legend.position = "bottom", legend.title = element_blank())+
  ylab("estimated annual performance detriment")+xlab("annual proportional performance detriment")+
  scale_color_viridis(discrete=TRUE)+geom_abline(slope=1, intercept=0)

#-----
#proportion days with 50% performance loss

#comparison plots
tol2d= cbind(tol.h, days_p50)
colnames(tol2d)[(ncol(tol.h)+1):(ncol(tol.h)+5)]=c('Perf.noAsym','Perf.dTopt','Perf.dSlope','Perf.aveAsym','Perf')
#drop fish
tol2d=tol2d[-which(tol2d$taxa=="fish"),] 

#drop unneeded columns
tol2sd= tol2d[,c("taxa","asym","Perf.noAsym",'Perf.dTopt','Perf.dSlope',"Perf.aveAsym","Perf")]

#to long format
tol.l <- melt(tol2sd, id=c("taxa","asym", "Perf"))

#scale to max
match1= match(tol.l$taxa, max.perf$taxa)
tol.l$Perf= tol.l$Perf/max.perf$Perf[match1]
tol.l$value= tol.l$value/max.perf$Perf[match1]

fig4c= ggplot(tol.l, aes(x=Perf,y=value, color=variable)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+theme(legend.position = "bottom", legend.title = element_blank())+
  ylab("estimated days with 50% performance loss")+xlab("days with 50% performance loss")+
  scale_color_viridis(discrete=TRUE)+geom_abline(slope=1, intercept=0)
#----
#Plot
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/figures/")
pdf("Figs4_TSM.pdf", height = 8, width = 8)

fig4a +fig4b +plot_annotation(tag_levels = 'a') +plot_layout(nrow=2) 

dev.off()

#----
#LATITUDINAL PLOT
#latitudinal figure for plankton
tol.p= tol2[which(tol2$taxa=="plankton"), c(1:12,15,17:21) ]

#to long format
tol.pl<- tol.p %>%
  gather("metric", "value", c("Perf","Perf.noAsym","Perf.aveAsym","Perf.dTopt","Perf.dSlope") ) #"days_p50"

#Scale
match1= match(tol.pl$taxa, max.perf$taxa)
tol.pl$value= tol.pl$value/max.perf$Perf[match1]
#set max to 0.5
tol.pl$value[tol.pl$value>0.5]=0.5

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
  xlab("absolute latitude (°)")+ylab("annual proportional performance detriment")

pdf("Figs5_TSMlat.pdf", height = 8, width = 8)
fig5a +fig5b +plot_annotation(tag_levels = 'a') +plot_layout(nrow=2) 
dev.off()

