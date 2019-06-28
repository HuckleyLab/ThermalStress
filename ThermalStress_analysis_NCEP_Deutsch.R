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

#convert longituide
convert.lon= function(r0) ifelse(r0 > 180, -360 + r0, r0)

#-------------------
#Load combined NCEP data
#FROM https://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NCEP+Reanalysis&Variable=Air+Temperature&group=0&submit=Search
#Summary: https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html
#surface temp

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/NCEP/")

ncep.nc= nc_open('air.sig995.2010.nc')
print(ncep.nc)

#extract metrics
ncep.temp=ncvar_get(ncep.nc,"air")
dim(ncep.temp) #dimensions lon, lat, time
#convert from K to C
ncep.temp= ncep.temp -273.15

#extract lon, lat, time
ncep.lons= ncvar_get(ncep.nc,"lon") #get info about long
dim(ncep.lons)
ncep.lats= ncvar_get(ncep.nc,"lat") #get info about latitude
dim(ncep.lats)
ncep.times= ncvar_get(ncep.nc,"time") #julian date, calendar day, since 0000-01-01
dim(ncep.times)
#change to dates
ncep.dates= as.POSIXct(ncep.times*3600,origin='1800-01-01 00:00') 
years= as.numeric(format(ncep.dates, "%Y"))
doy= as.numeric(format(ncep.dates, "%j"))
hours= as.numeric(format(ncep.dates, "%H"))

#daily min max
ncep.cell= cbind(years, doy, hours, ncep.temp[70,35,])
colnames(ncep.cell)[4]="temp"
ncep.cell= as.data.frame(ncep.cell)
#make day factor
ncep.cell$doy= as.factor(ncep.cell$doy)

#daily min and max
ncep.cmm= ncep.cell %>%
  group_by(doy) %>%
  summarise(min = min(temp), max= max(temp))

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

#=====================================
#LOAD BIOLOGICAL DATA

#Load Deutsch Data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/DeutschData/")
tol.dd= read.csv('Deutschetel.2008Insect.TPCdata.csv')

#partition marine and terrestrial
tol.dd$habitat= 'terrestrial'

#--------------------
#THERMAL STRESS ESTIMATES

ts= array(NA, dim= c(nrow(tol.dd), 365,4) )

#calculate degree days above Topt
for(spec.k in 1:nrow(tol.dd)){
  
  #find closest grid cell
  lon.ind= which.min(abs(ncep.lons.neg - tol.dd[spec.k, "Long" ]))
  lat.ind= which.min(abs(ncep.lats - tol.dd[spec.k, "Lat" ]))
  
  #extract data
  #daily min max
  ncep.cell= cbind(years, doy, hours, ncep.temp[lon.ind,lat.ind,])
  colnames(ncep.cell)[4]="temp"
  ncep.cell= as.data.frame(ncep.cell)
  #make day factor
  ncep.cell$doy= as.factor(ncep.cell$doy)
  
  #daily min and max
  ncep.cmm= ncep.cell %>%
    group_by(doy) %>%
    summarise(min = min(temp), max= max(temp))
  ncep.cmm= as.matrix(ncep.cmm)
  tmax.k= as.numeric( ncep.cmm[,3] )
  tmin.k= as.numeric( ncep.cmm[,2] )
  
  #daily safety margins
  ts[spec.k,,1]= tol.dd[spec.k,'CTmax']-tmax.k
  ts[spec.k,,2]= tol.dd[spec.k,'Topt']-tmax.k 
  
  #metabolic scaled thermal stress
  inds= which(tmax.k > tol.dd[spec.k,'Topt'])
  
  if(length(inds)>0) ts[spec.k,inds,3]= tmax.k[inds]- tol.dd[spec.k,'Topt']  
  
  #thermodynamic scale
  Topt.tt= thermo.temp(tol.dd[spec.k,'Topt'])
  tmax.k.tt= thermo.temp(tmax.k)
  
  inds= which(tmax.k.tt > Topt.tt)
  
  if(length(inds>0)) ts[spec.k,inds,4]= tmax.k.tt[inds]- Topt.tt  
  
} # end loop species

#------------
#TSM
tsm<-  array(NA, dim= c(nrow(tol.dd),11,4) )

#TSM daily: min, 10%, 50%, 90%
#CTmax-Tmax
tsm[,1:4,1]= t(apply(ts[,,1], MARGIN=1, FUN='quantile', probs=c(0,0.1,0.5,0.9)))
#Topt-Tmax
tsm[,1:4,2]= t(apply(ts[,,2], MARGIN=1, FUN='quantile', probs=c(0,0.1,0.5,0.9), na.rm=TRUE))

#lowest 7 day average
tsm[,5,1]= apply(rollmean(ts[,,1], 7, fill=NA), MARGIN=1, FUN='min', na.rm=TRUE)
tsm[,5,2]= apply(rollmean(ts[,,2], 7, fill=NA), MARGIN=1, FUN='min', na.rm=TRUE)
#lowest 14 day average
tsm[,6,1]= apply(rollmean(ts[,,1], 14, fill=NA), MARGIN=1, FUN='min', na.rm=TRUE)
tsm[,6,2]= apply(rollmean(ts[,,2], 14, fill=NA), MARGIN=1, FUN='min', na.rm=TRUE)

#count of days <5
count.sub5= function(x) sum(x<5)
tsm[,7,1]= apply(rollmean(ts[,,1], 14, fill=NA), MARGIN=1, FUN='count.sub5')
tsm[,7,2]= apply(rollmean(ts[,,2], 14, fill=NA), MARGIN=1, FUN='count.sub5')

#-----
#METABOLIC
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

tsm[,8,3]= apply(ts[,,3], MARGIN=1, FUN='maxc.overTopt')
tsm[,8,4]= apply(ts[,,4], MARGIN=1, FUN='maxc.overTopt')
#Replace -Inf with NA
tsm[which(is.infinite(tsm[,8,3])),8,3]=NA
tsm[which(is.infinite(tsm[,8,4])),8,4]=NA

#-----
#count
count=function(x) length(na.omit(x))

#aggregate
#sum, count, mean
tsm[,9,3]= apply(ts[,,3], MARGIN=1, FUN=sum, na.rm=T)
tsm[,9,4]= apply(ts[,,4], MARGIN=1, FUN=sum, na.rm=T)

tsm[,10,3]= apply(ts[,,3], MARGIN=1, FUN=count)
tsm[,10,4]= apply(ts[,,4], MARGIN=1, FUN=count)

tsm[,11,3]= apply(ts[,,3], MARGIN=1, FUN=mean, na.rm=T)
tsm[,11,4]= apply(ts[,,4], MARGIN=1, FUN=mean, na.rm=T)

#+++++++
# #AGGREGATE ACROSS YEARS
# #Calculate annual metrics
# ts.l= as.data.frame(t(rbind(years, ts[,,3]))) #linear
# ts.t= as.data.frame(t(rbind(years, ts[,,4]))) #thermodynamic
# 
# #aggregate for each year
# ts.l.sum= aggregate(ts.l, by=list(ts.l$years), FUN=sum, na.rm=T)
# ts.t.sum= aggregate(ts.t, by=list(ts.t$years), FUN=sum, na.rm=T)
# ts.l.count= aggregate(ts.l, by=list(ts.l$years), FUN=count)
# ts.t.count= aggregate(ts.t, by=list(ts.t$years), FUN=count)
# ts.l.mean= aggregate(ts.l, by=list(ts.l$years), FUN=mean, na.rm=T)
# ts.t.mean= aggregate(ts.t, by=list(ts.t$years), FUN=mean, na.rm=T)
# 
# #average across years
# tsm[,9,3]= colMeans(ts.l.sum[,3:ncol(ts.l.sum)], na.rm=T)
# tsm[,9,4]= colMeans(ts.t.sum[,3:ncol(ts.l.sum)], na.rm=T)
# tsm[,10,3]= colMeans(ts.l.count[,3:ncol(ts.l.sum)], na.rm=T)
# tsm[,10,4]= colMeans(ts.t.count[,3:ncol(ts.l.sum)], na.rm=T)
# tsm[,11,3]= colMeans(ts.l.mean[,3:ncol(ts.l.sum)], na.rm=T)
# tsm[,11,4]= colMeans(ts.t.mean[,3:ncol(ts.l.sum)], na.rm=T)
#+++++++

#add data
#dim 1: TSM CTmax
tol.ts= cbind(tol.dd, tsm[,1:7,1])
colnames(tol.ts)[12:18]=c('minTSM','TSM10p','TSM50p', 'TSM90p', 'low7d', 'low14d', 'count5' )
#dim 2: TSM Topt
tol.ts= cbind(tol.dd, tsm[,1:7,2])
colnames(tol.ts)[12:18]=c('minTSM','TSM10p','TSM50p', 'TSM90p', 'low7d', 'low14d', 'count5' )
#dim 3: Tmax - Topt
tol.ts= cbind(tol.dd, tsm[,8:11,3])
colnames(tol.ts)[12:15]=c('dTopt','sumI','countI','meanI')
#dim 4: thermodynamic Tmax - Topt
tol.ts= cbind(tol.dd, tsm[,8:11,4])
colnames(tol.ts)[12:15]=c('dTopt','sumI','countI','meanI')

#------------------
#PLOT
#TSM
plot.minTSM= ggplot(tol.ts, aes(x=abs(Lat),y=minTSM) ) +geom_point()+geom_smooth(method='loess',se=TRUE)
plot.TSM10p= ggplot(tol.ts, aes(x=abs(Lat),y=TSM10p) ) +geom_point()+geom_smooth(method='loess',se=TRUE)
plot.TSM50p= ggplot(tol.ts, aes(x=abs(Lat),y=TSM50p) ) +geom_point()+geom_smooth(method='loess',se=TRUE)
plot.TSM90p= ggplot(tol.ts, aes(x=abs(Lat),y=TSM90p) ) +geom_point()+geom_smooth(method='loess',se=TRUE)
#lowest 7 and 14 day average 
plot.low7d= ggplot(tol.ts, aes(x=abs(Lat),y=low7d) ) +geom_point()+geom_smooth(method='loess',se=TRUE)
plot.low14d= ggplot(tol.ts, aes(x=abs(Lat),y=low14d) ) +geom_point() +geom_smooth(method='loess',se=TRUE)
#count of days <5
plot.count5= ggplot(tol.ts, aes(x=abs(Lat),y=count5) ) +geom_point() +geom_smooth(method='loess',se=TRUE)

#Plot out
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/out/")
pdf("Figs_TSM_deutsch.pdf", height = 12, width = 12)
plot_grid(plot.minTSM, plot.TSM10p, plot.TSM50p, plot.TSM90p, labels = c("minTSM", "TSM10p", "TSM50p", "TSM90p"), ncol = 2)
dev.off()

pdf("Figs_TSMdur_deutsch.pdf", height = 12, width = 12)
plot_grid(plot.low7d, plot.low14d, plot.count5, labels = c("low7d", "low7d", "count5"), ncol = 2)
dev.off()

#---------
#Metabolic integration
#duration of days above Topt
plot.dTopt= ggplot(tol.ts, aes(x=abs(Lat),y=log(dTopt)) ) +geom_point()+geom_smooth(method='loess',se=TRUE)
#sum, count, and mean of intergration
plot.sumI= ggplot(tol.ts, aes(x=abs(Lat),y=log(sumI)) ) +geom_point()+geom_smooth(method='loess',se=TRUE)
plot.countI= ggplot(tol.ts, aes(x=abs(Lat),y=log(countI)) ) +geom_point() +geom_smooth(method='loess',se=TRUE)
plot.meanI= ggplot(tol.ts, aes(x=abs(Lat),y=log(meanI)) ) +geom_point() +geom_smooth(method='loess',se=TRUE)

#Plot out
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/out/")
pdf("Figs_MetIntegration_deutsch.pdf", height = 12, width = 12)
plot_grid(plot.dTopt, plot.sumI, plot.countI, plot.meanI, labels = c("dTopt", "sumI", "countI", "meanI"), ncol = 2)
dev.off()


