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

#convert longituide
convert.lon= function(r0) ifelse(r0 > 180, -360 + r0, r0)

#----------------------
#LOAD DATA

#LOAD TERRESTRIAL DATA
#https://www.metoffice.gov.uk/hadobs/hadghcnd/download.html

setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/ThermalStress/data/HadGHCND_TXTN_acts_1950-2014_15102015.nc/")

lst.nc= nc_open('HadGHCND_TXTN_acts_2001-2010_15102015.nc')
print(lst.nc)

#extract metrics
tmax=ncvar_get(lst.nc,"tmax")
dim(tmax) #dimensions lon, lat, time
tmin=ncvar_get(lst.nc,"tmin")
dim(tmin)

#extract lon, lat, time
lons= ncvar_get(lst.nc,"longitude") #get info about long
dim(lons)
lats= ncvar_get(lst.nc,"latitude") #get info about latitude
dim(lats)
times= ncvar_get(lst.nc,"time") #julian date, calendar day, since 0000-01-01
dim(times)
#change to dates
dates= as.Date(times, origin="0000-01-01")
years= as.numeric(format(dates, "%Y"))
doy= as.numeric(format(dates, "%j"))

#image.plot format longitude, latitude, matrix of data
image.plot(lons,lats[73:1],tmin[,73:1,1])
image.plot(lons,lats[73:1],tmax[,73:1,1])

#--------------------
#LOAD OCEAN DATA 
#https://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBListFiles.pl?did=132&tid=73559&vid=2423

setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/ThermalStress/data/SST/")

sst.nc= nc_open('sst.day.mean.2000.nc')
print(sst.nc)

#extract metrics
sst=ncvar_get(sst.nc,"sst")
dim(sst) #dimensions lon, lat, time

#extract lon, lat, time
s.lons= ncvar_get(sst.nc,"lon") #get info about long
dim(s.lons)
s.lats= ncvar_get(sst.nc,"lat") #get info about latitude
dim(s.lats)
s.times= ncvar_get(sst.nc,"time") #julian date, calendar day, since 0000-01-01
dim(s.times)
#change to dates
#dates= as.Date(s.times, origin="0000-01-01")

#image.plot format longitude, latitude, matrix of data
image.plot(s.lons,s.lats,sst[,,1])

#---------------------
#LOAD BIOLOGICAL DATA

setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")

tol.p= read.csv('Pinsky_dataset_1_hotwater.csv')
tol.gt= read.csv('GlobalTherm_upload_10_11_17.csv')

#check glob therm data
#species with CTmin and CTmax 
tol.gt.sub= tol.gt[ which(!is.na(tol.gt$Tmax) & !is.na(tol.gt$tmin) & tol.gt$max_metric=='ctmax' & tol.gt$min_metric=='ctmin'),]

#partition marine and terrestrial
tol.gt.sub$habitat= 'terrestrial'
tol.gt.sub$habitat[which(tol.gt.sub$Order %in% c('Perciformes','Decapoda','Cypriniformes','Cyprinodontiformes','Kurtiformes','Laminariales','Mugiliformes','Osmeriformes','Characiformes','Myliobatiformes','Salmoniformes') )]= 'marine'
#separate into terrestrial and marine, mostly terrestiral
tol.gt.ter= tol.gt.sub[tol.sub$habitat=='terrestrial',]
tol.gt.marine= tol.gt.sub[tol.sub$habitat=='marine',]

#extract values to find marine and terrestrial
#extract.pts <- cbind(tol.gt.sub$long_max-180,tol.gt.sub$lat_max)
#ext <- extract(ocean,extract.pts,method="bilinear")

#plot
setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/ThermalStress/data/SST/")
ocean <- brick('sst.day.mean.2000.nc')  
ocean= rotate(ocean[[1]])
plot(ocean)
points((tol.gt.ter$long_max),tol.gt.ter$lat_max)
points((tol.gt.marine$long_max),tol.gt.marine$lat_max, col='blue')

#------
#Match glob therm to hotwater data
tol.p$gen_spec= paste(tol.p$Genus, tol.p$Species, sep="_")
tol.gt$gen_spec= paste(tol.gt$Genus, tol.gt$Species, sep="_")
#match
match1= match(tol.p$gen_spec, tol.gt$gen_spec)

tol.p$tmin =tol.gt$tmin[match1]
tol.p$tmin_metric =tol.gt$min_metric[match1]

#species with CTmin and CTmax 
tol.sub= tol.p[ which(!is.na(tol.p$tmax) & !is.na(tol.p$tmin) & tol.p$tmax_metric=='crit' & tol.p$tmin_metric=='ctmin'),]
#Assume Topt is 70% of TTB
tol.sub$tmin = as.numeric( as.character(tol.sub$tmin))
tol.sub$Topt= tol.sub$tmin + (tol.sub$tmax -tol.sub$tmin)*0.7
#separate into terrestrial and marine, mostly terrestiral
tol.ter= tol.sub[tol.sub$habitat=='terrestrial',]
tol.marine= tol.sub[tol.sub$habitat=='marine',]
#check locations
plot(ocean)
points(tol.ter$lon,tol.ter$lat)
points(tol.marine$lon,tol.marine$lat, col="blue")

#------
#huey data
tol.h= read.csv('Hueyetal2009.csv')
#check locations
plot(ocean)
points(tol.h$Long,tol.h$Lat)

#===================================================
#Calculate thermal stress metrics

#Convert to thermodynamic scale
#E=0.757 #eV
#k= 1.38*10^-23 #J K^-1
#k=8.617* 10^-5 #eV K^-1
#function using temp in C
thermo.temp= function(t, E=0.757, k=8.617* 10^-5) exp(-E/(k*(t+273.15))) 

#----------
#TERRESTRIAL
ts= array(NA, dim= c(nrow(tol.h), length(dates),2) )

#calculate degree days above Topt
for(spec.k in 1:nrow(tol.h)){

#find closest grid cell
lons.neg= convert.lon(lons)
lon.ind= which.min(abs(lons.neg - tol.h[spec.k, "Long" ]))
lat.ind= which.min(abs(lats - tol.h[spec.k, "Lat" ]))

#extract data
tmin.k= tmin[lon.ind,lat.ind,]
tmax.k= tmax[lon.ind,lat.ind,]

#metabolic scaled thermal stress
Topt= tol.h[spec.k,'newTopt']
inds= which(tmax.k > Topt)

ts[spec.k,inds,1]= tmax.k[inds]- Topt  

#thermodynamic scale
Topt.tt= thermo.temp(tol.h[spec.k,'newTopt'])
tmax.k.tt= thermo.temp(tmax.k)

inds= which(tmax.k.tt > Topt.tt)

ts[spec.k,inds,2]= tmax.k.tt[inds]- Topt.tt  

} # end loop species

#------------
#Calaculate annual metrics

ts.l= as.data.frame(t(rbind(years, ts[,,1])))
ts.t= as.data.frame(t(rbind(years, ts[,,2])))
#count
count=function(x) length(na.omit(x))

ts.l.sum= aggregate(ts.l, by=list(ts.l$years), FUN=sum, na.rm=T)
ts.t.sum= aggregate(ts.t, by=list(ts.t$years), FUN=sum, na.rm=T)
ts.l.count= aggregate(ts.l, by=list(ts.l$years), FUN=count)
ts.t.count= aggregate(ts.t, by=list(ts.t$years), FUN=count)
ts.l.mean= aggregate(ts.l, by=list(ts.l$years), FUN=mean, na.rm=T)
ts.t.mean= aggregate(ts.t, by=list(ts.t$years), FUN=mean, na.rm=T)

#yearly
#ts.l.m= ts.l.agg[,3:ncol(ts.l.agg)]
#ts.t.m= ts.t.agg[,3:ncol(ts.t.agg)]
#row.names(ts.l.m)= ts.l.agg$Group.1
#row.names(ts.t.m)= ts.t.agg$Group.1
ts.l.sum.v= colMeans(ts.l.sum[,3:ncol(ts.l.sum)], na.rm=T)
ts.t.sum.v= colMeans(ts.t.sum[,3:ncol(ts.l.sum)], na.rm=T)
ts.l.count.v= colMeans(ts.l.count[,3:ncol(ts.l.sum)], na.rm=T)
ts.t.count.v= colMeans(ts.t.count[,3:ncol(ts.l.sum)], na.rm=T)
ts.l.mean.v= colMeans(ts.l.mean[,3:ncol(ts.l.sum)], na.rm=T)
ts.t.mean.v= colMeans(ts.t.mean[,3:ncol(ts.l.sum)], na.rm=T)

#add data
tol.h.ts= cbind(tol.h, ts.l.sum.v, ts.l.count.v, ts.l.mean.v, ts.t.sum.v, ts.t.count.v, ts.t.mean.v)

plot(tol.h.ts$Lat,tol.h.ts$ts.l.sum.v+1, log='y')
plot(tol.h.ts$Lat,tol.h.ts$ts.l.count.v+1, log='y')
plot(tol.h.ts$Lat,tol.h.ts$ts.l.mean.v)
plot(tol.h.ts$Lat,tol.h.ts$SafeZone)

#============================
#MARINE

tsm= array(NA, dim= c(nrow(tol.marine), length(s.times),2) )

#calculate degree days above Topt
for(spec.k in 1:nrow(tol.marine)){
  
  #find closest grid cell
  lon.ind= which.min(abs(s.lons - tol.h[spec.k, "Long" ]))
  lat.ind= which.min(abs(s.lats - tol.h[spec.k, "Lat" ]))
  
  #extract data
  tmin.k= tmin[lon.ind,lat.ind,]
  tmax.k= tmax[lon.ind,lat.ind,]
  
  #metabolic scaled thermal stress
  Topt= tol.h[spec.k,'newTopt']
  inds= which(tmax.k > Topt)
  
  ts[spec.k,inds,1]= tmax.k[inds]- Topt  
  
  #thermodynamic scale
  Topt.tt= thermo.temp(tol.h[spec.k,'newTopt'])
  tmax.k.tt= thermo.temp(tmax.k)
  
  inds= which(tmax.k.tt > Topt.tt)
  
  ts[spec.k,inds,2]= tmax.k.tt[inds]- Topt.tt  
  
} # end loop species



#%%%%%%%%%%%%%%
sst=ncvar_get(sst.nc,"sst")
dim(sst) #dimensions lon, lat, time

#extract lon, lat, time
s.lons= ncvar_get(sst.nc,"lon") #get info about long
dim(s.lons)
s.lats= ncvar_get(sst.nc,"lat") #get info about latitude
dim(s.lats)
s.times= ncvar_get(sst.nc,"time") #julian date, calendar day, since 0000-01-01
dim(s.times)
#change to dates
#dates= as.Date(s.times, origin="0000-01-01")

#image.plot format longitude, latitude, matrix of data
image.plot(s.lons,s.lats,sst[,,1])

