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

#huey data
tol.h= read.csv('Hueyetal2009.csv')

#===================================================
#Calculate thermal stress metrics
#terrestrial
spec.k=1

#find closest grid cell
lons.neg= lons-180
lon.ind= which.min(abs(lons.neg - tol.h[spec.k, "Long" ]))
lat.ind= which.min(abs(lats - tol.h[spec.k, "Lat" ]))

#extract data
tmin.k= tmin[lon.ind,lat.ind,]
tmax.k= tmax[lon.ind,lat.ind,]

#metabolic scaled thermal stress
Topt= tol.h[spec.k,'newTopt']
inds= which(tmax.k > Topt)
ts= rep(NA, length(dates))
ts[inds]= tmax.k[inds]- Topt  


