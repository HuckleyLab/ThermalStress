#LOAD CLIMDEX DATA
library(ncdf4)
library(fields)

setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/ThermalStress/data/climdex/")

#WDSI
wsdi.nc= nc_open('RawData_GHCNDEX_WSDI_1951-2019_ANN_from-90to90_from-180to180.nc')
print(wsdi.nc)
#TX90
tx90.nc<-nc_open('RawData_GHCNDEX_TX90p_1951-2019_ANN_from-90to90_from-180to180.nc')
print(tx90.nc)

#extract metrics
wsdi=ncvar_get(wsdi.nc,"WSDI")
dim(wsdi) #dimensions lon, lat, time
tx90=ncvar_get(tx90.nc,"TX90p")

#extract lon, lat, time
lons= ncvar_get(wsdi.nc,"lon") #get info about long
dim(lons)
lats= ncvar_get(wsdi.nc,"lat") #get info about latitude
dim(lats)
times= ncvar_get(wsdi.nc,"time") #julian date, calendar day, since 1800-01-01
times

#image.plot format longitude, latitude, matrix of data
image.plot(lons,lats,wsdi[,,1])
image.plot(lons,lats,tx90[,,1])
