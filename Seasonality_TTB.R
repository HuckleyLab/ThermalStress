


#------------------------------
#HADCRUT4
#land and sea

library(raster)
require(rgdal)
library(ncdf4)
library(maps)

setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/ThermalStress/data/HadCRUT4/")
cru<-nc_open('absolute.nc')
print(cru)

#extract temperature data
temp=ncvar_get(cru,"tem")
#lons
lons= ncvar_get(cru,"lon") #get info about long
#lats
lats= ncvar_get(cru,"lat") #get info about latitude

#month
temp1= temp[,,2]
#calculate SD across months
tempsd= apply(temp, c(1,2), sd)

#to raster
seas <- raster(t(tempsd), xmn=min(lons), xmx=max(lons), ymn=min(lats), ymx=max(lats), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

#plot
tempr <- flip(seas, direction='y')
plot(seas)
#add country outline
map('world', fill = FALSE, col = "grey", add=TRUE)

#---------------------
#relate thermal breath to latitude

#Load GlobTherm
setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
tol.gt= read.csv('GlobalTherm_upload_10_11_17.csv')
#species with CTmin and CTmax 
tol.gt= tol.gt[ which(tol.gt$max_metric=='ctmax' & tol.gt$min_metric=='ctmin'),]

#make tmin numeric
tol.gt$tmin= as.numeric(as.character(tol.gt$tmin))
#Thermal tolerance breadth
tol.gt$TTB= tol.gt$Tmax - tol.gt$tmin

#extract SD
tol.gt$SD= extract(seas, tol.gt[,c('long_max','lat_max')] )

#partition marine and terrestrial
tol.gt$habitat= 'terrestrial'
tol.gt$habitat[which(tol.gt$Order %in% c('Perciformes','Decapoda','Cypriniformes','Cyprinodontiformes','Kurtiformes','Laminariales','Mugiliformes','Osmeriformes','Characiformes','Myliobatiformes','Salmoniformes') )]= 'marine'
##separate into terrestrial and marine, mostly terrestiral
#tol.gt.ter= tol.gt[tol.gt$habitat=='terrestrial',]
#tol.gt.marine= tol.gt[tol.gt$habitat=='marine',]

#============================

#-------------------------
#DATA from Buckley et al. MOVEMENT ANALYSIS

#Seasonality data
setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/ThermalStress/data/CRU_Movement/")

tseas<- read.csv("TempSeasonality3.csv") 

#T SEASONALITY
#Make spatial points data frame
xy.sp= cbind(tseas$lon, tseas$lat)
xy.cc= coordinates(xy.sp)
bbox(xy.sp)

#Make Spatial pixels data frame
#grd1 <- SpatialPixelsDataFrame(points=xy.sp, data = dat50[,10:13], tolerance=0.1, proj4string=CRS("+proj=longlat +proj=lcc"))
grd1 <- SpatialPixelsDataFrame(points=xy.sp, data = tseas[,4:5], tolerance=0.1)

#Plot SD*100
spplot(grd1, "SD")

#extract values
sdr= raster(grd1, values=TRUE)
plot(sdr)

#extract SD
tol.gt$SDm= extract(sdr, tol.gt[,c('long_max','lat_max')] )

#--------------------------
#Plot and fit models

#plot relationship with SD
ggplot(tol.gt, aes(SD, TTB, color=Class, shape=habitat)) +geom_point()+facet_wrap(~habitat)
#plot relationship with SD from movement analysis
ggplot(tol.gt, aes(SDm, TTB, color=Class)) +geom_point()
#with absolute latitude
ggplot(tol.gt, aes(abs(lat_max), TTB, color=Class, shape=habitat)) +geom_point()
#with absolute latitude
ggplot(tol.gt, aes(abs(lat_max), TTB, color=Class, shape=habitat)) +geom_point()
#plot terrestrial or marine
ggplot(tol.gt, aes(abs(lat_max), TTB)) +geom_point()+facet_wrap(~habitat)+geom_smooth(method='lm',se=FALSE)

#fit models
mod1= lm(TTB~SD*habitat, data=tol.gt)
mod1= lm(TTB~SDm, data=tol.gt)

mod1= lm(TTB~abs(lat_max)*habitat, data=tol.gt)
mod1= lm(TTB~abs(lat_max), data=tol.gt[tol.gt$habitat=='terrestrial',])
mod1= lm(TTB~abs(lat_max), data=tol.gt[tol.gt$habitat=='marine',])
summary(mod1)

#======================
#Fit Huey data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
tol.h= read.csv('Hueyetal2009.csv', na.strings ='-9999')

#estimate TTB
tol.h$TTB= tol.h$CTmax - tol.h$CTmin
#extract SD
tol.h$SDm= extract(sdr, tol.h[,c('Long','Lat')] )
#estimate Topt as percent
tol.h$ToptPer= (tol.h$newTopt-tol.h$CTmin)/(tol.h$CTmax-tol.h$CTmin)

#--------------------------
#Plot and fit models

#plot relationship with SD from movement analysis
ggplot(tol.h, aes(SDm, TTB, color=Family)) +geom_point()
#with absolute latitude
ggplot(tol.h, aes(AbsLat, TTB, color=Family)) +geom_point()

#plot position of Topt
ggplot(tol.h, aes(AbsLat, newTopt, color=Family)) +geom_point()
ggplot(tol.h, aes(AbsLat, ToptPer, color=Family)) +geom_point()

#Topt as percent of TTB
mod1=lm(ToptPer ~AbsLat , dat=tol.h)

#-------------------------
#Fit TPC

#Performance Curve Function from Deutsch et al. 2008
TPC= function(T,topt,sigma, ctmax){
  F=rep(NA, length(T))
  ind=which(T<=topt)
  F[ind]= exp(-((T[ind]-topt)/(2*sigma))^2) 
  ind=which(T>topt)
  F[ind]= 1- ((T[ind]-topt)/(topt-ctmax))^2
  return(F)
}

#==========================================
#SUMMARY FUNCTION

#TTB, based on GlobTherm
#terrestrial
TTB.terr= function(AbsLat) 29.15383 + 0.19897 * AbsLat
#marine
TTB.mar= function(AbsLat) 0.6945813 + 0.0020061 * AbsLat
#TTB= 26.25588 + 0.09722 * AbsLat

#Topt as percent of TTB, based on Huey data
ToptPer= function(AbsLat) 0.6945813 + 0.0020061 * AbsLat




