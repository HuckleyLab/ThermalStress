#Load packages
#library(rerddap)
library(ncdf4)
#library(heatwaveR)
#library(RNCEP)
library(lubridate) #date and time manipulation
library(tidyverse) #data manipulation and visualization
library(RColorBrewer) #color schemes
#library(sf) #to import a spatial object and to work with geom_sf in ggplot2
library(raster)
library(zoo)
library(cowplot)
library(patchwork)
library(viridis)
library(reshape2)
library(tidyr)
library(plyr)
library(dplyr)

#convert longitude
convert.lon= function(r0) ifelse(r0 > 180, -360 + r0, r0)

#----------------------
#LOAD ERM5, check whether surface
library(raster)
library(ncdf4)

#https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview
#API_User = "78176",
#API_Key = "062c3e77-bcc8-4c56-8e72-4872e7a92be6"

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/ERA5/")

#July 2015 data, load others
era.z1= nc_open('ERA5_Zone1.nc')
era.z2= nc_open('ERA5_Zone2.nc')
era.z3= nc_open('ERA5_Zone3.nc')
#load extra points
era.z1.p1= nc_open('ERA5_Zone1_pt1.nc') #-64.78 lon:-64.06
era.z1.p2= nc_open('ERA5_Zone1_pt2.nc')  #lat:76.28 lon:-74.75

print(era.z1.p2)

#extract lon, lat, time
lons.z1= ncvar_get(era.z1,"longitude") #get info about long
lats.z1= ncvar_get(era.z1,"latitude") #get info about latitude
lons.z2= ncvar_get(era.z2,"longitude") #get info about long
lats.z2= ncvar_get(era.z2,"latitude") #get info about latitude
lons.z3= ncvar_get(era.z3,"longitude") #get info about long
lats.z3= ncvar_get(era.z3,"latitude") #get info about latitude

#extract lon, lat, time
times= ncvar_get(era.z1,"time") #julian date, calendar day, since 0000-01-01
#times.z2= ncvar_get(era.z3,"time")
#times.z3= ncvar_get(era.z3,"time")
#all same

#units: hours since 1900-01-01 00:00:00.0
#change to dates

era.dates= as.POSIXct(times*3600,origin='1900-01-01 00:00') 
year= as.numeric(format(era.dates, "%Y"))
doy= as.numeric(format(era.dates, "%j"))
hours= as.numeric(format(era.dates, "%H"))
day= as.numeric(format(era.dates, "%d"))
months= month(era.dates)
#need to add a few hours?

#extract skt
#skt.z1=ncvar_get(era.z1,"skt",start=c(1,1,1), count=c(1,1,43824))
#skt.z2=ncvar_get(era.z2,"skt")
#skt.z3=ncvar_get(era.z3,"skt")
skt.z1.p1=ncvar_get(era.z1.p1,"skt")
skt.z1.p2=ncvar_get(era.z1.p2,"skt")

#use raster
skt.z1= brick('ERA5_Zone1.nc', var="skt")
skt.z2= brick('ERA5_Zone2.nc', var="skt")
skt.z3= brick('ERA5_Zone3.nc', var="skt")

#figure out data subsets
#1st subset: -158.1 to -49.25; -1.4, 52,26
# + 2 points: lat: -64.78 lon:-64.06 and lat:76.28 lon:-74.75
#2nd subset: -36.25 to 74.45; -58, 80
#3rd subset: 108 to 174.9; -77.5 to 55.65

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
tol.h=read.csv("tpcs_wSource_Mar2021.csv")  #  #("tpcs.csv")

#species with lat lon
tol.h= subset(tol.h, !is.na(tol.h$lat) & !is.na(tol.h$lon) )
#drop data without all metrics
tol.h= subset(tol.h, !is.na(tol.h$CTmin) & !is.na(tol.h$Topt) & !is.na(tol.h$CTmax) )

#add asym
tol.h$asym= (2*tol.h$Topt-tol.h$CTmax - tol.h$CTmin)/(tol.h$CTmax-tol.h$CTmin )

plot(tol.h$lon, tol.h$lat)

tol.h$taxa= as.character(tol.h$taxa)

#change names
tol.h$taxa[tol.h$taxa=="ants"]="ant performance"
tol.h$taxa[tol.h$taxa=="lizards"]="lizard performance"
tol.h$taxa[tol.h$taxa=="fish"]="fish growth"
tol.h$taxa[tol.h$taxa=="insects"]="insect fitness"
tol.h$taxa[tol.h$taxa=="plankton"]="plankton fitness"

taxas= c("ant performance", "lizard performance", "fish growth","insect fitness", "plankton fitness")

#===================================================
#Estimate thermal stress with increasing information

#Deutsch et al. TPC
#Performance Curve Function from Deutsch et al. 2008
tpc.plot= function(T,Topt,CTmin, CTmax){
  F=T
  F[]=NA
  sigma= (Topt-CTmin)/4
  F[T<=Topt & !is.na(T)]= exp(-((T[T<=Topt & !is.na(T)]-Topt)/(2*sigma))^2) 
  F[T>Topt & !is.na(T)]= 1- ((T[T>Topt & !is.na(T)]-Topt)/(Topt-CTmax))^2
  #set negetative to zero
  F[F<0]<-0
  
  return(F)
}

#fit asymetry vs Topt by taxa
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

years=2015:2019

ts= array(0, dim= c(length(years), nrow(tol.h), 18) )

#counts of Topt exceedences
ts.exceed= array(NA, dim= c(length(years), nrow(tol.h), 2) )

#calculate degree days above Topt
for(spec.k in 1:nrow(tol.h)){
  print(spec.k)
  
  #find zone
  zone=2
  if(tol.h[spec.k, "lon" ]< -40) zone=1
  if(tol.h[spec.k, "lon" ]>80) zone=3
  
  if(zone==1 & tol.h[spec.k, "lat"]< -50) zone=1.1
  if(zone==1 & tol.h[spec.k, "lat"]> 70) zone=1.2
  
  
  if(zone==1){
    #find closest grid cell
    lon.ind= which.min(abs(lons.z1 - tol.h[spec.k, "lon" ]))
    lat.ind= which.min(abs(lats.z1 - tol.h[spec.k, "lat" ]))
    
    tmax.k= skt.z1[lat.ind,lon.ind,]
    tmax.k.yrs= as.vector(tmax.k[1,])-273.15
  }
  
  if(zone==2){
    #find closest grid cell
    lon.ind= which.min(abs(lons.z2 - tol.h[spec.k, "lon" ]))
    lat.ind= which.min(abs(lats.z2 - tol.h[spec.k, "lat" ]))
    
    tmax.k= skt.z2[lat.ind,lon.ind,]
    tmax.k.yrs= as.vector(tmax.k[1,])-273.15
  }
  
  if(zone==3){
    #find closest grid cell
    lon.ind= which.min(abs(lons.z3 - tol.h[spec.k, "lon" ]))
    lat.ind= which.min(abs(lats.z3 - tol.h[spec.k, "lat" ]))
    
    tmax.k= skt.z3[lat.ind,lon.ind,]
    tmax.k.yrs= as.vector(tmax.k[1,])-273.15
  }
  
  if(zone==1.1)tmax.k.yrs= skt.z1.p1-273.15
  if(zone==1.2)tmax.k.yrs= skt.z1.p2-273.15
  
  #bind times
  tmat= as.data.frame(cbind(year,months, doy,hours,tmax.k.yrs))
  #drop 2014 data
  tmat= subset(tmat, tmat$year>2014)
  
  #max temp
  #hourly data
  tmax.hr <- tmat %>% group_by(year)  %>% summarise(tmax=max(tmax.k.yrs))
  #daily data
  tmat1 <- tmat %>% group_by(year,doy) %>% summarise(tmean=mean(tmax.k.yrs))
  tmax.day<- tmat1 %>% group_by(year)  %>% summarise(tmax=max(tmean))
  #monthly data, monthly mean of daily max
  tmat1 <- tmat %>% group_by(year,doy,months) %>% summarise(tmean=max(tmax.k.yrs)) #daily max
  tmat1 <- tmat1 %>% group_by(year,months)  %>% summarise(tmean=mean(tmean))
  tmax.month<- tmat1 %>% group_by(year)  %>% summarise(tmax=max(tmean))
  
  
  #thermal safety margins
  #hourly
  tsm= tol.h[spec.k,'CTmax']-tmax.hr[,"tmax"]
  ts[,spec.k,1]= tsm[,1]
  #TSM days
  tsm= tol.h[spec.k,'CTmax']-tmax.day[,"tmax"]
  ts[,spec.k,2]= tsm[,1] 
  #TSM monthly
  tsm= tol.h[spec.k,'CTmax']-tmax.month[,"tmax"] 
  ts[,spec.k,8]= tsm[,1]
  
  #loop years
  for(year.k in 1:length(years)){
    
    #extract yearly data
    tmax.k= tmax.k.yrs[which(year==years[year.k])]
    
    #find Topt exceedences
    inds.s= which(tmax.k > tol.h[spec.k,'Topt.noasym'])
    inds= which(tmax.k > tol.h[spec.k,'Topt'])
    
    #ts[year.k,spec.k,inds,2]= tmax.k[inds]- tol.h[spec.k,'Topt']  
    
    #find Topt exceedences
    inds.s= which(tmax.k > tol.h[spec.k,'Topt.noasym'])
    inds= which(tmax.k > tol.h[spec.k,'Topt'])
    
    #count exceed Topt
    ts.exceed[year.k,spec.k,1]= length(inds.s)
    ts.exceed[year.k,spec.k,2]= length(inds)
    
    #performance detriment
    if(length(inds.s)>0) {
      
      #quadratic
      #Topt no asymetry
      perf= sapply(tmax.k[inds.s],FUN=quad, Topt=tol.h[spec.k,'Topt.noasym'], CTmax=tol.h[spec.k,'CTmax'])
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt.noasym']]=1
      ts[year.k,spec.k,3]= sum(1-perf)/length(tmax.k)
      #shift Topt
      Topt.shift= tol.h[spec.k,'Topt']-tol.h[spec.k,'Topt.noasym']
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt']]=1
      perf= sapply(tmax.k[inds.s],FUN=quad, Topt=tol.h[spec.k,'Topt'], CTmax=(tol.h[spec.k,'CTmax']+Topt.shift))
      ts[year.k,spec.k,4]= sum(1-perf)/length(tmax.k)
      #shift slope
      perf= sapply(tmax.k[inds.s],FUN=quad, Topt=tol.h[spec.k,'Topt.noasym'], CTmax=(tol.h[spec.k,'CTmax']-Topt.shift))
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt.noasym']]=1
      ts[year.k,spec.k,5]=  sum(1-perf)/length(tmax.k)
      #Topt average assymetry
      perf= sapply(tmax.k[inds.s],FUN=quad, Topt=tol.h[spec.k,'Topt.aveasym'], CTmax=tol.h[spec.k,'CTmax'])
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt.aveasym']]=1
      ts[year.k,spec.k,6]=  sum(1-perf)/length(tmax.k)
      #actual Topt
      perf= sapply(tmax.k[inds.s],FUN=quad, Topt=tol.h[spec.k,'Topt'], CTmax=tol.h[spec.k,'CTmax'])
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt']]=1
      ts[year.k,spec.k,7]=  sum(1-perf)/length(tmax.k)
      
      #linear
      perf= sapply(tmax.k[inds.s],FUN=lin, Topt=tol.h[spec.k,'Topt.noasym'], CTmax=tol.h[spec.k,'CTmax'])
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt.noasym']]=1
      ts[year.k,spec.k,9]= sum(1-perf)/length(tmax.k)
      #shift Topt
      Topt.shift= tol.h[spec.k,'Topt']-tol.h[spec.k,'Topt.noasym']
      perf= sapply(tmax.k[inds.s],FUN=quad, Topt=tol.h[spec.k,'Topt'], CTmax=(tol.h[spec.k,'CTmax']+Topt.shift))
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt']]=1
      ts[year.k,spec.k,10]= sum(1-perf)/length(tmax.k)
      #shift slope
      perf= sapply(tmax.k[inds.s],FUN=lin, Topt=tol.h[spec.k,'Topt.noasym'], CTmax=(tol.h[spec.k,'CTmax']-Topt.shift))
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt.noasym']]=1
      ts[year.k,spec.k,11]=  sum(1-perf)/length(tmax.k)
      #Topt average assymetry
      perf= sapply(tmax.k[inds.s],FUN=lin, Topt=tol.h[spec.k,'Topt.aveasym'], CTmax=tol.h[spec.k,'CTmax'])
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt.aveasym']]=1
      ts[year.k,spec.k,12]=  sum(1-perf)/length(tmax.k)
      #actual Topt
      perf= sapply(tmax.k[inds.s],FUN=lin, Topt=tol.h[spec.k,'Topt'], CTmax=tol.h[spec.k,'CTmax'])
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt']]=1
      ts[year.k,spec.k,13]=  sum(1-perf)/length(tmax.k)
      
      #gaussian
      perf= sapply(tmax.k[inds.s],FUN=gaus, Topt=tol.h[spec.k,'Topt.noasym'], CTmax=tol.h[spec.k,'CTmax'])
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt.noasym']]=1
      ts[year.k,spec.k,14]= sum(1-perf)/length(tmax.k)
      #shift Topt
      Topt.shift= tol.h[spec.k,'Topt']-tol.h[spec.k,'Topt.noasym']
      perf= sapply(tmax.k[inds.s],FUN=gaus, Topt=tol.h[spec.k,'Topt'], CTmax=(tol.h[spec.k,'CTmax']+Topt.shift))
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt']]=1
      ts[year.k,spec.k,15]= sum(1-perf)/length(tmax.k)
      #shift slope
      perf= sapply(tmax.k[inds.s],FUN=gaus, Topt=tol.h[spec.k,'Topt.noasym'], CTmax=(tol.h[spec.k,'CTmax']-Topt.shift))
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt.noasym']]=1
      ts[year.k,spec.k,16]=  sum(1-perf)/length(tmax.k)
      #Topt average assymetry
      perf= sapply(tmax.k[inds.s],FUN=gaus, Topt=tol.h[spec.k,'Topt.aveasym'], CTmax=tol.h[spec.k,'CTmax'])
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt.aveasym']]=1
      ts[year.k,spec.k,17]=  sum(1-perf)/length(tmax.k)
      #actual Topt
      perf= sapply(tmax.k[inds.s],FUN=gaus, Topt=tol.h[spec.k,'Topt'], CTmax=tol.h[spec.k,'CTmax'])
      perf[tmax.k[inds.s]<tol.h[spec.k,'Topt']]=1
      ts[year.k,spec.k,18]=  sum(1-perf)/length(tmax.k)
      
    } #end check length
    
  } #loop years
  
  # #thermodynamic scale
  # Topt.tt= thermo.temp(tol.h[spec.k,'topt'])
  # tmax.k.tt= thermo.temp(tmax.k)
  # 
  # inds= which(tmax.k.tt > Topt.tt)
  # 
  # if(length(inds>0)) ts[year.k,spec.k,inds,4]= tmax.k.tt[inds]- Topt.tt  
  
} # end loop species

#save output
# Save an object to a file
#saveRDS(ts, file = "ts.rds")
#saveRDS(ts.exceed, file = "tsexceed.rds")

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
ts=readRDS("ts.rds")

#------------
# #AGGREGATE ACROSS YEARS
tsm.yrs= apply(ts, MARGIN=c(2,3), FUN=mean, na.rm=T)

colnames(tsm.yrs)=c('TSMhr',"TSMday",'Perf.noAsym','Perf.dTopt','Perf.dSlope','Perf.aveAsym','Perf',"TSMmonth",
                    'Perf.noAsym.lin','Perf.dTopt.lin','Perf.dSlope.lin','Perf.aveAsym.lin','Perf.lin',
                    'Perf.noAsym.gaus','Perf.dTopt.gaus','Perf.dSlope.gaus','Perf.aveAsym.gaus','Perf.gaus')

#----------------------
#PLOT
#comparison plots
tol2= cbind(tol.h, tsm.yrs)
#drop fish
tol2=tol2[-which(tol2$taxa=="fish growth"),] 

#drop unneeded columns
tol2s= tol2[,c("taxa","asym","TSMhr","Perf.noAsym",'Perf.dTopt','Perf.dSlope',"Perf.aveAsym","Perf","TSMday","TSMmonth", 
               'Perf.noAsym.lin','Perf.dTopt.lin','Perf.dSlope.lin','Perf.aveAsym.lin','Perf.lin',
               'Perf.noAsym.gaus','Perf.dTopt.gaus','Perf.dSlope.gaus','Perf.aveAsym.gaus','Perf.gaus')]
#change names
#names(tol2s)[5:6]=c("without asymetry","fitted asymetry")

#-----
#compare to TSM

fig4a= ggplot(tol2, aes(x=Perf.lin,y=TSMday, color=asym)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+scale_color_viridis(name="asymmetry")+ theme(legend.position = "bottom")+
  xlab("CPD (proportion)")+ylab("TSM (°C)")+
  ylim(-10,15)

fig4a.gaus= ggplot(tol2, aes(x=Perf.gaus,y=TSMday, color=asym)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+scale_color_viridis(name="asymmetry")+ theme(legend.position = "bottom")+
  xlab("CPD (proportion)")+ylab("TSM (°C)")+
  ylim(-10,15)

fig4a.quad= ggplot(tol2, aes(x=Perf,y=TSMday, color=asym)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+scale_color_viridis(name="asymmetry")+ theme(legend.position = "bottom")+
  xlab("CPD (proportion)")+ylab("TSM (°C)")+
  ylim(-10,15)

#-----
#compare TSM
tol.tsm= tol2s[,c("taxa","asym","TSMhr","TSMday","TSMmonth")]
tol.tsm <- melt(tol.tsm, id=c("taxa","asym","TSMday"))

fig4.tsm= ggplot(tol.tsm, aes(x=TSMday,y=value, color=variable)) +geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+ theme(legend.position = "bottom")+geom_smooth(method="lm") +geom_abline(col="gray")+
  xlab("TSM (daily)")+ylab("TSM (°C)")
#+ylim(-10,15)

#-----
#Compare performance estimates

#compare linear, quad, gaus
tol.p= tol2s[,c("taxa","asym","Perf","Perf.lin","Perf.gaus")]
tol.pl= melt(tol.p, id=c("taxa","asym","Perf.lin"))

fig4p1= ggplot(tol.pl, aes(x=(Perf.lin),y=(value), color=variable)) +geom_point()+
  facet_wrap(~taxa, nrow=1) +
  theme_bw()+ theme(legend.position = "bottom", legend.title = element_blank())+
  ylab("CPD for observed TPC")+xlab("CPD for observed TPC")+
  scale_color_viridis(discrete=TRUE)+geom_abline(slope=1, intercept=0)

#compare linear, quadratic, guassian
#to long format
tol.l <- melt(tol2s, id=c("taxa","asym","TSMhr","TSMday","TSMmonth", "Perf","Perf.lin","Perf.gaus"))

#gather Perf, Perf type, qud / lin/ gaus
tol.l$perf="perf"
tol.l$perf[grepl("noAsym", tol.l$variable)]="no asymmetry"
tol.l$perf[grepl("dTopt", tol.l$variable)]="omit slope"
tol.l$perf[grepl("dSlope", tol.l$variable)]="omit Topt shift"
tol.l$perf[grepl("aveAsym", tol.l$variable)]="taxa asymmetry"
#drop no asymmetry
tol.l= tol.l[-which(tol.l$perf=="no asymmetry"),]

#tpc type
tol.l$tpc="quadratic"
tol.l$tpc[grepl("lin", tol.l$variable)]="linear"
tol.l$tpc[grepl("gaus", tol.l$variable)]="gaussian"

#linear
tol.ls= subset(tol.l, tol.l$tpc=="linear")
  
fig4b= ggplot(tol.ls, aes(x=Perf.lin, y=value, color=perf))+geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+ theme(legend.position = "bottom", legend.title = element_blank())+
  ylab("CPD for altered TPC (proportion)")+xlab("CPD for observed TPC (proportion)")+
  scale_color_viridis(discrete=TRUE)+geom_abline(slope=1, intercept=0)

#gaussian
tol.ls= subset(tol.l, tol.l$tpc=="gaussian")

fig4b.gaus= ggplot(tol.ls, aes(x=Perf.lin, y=value, color=perf))+geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+ theme(legend.position = "bottom", legend.title = element_blank())+
  ylab("CPD")+xlab("CPD for observed TPC")+
  scale_color_viridis(discrete=TRUE)+geom_abline(slope=1, intercept=0)

#quadratic
tol.ls= subset(tol.l, tol.l$tpc=="quadratic")

fig4b.quad= ggplot(tol.ls, aes(x=Perf.lin, y=value, color=perf))+geom_point()+facet_wrap(~taxa, nrow=1) +
  theme_bw()+ theme(legend.position = "bottom", legend.title = element_blank())+
  ylab("CPD")+xlab("CPD for observed TPC")+
  scale_color_viridis(discrete=TRUE)+geom_abline(slope=1, intercept=0)

#-----
#Plot
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/figures/")
pdf("Figs4_TSM.pdf", height = 8, width = 8)
fig4a +fig4b +plot_annotation(tag_levels = 'A') +plot_layout(nrow=2) 
dev.off()

#Supplementary plots
#compare to TSM
pdf("Figs4_TSM_spp.pdf", height = 8, width = 8)
fig4a.gaus +fig4a.quad +plot_annotation(tag_levels = 'A') +plot_layout(nrow=2) 
dev.off()

pdf("Figs4Sup_TSM.pdf", height = 8, width = 8)
fig4.tsm + fig4b.gaus + fig4b.quad +plot_annotation(tag_levels = 'A') +plot_layout(nrow=3)
dev.off()

#-----
#stats
#fig 4a
mod1= lm(minTSM~poly(PerfS,2)*asym, data=tol2[tol2$taxa=="insects",])
mod1= lm(minTSM~poly(PerfS,2)*asym, data=tol2[tol2$taxa=="lizards",])
mod1= lm(minTSM~poly(PerfS,2)*asym, data=tol2[tol2$taxa=="plankton",])

#combine taxa
mod1= lm(minTSM~poly(PerfS,2)*asym*taxa, data=tol2)

library(nlme)
mod1= lme(minTSM~poly(PerfS,2)*asym, random =~1|taxa, data=tol2)

#fig 4b
mod1= lme(value~Perf-1, random =~1|taxa, data=tol.l[tol.l$metric.lab %in% c("omit Topt shift"),])
mod1= lme(value~Perf-1, random =~1|taxa, data=tol.l[tol.l$metric.lab %in% c("taxa asymmetry"),])
mod1= lme(value~Perf-1, random =~1|taxa, data=tol.l[tol.l$metric.lab %in% c("omit slope"),])

#95% CI
#Topt: 1.23 + 0.33
#taxa asym: 1.003785 +0.044
#omit slope: 0.73 + 0.038

summary(mod1)
anova(mod1)

#----
#LATITUDINAL PLOT
#latitudinal figure for plankton
tol.p= tol2[which(tol2$taxa %in% c("plankton fitness")), c(1:13,17,24:28) ] #"insects","lizards",

#to long format
tol.pl<- tol.p %>%
  gather("metric", "value", c("Perf.lin","Perf.noAsym.lin","Perf.aveAsym.lin","Perf.dTopt.lin","Perf.dSlope.lin") ) #"days_p50"

#make labels
tol.pl$metric.lab<-NA
tol.pl$metric.lab[tol.pl$metric=="TSMday"]<- "annual minimum of daily TSM"
tol.pl$metric.lab[tol.pl$metric=="Perf.lin"]<- "observed" #"log annual performance detriment"
tol.pl$metric.lab[tol.pl$metric=="Perf.noAsym.lin"]<- "without asymetry"
tol.pl$metric.lab[tol.pl$metric=="Perf.aveAsym.lin"]<- "average asymetry"
tol.pl$metric.lab= factor(tol.pl$metric.lab, levels=c("annual minimum of daily TSM","observed","without asymetry","average asymetry"))

#TSM plot
fig5a= ggplot(tol.p, aes(x=abs(lat),y=TSMday, color=asym) ) +
  # facet_wrap(~taxa, nrow=1)+
  geom_point() +geom_smooth(method='loess',se=TRUE) +
  theme_bw()+scale_color_viridis(name="asymmetry")+ theme(legend.position = "bottom",legend.key.width = unit(2, "cm"))+
  xlab("absolute latitude (°)")+ylab("TSM (°C)")+
  ylim(-10,15) +geom_vline(xintercept=23.55)+geom_vline(xintercept=66.6)

#Perf plot
tol.p2= tol.pl[which(tol.pl$metric %in% c("Perf.lin","Perf.noAsym.lin")),]

fig5b= ggplot(tol.p2, aes(x=abs(lat),y=value, color=metric.lab) ) +
  geom_point()+geom_smooth(method='loess',se=TRUE) +
  # facet_wrap(~taxa, nrow=1)+
  theme_bw()+scale_color_viridis(name="", discrete=TRUE)+ theme(legend.position = "bottom",legend.key.width = unit(2, "cm"))+
  xlab("absolute latitude (°)")+ylab("CPD (normalized)")+
  geom_vline(xintercept=23.55)+geom_vline(xintercept=66.6)

pdf("Figs5_TSMlat.pdf", height = 8, width = 8)
fig5a +fig5b +plot_annotation(tag_levels = 'A') +plot_layout(nrow=2) 
dev.off()
