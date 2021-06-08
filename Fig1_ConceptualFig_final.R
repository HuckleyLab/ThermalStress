library(ggplot2)
library(cowplot)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(cowplot)
library(viridis)
library(patchwork)
library(latex2exp)
library(ggnewscale)

### TSM defs
#Sinclair et al. 2016, CTmax- maximum Tb (or Topt?)

#Deutsch et al. 2008, WT= CTmax-Thab, TSM= Topt - Thab
#Although these heuristic indicators are defined by using annual mean temperature (Thab)

#Kingsolver et al 2013, 
#How close maximum summer temperatures are to the critical thermal maximum of a species – the thermal buffer
#TSM= Topt-Thab
#monthly mean surface air temperatures
#B = CTmax−max(Ta), where max(Ta) is the hottest mean monthly temperature in a given year

#Sunday et al 2014, CTmax > Te,max
#Maximum air temperatures [highest monthly mean of daily maximum air temperature (Ta,max)] 

#Pincebourde and Casas, WT= CTmax-Te
#99th percentile of air temperature distribution

#Pinsky et al 2019, TSM= CTmax-Te
#TSM: acute upper thermal limit of a species and the extreme hot hourly body temperature

#ADD: 
#mean annual
#hourly

#max of monthly temperatures
#max monthly mean of daily air temperatures

#ANALYSIS:
#CPD based on hourly temps
#Compare to: hourly, monthly, annual, highest monthly mean of daily maximum air temperature

#Find new dataset?
#Currently use: NCEP Reanalysis 1 project provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA (www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html). Data are derived from a forecast model incorporating empirical data assimilation and are available as a 2.5° latitude x 2.5° longitude global grid. 
#CPD: currently use Tdaily max
#Use ERA5?

####
#performance detriment functions
#set 0 min

#quadratic
quad=function(x, Topt, CTmax) max( 1-((x-Topt)/(Topt-CTmax))^2, 0)
#linear
lin=function(x, Topt, CTmax)  max( 1-((x-Topt)/(CTmax-Topt)), 0)

#gaussian
#2 sigma: 95%, 3 sigma 99.7%, Deutsch uses 4
gaus=function(x, Topt, CTmax, a=1) {
  b=Topt
  c= (CTmax-Topt)/3 #for 95%
  max( a*exp(-(x-b)^2/(2*c^2)), 0)
} 
  
####
#Need to run Fig4_TSManalysis_ERM5.R first to load envi data

#load data
setwd("./data/")
tpc=read.csv("tpcs.csv")

#drop data without all metrics
tpc= subset(tpc, !is.na(tpc$CTmin) & !is.na(tpc$Topt) & !is.na(tpc$CTmax) )
#change photosynthesis to plants
tpc$taxa[tpc$taxa=="photosynthesis"]="plants"

taxas= c("insects","lizards","plankton","fish","plants","ants")

#setwd for figures
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/figures/")
#------

#calculate asymetry
#tpc$asym= (tpc$CTmax - tpc$Topt)/(tpc$Topt- tpc$CTmin )
#OR INVERSE: tpc$asym= (tpc$Topt- tpc$CTmin )/(tpc$CTmax - tpc$Topt)
tpc$asym= (2*tpc$Topt-tpc$CTmax - tpc$CTmin)/(tpc$CTmax-tpc$CTmin )
#check relationship of assymetry metrics
#plot(tpc$asym,tpc$asym2, ylab="Martin Huey assymetry", xlab="Deutsch et al assymetry")

#calculate breadth
tpc$breadth= tpc$CTmax -tpc$CTmin
#calculate declining breadth
tpc$CTmax.Topt.breadth= tpc$CTmax - tpc$Topt

#Convert to thermodynamic scale
#E=0.757 #eV
#k= 1.38*10^-23 #J K^-1
#k=8.617* 10^-5 #eV K^-1
#function using temp in C
thermo.temp= function(t, E=0.757, k=8.617* 10^-5) exp(-E/(k*(t+273.15))) 

#================
#Plots
#PLOT TPCs

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

tpc.mat= function(tpc.dat){
  T=-5:60
  Topt= tpc.dat[1]
  CTmin= tpc.dat[2]
  CTmax= tpc.dat[3]
  F=T
  F[]=NA
  sigma= (Topt-CTmin)/4
  F[T<=Topt & !is.na(T)]= exp(-((T[T<=Topt & !is.na(T)]-Topt)/(2*sigma))^2) 
  F[T>Topt & !is.na(T)]= 1- ((T[T>Topt & !is.na(T)]-Topt)/(Topt-CTmax))^2
  #set negative to zero
  F[F<0]<-0
  
  return(F)
}

temps=-5:60

#================
#Fig 1 Conceptual

#make TPC matrix
solve.asym= function(asym, CTmin=0, CTmax=40){
Topt= (asym*(CTmax-CTmin)+CTmin+CTmax)/2  
return(Topt)
}

tpc0=as.data.frame(rbind( c(solve.asym(0.0), 0, 40),c(solve.asym(0.25), 0, 40), c(solve.asym(0.5), 0, 40)) )
tpc0$scen="shift in asymmetry (Topt shift)"
tpc0.2= as.data.frame(rbind( c(20, -5, 35),c(25, 0, 40), c(30, 5, 45)) )
tpc0.2$scen="shift in TPC (CTmin, Topt, CTmax shift)"
tpc0= rbind(tpc0, tpc0.2)

colnames(tpc0)[1:3]=c("Topt","CTmin","CTmax")
tpc0$asym= (2*tpc0$Topt-tpc0$CTmax - tpc0$CTmin)/(tpc0$CTmax-tpc0$CTmin )  
tpc0$curve= c(0.0, 0.25, 0.5, 35, 40, 45)

#estimate performances
out=t(apply(tpc0[,1:3], MARGIN=1, FUN=tpc.mat))

colnames(out)= temps
tpc.pred= cbind(tpc0, out)

#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 7:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))
tpc.l$tpc="quadratic"

#add linear and gaussian
tpc.l.lin=tpc.l[which(tpc.l$Topt==25 & tpc.l$CTmax==40 & tpc.l$scen=="shift in asymmetry (Topt shift)" & tpc.l$temperature>24),]
tpc.l.lin$performance= lin(tpc.l.lin$temperature, 25, 40)
tpc.l.lin$tpc="linear"

tpc.l.gaus=tpc.l[which(tpc.l$Topt==25 & tpc.l$CTmax==40 & tpc.l$scen=="shift in asymmetry (Topt shift)"& tpc.l$temperature>24),]
tpc.l.gaus$performance= gaus(tpc.l.lin$temperature, 25, 40)
tpc.l.gaus$tpc="gaussian"

tpc.l= rbind(tpc.l, tpc.l.lin, tpc.l.gaus)
tpc.l$curve2= paste(tpc.l$curve, tpc.l$tpc,sep="_")
#-------

#make data frame  
tpc.lab= as.data.frame(matrix(NA, 6, 3))
names(tpc.lab)=c("temperature","performance","lab")
tpc.lab[,1]=c(0, 43, tpc0[1,1]-3, 9, 5, 4)
tpc.lab[,2]=c(0.05, 0.05, 1, 0.9, 0.8, 0.6)
tpc.lab[,3]=c("CTmin","CTmax","Topt","","","" )
tpc.lab$scen= "shift in asymmetry (Topt shift)"

#plot
tpc.l$asym= factor(tpc.l$asym)
fig0a= ggplot(tpc.l)+aes(x=temperature, y = performance)+ facet_wrap(~scen)+
  geom_line(aes(color=asym, group=curve2))+
  theme_bw(base_size=14)+theme(legend.position="none")+scale_color_viridis(discrete=TRUE, name="asymmetry")+
  xlim(-2,45)+ylim(-0.01,1)+
  ylab("relative performance")+xlab("body temperature (°C)")+
  theme(plot.margin = unit(c(1,0.5,0,0.5), "lines"))

#add performance
fig0a= fig0a +
  geom_segment(data=tpc0, aes(x = 35, y = 1, xend = 35, yend = 0), lty="dashed", size=1)+ #draw performance detriment
  geom_segment(data=tpc0, aes(x = 35, y = 1, xend = 35, yend = tpc.plot(35,tpc0[2,1],tpc0[2,2],tpc0[2,3])), colour = "blue", size=1.5, arrow = arrow(length = unit(0.25, "cm")) ) #draw performance detriment

#add TSM
fig0a= fig0a +
  geom_segment(data=tpc0, aes(x = 35, y = 0, xend = 40, yend = 0), colour = "red", size=1, arrow = arrow(length = unit(0.15, "cm")))+  #draw TSM
  geom_segment(data=tpc0, aes(x = 40, y = 0, xend = 35, yend = 0), colour = "red", size=1, arrow = arrow(length = unit(0.15, "cm")))
  
#add text
fig0a= fig0a +
  geom_text(data=tpc0, x=37, y=0.9, label="PD", color="blue", angle=90)+
  geom_text(data=tpc0,x=33, y=0.05, label="TSM", color="red")
  #geom_segment(data=tpc0, aes(x = 35, y = 1, xend = 35, yend = 0),arrow = arrow(length = unit(0.5, "cm")) )

#add equations
fig0a= fig0a + facet_wrap(~scen) +
  geom_text(data=tpc.lab, aes(x=temperature, y=performance, label=lab))+
  geom_hline(yintercept=0)
 # geom_text(data=tpc0,x=9, y=0.9, label=TeX("annual performance detriment="), size=4)+
 # geom_text(data=tpc0,x=5, y=0.8, label=TeX("$\\sum_{d=1}^{365} 1-P(T_d)$ if $T_d>Topt$"), size=4)+
 # geom_text(data=tpc0,x=4, y=0.6, label=TeX("TSM=$min_{d \\in \\lbrack 1,365 \\rbrack} CT_{max}-T_d$"), size=4)         
#----
# metrics for range of asymmetry

asyms= seq(0, 0.5, 0.05)

#make vectors to store data
tsm=matrix(NA, nrow=length(asyms), ncol=2 )
pd=matrix(NA, nrow=length(asyms), ncol=2 )
pd.lin=matrix(NA, nrow=length(asyms), ncol=2 )
pd.gaus=matrix(NA, nrow=length(asyms), ncol=2 )

#make tpc parameters
tpc.p=array(NA, dim=c(length(asyms),3,2))
#vary assym
tpc.p[,1,1]=0
tpc.p[,2,1]=solve.asym(asyms)
tpc.p[,3,1]=40
#vary CTmin, CTmax
tpc.p[,1,2]=solve.asym(asyms)-25
tpc.p[,2,2]=solve.asym(asyms)
tpc.p[,3,2]=solve.asym(asyms)+15

#=================
#EXAMINE TIMESCALE
library("accelerometry")

#Aggregate data to different resolutions as in NCC
#5,10,30,60
#TSM and CPD across resolutions

#Timescale: hourly estimates? Plot similar to Buckley NCC and also calculate CPD for various locations.
#For CPD: keep Topt and CTmax constant

#load climate data
#Los Alamos surface temp, 5 minute interval
#ftp://ftp.ncdc.noaa.gov/pub/data/uscrn/products/subhourly01/2019/
#data from 2013-2019 available, use 2019

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/USCRN/")

# 4    LST_DATE                       YYYYMMDD
# 5    LST_TIME                       HHmm
# 9    AIR_TEMPERATURE                Celsius
# 13   SURFACE_TEMPERATURE            Celsius
# 14   ST_TYPE                        X
# 15   ST_FLAG                        X

#clim= read.table("CRNS0101-05-2019-NM_Los_Alamos_13_W.txt", na.strings = "-9999.0")
#clim= read.table("CRNS0101-05-2019-NM_Las_Cruces_20_N.txt", na.strings = "-9999.0")
#clim=clim[,c(4,5,9,13)]
#names(clim)<- c("date","time","Tair","Tsurf")

#ERA5 for Santa Fe
#find closest grid cell
lon.ind= which.min(abs(lons.z1 - (-105.94)))
lat.ind= which.min(abs(lats.z1 - 35.69))

tmax.k= skt.z1[lat.ind,lon.ind,]
tmax.k.yrs= as.vector(tmax.k[1,])-273.15

#link time
#bind times
tmat= as.data.frame(cbind(year,months, doy,day,hours,tmax.k.yrs))
#drop 2014 data
tmat= subset(tmat, tmat$year>2014)

#hourly bin
tmat$hr.bin[tmat$hours %in% 0:5]=1
tmat$hr.bin[tmat$hours %in% 6:11]=2
tmat$hr.bin[tmat$hours %in% 12:17]=3
tmat$hr.bin[tmat$hours %in% 18:23]=4
tmat$hr.bin= as.numeric(tmat$hr.bin)

#weekly and 2 week
tmat$wk.bin[tmat$day %in% 1:7]=1
tmat$wk.bin[tmat$day %in% 8:14]=2
tmat$wk.bin[tmat$day %in% 15:21]=3
tmat$wk.bin[tmat$day %in% 21:31]=4
tmat$wk2.bin=1
tmat$wk2.bin[tmat$wk.bin %in% 3:4]=2

#hourly, 6 hourly, daily, weekly, monthly, quarterly, annually
clim.hr <- tmax.k.yrs
clim.6hr <- tmat %>% group_by(year,doy,hr.bin) %>% summarise(tmean=mean(tmax.k.yrs))
clim.day <- tmat %>% group_by(year,doy) %>% summarise(tmean=mean(tmax.k.yrs))
clim.week <- tmat %>% group_by(year,months, wk.bin) %>% summarise(tmean=mean(tmax.k.yrs))
clim.month <- tmat %>% group_by(year,months) %>% summarise(tmean=mean(tmax.k.yrs))
clim.ann <- tmat %>% group_by(year) %>% summarise(tmean=mean(tmax.k.yrs))
#monthly data, monthly mean of daily max
tmat1 <- tmat %>% group_by(year,doy,months) %>% summarise(tmean=max(tmax.k.yrs)) #daily max
clim.momdm <- tmat1 %>% group_by(year,months)  %>% summarise(tmean=mean(tmean))

#--------------

CTmax=40
Topt=25

#make data array
tr= array(data=NA, dim= c(7,4,1) ) #dims are time, metric, mean and max

#estimate metrics
for(time.k in 1:7){
  
    if(time.k==1)temps= clim.hr
    if(time.k==2)temps= clim.6hr$tmean
    if(time.k==3)temps= clim.day$tmean
    if(time.k==4)temps= clim.week$tmean
    if(time.k==5)temps= clim.month$tmean
    if(time.k==6)temps= clim.ann$tmean
    if(time.k==7)temps= clim.momdm$tmean

  #TSM
   tr[time.k, 1, 1]= CTmax - max(temps, na.rm=TRUE)
   #CPDs
   inds= which(temps > Topt)
   
   if(time.k<7){
   pd1= apply(as.data.frame(temps[inds]), FUN="quad", MARGIN=1, Topt=Topt, CTmax=CTmax)
   tr[time.k, 2, 1]= sum(1- pd1)/length(temps)
   pd1= apply(as.data.frame(temps[inds]), FUN="lin", MARGIN=1, Topt=Topt, CTmax=CTmax)
   tr[time.k, 3, 1]= sum(1- pd1)/length(temps)
   pd1= apply(as.data.frame(temps[inds]), FUN="gaus", MARGIN=1, Topt=Topt, CTmax=CTmax)
   tr[time.k, 4, 1]= sum(1- pd1)/length(temps)
   }
   
  } #end loop time


#flatten array
tr.l=as.data.frame(tr[,,1])
colnames(tr.l)=c("tsm","quadratic","linear","gaussian")
tr.l$agg= "mean"
tr.l$time= c("hr","6hr", "day","week","month","annual","Td,max") 
tr.l$hours= 1:7   #c(0.083, 1, 24, 24*30) 

#to long format
tr.l<- tr.l %>%
  gather("metric", "value", 1:4)
#separate TSM
tr.l$met="CPD"
tr.l$met[tr.l$metric=="tsm"]="TSM"

tr.l$time= factor(tr.l$time, levels=c("hr","6hr", "day","week","month","annual","Td,max"))

#tr.l$metric[tr.l$metric=="tsm"]="NA"

#-------
#Fig. Plot metrics as a function of exposure time

fig0c= ggplot(tr.l)+aes(x=time, y=value, shape=metric)+geom_line(aes(group=metric))+geom_point()+
  facet_grid(met~., scales="free_y", switch="y")+ 
  xlab("temporal aggregation of temperature data")+
  theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1))+
  theme(legend.position="bottom")+
  guides(color=guide_legend(title="TPC type"), shape=guide_legend(title="TPC type"))+
 scale_shape_discrete(breaks = c("gaussian", "linear", "quadratic"), name="TPC type")+
  theme(legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"))+
  theme(plot.margin = unit(c(1,0.5,0,0.5), "lines"))
#+scale_color_viridis(name="times scale") + 

#================================
#calculation across curves
tmax.k=clim.hr

#calc tsm and pd
for(k in 1:length(asyms)){
  
  #TSM 
  tsm[k,1]= tpc.p[k,3,1]-max(clim.day$tmean)
  tsm[k,2]= tpc.p[k,3,2]-max(clim.day$tmean)
  
  #perform det
  inds= which(tmax.k > tpc.p[k,2,1])
  perf= sapply(tmax.k[inds],FUN=quad, Topt=tpc.p[k,2,1], CTmax=tpc.p[k,3,1])
  pd[k,1]= sum(1- perf)/length(tmax.k)
  perf= sapply(tmax.k[inds],FUN=lin, Topt=tpc.p[k,2,1], CTmax=tpc.p[k,3,1])
  pd.lin[k,1]= sum(1- perf)/length(tmax.k)
  perf= sapply(tmax.k[inds],FUN=gaus, Topt=tpc.p[k,2,1], CTmax=tpc.p[k,3,1])
  pd.gaus[k,1]= sum(1- perf)/length(tmax.k)
  
  inds= which(tmax.k > tpc.p[k,2,2])
  perf= sapply(tmax.k[inds],FUN=quad, Topt=tpc.p[k,2,2], CTmax=tpc.p[k,3,2])
  pd[k,2]= sum(1- perf)/length(tmax.k)
  perf= sapply(tmax.k[inds],FUN=lin, Topt=tpc.p[k,2,2], CTmax=tpc.p[k,3,2])
  pd.lin[k,2]= sum(1- perf)/length(tmax.k)
  perf= sapply(tmax.k[inds],FUN=gaus, Topt=tpc.p[k,2,2], CTmax=tpc.p[k,3,2])
  pd.gaus[k,2]= sum(1- perf)/length(tmax.k)
  
}# loop asymmetry

##normalize pd
#pd=pd/max(pd)

#est asym of shift in TPC
ctmin= tpc.p[,1,2]
topt= tpc.p[,2,2]
ctmax= tpc.p[,3,2]
asym.shift= (2*topt-ctmax - ctmin)/(ctmax-ctmin )

#plot
tpc1= as.data.frame(cbind(tpc.p[,2,1], tsm[,1], 1:length(asyms),asyms ))
tpc1$tpc= NA
tpc1$scen="shift in asymmetry (Topt shift)"
tpc1$var="TSM"
names(tpc1)[4]="asym"

tpc2= as.data.frame(cbind(tpc.p[,2,2], tsm[,2], 1:length(asyms),asym.shift ))
tpc2$tpc= NA
tpc2$scen="shift in TPC (CTmin, Topt, CTmax shift)"
tpc2$var="TSM"
names(tpc2)[4]="asym"

tpc3= as.data.frame(cbind(tpc.p[,2,1], pd[,1], 1:length(asyms),asyms ))
tpc3$tpc= "quadratic"
tpc3l= as.data.frame(cbind(tpc.p[,2,1], pd.lin[,1], 1:length(asyms),asyms ))
tpc3l$tpc= "linear"
tpc3g= as.data.frame(cbind(tpc.p[,2,1], pd.gaus[,1], 1:length(asyms),asyms ))
tpc3g$tpc= "gaussian"
tpc3= rbind(tpc3, tpc3l, tpc3g)

tpc3$scen="shift in asymmetry (Topt shift)"
tpc3$var="CPD"
names(tpc3)[4]="asym"

tpc4= as.data.frame(cbind(tpc.p[,2,2], pd[,2], 1:length(asyms),asym.shift ))
tpc4$tpc= "quadratic"
tpc4l= as.data.frame(cbind(tpc.p[,2,2], pd.lin[,2], 1:length(asyms),asym.shift ))
tpc4l$tpc= "linear"
tpc4g= as.data.frame(cbind(tpc.p[,2,2], pd.gaus[,2], 1:length(asyms),asym.shift ))
tpc4g$tpc= "gaussian"
tpc4= rbind(tpc4, tpc4l, tpc4g)

tpc4$scen="shift in TPC (CTmin, Topt, CTmax shift)"
tpc4$var="CPD"
names(tpc4)[4]="asym"

#tpc4$scen= factor(tpc4$scen, levels=c("TSM","performance detriment"), ordered=TRUE)

tpc.pl= rbind(tpc1, tpc2, tpc3, tpc4)
names(tpc.pl)[1:3]=c("Topt","metric","k")

#plot
fig0b= ggplot(tpc.pl)+aes(x=Topt, y = metric, color=asym, shape=tpc)+geom_point()+
  facet_grid(var~scen, scales="free_y", switch="y")+
  theme_bw(base_size=14)+scale_color_viridis(name="asymmetry") + 
  theme(legend.position="bottom", legend.key.width=unit(2,"cm"),legend.box = "vertical")+geom_line()+guides(shape=FALSE)+
  theme(legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"))+
  theme(plot.margin = unit(c(1,0.5,0,0.5), "lines"))
  #+scale_shape_discrete(breaks = c("gaussian", "linear", "quadratic"), name="TPC type")  

#----
combined <- fig0a +fig0b +fig0c +plot_annotation(tag_levels = 'A') +plot_layout(nrow=3)+plot_layout(heights = c(1.6, 1,1))

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/figures/")
pdf("Fig0.pdf", height = 10, width = 8)
combined
dev.off()

#=====================
#analyze metrics

#vary CTmax by error
# estimate TSM
# estimate CPD
# calcuate SE across estimates

#load climate data
#Los Alamos surface temp, 5 minute interval
#ftp://ftp.ncdc.noaa.gov/pub/data/uscrn/products/subhourly01/2019/
#data from 2013-2019 available, use 2019

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/USCRN/")

clim= read.table("CRNS0101-05-2019-NM_Los_Alamos_13_W.txt", na.strings = "-9999.0")
#clim= read.table("CRNS0101-05-2019-NM_Las_Cruces_20_N.txt", na.strings = "-9999.0")
clim=clim[,c(4,5,9,13)]
names(clim)<- c("date","time","Tair","Tsurf")

#make tpc parameters
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }
#se ~1C, sample size ~20, se= sd/sqrt(n)
#sd=1*sqrt(20)~4.5

tpc.p=array(NA, dim=c(1000,3,4.5))
#vary CTmax
tpc.p[,1,1]=0
tpc.p[,2,1]=25
tpc.p[,3,1]= as.vector(rnorm2(1000,40,4.5))

#vary Topt and CTmax
tpc.p[,1,2]=0
tpc.p[,2,2]=as.vector(rnorm2(1000,25,4.5))
tpc.p[,3,2]= as.vector(rnorm2(1000,40,4.5))

#--------------------
#calculation across curves
tmax.k=na.omit(clim$Tsurf)

#data structures
#make vectors to store data
tsm=matrix(NA, nrow=1000, ncol=2 )
pd=matrix(NA, nrow=1000, ncol=2 )
pd.lin=matrix(NA, nrow=1000, ncol=2 )
pd.gaus=matrix(NA, nrow=1000, ncol=2 )

#calc tsm and pd
for(k in 1:1000){
  
  #TSM 
  tsm[k,1]= tpc.p[k,3,1]-max(tmax.k)
  tsm[k,2]= tpc.p[k,3,2]-max(tmax.k)
  
  #perform det
  inds= which(tmax.k > tpc.p[k,2,1])
  perf= sapply(tmax.k[inds],FUN=quad, Topt=tpc.p[k,2,1], CTmax=tpc.p[k,3,1])
  pd[k,1]= sum(1- perf)/length(tmax.k)
  perf= sapply(tmax.k[inds],FUN=lin, Topt=tpc.p[k,2,1], CTmax=tpc.p[k,3,1])
  pd.lin[k,1]= sum(1- perf)/length(tmax.k)
  perf= sapply(tmax.k[inds],FUN=gaus, Topt=tpc.p[k,2,1], CTmax=tpc.p[k,3,1])
  pd.gaus[k,1]= sum(1- perf)/length(tmax.k)
  
  inds= which(tmax.k > tpc.p[k,2,2])
  perf= sapply(tmax.k[inds],FUN=quad, Topt=tpc.p[k,2,2], CTmax=tpc.p[k,3,2])
  pd[k,2]= sum(1- perf)/length(tmax.k)
  perf= sapply(tmax.k[inds],FUN=lin, Topt=tpc.p[k,2,2], CTmax=tpc.p[k,3,2])
  pd.lin[k,2]= sum(1- perf)/length(tmax.k)
  perf= sapply(tmax.k[inds],FUN=gaus, Topt=tpc.p[k,2,2], CTmax=tpc.p[k,3,2])
  pd.gaus[k,2]= sum(1- perf)/length(tmax.k)
  
}# loop asymmetry

#actual values
tsm= -2.7, 0.027

#summarize variance
#paper uses linear
sd1= c(sd(tsm[,1]), sd(pd.lin[,1]) )
m1= c(mean(tsm[,1]), mean(pd.lin[,1]) )
sd1/m1
sd1/c(-2.7, 0.027)

sd2= c(sd(tsm[,2]), sd(pd.lin[,2]) )
m2= c(mean(tsm[,2]), mean(pd.lin[,2]) )
sd2/m2
sd2/c(-2.7, 0.027)

#compare uncertainty Topt, CTmax
#fish aerobic metabolic scope: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2435.12618 (Payne et al. 2015)
#CTmax: https://cdnsciencepub.com/doi/pdf/10.1139/cjz-2012-0300 (Chen et al.)
#CTmax varied with incubation temp and whether by time or to common body mass
#use 90 days posthatch at standard temperature

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/FishCTvar/")
fits= read.csv("ToptTcritUncertainty.csv")
cts= read.csv("SockeyeCTs.csv")


