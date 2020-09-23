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

#load data
#setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
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
  #set negetative to zero
  F[F<0]<-0
  
  return(F)
}

temps=-5:60

#===============
#Run ToptCTmax_analysis_weather.R to load envi data

#================
#Fig 0 Conceptual

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

out=t(apply(tpc0[,1:3], MARGIN=1, FUN=tpc.mat))
colnames(tpc0)[1:3]=c("Topt","CTmin","CTmax")
tpc0$asym= (2*tpc0$Topt-tpc0$CTmax - tpc0$CTmin)/(tpc0$CTmax-tpc0$CTmin )  
tpc0$curve= c(0.0, 0.25, 0.5, 35, 40, 45)

colnames(out)= temps
tpc.pred= cbind(tpc0, out)

#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 7:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))

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
  geom_line(aes(color=asym, group=curve))+
  theme_bw(base_size=14)+theme(legend.position="none")+scale_color_viridis(discrete=TRUE, name="asymmetry")+
  xlim(-2,45)+ylim(-0.01,1)+
  ylab("relative performance")+xlab("body temperature (°C)")

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

#santa fe, nm
lat= 35.6870
lon= -105.9378
year.k=1

#find closest grid cell
lon.ind= which.min(abs(ncep.lons.neg - lon))
lat.ind= which.min(abs(ncep.lats - lat))

#extract data
#daily min max
ncep.cell= cbind(year, doy, hours, ncep.temp.yrs[year.k,lon.ind,lat.ind,])
colnames(ncep.cell)[4]="temp"
ncep.cell= as.data.frame(ncep.cell)
#make day factor
ncep.cell$doy= as.factor(ncep.cell$doy)

# #daily min and max
tmax.k= tapply(ncep.cell$temp, ncep.cell$doy, max)
tmin.k= tapply(ncep.cell$temp, ncep.cell$doy, min)

hist(tmax.k)

#----
#calc tsm and pd
for(k in 1:length(asyms)){

#TSM 
tsm[k,1]= tpc.p[k,3,1]-max(tmax.k)
tsm[k,2]= tpc.p[k,3,2]-max(tmax.k)

#perform det
inds= which(tmax.k > tpc.p[k,2,1])
pd[k,1]= sum(1- tpc.plot(tmax.k[inds], tpc.p[k,1,1], tpc.p[k,2,1], tpc.p[k,3,1]))  
inds= which(tmax.k > tpc.p[k,2,2])
pd[k,2]= sum(1- tpc.plot(tmax.k[inds], tpc.p[k,1,2], tpc.p[k,2,2], tpc.p[k,3,2]))  

}# loop asymmetry

#normalize pd
pd=pd/max(pd)

#est asym of shift in TPC
ctmin= tpc.p[,1,2]
topt= tpc.p[,2,2]
ctmax= tpc.p[,3,2]
asym.shift= (2*topt-ctmax - ctmin)/(ctmax-ctmin )

#plot
tpc1= as.data.frame(cbind(tpc.p[,2,1], tsm[,1], 1:length(asyms),asyms ))
tpc1$scen="shift in asymmetry (Topt shift)"
tpc1$var="TSM"
names(tpc1)[4]="asym"

tpc2= as.data.frame(cbind(tpc.p[,2,2], tsm[,2], 1:length(asyms),asym.shift ))
tpc2$scen="shift in TPC (CTmin, Topt, CTmax shift)"
tpc2$var="TSM"
names(tpc2)[4]="asym"

tpc3= as.data.frame(cbind(tpc.p[,2,1], pd[,1], 1:length(asyms),asyms ))
tpc3$scen="shift in asymmetry (Topt shift)"
tpc3$var="CPD"
names(tpc3)[4]="asym"

tpc4= as.data.frame(cbind(tpc.p[,2,2], pd[,2], 1:length(asyms),asym.shift ))
tpc4$scen="shift in TPC (CTmin, Topt, CTmax shift)"
tpc4$var="CPD"
names(tpc4)[4]="asym"

#tpc4$scen= factor(tpc4$scen, levels=c("TSM","performance detriment"), ordered=TRUE)

tpc.pl= rbind(tpc1, tpc2, tpc3, tpc4)
names(tpc.pl)[1:3]=c("Topt","metric","k")

#plot
fig0b= ggplot(tpc.pl)+aes(x=Topt, y = metric, color=asym)+geom_point()+
  facet_grid(var~scen, scales="free_y", switch="y")+
  theme_bw(base_size=14)+scale_color_viridis(name="asymmetry") + 
  theme(legend.position="bottom", legend.key.width=unit(2,"cm"))+geom_line()

#----
combined <- fig0a +fig0b + plot_annotation(tag_levels = 'A') +plot_layout(nrow=2)+plot_layout(heights = c(1.6, 1))

pdf("Fig0.pdf", height = 8, width = 10)
combined
dev.off()

