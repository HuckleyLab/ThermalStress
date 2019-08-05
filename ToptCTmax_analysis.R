library(ggplot2)

#load data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
tpc=read.csv("tpcs.csv")
#drop data without all metrics
tpc= subset(tpc, !is.na(tpc$CTmin) & !is.na(tpc$Topt) & !is.na(tpc$CTmax) )

#calculate asymetry
tpc$asym= (tpc$CTmax - tpc$Topt)/(tpc$Topt- tpc$CTmin )
tpc$asym2= (2*tpc$Topt-tpc$CTmax - tpc$CTmin)/(tpc$CTmax-tpc$CTmin )
#check relationship of assymetry metrics
plot(tpc$asym,tpc$asym2)

#calculate declining breadth
tpc$CTmax.Topt.breadth= tpc$CTmax - tpc$Topt

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

temps=-5:50

par(mfrow=c(2,2))
taxas= c("insects","lizards","lizards_Tp","plankton")
for(taxa in 1:4){
  
tpc.sub= tpc[which(tpc$taxa==taxas[taxa]),]
#add a column of color values based on assymetry values
tpc.sub$col <- rainbow(20)[as.numeric(cut(tpc.sub$asym,breaks = 20))]

plot(temps, tpc.plot(temps, Topt=tpc.sub[1,"Topt"] , CTmin=tpc.sub[1,"CTmin"], CTmax=tpc.sub[1,"CTmax"]), type="l", xlab="temperature",ylab="performance", main=taxas[taxa])
for(r in 2:nrow(tpc.sub)) points(temps, tpc.plot(temps, Topt=tpc.sub[r,"Topt"] , CTmin=tpc.sub[r,"CTmin"], CTmax=tpc.sub[r,"CTmax"]), type="l", col=tpc.sub[r,"col"]) 
}

#-----------
#PLOT RELATIONSHIPS
#assymetry vs Topt
ggplot(tpc) + aes(x=Topt, y = asym, color=habitat)+geom_point()+facet_wrap(~taxa) #, scales="free"

#TPC vs assymetry
ggplot(tpc) + aes(x=asym, y = CTmax.Topt.breadth, color=habitat)+geom_point()+facet_wrap(~taxa)
#+geom_smooth(method="lm", se=FALSE)

