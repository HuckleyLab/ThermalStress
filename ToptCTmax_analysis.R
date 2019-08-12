library(ggplot2)

#load data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
tpc=read.csv("tpcs.csv")
#drop data without all metrics
tpc= subset(tpc, !is.na(tpc$CTmin) & !is.na(tpc$Topt) & !is.na(tpc$CTmax) )

#calculate asymetry
tpc$asym= (tpc$CTmax - tpc$Topt)/(tpc$Topt- tpc$CTmin )
#OR INVERSE: tpc$asym= (tpc$Topt- tpc$CTmin )/(tpc$CTmax - tpc$Topt)
tpc$asym2= (2*tpc$Topt-tpc$CTmax - tpc$CTmin)/(tpc$CTmax-tpc$CTmin )
#check relationship of assymetry metrics
plot(tpc$asym,tpc$asym2, ylab="Martin Huey assymetry", xlab="Deutsch et al assymetry")

#calculate declining breadth
tpc$CTmax.Topt.breadth= tpc$CTmax - tpc$Topt

#Convert to thermodynamic scale
#E=0.757 #eV
#k= 1.38*10^-23 #J K^-1
#k=8.617* 10^-5 #eV K^-1
#function using temp in C
thermo.temp= function(t, E=0.757, k=8.617* 10^-5) exp(-E/(k*(t+273.15))) 

tpc$CTmax.Topt.breadth.thermo= thermo.temp(tpc$CTmax) - thermo.temp(tpc$Topt)
tpc$asym.thermo= (thermo.temp(tpc$CTmax) - thermo.temp(tpc$Topt))/(thermo.temp(tpc$Topt)- thermo.temp(tpc$CTmin) )

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

library(viridis)
temps=-5:50

par(mfrow=c(3,2))
taxas= c("lizards","lizards_Tp","insects","flies","plankton")
for(taxa in 1:5){
  
tpc.sub= tpc[which(tpc$taxa==taxas[taxa]),]
#add a column of color values based on assymetry values
tpc.sub$col <- viridis(20)[as.numeric(cut(tpc.sub$asym,breaks = quantile(tpc.sub$asym, probs = seq(0, 1, 0.05)) ))]

plot(temps, tpc.plot(temps, Topt=tpc.sub[1,"Topt"] , CTmin=tpc.sub[1,"CTmin"], CTmax=tpc.sub[1,"CTmax"]), type="l", xlab="temperature",ylab="performance", main=taxas[taxa])
for(r in 2:nrow(tpc.sub)) points(temps, tpc.plot(temps, Topt=tpc.sub[r,"Topt"] , CTmin=tpc.sub[r,"CTmin"], CTmax=tpc.sub[r,"CTmax"]), type="l", col=tpc.sub[r,"col"]) 
}

#-----------
#PLOT RELATIONSHIPS
#assymetry vs Topt
ggplot(tpc) + aes(x=Topt, y = asym, color=habitat)+geom_point()+facet_wrap(~taxa)+ylim(0,2) #, scales="free"

#TPC vs assymetry
ggplot(tpc) + aes(x=asym, y = CTmax.Topt.breadth, color=habitat)+geom_point()+facet_wrap(~taxa)+xlim(0,2)
#+geom_smooth(method="lm", se=FALSE)

ggplot(tpc) + aes(x=Topt, y = CTmax.Topt.breadth, color=habitat)+geom_point()+facet_wrap(~taxa)
ggplot(tpc) + aes(x=CTmax, y = CTmax.Topt.breadth, color=habitat)+geom_point()+facet_wrap(~taxa)

#thermodynamic plots
#assymetry vs Topt
ggplot(tpc) + aes(x=thermo.temp(Topt), y = asym.thermo, color=habitat)+geom_point()+facet_wrap(~taxa) #, scales="free"
#TPC vs assymetry
ggplot(tpc) + aes(x=asym.thermo, y = CTmax.Topt.breadth.thermo, color=habitat)+geom_point()+facet_wrap(~taxa)

#-------
#JK suggestions:
#PCA analysis as in Knies et al. 
#Null analysis of asymmetry relationship based on TPC constraints

taxas= c("lizards","lizards_Tp","insects","flies","plankton")
tpc.sub= tpc[tpc$taxa==taxas[5],]

#plot segments
plot(tpc.sub$Topt,tpc.sub$CTmax-tpc.sub$Topt)
plot(tpc.sub$Topt,tpc.sub$Topt-tpc.sub$CTmin)

#---
pc=princomp(tpc.sub[,c("CTmin","Topt","CTmax")])
tpc.sub= cbind(tpc.sub, pc$scores)

pc$loadings
#pc1: all parameters increase together, shift in asymmetry
#pc2: as CTmin increases, Topt and CTmax decrease, narrowing
#pc3: gets steeper

#plot along pca axes
p=tpc.sub[,c("CTmin","Topt","CTmax")]
p.mean= colMeans(p)

#pc1
#a_i= mean(a)+score*loading
p.pc1=p
p.pc1[,1]=p.mean[1]+pc$scores[,1]*pc$loadings[1,1]
p.pc1[,2]=p.mean[2]+pc$scores[,1]*pc$loadings[1,2]
p.pc1[,3]=p.mean[3]+pc$scores[,1]*pc$loadings[1,3]

#estimate asmmetry
p.pc1$asym= (p.pc1[,3] - p.pc1[,2])/(p.pc1[,2]- p.pc1[,1])
#p.pc1$asym= (2*p.pc1[,2]-p.pc1[,3] - p.pc1[,1])/(p.pc1[,3]-p.pc1[,1] )
plot(p.pc1[,2], p.pc1$asym)
plot(p.pc1$asym, p.pc1[,3] - p.pc1[,2])
#first pc captures asymmetry
#add a column of color values based on assymetry values
p.pc1$col <- viridis(20)[as.numeric(cut(tpc.sub$asym,breaks = quantile(p.pc1$asym, probs = seq(0, 1, 0.05)) ))]

#plot
plot(temps, tpc.plot(temps, Topt=p.pc1[1,2] , CTmin=p.pc1[1,1], CTmax=p.pc1[1,3]), type="l", xlab="temperature",ylab="performance")
for(r in 2:nrow(p.pc1)) points(temps, tpc.plot(temps, Topt=p.pc1[r,2] , CTmin=p.pc1[r,1], CTmax=p.pc1[r,3]), type="l", col=p.pc1[r,"col"]) 

#pc2
p.pc2=p
p.pc2[,1]=p.mean[1]+pc$scores[,2]*pc$loadings[2,1]
p.pc2[,2]=p.mean[2]+pc$scores[,2]*pc$loadings[2,2]
p.pc2[,3]=p.mean[3]+pc$scores[,2]*pc$loadings[2,3]

#plot
plot(temps, tpc.plot(temps, Topt=p.pc2[1,2] , CTmin=p.pc2[1,1], CTmax=p.pc2[1,3]), type="l", xlab="temperature",ylab="performance")
for(r in 2:nrow(p.pc2)) points(temps, tpc.plot(temps, Topt=p.pc2[r,2] , CTmin=p.pc2[r,1], CTmax=p.pc2[r,3]), type="l") 

#pc3
p.pc3=p
p.pc3[,1]=p.mean[1]+pc$scores[,3]*pc$loadings[3,1]
p.pc3[,2]=p.mean[2]+pc$scores[,3]*pc$loadings[3,2]
p.pc3[,3]=p.mean[3]+pc$scores[,3]*pc$loadings[3,3]

#plot
plot(temps, tpc.plot(temps, Topt=p.pc3[1,2] , CTmin=p.pc3[1,1], CTmax=p.pc3[1,3]), type="l", xlab="temperature",ylab="performance")
for(r in 2:nrow(p.pc3)) points(temps, tpc.plot(temps, Topt=p.pc3[r,2] , CTmin=p.pc3[r,1], CTmax=p.pc3[r,3]), type="l") 

#--------------------
#variances
var(tpc.sub$CTmin)
var(tpc.sub$Topt)
var(tpc.sub$CTmax)

#tpc.m= melt(tpc.sub, id.vars=c("species","genus","CTmin","CTmax","Topt","asym") , measure.vars=c("Comp.1","Comp.2","Comp.3"))

ggplot(data=tpc.sub, aes(x=Comp.1, y=Comp.2, color=Topt))+geom_point()
ggplot(data=tpc.sub, aes(x=Comp.1, y=Comp.3, color=Topt))+geom_point()
ggplot(data=tpc.sub, aes(x=Comp.1, y=Comp.2, color=asym))+geom_point()

# P matrix
#https://www.biorxiv.org/content/biorxiv/early/2015/09/11/026518.full.pdf
#install.packages("evolqg")
#library(evolqg)

tpc.sub$gen_spec=paste(tpc.sub$genus,tpc.sub$species,sep="_")
tpc.lm= lm(as.matrix(tpc.sub[,c("CTmin","CTmax","Topt")])~tpc.sub[,"gen_spec"])
cov.matrix <- CalculateMatrix(tpc.lm)
#To obtain a correlation matrix, use:
cor.matrix <- cov2cor(cov.matrix)

MeanMatrixStatistics(cov.matrix)

#https://rdrr.io/cran/QGglmm/man/QGvcov.html





