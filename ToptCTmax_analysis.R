library(ggplot2)
library(cowplot)
library(dplyr)

#load data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
tpc=read.csv("tpcs.csv")
#drop data without all metrics
tpc= subset(tpc, !is.na(tpc$CTmin) & !is.na(tpc$Topt) & !is.na(tpc$CTmax) )
#drop sea urchin data
tpc= subset(tpc, tpc$taxa!="Sea urchins" & tpc$taxa!="Isopod" & tpc$taxa!="Bonefish")

#setwd for figures
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/figures/")
#------

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

#combine fish and lizards
tpc[which(tpc$taxa %in% c("Charr","Trout","Salmon","Bonefish")),"taxa"]<-"fish"
tpc[which(tpc$taxa=="Australian lizards"),"taxa"]<-"lizards"

#taxas= c("lizards","lizards_Tp","insects","flies","plankton","Sea urchins", "Australian lizards", "Charr","Trout","Salmon") #"Bonefish"
taxas= c("flies","insects","lizards","lizards_Tp","plankton","fish")

pdf("Fig1a_TPCs.pdf", height = 12, width = 4)
par(mfrow=c(6,1), cex=1.1, lwd=1, mar=c(3,3,0,0), mgp=c(1.3, 0.5, 0), oma=c(0,2,0,0), bty="l", cex.lab=1.2)

for(taxa in 1:length(taxas) ){
tpc.sub= tpc[which(tpc$taxa==taxas[taxa]),]
#add a column of color values based on assymetry values
tpc.sub$col<-viridis(1)
try(tpc.sub$col <- viridis(20)[as.numeric(cut(tpc.sub$asym,breaks = quantile(tpc.sub$asym, probs = seq(0, 1, 0.05)) ))])

if(taxa %in% c(1,2,4,5)) plot(temps, tpc.plot(temps, Topt=tpc.sub[1,"Topt"] , CTmin=tpc.sub[1,"CTmin"], CTmax=tpc.sub[1,"CTmax"]), type="l", xlab="",ylab="", main="") #taxas[taxa]
if(taxa==3) plot(temps, tpc.plot(temps, Topt=tpc.sub[1,"Topt"] , CTmin=tpc.sub[1,"CTmin"], CTmax=tpc.sub[1,"CTmax"]), type="l", xlab="",ylab="performance", main="")
if(taxa==6) plot(temps, tpc.plot(temps, Topt=tpc.sub[1,"Topt"] , CTmin=tpc.sub[1,"CTmin"], CTmax=tpc.sub[1,"CTmax"]), type="l", xlab="temperature (C)",ylab="", main="")

for(r in 2:nrow(tpc.sub)) points(temps, tpc.plot(temps, Topt=tpc.sub[r,"Topt"] , CTmin=tpc.sub[r,"CTmin"], CTmax=tpc.sub[r,"CTmax"]), type="l", col=tpc.sub[r,"col"]) 
}
dev.off()

#-----------
#PLOT RELATIONSHIPS

tpc$taxa= factor( tpc$taxa, levels=c("flies","insects","lizards","lizards_Tp","plankton","fish"))

#Fig 1b: warm cool and warm sections of tpcs
#plot segments
tpc.cold= tpc
tpc.cold$section<- "cold" 
tpc.cold$breadth<- tpc.cold$Topt-tpc.cold$CTmin
tpc.warm= tpc
tpc.warm$section<- "warm" 
tpc.warm$breadth<- tpc.cold$CTmax-tpc.cold$Topt
tpc2= rbind(tpc.cold, tpc.warm)

fig1b= ggplot(tpc2) + aes(x=Topt, y = breadth, color=habitat, group=taxa)+geom_point()+facet_grid(taxa~section)+theme(legend.position="bottom")+ geom_smooth(method="lm")

#Shuffling 
#Need to bootstrap
tpc2<- tpc %>% group_by(taxa) %>% mutate(CTmin=sample(CTmin), CTmax=sample(CTmax), Topt=sample(Topt) )
tpc2.cold= tpc2
tpc2.cold$section<- "cold" 
tpc2.cold$breadth<- tpc2.cold$Topt-tpc2.cold$CTmin
tpc2.warm= tpc2
tpc2.warm$section<- "warm" 
tpc2.warm$breadth<- tpc2.cold$CTmax-tpc2.cold$Topt
tpc2r= rbind(tpc2.cold, tpc2.warm)

#add random to initial plot
fig1b= fig1b + geom_smooth(data=tpc2r,aes(x=Topt, y = breadth, group=taxa), method="lm", lty="dashed")

pdf("Fig1b_TPCbreadth.pdf", height = 10, width = 6)
fig1b
dev.off()

#-----------------------------------
#assymetry vs Topt
fig2a= ggplot(tpc) + aes(x=Topt, y = asym2, color=habitat, group=taxa)+geom_point()+ylim(0,1)+ylab("assymetry")+facet_grid(taxa~1)+theme(legend.position="bottom")+geom_smooth(method="lm")

#TPC vs assymetry
fig2b= ggplot(tpc) + aes(x=asym2, y = CTmax.Topt.breadth, color=habitat, group=taxa)+geom_point()+xlim(0,1)+xlab("assymetry")+ylab("warm side breadth (CTmax-Topt)")+facet_grid(taxa~1)+theme(legend.position="bottom")+geom_smooth(method="lm")
#+geom_smooth(method="lm", se=FALSE)

#Shuffling 
#Need to bootstrap
tpc2<- tpc %>% group_by(taxa) %>% mutate(CTmin=sample(CTmin), CTmax=sample(CTmax), Topt=sample(Topt) )
#calculate asymetry
tpc2$asym= (tpc2$CTmax - tpc2$Topt)/(tpc2$Topt- tpc2$CTmin )
tpc2$asym2= (2*tpc2$Topt-tpc2$CTmax - tpc2$CTmin)/(tpc2$CTmax-tpc2$CTmin )
#estimate breadth
tpc2$CTmax.Topt.breadth= tpc2$CTmax - tpc$Topt

#add random to initial plot
#assymetry vs Topt
#fig2a= ggplot(tpc2) + aes(x=Topt, y = asym2, color=habitat, group=taxa)+geom_point()+ylim(0,1)+ylab("assymetry")+facet_grid(taxa~1)+theme(legend.position="bottom")+geom_smooth(method="lm")
fig2a= fig2a + geom_smooth(data=tpc2,aes(x=Topt, y = asym2, color=habitat, group=taxa), method="lm", lty="dashed")

#TPC vs assymetry
#fig2b= ggplot(tpc2) + aes(x=asym2, y = CTmax.Topt.breadth, color=habitat, group=taxa)+geom_point()+xlim(0,1)+xlab("assymetry")+ylab("warm side breadth (CTmax-Topt)")+facet_grid(taxa~1)+theme(legend.position="bottom")+geom_smooth(method="lm")
fig2b= fig2b + geom_smooth(data=tpc2,aes(x=asym2, y = CTmax.Topt.breadth, color=habitat, group=taxa), method="lm", lty="dashed")

pdf("Fig2_Assym.pdf", height = 10, width = 8)
plot_grid(fig2a, fig2b, nrow=1)
dev.off()

#-----
#thermodynamic plots
tpc$asym.thermo<- (2*thermo.temp(tpc$Topt)-thermo.temp(tpc$CTmax) - thermo.temp(tpc$CTmin) )/(thermo.temp(tpc$CTmax)-thermo.temp(tpc$CTmin) )
tpc$CTmax.Topt.breadth.thermo= thermo.temp(tpc$CTmax) - thermo.temp(tpc$Topt)

#assymetry vs Topt
fig2at= ggplot(tpc) + aes(x=thermo.temp(Topt), y = asym.thermo, color=habitat, group=taxa)+geom_point()+ylab("assymetry")+xlab("thermodynamic Topt")+facet_grid(taxa~1)+theme(legend.position="bottom")+geom_smooth(method="lm")
#TPC vs assymetry
fig2bt= ggplot(tpc) + aes(x=asym.thermo, y = CTmax.Topt.breadth.thermo, color=habitat, group=taxa)+geom_point()+xlab("assymetry")+ylab("thermodynamic breadth (CTmax-Topt)")+facet_grid(taxa~1)+theme(legend.position="bottom")+geom_smooth(method="lm")

pdf("FigSX_ThermoAssym.pdf", height = 10, width = 8)
plot_grid(fig2at, fig2bt, nrow=1)
dev.off()

#=================
#JK suggestions:
#Null analysis of asymmetry relationship based on TPC constraints

#PCA analysis as in Knies et al. 
pdf("Fig3_PCAs.pdf", height = 10, width = 10)

par(mfrow=c(6,4), cex=1.1, lwd=1, mar=c(3,3,1,0), mgp=c(1.3, 0.5, 0), oma=c(0,0,0,0), bty="l", cex.lab=1.2)
#taxas= c("flies","insects","lizards","lizards_Tp","plankton")

for(taxa in 1:6){

  tpc.sub= tpc[which(tpc$taxa==taxas[taxa]),]
  #add a column of color values based on assymetry values
   
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
#p.pc1$asym= (p.pc1[,3] - p.pc1[,2])/(p.pc1[,2]- p.pc1[,1])
p.pc1$asym= (2*p.pc1[,2]-p.pc1[,3] - p.pc1[,1])/(p.pc1[,3]-p.pc1[,1] )
plot(p.pc1[,2], p.pc1$asym, ylim=c(-0.2,1), xlab="PC1", ylab="assymetry")
#plot(p.pc1$asym, p.pc1[,3] - p.pc1[,2])
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
plot(temps, tpc.plot(temps, Topt=p.pc2[1,2] , CTmin=p.pc2[1,1], CTmax=p.pc2[1,3]), type="l", xlab="temperature",ylab="performance", main=taxas[taxa])
for(r in 2:nrow(p.pc2)) points(temps, tpc.plot(temps, Topt=p.pc2[r,2] , CTmin=p.pc2[r,1], CTmax=p.pc2[r,3]), type="l") 

#pc3
p.pc3=p
p.pc3[,1]=p.mean[1]+pc$scores[,3]*pc$loadings[3,1]
p.pc3[,2]=p.mean[2]+pc$scores[,3]*pc$loadings[3,2]
p.pc3[,3]=p.mean[3]+pc$scores[,3]*pc$loadings[3,3]

#plot
plot(temps, tpc.plot(temps, Topt=p.pc3[1,2] , CTmin=p.pc3[1,1], CTmax=p.pc3[1,3]), type="l", xlab="temperature",ylab="performance")
for(r in 2:nrow(p.pc3)) points(temps, tpc.plot(temps, Topt=p.pc3[r,2] , CTmin=p.pc3[r,1], CTmax=p.pc3[r,3]), type="l") 

} #end loop taxa

dev.off()

#--------------------
#variances
tpc.sub= tpc[which(tpc$taxa==taxas[5]),]

var(tpc.sub$CTmin)
var(tpc.sub$Topt)
var(tpc.sub$CTmax)

# P matrix
#https://www.biorxiv.org/content/biorxiv/early/2015/09/11/026518.full.pdf
#install.packages("evolqg")
library(evolqg)

tpc.sub$gen_spec=paste(tpc.sub$genus,tpc.sub$species,sep="_")
tpc.lm= lm(as.matrix(tpc.sub[,c("CTmin","CTmax","Topt")])~tpc.sub[,"gen_spec"])
cov.matrix <- CalculateMatrix(tpc.lm)
#To obtain a correlation matrix, use:
cor.matrix <- cov2cor(cov.matrix)

MeanMatrixStatistics(cov.matrix)

#https://rdrr.io/cran/QGglmm/man/QGvcov.html

#VISUALIZATION scores*loadings
#PC1 scores plot along pc space
#PC2 
#how does assymetry vary along axis

#cov: opt and max more strongly correlation
#where is Topt
#selection gradient CTmax
#selection on tpc, 
#plankton function valued, curves vs parameters selection, selection on curves not parameters, e.g., Ann Rev Gomulkiwitz, Kingsolver

#--------------------------

