library(ggplot2)
library(cowplot)
library(dplyr)
library(reshape2)
library(tidyr)
library(cowplot)

#load data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
tpc=read.csv("tpcs.csv")
#drop data without all metrics
tpc= subset(tpc, !is.na(tpc$CTmin) & !is.na(tpc$Topt) & !is.na(tpc$CTmax) )

taxas= c("insects","lizards","plankton","fish","photosynthesis","ants")

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
  T=-5:50
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

library(viridis)
temps=-5:50

#ggplot
out=t(apply(tpc[,c("Topt","CTmin","CTmax")], MARGIN=1, FUN=tpc.mat))
colnames(out)= temps
tpc.pred= cbind(tpc, out)
#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 15:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))

#plot
fig1a= ggplot(tpc.l)+aes(x=temperature, y = performance, color=asym, group=X)+facet_grid(taxa~1)+geom_line()+
  theme_bw()+theme(legend.position="bottom")+scale_color_viridis(name="asymmetry")+
  ylab("relative performance")+xlab("temperature (°C)")

#-----------------------------------
#Topt vs
#assymetry
fig1b= ggplot(tpc) + aes(x=Topt, y = asym, color=asym, group=taxa)+geom_point()+ylab("asymmetry")+facet_grid(taxa~1)+
  ylim(-0.5,0.75)+
  theme_bw()+theme(legend.position="bottom")+geom_smooth(method="lm")+scale_color_viridis(name="asymmetry")+
  xlab("thermal optima, Topt (°C)")

#CTmax
fig1c= ggplot(data=tpc, aes(x=Topt, y = CTmax, color=asym, group=taxa))+geom_point()+facet_grid(taxa~1,scales="free_y")+ geom_smooth(method="lm")+
  theme_bw()+theme(legend.position="bottom")+scale_color_viridis(name="asymmetry")+
  ylab("CTmax (°C)")+xlab("Topt (°C)")

#breadth
fig1d= ggplot(data=tpc, aes(x=Topt, y = breadth, color=asym, group=taxa))+geom_point()+facet_grid(taxa~1,scales="free_y")+ geom_smooth(method="lm")+
  theme_bw()+theme(legend.position="bottom")+scale_color_viridis(name="asymmetry")+
  ylab("breath (°C)")+xlab("Topt (°C)")#+ylim(0,50)

#----
#Null analysis
#mean CTmin CTmax
tpc.agg <-aggregate(tpc, by=list(tpc$taxa), FUN=mean, na.rm=TRUE)
tpc.agg$taxa= tpc.agg$Group.1
match1= match(tpc$taxa, tpc.agg$taxa)
tpc$CTmin.mean= tpc$CTmin[match1]
tpc$CTmax.mean= tpc$CTmax[match1]
tpc$asym.null= (2*tpc$Topt-tpc$CTmax.mean - tpc$CTmin.mean)/(tpc$CTmax.mean-tpc$CTmin.mean )
  
#add null to asymetry plot
#assymetry vs Topt
#fig1b= fig1b + geom_smooth(data=tpc,aes(x=Topt, y = asym.null, group=taxa), method="lm", se=FALSE, lty="dashed", color="blue")

#-----
pdf("Fig1.pdf", height = 10, width = 12)
plot_grid(fig1a, fig1b, fig1c, fig1d, labels = c('A', 'B','C','D'), rel_widths = c(1.3,1,1,1), nrow=1)
dev.off()

#-----
#thermodynamic plots
tpc$asym.thermo<- (2*thermo.temp(tpc$Topt)-thermo.temp(tpc$CTmax) - thermo.temp(tpc$CTmin) )/(thermo.temp(tpc$CTmax)-thermo.temp(tpc$CTmin) )
tpc$CTmax.Topt.breadth.thermo= thermo.temp(tpc$CTmax) - thermo.temp(tpc$Topt)

#assymetry vs Topt
fig2at= ggplot(tpc) + aes(x=thermo.temp(Topt), y = asym.thermo, group=taxa)+geom_point()+ylab("asymmetry")+facet_grid(taxa~1)+
  ylim(-0.5,0.75)+
  theme_bw()+theme(legend.position="bottom")+geom_smooth(method="lm", color="black")+xlab("thermodynamic Topt (°C)")

#TPC vs assymetry
fig2bt= ggplot(tpc) + aes(x=asym.thermo, y = CTmax.Topt.breadth.thermo, group=taxa)+geom_point()+xlab("asymmetry")+ylab("thermodynamic breadth (CTmax-Topt, °C)")+facet_grid(taxa~1)+
  theme_bw()+theme(legend.position="bottom")+geom_smooth(method="lm", color="black")
#+geom_smooth(method="lm", se=FALSE)

pdf("FigSX_ThermoAssym.pdf", height = 10, width = 8)
plot_grid(fig2at, fig2bt, nrow=1)
dev.off()

#=================
#Null analysis of asymmetry relationship based on TPC constraints

ps= matrix(NA, nrow=nrow(tpc), ncol=10)
pc.var= as.data.frame(matrix(NA, nrow=length(taxas), ncol=4))
pc.var[,1]= taxas

for(taxa in 1:length(taxas)){
  
  inds= which(tpc$taxa==taxas[taxa])
  tpc.sub= tpc[inds,]
  #add a column of color values based on assymetry values
  
  pc=princomp(tpc.sub[,c("CTmin","Topt","CTmax")], cor=FALSE)
  tpc.sub= cbind(tpc.sub, pc$scores)
  
  pc$loadings
  #pc1: all parameters increase together, shift in asymmetry
  #pc2: as CTmin increases, Topt and CTmax decrease, narrowing
  #pc3: gets steeper
  
  #variances
  pc.var[taxa,2:4] <- pc$sdev^2/sum(pc$sdev^2)
  
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
  
  #pc2
  p.pc2=p
  p.pc2[,1]=p.mean[1]+pc$scores[,2]*pc$loadings[2,1]
  p.pc2[,2]=p.mean[2]+pc$scores[,2]*pc$loadings[2,2]
  p.pc2[,3]=p.mean[3]+pc$scores[,2]*pc$loadings[2,3]
  
  #pc3
  p.pc3=p
  p.pc3[,1]=p.mean[1]+pc$scores[,3]*pc$loadings[3,1]
  p.pc3[,2]=p.mean[2]+pc$scores[,3]*pc$loadings[3,2]
  p.pc3[,3]=p.mean[3]+pc$scores[,3]*pc$loadings[3,3]

  #combine data
  add=cbind(p.pc1, p.pc2, p.pc3)
  ps[inds,]=as.matrix(add)
    
} #end loop taxa

#add names
colnames(ps)=c('pc1.CTmin','pc1.Topt','pc1.CTmax','pc1.asym','pc2.CTmin','pc2.Topt','pc2.CTmax','pc3.CTmin','pc3.Topt','pc3.CTmax')
#combine
tpc.pca= cbind(tpc,ps)

colnames(pc.var)=c("taxa","pc1.var","pc2.var", "pc3.var")

#------
#plots
#assymetry
tpc.pca$lab="PC1"
fig3a= ggplot(tpc.pca) + aes(x=pc1.Topt, y = pc1.asym)+geom_point(size=2)+ylab("asymmetry")+xlab("Topt (°C)")+facet_grid(taxa~lab)+
  theme_bw()+theme(legend.position="bottom")+scale_color_viridis(name="asymmetry")

pc.var$label=paste("PC1=",round(pc.var$pc1.var,2),"PC2=",round(pc.var$pc2.var,2),"PC3=",round(pc.var$pc3.var,2), sep=" ")
fig3a= fig3a +geom_text(data= pc.var, mapping=aes(x=12, y=0.9,label=label))

#pc1 plots
out=t(apply(tpc.pca[,c("pc1.Topt","pc1.CTmin","pc1.CTmax","var.pc1")], MARGIN=1, FUN=tpc.mat))
colnames(out)= temps
tpc.pred= cbind(tpc, out)
#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 15:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))
tpc.l$lab="PC1"

#plot
fig3b= ggplot(tpc.l)+aes(x=temperature, y = performance, color=asym, group=X)+facet_grid(taxa~lab)+geom_line()+
  theme_bw()+theme(legend.position="bottom")+scale_color_viridis(name="asymmetry")+ylim(0,1)+xlab("temperature (°C)")

#pc2 plots
out=t(apply(tpc.pca[,c("pc2.Topt","pc2.CTmin","pc2.CTmax")], MARGIN=1, FUN=tpc.mat))
colnames(out)= temps
tpc.pred= cbind(tpc, out)
#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 15:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))
tpc.l$lab="PC2"

#plot
fig3c= ggplot(tpc.l)+aes(x=temperature, y = performance, color=asym, group=X)+facet_grid(taxa~lab)+geom_line()+
  theme_bw()+theme(legend.position="bottom")+scale_color_viridis(name="asymmetry")+ylim(0,1)+xlab("temperature (°C)")

#pc3 plots
out=t(apply(tpc.pca[,c("pc3.Topt","pc3.CTmin","pc3.CTmax")], MARGIN=1, FUN=tpc.mat))
colnames(out)= temps
tpc.pred= cbind(tpc, out)
#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 15:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))
tpc.l$lab="PC3"

#plot
fig3d= ggplot(tpc.l)+aes(x=temperature, y = performance, color=asym, group=X)+facet_grid(taxa~lab)+geom_line()+
  theme_bw()+theme(legend.position="bottom")+scale_color_viridis(name="asymmetry")+ylim(0,1)+xlab("temperature (°C)")

#PCA analysis as in Knies et al.
pdf("Fig3_PCAs.pdf", height = 10, width = 10)
plot_grid(fig3a, fig3b, fig3c, fig3d, nrow=1)
dev.off()

#================
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

#===============================
#DROP


#Fig 1b: warm cool and warm sections of tpcs
#plot segments
tpc.cold= tpc
tpc.cold$section<- "cold" 
tpc.cold$breadth<- tpc.cold$Topt-tpc.cold$CTmin
tpc.warm= tpc
tpc.warm$section<- "warm" 
tpc.warm$breadth<- tpc.cold$CTmax-tpc.cold$Topt
tpc2= rbind(tpc.cold, tpc.warm)

fig1b= ggplot(tpc2) + aes(x=Topt, y = breadth, color=asym, group=taxa)+geom_point()+facet_grid(taxa~section,scales="free_y")+ geom_smooth(method="lm")+
  theme_bw()+theme(legend.position="bottom")+scale_color_viridis(name="asymmetry")+
  ylab("breath of side (°C)")+xlab("thermal optima, Topt (°C)")#+ylim(0,50)

#Shuffling
#NEEDS FIXING SINCE CURRENTLY SHUFFLES ALL
bootstrap_tpc<- function(CTmins, CTmaxs, Topts, runs){
  for(k in 1:runs){
    CTmins1= sample(CTmins)
    CTmaxs1= sample(CTmaxs)
    cold.breadth<- Topts-CTmins1
    warm.breadth<- CTmaxs1-Topts
    asym= (2*Topts-CTmaxs1 - CTmins1)/(CTmaxs-CTmins1 )
    if(k==1){cb= cold.breadth; wb=warm.breadth; as= asym} 
    if(k>1){cb= cbind(cb,cold.breadth); wb= cbind(wb,warm.breadth); as=cbind(as,asym)} 
  }
  cm=rowMeans(cb); wm= rowMeans(wb); am=rowMeans(as)
  cse=apply(cb, MARGIN=1, FUN=function(x, runs){sd(x)/sqrt(runs)}, runs=runs)
  wse=apply(wb, MARGIN=1, FUN=function(x, runs){sd(x)/sqrt(runs)}, runs=runs)
  ase=apply(as, MARGIN=1, FUN=function(x, runs){sd(x)/sqrt(runs)}, runs=runs)
  
  return(cbind(cm,wm,cse,wse,am,ase))
}

boot= bootstrap_tpc(CTmins=tpc$CTmin, CTmaxs=tpc$CTmax, Topts=tpc$Topt,runs=100)
colnames(boot)= c("cold","warm","cold.se","warm.se","asym","asym.se")
tpc.boot= cbind(tpc, boot[,1:4])

#to long format
tpc.r<- tpc.boot %>%
  gather("section", "breadth", 16:17)
#drop negative breaths
tpc.r= tpc.r[which(tpc.r$breadth>0),]

#add random to initial plot
fig1b= fig1b + geom_smooth(data=tpc.r,aes(x=Topt, y = breadth, group=taxa), method="lm", se=FALSE, lty="dashed")

#-----
#thermodynamic plots
tpc$asym.thermo<- (2*thermo.temp(tpc$Topt)-thermo.temp(tpc$CTmax) - thermo.temp(tpc$CTmin) )/(thermo.temp(tpc$CTmax)-thermo.temp(tpc$CTmin) )
tpc$CTmax.Topt.breadth.thermo= thermo.temp(tpc$CTmax) - thermo.temp(tpc$Topt)

#assymetry vs Topt
fig2at= ggplot(tpc) + aes(x=thermo.temp(Topt), y = asym.thermo, group=taxa)+geom_point()+ylab("asymmetry")+facet_grid(taxa~1)+
  ylim(-0.5,0.75)+
  theme_bw()+theme(legend.position="bottom")+geom_smooth(method="lm", color="black")+xlab("thermodynamic Topt (°C)")

#TPC vs assymetry
fig2bt= ggplot(tpc) + aes(x=asym.thermo, y = CTmax.Topt.breadth.thermo, group=taxa)+geom_point()+xlab("asymmetry")+ylab("thermodynamic breadth (CTmax-Topt, °C)")+facet_grid(taxa~1)+
  theme_bw()+theme(legend.position="bottom")+geom_smooth(method="lm", color="black")
#+geom_smooth(method="lm", se=FALSE)

pdf("FigSX_ThermoAssym.pdf", height = 10, width = 8)
plot_grid(fig2at, fig2bt, nrow=1)
dev.off()