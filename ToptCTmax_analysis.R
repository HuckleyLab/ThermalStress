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

#================
#Fig 0 Conceptual

#make TPC matrix
solve.asym= function(asym, CTmin=0, CTmax=40){
Topt= (asym*(CTmax-CTmin)+CTmin+CTmax)/2  
return(Topt)
}

tpc0=as.data.frame(rbind( c(solve.asym(0.0), 0, 40), c(solve.asym(0.5), 0, 40)) )

out=t(apply(tpc0, MARGIN=1, FUN=tpc.mat))
tpc0$asym= c(0.0, 0.5)
colnames(tpc0)[1:3]=c("Topt","CTmin","CTmax")

colnames(out)= temps
tpc.pred= cbind(tpc0, out)

#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 5:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))

#plot
tpc.l$asym= factor(tpc.l$asym)
fig0a= ggplot(tpc.l)+aes(x=temperature, y = performance, group=asym)+geom_line(aes(color=asym))+
  theme_bw(base_size=14)+theme(legend.position="bottom")+scale_color_viridis(discrete=TRUE, name="asymmetry")+
  xlim(-2,45)+
  ylab("relative performance")+xlab("body temperature (°C)")+
  geom_segment(data=tpc0, aes(x = 35, y = 1, xend = 35, yend = 0), colour = "blue", size=1, linetype = "dashed")+ #draw performance detriment
  geom_segment(data=tpc0, aes(x = 35, y = 0, xend = 40, yend = 0), colour = "red", size=2)+  #draw TSM
  geom_segment(data=tpc0, aes(x = 35, y = 1, xend = 35, yend = tpc.plot(35,tpc0[1,1],tpc0[1,2],tpc0[1,3]) ), colour = "blue", size=2)+ #draw performance detriment
  geom_segment(data=tpc0, aes(x = 35.3, y = 1, xend = 35.3, yend = tpc.plot(35,tpc0[2,1],tpc0[2,2],tpc0[2,3]) ), colour = "blue", size=2)+ #draw performance detriment
  geom_text(data=tpc0, x=0, y=0.05, label="CTmin")+
  geom_text(data=tpc0, x=43, y=0.05, label="CTmax")+
  geom_text(data=tpc0, x=tpc0[1,1]-3, y=1, label="Topt")+
  geom_text(data=tpc0, x=35, y=0.85, label="performance detriment")+
  geom_text(data=tpc0,x=33, y=0, label="TSM")+
  geom_text(data=tpc0,x=33, y=0, label="TSM")

#add equations
fig0a= fig0a + 
  geom_text(data=tpc0,x=9, y=0.9, label="annual performance detriment=", size=4)+
  geom_text(data=tpc0,x=5, y=0.8, label=TeX("$\\sum_{d=1}^{365} 1-P(T_d)$ if $T_d>Topt$"), size=4)+
  geom_text(data=tpc0,x=4, y=0.6, label=TeX("TSM=$min_{d \\in \\lbrack 1,365 \\rbrack} T_d$"), size=4)

#----
#plot other curves

tpc0=as.data.frame(rbind( c(solve.asym(0.5), 0, 40), c(solve.asym(0.5), 0, 40+10), c(solve.asym(0.0), 0, 40-10) ))

out=t(apply(tpc0, MARGIN=1, FUN=tpc.mat))
tpc0$lab= c("observed", "omit slope", "omit Topt shift")
colnames(tpc0)[1:3]=c("Topt","CTmin","CTmax")

colnames(out)= temps
tpc.pred= cbind(tpc0, out)

#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 5:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))

#plot
tpc.l$lab= factor(tpc.l$lab, levels=c("omit Topt shift","omit slope","observed"))
fig0b= ggplot(tpc.l)+aes(x=temperature, y = performance, group=lab)+geom_line(aes(color=lab))+
  theme_bw(base_size=14)+theme(legend.position="bottom")+scale_color_viridis(discrete=TRUE, name="")+
  xlim(-2,45)+
  ylab("relative performance")+xlab("body temperature (°C)")+
  geom_segment(data=tpc0, aes(x = 35, y = 1, xend = 35, yend = 0), colour = "blue", linetype = "dashed") #draw performance detriment

#----
combined <- fig0a +fig0b + plot_annotation(tag_levels = 'a') +plot_layout(nrow=1)

pdf("Fig0.pdf", height = 8, width = 12)
combined
dev.off()

#================
#FIGURE 1
out=t(apply(tpc[,c("Topt","CTmin","CTmax")], MARGIN=1, FUN=tpc.mat))
colnames(out)= temps
tpc.pred= cbind(tpc, out)
#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 15:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))

#plot
fig1a= ggplot(tpc.l)+aes(x=temperature, y = performance, color=asym, group=X)+facet_grid(taxa~.)+geom_line()+
  theme_bw()+theme(legend.position="bottom")+scale_color_viridis(name="asymmetry")+
  xlim(-5,55)+ylim(0,1)+
  ylab("relative performance")+xlab("temperature (°C)")

#-----------------------------------
# #TRY BUT ABANDON PLOTTING TOGETHER
# #Topt vs plots
tpc.plot=tpc[,c("taxa","asym","Topt","CTmax","CTmin","breadth","CTmax.Topt.breadth") ]
tpc.plot$asymmetry= tpc.plot$asym
names(tpc.plot)[7]="warm side breadth"
 
#to long format
tpc.plot <- melt(tpc.plot, id=c("taxa","asym","Topt"))
tpc.plot$variable= factor(tpc.plot$variable, levels=c("CTmax","CTmin","breadth","warm side breadth","asymmetry") )
 
# #plot
# ggplot(tpc.plot) + aes(x=Topt, y = value, color=asym, group=taxa)+geom_point()+
#   facet_grid(taxa~variable, scales="free")+
#   theme_bw()+theme(legend.position="bottom")+geom_smooth(method="lm")+scale_color_viridis(name="asymmetry")+
#   xlab("thermal optima, Topt (°C)")
#-----------------------------------

#assymetry
fig1b= ggplot(tpc) + aes(x=Topt, y = asym, color=asym, group=taxa)+geom_point()+ylab("asymmetry")+
  facet_grid(taxa~.)+
  ylim(-0.5,0.75)+
  theme_bw()+
  theme(legend.position="bottom",strip.background = element_blank(), strip.text = element_blank())+
  geom_smooth(method="lm", color="black")+scale_color_viridis(name="asymmetry")+
  xlab("Topt (°C)")+
  guides(color = FALSE)

#combine CTmin and CTmax
tpc.plot2= subset(tpc.plot, tpc.plot$variable %in% c("CTmax","CTmin") )
fig1c= ggplot(data=tpc.plot2, aes(x=Topt, y = value, shape=variable))+geom_point(aes(color=asym))+
  facet_grid(taxa~.)+ 
  scale_color_viridis(name="asymmetry")+
  new_scale_color() +
  geom_smooth(method="lm", aes(color=variable))+
  theme_bw()+theme(legend.position="bottom",strip.background = element_blank(), strip.text = element_blank())+
  ylab("CTmin and CTmax (°C)")+xlab("Topt (°C)")+
  guides(color = FALSE, shape=FALSE)

#combine breadth and warm side breadth
tpc.plot2= subset(tpc.plot, tpc.plot$variable %in% c("breadth", "warm side breadth") )
fig1d= ggplot(data=tpc.plot2, aes(x=Topt, y = value, shape=variable))+geom_point(aes(color=asym))+
  facet_grid(taxa~.)+ 
  scale_color_viridis(name="asymmetry")+
  new_scale_color() +
  geom_smooth(method="lm", aes(color=variable))+
  theme_bw()+theme(legend.position="bottom",strip.background = element_blank(), strip.text = element_blank())+
  ylab("breadth (CTmax-CTmin and CTmax-Topt, °C)")+xlab("Topt (°C)")+
  guides(color = FALSE, shape=FALSE)

#----
#Null analysis
#mean CTmin CTmax
tpc.agg <-aggregate(tpc, by=list(tpc$taxa), FUN=mean, na.rm=TRUE)
tpc.agg$taxa= tpc.agg$Group.1
match1= match(tpc$taxa, tpc.agg$taxa)
tpc$CTmin.mean= tpc.agg$CTmin[match1]
tpc$CTmax.mean= tpc.agg$CTmax[match1]
tpc$asym.null= (2*tpc$Topt-tpc$CTmax.mean - tpc$CTmin.mean)/(tpc$CTmax.mean-tpc$CTmin.mean )
  
#add null to asymetry plot
#assymetry vs Topt
#fig1b= fig1b + geom_smooth(data=tpc,aes(x=Topt, y = asym.null, group=taxa), method="lm", se=FALSE, lty="dashed", color="blue")

#-----
combined <- fig1a +fig1b +fig1c +fig1d + plot_annotation(tag_levels = 'a') +plot_layout(nrow=1, guides = "collect") & theme(legend.position = "bottom") 
#+fig1e

pdf("Fig1.pdf", height = 8, width = 12)
combined
dev.off()

#-----
#thermodynamic plots
tpc$asym.thermo<- (2*thermo.temp(tpc$Topt)-thermo.temp(tpc$CTmax) - thermo.temp(tpc$CTmin) )/(thermo.temp(tpc$CTmax)-thermo.temp(tpc$CTmin) )
tpc$CTmax.Topt.breadth.thermo= thermo.temp(tpc$CTmax) - thermo.temp(tpc$Topt)

#assymetry vs Topt
fig2at= ggplot(tpc) + aes(x=thermo.temp(Topt), y = asym.thermo, group=taxa)+geom_point()+ylab("asymmetry")+facet_grid(taxa~1)+
  ylim(-0.5,0.75)+
  theme_bw()+theme(legend.position="bottom")+geom_smooth(method="lm", color="black")+xlab("thermodynamic Topt (°C)")

pdf("FigSX_ThermoAssym.pdf", height = 10, width = 4)
fig2at
dev.off()

#------------
#Fig 1 stats by taxa

#MODELS:
#Asymm~Topt
#CTmax~Topt
#breadth~Topt
#warm side breadth ~Topt

for(var.k in 1:5){

  if(var.k==1) models <- dlply(tpc, "taxa", function(df) lm(asym ~ Topt, data = df))
  if(var.k==2) models <- dlply(tpc, "taxa", function(df) lm(CTmin ~ Topt, data = df))
  if(var.k==3) models <- dlply(tpc, "taxa", function(df) lm(CTmax ~ Topt, data = df))
  if(var.k==4) models <- dlply(tpc, "taxa", function(df) lm(breadth ~ Topt, data = df))
  if(var.k==5) models <- dlply(tpc, "taxa", function(df) lm(CTmax.Topt.breadth ~ Topt, data = df))

# Print the summary of each model
#l_ply(models, summary, .print = TRUE)

# Apply coef to each model and return a data frame
# see https://stats.stackexchange.com/questions/5135/interpretation-of-rs-lm-output
coefs= ldply(models, coef)
ses= ldply(models, function(mod) sqrt(diag(vcov(mod))) )
ts= ldply(models, function(mod) coef(mod) / sqrt(diag(vcov(mod))) )
ps= ldply(models, function(mod) 2 * pt(abs(coef(mod) / sqrt(diag(vcov(mod)))), df = df.residual(mod), lower.tail = FALSE)  )

#combine slope data
asymm.mod= cbind(coefs$Topt, ses$Topt, ts$Topt, ps$Topt)

if(var.k==1) slope.mod= asymm.mod
if(var.k>1) slope.mod= cbind(slope.mod, asymm.mod)
} #end loop var.k

colnames(slope.mod)= rep(c("slope","se","t value","p"),5)
#round
slope.mod= round(slope.mod, 2)
#add taxa
rownames(slope.mod)= coefs$taxa

#write out
write.csv(slope.mod, "Fig1stats.csv", row.names = TRUE)

#==============================================
#Null analysis of asymmetry relationship based on TPC constraints

ps= matrix(NA, nrow=nrow(tpc), ncol=11)
pc.var= as.data.frame(matrix(NA, nrow=length(taxas), ncol=4))
pc.var[,1]= taxas
#loadings
pc.load= array(NA, dim=c(3,3,length(taxas) ))

tpc$CTmin.Topt.breadth= tpc$Topt -tpc$CTmin

#pick PC variables
pc.vars=1
# 1 is Topt, CTmax, breadth
# 2 is Topt, cool breadth, warm breadth
# 3 is Topt, CTmin, CTmax
# 4 is CTmin, Topt, CTmax

#for(pc.vars in 1:1){

for(taxa in 1:length(taxas)){
  
  inds= which(tpc$taxa==taxas[taxa])
  tpc.sub= tpc[inds,]
  #add a column of color values based on assymetry values
  
  # if(pc.vars==1) pc=prcomp(tpc.sub[,c("Topt","CTmax","breadth")], center = TRUE,scale. = FALSE)
  # if(pc.vars==2) pc=prcomp(tpc.sub[,c("Topt","CTmin.Topt.breadth","CTmax.Topt.breadth")], center = TRUE,scale. = FALSE)
  # if(pc.vars==3) pc=prcomp(tpc.sub[,c("Topt", "CTmin","CTmax")], center = TRUE,scale. = FALSE)
  # if(pc.vars==4) pc=prcomp(tpc.sub[,c("CTmin","Topt","CTmax")], center = TRUE,scale. = FALSE)
  if(pc.vars==1) pc=princomp(tpc.sub[,c("Topt","CTmax","breadth")], cor=FALSE, fix_sign=TRUE)
  if(pc.vars==2) pc=princomp(tpc.sub[,c("Topt","CTmin.Topt.breadth","CTmax.Topt.breadth")], cor=FALSE, fix_sign=TRUE)
  if(pc.vars==3) pc=princomp(tpc.sub[,c("Topt", "CTmin","CTmax")], cor=FALSE, fix_sign=TRUE)
  if(pc.vars==4) pc=princomp(tpc.sub[,c("CTmin","Topt","CTmax")], cor=FALSE, fix_sign=TRUE)
  
  #tpc.sub= cbind(tpc.sub, pc$x)
  tpc.sub= cbind(tpc.sub, pc$scores)
  
  #pc.load[,,taxa]=pc$rotation
  pc.load[,,taxa]=pc$loadings
  
  #variances
  pc.var[taxa,2:4] <- pc$sdev^2/sum(pc$sdev^2)
  
  #plot along pca axes
  if(pc.vars==1) p=tpc.sub[,c("Topt","CTmax","breadth")]
  if(pc.vars==2) p=tpc.sub[,c("Topt","CTmin.Topt.breadth","CTmax.Topt.breadth")]
  if(pc.vars==3) p=tpc.sub[,c("Topt", "CTmin","CTmax")]
  if(pc.vars==4) p=tpc.sub[,c("CTmin","Topt", "CTmax")]
  
  p.mean= colMeans(p)
  
  #pc1
  #a_i= mean(a)+score*loading
  p.pc1=p
  # p.pc1[,1]=p.mean[1]+pc$x[,1]*pc$rotation[1,1]
  # p.pc1[,2]=p.mean[2]+pc$x[,1]*pc$rotation[1,2]
  # p.pc1[,3]=p.mean[3]+pc$x[,1]*pc$rotation[1,3]
  p.pc1[,1]=p.mean[1]+pc$scores[,1]*pc$loadings[1,1]
  p.pc1[,2]=p.mean[2]+pc$scores[,1]*pc$loadings[1,2]
  p.pc1[,3]=p.mean[3]+pc$scores[,1]*pc$loadings[1,3]
  
  #estimate asmmetry
  p.Topt= p.pc1[,1]
  if(pc.vars==1) {p.CTmin= p.pc1[,2]-p.pc1[,3]; p.CTmax= p.pc1[,2]}
  if(pc.vars==2) {p.CTmin= p.pc1[,1]-p.pc1[,2]; p.CTmax= p.pc1[,1]+p.pc1[,3]}
  if(pc.vars==3) {p.CTmin= p.pc1[,2]; p.CTmax= p.pc1[,3]}
  if(pc.vars==4) {p.Topt= p.pc1[,2]; p.CTmin= p.pc1[,1]; p.CTmax= p.pc1[,3]}
  
  p.pc1$asym= (2*p.Topt -p.CTmax -p.CTmin)/(p.CTmax -p.CTmin )
  
  #pc2
  p.pc2=p
  # p.pc2[,1]=p.mean[1]+pc$x[,2]*pc$rotation[2,1]
  # p.pc2[,2]=p.mean[2]+pc$x[,2]*pc$rotation[2,2]
  # p.pc2[,3]=p.mean[3]+pc$x[,2]*pc$rotation[2,3]
  p.pc2[,1]=p.mean[1]+pc$scores[,2]*pc$loadings[2,1]
  p.pc2[,2]=p.mean[2]+pc$scores[,2]*pc$loadings[2,2]
  p.pc2[,3]=p.mean[3]+pc$scores[,2]*pc$loadings[2,3]
  
  #estimate asmmetry
  p.Topt= p.pc1[,1]
  if(pc.vars==1) {p.CTmin= p.pc2[,2]-p.pc2[,3]; p.CTmax= p.pc2[,2]}
  p.pc2$asym= (2*p.Topt -p.CTmax -p.CTmin)/(p.CTmax -p.CTmin )
  
  #pc3
  p.pc3=p
  # p.pc3[,1]=p.mean[1]+pc$x[,3]*pc$rotation[3,1]
  # p.pc3[,2]=p.mean[2]+pc$x[,3]*pc$rotation[3,2]
  # p.pc3[,3]=p.mean[3]+pc$x[,3]*pc$rotation[3,3]
  p.pc3[,1]=p.mean[1]+pc$scores[,3]*pc$loadings[3,1]
  p.pc3[,2]=p.mean[2]+pc$scores[,3]*pc$loadings[3,2]
  p.pc3[,3]=p.mean[3]+pc$scores[,3]*pc$loadings[3,3]

  #combine data
  add=cbind(p.pc1, p.pc2, p.pc3)
  ps[inds,]=as.matrix(add)
    
} #end loop taxa

#add names
if(pc.vars==1) colnames(ps)=c('pc1.Topt','pc1.CTmax','pc1.breadth','pc1.asym','pc2.Topt','pc2.CTmax','pc2.breadth','pc2.asym','pc3.Topt','pc3.CTmax','pc3.breadth')
if(pc.vars==2) colnames(ps)=c('pc1.Topt','pc1.coolb', 'pc1.warmb','pc1.asym','pc2.Topt','pc2.coolb', 'pc2.warmb','pc3.Topt','pc3.coolb', 'pc3.warmb')
if(pc.vars==3) colnames(ps)=c('pc1.Topt','pc1.CTmin', 'pc1.CTmax','pc1.asym','pc2.Topt','pc2.CTmin', 'pc2.CTmax','pc3.Topt','pc3.CTmin', 'pc3.CTmax')
if(pc.vars==4) colnames(ps)=c('pc1.CTmin','pc1.Topt', 'pc1.CTmax','pc1.asym','pc2.CTmin','pc2.Topt', 'pc2.CTmax','pc3.CTmin','pc3.Topt', 'pc3.CTmax')

#combine
tpc.pca= cbind(tpc,ps)

colnames(pc.var)=c("taxa","pc1.var","pc2.var", "pc3.var")
pc.var$asym=0

#estimate CTmin and CTmax
if(pc.vars==1){tpc.pca$pc1.CTmin= tpc.pca$pc1.CTmax - tpc.pca$pc1.breadth
tpc.pca$pc2.CTmin= tpc.pca$pc2.CTmax - tpc.pca$pc2.breadth
tpc.pca$pc3.CTmin= tpc.pca$pc3.CTmax - tpc.pca$pc3.breadth}

if(pc.vars==2){tpc.pca$pc1.CTmin= tpc.pca$pc1.Topt - tpc.pca$pc1.coolb
tpc.pca$pc1.CTmax= tpc.pca$pc1.Topt + tpc.pca$pc1.warmb
tpc.pca$pc2.CTmin= tpc.pca$pc2.Topt - tpc.pca$pc2.coolb
tpc.pca$pc2.CTmax= tpc.pca$pc2.Topt + tpc.pca$pc2.warmb
tpc.pca$pc3.CTmin= tpc.pca$pc3.Topt - tpc.pca$pc3.coolb
tpc.pca$pc3.CTmax= tpc.pca$pc3.Topt + tpc.pca$pc3.warmb}

#------
#plots
#asymmetry
tpc.a1= tpc.pca[,c("pc1.Topt","pc1.asym", "asym","taxa")]
tpc.a1$pc="PC1"
names(tpc.a1)[1:2]=c("Topt","pc.asym")
tpc.a2= tpc.pca[,c("pc2.Topt","pc2.asym", "asym","taxa")]
tpc.a2$pc="PC2"
names(tpc.a2)[1:2]=c("Topt","pc.asym")
tpc.a= rbind(tpc.a1, tpc.a2)
tpc.a$lab2="PC1:blue, PC2:yellow"

fig3a= ggplot(tpc.a) + aes(x=Topt, y = pc.asym, color=pc)+geom_point(size=2, alpha=0.7)+ylab("asymmetry")+
  xlab("Topt (°C)")+facet_grid(taxa~lab2)+
  theme_bw()+theme(legend.position="bottom")+ylim(-1,1) + scale_color_viridis(discrete=TRUE)+guides(color=FALSE)

pc.var$pc=""
pc.var$label=paste(" PC1=",round(pc.var$pc1.var,2),"\n PC2=",round(pc.var$pc2.var,2), sep=" ") #"\n PC3=",round(pc.var$pc3.var,2),
fig3a= fig3a +geom_text(data= pc.var, mapping=aes(x=10, y=0.6,label=label))

#pc plots
#pc1 plots
out=t(apply(tpc.pca[,c("pc1.Topt","pc1.CTmin","pc1.CTmax")], MARGIN=1, FUN=tpc.mat))
colnames(out)= temps
tpc.pred= cbind(tpc, out)
#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 18:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))
tpc.l$lab="PC1"
tpc.pc1=tpc.l

#pc2 plots
out=t(apply(tpc.pca[,c("pc2.Topt","pc2.CTmin","pc2.CTmax")], MARGIN=1, FUN=tpc.mat))
colnames(out)= temps
tpc.pred= cbind(tpc, out)
#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 18:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))
tpc.l$lab="PC2"
tpc.pc2=tpc.l

#pc3 plots
out=t(apply(tpc.pca[,c("pc3.Topt","pc3.CTmin","pc3.CTmax")], MARGIN=1, FUN=tpc.mat))
colnames(out)= temps
tpc.pred= cbind(tpc, out)
#to long format
tpc.l<- tpc.pred %>%
  gather("temperature", "performance", 18:ncol(tpc.pred))
tpc.l$temperature= as.numeric(as.character(tpc.l$temperature))
tpc.l$lab="PC3"
tpc.pc3=tpc.l

#combine
tpc.l=rbind(tpc.pc1, tpc.pc2)

#plot
fig3b= ggplot(tpc.l)+aes(x=temperature, y = performance, color=asym, group=X)+facet_grid(taxa~lab)+geom_line()+
  theme_bw()+theme(legend.position="bottom")+scale_color_viridis(name="asymmetry")+ylim(0,1)+xlab("temperature (°C)")

#labels
pc.lab= c("ToptCTmaxBreadth","ToptCoolbWarmb","ToptCTminCTmax","CTminToptCTmax")

#PCA analysis as in Knies et al.
pdf(paste("Fig3_PCAs_",pc.lab[pc.vars],".pdf",sep=""), height = 10, width = 8)
fig3a +fig3b +plot_annotation(tag_levels = 'a') +plot_layout(nrow=1, guides = "collect", widths = c(0.35, 1)) & theme(legend.position = "bottom") 
dev.off()

#------
#plot write out loadings
pc.out=rbind(pc.load[,,1],pc.load[,,2],pc.load[,,3],pc.load[,,4],pc.load[,,5],pc.load[,,6])
pc.out=as.data.frame(round(pc.out,2))
colnames(pc.out)=c("PC1","PC2","PC3")
pc.out$var=rep(c("Topt","CTmax","breadth"),6)
pc.out$var= ordered(pc.out$var, levels=c("Topt","CTmax","breadth") )
pc.out$taxa=rep(taxas, each=3)

#write.csv(pc.out, "pcload_ToptCTmaxBreadth.csv")

#to long format
pc.plot <- melt(pc.out, id=c("var","taxa"))

#plot
pdf(paste("LoadingsPCAs_",pc.lab[pc.vars],".pdf",sep=""), height = 8, width = 6)
ggplot(pc.plot)+aes(x=var, y = value)+facet_grid(taxa~variable)+geom_col()+
  theme_bw()+theme(legend.position="bottom")+xlab("parameter")+ylab("loading")+geom_hline(yintercept=0)
dev.off()

#}#end PC vars loop


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
