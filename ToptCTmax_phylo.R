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

#--------------
#estimate trait evolution across phylogeny

#http://treethinkers.org/wp-content/uploads/2013/01/Bodega2013.pdf
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12977
library(phytools)

#lizard phylogeny
#https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-13-93

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/lizard_phy/")
phy<-read.newick("pyronetal2013.txt")

#match to tips
#LIZARDS- Huey
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
liz= read.csv('Hueyetal2009.csv', na.strings ='-9999')
#match
match1= match(liz$Species, phy$tip.label)
phy.prune<-drop.tip(phy,phy$tip.label[-na.omit(match1)])

#order traits
liz.match=liz[match(phy.prune$tip.label, liz$Species),c(4,17:19)]
liz.match= na.omit(liz.match)
liz.match$Species= as.character(liz.match$Species)
rownames(liz.match)= liz.match$Species

#fit
#http://www.eve.ucdavis.edu/~wainwrightlab/Roi/Site/Teaching_files/Intro2Phylo_S5.R
#http://www.phytools.org/Bariloche2016/ch/3/Cont-char-models-solution.html
library(geiger)
#sigmasq is evolutionary rate parameters
fitBM<- fitContinuous(phy.prune, liz.match[,2:4], model="BM")
fitOU<- fitContinuous(phy.prune, liz.match[,2:4], model="OU")
fitEB<- fitContinuous(phy.prune, liz.match[,2:4], model="EB")

aic.vals<-setNames(c(fitBM$newTopt$opt$aicc,fitOU$newTopt$opt$aicc,fitEB$newTopt$opt$aicc),
                   c("BM","OU","EB"))
aic.w(aic.vals)

##Partition variance to phylogeney and environment
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/code/")
#source functions
source("brutal.r")
source("distance2_25Aug2010.r")
source("pglm.r")
library(spdep)

#compile lizard data
phy.prune<-drop.tip(phy.prune, phy.prune$tip.label[-match(liz.match$Species, phy.prune$tip.label) ])
liz.match= liz[match(phy.prune$tip.label, liz$Species),]
liz.match$Species= as.character(liz.match$Species)

#APPROACH 1: uses distance2
fm <- cpglm(newTopt ~ Lat - 1, liz.match, phy.prune)
fm1 <- cpglm.lambda(newTopt ~ Lat - 1, liz.match, phy.prune, 1) #sets lambda to 1
fm2 <- pglm( newTopt ~ Lat, liz.match, vcv.phylo(phy.prune) )

fm.estlambda<- cpglm.estlambda(newTopt ~ Lat, liz.match, phy.prune) #PGLM for estimation of lambda
fm3 <-  cpglm.spatial( newTopt ~ Lat, liz.match, phy.prune, liz.match$Lat, liz.match$Long, 1, 0) #function( formula, data, phylo, lat, lon, lambda, phi) 

#Fit phi and lambda
fit.temp <- optimise.cpglm(newTopt ~  -1, liz.match, phy.prune, liz.match$Lat, liz.match$Long) #estimate phy/lambda
fitphi.temp <- optimisephi.cpglm( newTopt ~  -1, liz.match, phy.prune, liz.match$Lat, liz.match$Long) #set lambda to zero to estimate phi singly
#may need to manually step through to get to work

#APPROACH 2: uses brutual
Vliz <- vcv.phylo(phy.prune)
#Dliz <- dist.mat( tb$Lat, tb$Lon, tb.back$Genus_sp)
Dliz <- spDists(as.matrix(liz.match[,c(13,11)]), longlat = TRUE);
rownames(Dliz)<-liz.match$Species

fm <- pglmSpatial( newTopt ~ 1, liz.match, Vliz, Dliz)

fm <- pglm( newTopt ~ 1, liz.match, Vliz)
fm1 <- pglmEstLambda( newTopt ~ Lat, liz.match, Vliz)
fm2 <- pglmSpatialFit( newTopt ~ Lat, liz.match, Vliz, Dliz)


