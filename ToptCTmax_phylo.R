library(ggplot2)

#load data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
tpc=read.csv("tpcs_wSource_Mar2021.csv")

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
#load phylogenies

#http://treethinkers.org/wp-content/uploads/2013/01/Bodega2013.pdf
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12977
library(phytools)

#Lizards
#https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-13-93

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/phylogeny/")
phy<-read.newick("pyronetal2013.txt")

#match to tips
liz= tpc[tpc$taxa=="lizards",]
#match
match1= match(liz$species, phy$tip.label)
liz.prune<-drop.tip(phy,phy$tip.label[-na.omit(match1)])
liz.match=liz[match(liz.prune$tip.label, liz$species),]
#liz.match= na.omit(liz.match)

#----
#ant newick
#https://www.antwiki.org/wiki/Phylogeny_of_Formicidae
phy= read.newick("ant.txt")

#match to tips
ant= tpc[tpc$taxa=="ants",]
#match
match1= match(ant$genus, phy$tip.label)
ant.prune<-drop.tip(phy,phy$tip.label[-na.omit(match1)])
ant.match=ant[match(ant.prune$tip.label, ant$genus),]
#use genus name for species for code simplicity
ant.match$species= ant.match$genus

----
#insects

#match to tips
ins= tpc[tpc$taxa=="insects",]
ins$genspec= paste(ins$genus, ins$species, sep=" ")

#https://phylot.biobyte.de/
#approaches:
#https://taylorreiter.github.io/2017-07-28-Taxonomy-from-Species-Name-in-R/

library(taxize)
library(rotl)
library(ape)

ins_class <- classification(unique(ins$genspec), db = "ncbi")
ins_tree <- class2tree(ins_class, check = TRUE)
plot(ins_tree)
ins_phy= ins_tree$phylo

# specs <- tnrs_match_names(ins$genspec)
# specs <- specs[!is.na(specs$ott_id),]
# specs <- specs[-which(specs$ott_id==3360350),]
# specs <- specs[-which(specs$ott_id==3359321),]
# 
# tree <- tol_induced_subtree(ott_ids = specs$ott_id)
# plot(tree)

#---
#match
match1= match(ins$genspec, ins_phy$tip.label)

ins.prune<-drop.tip(ins_phy,ins_phy$tip.label[-na.omit(match1)])
ins.match=ins[match(ins.prune$tip.label, ins$genspec),]

#----
#phytoplankton
#https://www.ibi.vu.nl/programs/phylopars/

#match to tips
plank= tpc[tpc$taxa=="plankton",]
plank$genspec= paste(plank$genus, plank$species, sep=" ")

plank_class <- classification(unique(plank$genspec), db = "ncbi")
plank_tree <- class2tree(plank_class, check = TRUE)
plot(plank_tree)
plank_phy= plank_tree$phylo

#---
#match
match1= match(plank$genspec, plank_phy$tip.label)

plank.prune<-drop.tip(plank_phy,plank_phy$tip.label[-na.omit(match1)])
plank.match=plank[match(plank.prune$tip.label, plank$genspec),]

#----
#fish

#match to tips
fish= tpc[tpc$taxa=="fish",]
fish$genspec= paste(fish$genus, fish$species, sep=" ")
#only three species

#----------------
#phyr
#https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13471
#https://daijiang.github.io/phyr/articles/phyr_example_empirical.html
#Zheng, L., A. R. Ives, T. Garland, B. R. Larget, Y. Yu, and K. F. Cao. 2009. New multivariate tests for phylogenetic signal and trait correlations applied to ecophysiological phenotypes of nine Manglietia species. Functional Ecology 23:1059–1069.

library(phyr)
library(ape)
library(dplyr)

phy= liz.prune
tmat= liz.match

#lizards
topt.asym=cor_phylo(variates = ~ Topt + asym, species = ~ species, phy = phy,
  data = tmat, boot = 100)
topt.ctmin=cor_phylo(variates = ~ Topt + CTmin, species = ~ species, phy = phy,
                        data = tmat, boot = 100)
topt.ctmax=cor_phylo(variates = ~ Topt + CTmax, species = ~ species, phy = phy,
                        data = tmat, boot = 100)
topt.tol=cor_phylo(variates = ~ Topt + CTmax.Topt.breadth, species = ~ species, phy = phy,
                        data = tmat, boot = 100)

#B matrix
topt.asym$B

#--------
library(phylolm)
#phyloglm

phy= plank.prune
tmat= plank.match
rownames(tmat)= tmat$genspec

fit = phyloglm(asym~Topt,phy=phy,data=tmat,boot=50)
summary(fit)
coef(fit)
vcov(fit)

set.seed(123456)
tre = phy
x = rTrait(n=1,phy=phy)
X = cbind(rep(1,50),x)
y = rbinTrait(n=1,phy=tre, beta=c(-1,0.5), alpha=1 ,X=X)
dat = data.frame(trait01 = y, predictor = x)
fit = phyloglm(trait01~predictor,phy=tre,data=dat,boot=100)

#====================
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

