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

#----
#phytoplankton
#https://www.ibi.vu.nl/programs/phylopars/

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
  data = tmat, boot = 200)
topt.ctmin=cor_phylo(variates = ~ Topt + CTmin, species = ~ species, phy = phy,
                        data = tmat, boot = 200)
topt.ctmax=cor_phylo(variates = ~ Topt + CTmax, species = ~ species, phy = phy,
                        data = tmat, boot = 200)
topt.tol=cor_phylo(variates = ~ Topt + CTmax.Topt.breadth, species = ~ species, phy = phy,
                        data = tmat, boot = 200)





#B matrix
