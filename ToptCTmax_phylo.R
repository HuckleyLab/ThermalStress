library(ggplot2)
library(taxize)
library(rotl)
library(ape)
library(nlme)
library(geiger)
library(phytools)

#load data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
tpc=read.csv("tpcs_wSource_Mar2021.csv")

#drop data without all metrics
tpc= subset(tpc, !is.na(tpc$CTmin) & !is.na(tpc$Topt) & !is.na(tpc$CTmax) )

#calculate asymetry
#tpc$asym= (tpc$CTmax - tpc$Topt)/(tpc$Topt- tpc$CTmin )
#OR INVERSE: tpc$asym= (tpc$Topt- tpc$CTmin )/(tpc$CTmax - tpc$Topt)
tpc$asym= (2*tpc$Topt-tpc$CTmax - tpc$CTmin)/(tpc$CTmax-tpc$CTmin )

#calculate breadth
tpc$breadth= tpc$CTmax -tpc$CTmin
#calculate declining breadth
tpc$CTmax.Topt.breadth= tpc$CTmax - tpc$Topt

#--------------
#load phylogenies

#http://treethinkers.org/wp-content/uploads/2013/01/Bodega2013.pdf
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12977

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

ins_class <- classification(unique(ins$genspec), db = "ncbi")
ins_tree <- class2tree(ins_class, check = TRUE)
plot(ins_tree)
ins_phy= ins_tree$phylo

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
#http://www.phytools.org/Cordoba2017/ex/4/PGLS.html

for(tax in 1:4){

#ant, liz, insects, plankton
if(tax==1){  
phy=ant.prune
tmat=ant.match
rownames(tmat)= tmat$genus
}
  if(tax==2){ 
phy= liz.prune
tmat= liz.match
rownames(tmat)= tmat$species
  }
  if(tax==3){
phy= ins.prune
tmat= ins.match
rownames(tmat)= tmat$genspec
  }
  if(tax==4){
phy= plank.prune
tmat= plank.match
rownames(tmat)= tmat$genspec
}
#---
#name.check(phy,tmat)

#ape approach
#mod1<-gls(asym~Topt, data=tmat, correlation=corBrownian(1, phy))
mod.asym<-gls(asym~Topt, data=tmat, correlation=corPagel(0.5,phy))
mod.ctmin<-gls(CTmin~Topt, data=tmat, correlation=corPagel(1,phy))
mod.ctmax<-gls(CTmax~Topt, data=tmat, correlation=corPagel(1,phy))
mod.breadth<-gls(breadth~Topt, data=tmat, correlation=corPagel(1,phy))
mod.wbr<-gls(CTmax.Topt.breadth~Topt, data=tmat, correlation=corPagel(1,phy))

coefs= c(summary(mod.asym)$tTable[2,],summary(mod.ctmin)$tTable[2,],summary(mod.ctmax)$tTable[2,],summary(mod.breadth)$tTable[2,],summary(mod.wbr)$tTable[2,])
lambdas= c(intervals(mod.asym)$corStruct[2],intervals(mod.ctmin)$corStruct[2],intervals(mod.ctmax)$corStruct[2],intervals(mod.breadth)$corStruct[2],intervals(mod.wbr)$corStruct[2])

if(tax==1){coefs.all=coefs; lambdas.all=lambdas}
if(tax>1){
  coefs.all= rbind(coefs.all, coefs)
  lambdas.all= rbind(lambdas.all, lambdas)
}

} #end loop taxa

colnames(coefs.all)= rep(c("slope","se","t value","p"),5)
#round
coefs.all= signif(coefs.all, 2)
#add taxa
rownames(coefs.all)= taxas[c(1:2,4:5)]

#round
lambdas.all= signif(lambdas.all, 2)
#add taxa
rownames(lambdas.all)= taxas[c(1:2,4:5)]

#write out
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/figures/")
write.csv(coefs.all, "Fig2stats_phylo.csv", row.names = TRUE)
write.csv(lambdas.all, "Fig2stats_lambdas.csv", row.names = TRUE)

