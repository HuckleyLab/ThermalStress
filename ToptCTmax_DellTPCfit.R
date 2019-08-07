# Code to fit thermal reaction norms (also known as thermal performance curves) and estimate 5 temperature traits (Topt, Tmax, Tmin, temperature niche width, maximum growth rate)
# See the following papers for details about these traits and the thermal reaction norm function: 
# Thomas, M. K., C. T. Kremer, C. A. Klausmeier, and E. Litchman. 2012. A global pattern of thermal adaptation in marine phytoplankton. Science 338:1085?1088.
# Thomas, M. K., C. T. Kremer, and E. Litchman. 2016. Environment and evolutionary history determine the global biogeography of phytoplankton temperature traits. 
# Global Ecology and Biogeography 25:75?86.
#
# Script originally developed by Colin T. Kremer, but tinkered with by Mridul K. Thomas (the latter claims responsibility for all errors)
# For the latest version of this script, please visit http://mridulkthomas.weebly.com/data--code.html

# Check if required packages are installed and install them if unavailable
list.of.packages <- c('bbmle', 'rootSolve')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
library(bbmle)
library(rootSolve)

### Load data. Use your dataset here. Make sure you the data is formatted as this example file is, with columns titled 'curve.id', 'temperature', and 'growth.rate'
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/")
dat.full<-read.csv('Delletal2013_forfitting.csv')
#adjust names
dat.full$growth.rate=dat.full$TraitValueSI
str(dat.full)

##cut to unimodal
#dat.full=dat.full[which(dat.full$CitationID %in% c(61,91,3,36,6,101,49,136,176)),]
#dat.full$init.curve.id= dat.full$curve.id

ids=unique(dat.full$curve.id)
dat.full$curve.id= match(dat.full$curve.id, ids)

# Thermal reaction norm function 
# Modified from Norberg (2004). See Thomas et al. (2012) for explanation.

nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

# Create new vector of unique curve.id values
curve.id.list<-unique(dat.full$curve.id)	

# Create empty vectors to populate with parameter values and trait estimates
z.list<-rep(NA, length(curve.id.list))				#Parameter 'z'
w.list<-rep(NA, length(curve.id.list))				#Parameter 'w', which is the temperature niche width
a.list<-rep(NA, length(curve.id.list))				#Parameter 'a'
b.list<-rep(NA, length(curve.id.list))				#Parameter 'b'
topt.list<-rep(NA, length(curve.id.list))			#Topt, Optimum temperature for growth
maxgrowth.list<-rep(NA, length(curve.id.list))		#Maximum growth rate (i.e. growth rate at Topt)
rsqr.list<-rep(NA, length(curve.id.list))			#R^2 values for all fits
s.list<-rep(NA, length(curve.id.list))				#Error
n.list<-rep(NA, length(curve.id.list))				#Number of growth rate measurements used in the curve

# Loop through all curve.id.list values to estimate parameters for all curves

#change directory
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/BioTraits/")

for(i in 1:length(curve.id.list)){
  print(i)
  
  # Take a subset of the data corressponding to the ith curve.id.list value
  dat<-subset(dat.full,dat.full$curve.id==curve.id.list[i])
  
  # guess starting values for parameters 'z' and 'w'
  z.guess<-mean(dat$temperature[dat$growth.rate==max(dat$growth.rate)])		#starting estimates for 'z'
  w.guess<-diff(range(dat$temperature))										#starting estimates for niche width
  
  ## This loop fits the model using a range of different starting guesses. We choose the best one using AIC. This helps find good solutions even if there are
  # convergence problems.
  # Starting estimates for parameters 'a' and 'b' use a plausible range but with broadly spaced estimates to speed up fitting. 
  avals<-seq(-0.2,1.2,0.1)		
  bvals<-seq(-0.2,0.3,0.05)
  mod.list<-list()
  AIC.list<-c()
  
  for(ia in 1:length(avals)){
    for(ib in 1:length(bvals)){
      a.guess<-avals[ia]
      b.guess<-bvals[ib]
      res2<-try(fit<-mle2(dat$growth.rate~dnorm(mean=nbcurve(dat$temperature,z=z,w=w,a=a,b=b),sd=s),start=list(z=z.guess,w=w.guess,a=a.guess,b=b.guess,s=0.3),
                          skip.hessian=TRUE,data=dat))
      if(class(res2)!="try-error"){
        mod.list<-append(mod.list,fit)
        AIC.list<-append(AIC.list,AIC(fit))
      }
    }
  }
  
  # Identify the best model from the list and save coefficients and R^2 values
  if(!is.null(AIC.list)){
    bestmodind<-which(AIC.list==min(AIC.list))
    if(length(bestmodind)>1){
      bestmodind<-sample(bestmodind,1)
    }
    bestmod<-mod.list[[bestmodind]]
    cfs<-coef(bestmod)
    expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
    rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)
  }
  
  # If the quick fit yielded poor results (low R^2), try a more thorough search through parameter space
  if(rsqr<0.95){
    avals<-seq(-0.2,1.2,0.02)
    bvals<-seq(-0.2,0.3,0.02)
    mod.list<-list()
    AIC.list<-c()
    for(ia in 1:length(avals)){
      for(ib in 1:length(bvals)){
        a.guess<-avals[ia]
        b.guess<-bvals[ib]
        res2<-try(fit<-mle2(dat$growth.rate~dnorm(mean=nbcurve(dat$temperature,z=z,w=w,a=a,b=b),sd=s),start=list(z=z.guess,w=w.guess,a=a.guess,b=b.guess,s=0.3),
                            skip.hessian=TRUE,data=dat))
        if(class(res2)!="try-error"){
          mod.list<-append(mod.list,fit)
          AIC.list<-append(AIC.list,AIC(fit))
        }
      }
    }
    # Identify the best model from the list and save coefficients and R^2 values
    bestmodind<-which(AIC.list==min(AIC.list))
    if(length(bestmodind)>1){
      bestmodind<-sample(bestmodind,1)
    }
    
    bestmod<-mod.list[[bestmodind]]
    cfs<-coef(bestmod)
    expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
    rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)
  }
  #Save .png plot with fitted curve. File is saved with the curve.id.list value as the name
  png(paste(curve.id.list[i],'.png',sep=''))
  plot(dat$growth.rate~dat$temperature,ylim=c(pmin(0,min(dat$growth.rate)),max(dat$growth.rate)+(0.2)*max(dat$growth.rate)),main=curve.id.list[i],
       xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
  curve(nbcurve(x,cfs[1],cfs[2],cfs[3],cfs[4]),col='red', lwd=2,add=TRUE)
  dev.off()
  
  # Use the curve fit to find Topt and the estimated maximum growth rate (i.e. growth rate at Topt)
  grfunc<-function(x){
    -nbcurve(x,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
  }
  optinfo<-optim(c(x=cfs[[1]]),grfunc)
  opt<-optinfo$par[[1]]
  maxgrowth<- -optinfo$value
  
  #stash results		
  rsqr.list[i]<-rsqr
  z.list[i]<-cfs[[1]]
  w.list[i]<-cfs[[2]]
  a.list[i]<-cfs[[3]]
  b.list[i]<-cfs[[4]]
  s.list[i]<-cfs[[5]]
  topt.list[i]<-opt
  maxgrowth.list[i]<-maxgrowth
  n.list[i]<-length(dat$temperature)
}

fits<-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list)

# There are 2 ways of estimating Tmax and Tmin. Both work in most cases, but fail for rare parameter combinations. I favour one method for simplicity, 
# but I include code for the second as well (commented out) in case this gives poor results. 
# Please do not trust these values without looking at your plots to make sure that they make sense! 
# Additionally, if these estimates are well outside the range of your measurements (>5 degrees C is my rule of thumb), I would not trust them either. 

# Method 1:

fits$tmax<-fits$z.list+(fits$w.list/2)
fits$tmin<-fits$z.list-(fits$w.list/2)

#Method 2. 
# for (i in 1:length(curve.id.list)){
# print(i)
# nbcurve.tmax<-function(x){
# nb<-nbcurve(x,fits$z.list[i],fits$w.list[i],fits$a.list[i],fits$b.list[i])
# nb
# }
# fits$tmax[i]<-uniroot.all(nbcurve.tmax,c(fits$topt.list[i],150))[1]
# fits$tmin[i]<-uniroot.all(nbcurve.tmax,c(0,fits$topt.list[i]))[1] 
# }

write.csv(fits,'BioTraitsTPCfits.csv')


