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

#=========================
#Fit Rezende data 
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/Rezende")
rez.photo= read.csv("RezendeTableA1.csv")
rez.fit= read.csv("RezendeTableA3.csv")

### Load data. Use your dataset here. Make sure you the data is formatted as this example file is, with columns titled 'curve.id', 'temperature', and 'growth.rate'
dat.full= rez.fit
#adjust names
colnames(dat.full)[2:4]=c('curve.id', 'temperature', 'growth.rate')

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
mint.list<-rep(NA, length(curve.id.list))	  #min measurement T
maxt.list<-rep(NA, length(curve.id.list)) #max measurement T

# Loop through all curve.id.list values to estimate parameters for all curves

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
  mint.list[i]<- min(dat$temperature)
  maxt.list[i]<- max(dat$temperature)
}

fits<-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list,mint.list,maxt.list)

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

write.csv(fits,'RezFitFits_Mar2021.csv')

#==================
#Add species list

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/Rezende")
#photosynthesis
rez.photo.fits= read.csv("RezPhotoFits_Mar2021.csv")
rez.photo.fits$species= rez.photo[match(rez.photo.fits$curve.id.list, rez.photo$ID),"Species"]
write.csv(rez.photo.fits,'RezPhotoFits_Mar2021.csv')

#insect fitness
rez.fit.fits= read.csv("RezFitFits_Mar2021.csv")
rez.fit.fits$species= rez.fit[match(rez.fit.fits$curve.id.list, rez.fit$ID),"Species"]
write.csv(rez.fit.fits,'RezFitFits_Mar2021.csv')

#==========================
#Combine insect fitness data
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/Rezende")

rez= read.csv("RezFitFits_Mar2021.csv")

#add reference
rez.a3= read.csv("RezendeTableA3.csv")
match1= match(rez$curve.id.list, rez.a3$ID)
rez$reference= rez.a3$Reference[match1]

#check distance between estimates and measurements
plot(rez$tmax-rez$maxt.list, rez$tmax)
abline(v=7)
plot(-rez$tmin+rez$mint.list, rez$tmin)
abline(v=10)

#cut tpcs with >x degrees between last temperature and CTmin or CTmax estimate
rez= rez[-which((rez$tmax-rez$maxt.list)>7), ]
rez= rez[-which((-rez$tmin+rez$mint.list)>10), ]

#cut outliers
#rez=rez[-c(18,33),]

#read Deutsch et al.
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/DeutschData/")
ins= read.csv('InsectFit_Deutsch.csv')
ins$spgen=paste(ins$gen, ins$sp, sep=" ")

#find duplicate species
match1= match(ins$spgen, rez$species)
matched= which(!is.na(match1))

ins[!is.na(match1),"Reference"]
rez[na.omit(match1),"reference"]

ins[matched[c(3,4)],]
rez[na.omit(match1)[c(3,4)],]

#add rez data
rez.add= na.omit(match1)
rez.drop= rez.add[-c(3,4)]
rez.add= rez[-rez.drop,]
rez.add$Lat=NA
rez.add$Long=NA
rez.add$Order=NA
rez.add$Location=NA
rez.add= rez.add[,c("Lat","Long","Order","species","Location","reference","tmin","topt.list","tmax")]

#combine
ins.add= ins[,c("Lat","Long","Order","spgen","Location","Reference","Ctmin","Topt","CTmax")]
ins.comb= rbind(ins.add, setNames(rez.add, names(ins.add)))

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/DeutschData/")
write.csv(ins.comb, 'InsectFit_DeutschRezende_Mar2021.csv')

#=======================
#Add Rezende data 
setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/ThermalStress/data/CTlimits/Rezende")
rez= read.csv("RezendeAppendixC.csv")

#library(rTPC), rezende_2019

tpc.rezende<- function (temp, q10, C, Tth, d) 
{ temp=as.numeric(temp)
p= rep(NA, length(temp))
inds= which(temp<=Tth)
p[inds]<- C*2.71828^(temp[inds]*log(q10)/10)
inds= which(temp>Tth)
p[inds]<- C*2.71828^(temp[inds]*log(q10)/10)*(1-d*(temp[inds]-Tth)^2)
p[p<0]<-0
return(p)
}

ctmin.rezende<- function(tpc){
  tpc=as.numeric(tpc)
  temp= seq(tpc[5]-50,tpc[5],0.2)
  ps=tpc.rezende(temp, q10=tpc[1], C=tpc[2], Tth=tpc[3], d=tpc[4])
  temp[which.max(ps<0.2*tpc[6])]
}

ctmax.rezende<- function(tpc){
  tpc=as.numeric(tpc)
  temp= seq(tpc[5],tpc[5]+20,0.2)
  ps=tpc.rezende(temp, q10=tpc[1], C=tpc[2], Tth=tpc[3], d=tpc[4])
  temp[which.max(ps<0.2*tpc[6])]
}

CTmax= apply(rez[,c("Q10","C","Tth","d","Topt","Pmax")], FUN=ctmax.rezende, MARGIN=1)
plot(CTmax, rez$Ctmax)
abline(a=0, b=1)

ind=10
plot(1:70, tpc.rezende(1:70,q10=rez[ind,"Q10"], C=rez[ind,"C"], Tth=rez[ind,"Tth"], d=rez[ind,"d"]), type="l")
points(ctmax.rezende(tpc=rez[ind,c("Q10","C","Tth","d","Topt","Pmax")]),0)
points(rez[ind,"Ctmax"],0,pch="*")

rez$Ctmin= apply(rez[,c("Q10","C","Tth","d","Topt","Pmax")], FUN=ctmin.rezende, MARGIN=1)
rez$Genus=NA
rez$family=NA
rez$lat=NA

#add data
tpc2= rez[,c("Species","Genus","family","Ctmin","Ctmax","Topt")] 
tpc2$habitat="terrestrial"
tpc2$lat=NA
tpc2$lon= NA
tpc2$taxa= rez$Type
#bind
tpc= rbind(tpc, setNames(tpc2, names(tpc)))

#===============
#Using R packages

#try rTPC package
#https://github.com/padpadpadpad/rTPC

nls_multstart(rate ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
              data = dat.tpc,
              iter = 500,
              start_lower = get_start_vals(dat.tpc$temp, dat.tpc$rate, model_name = 'lactin2_1995') - 2,
              start_upper = get_start_vals(dat.tpc$temp, dat.tpc$rate, model_name = 'lactin2_1995') + 2,
              supp_errors = 'Y')

dat.tpc$process="perf"
d_1=dat.tpc


d_models <- group_by(dat.tpc, process) %>%
  nest() %>%
  mutate(., thomas_2012 = map(data, ~nls_multstart(rate ~ thomas_2012(temp = temp, a, b, c, topt),
                                                   data = .x,
                                                   iter = 500,
                                                   start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') - 1,
                                                   start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') + 2,
                                                   supp_errors = 'Y',
                                                   lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2012'))),
         gaussian = map(data, ~nls_multstart(rate ~ gaussian_1987(temp = temp, rmax, topt, a),
                                             data = .x,
                                             iter = 500,
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 2,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 2,
                                             supp_errors = 'Y')),
         quadratic = map(data, ~nls_multstart(rate ~ quadratic_2008(temp = temp, a, b, c),
                                              data = .x,
                                              iter = 500,
                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 1,
                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 1,
                                              supp_errors = 'Y')),
         weibull = map(data, ~nls_multstart(rate ~ weibull_1995(temp = temp, a, topt, b, c),
                                            data = .x,
                                            iter = 1000,
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') -2,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') + 2,
                                            supp_errors = 'Y')),
         rezende = map(data, ~nls_multstart(rate ~ rezende_2019(temp = temp, a, q10, b, c),
                                            data = .x,
                                            iter = 500,
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 0.8,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 1.2,
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                            supp_errors = 'Y')),
         beta = map(data, ~nls_multstart(rate ~ beta_2012(temp = temp, a, b, c, d, e),
                                         data = .x,
                                         iter = 500,
                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') -10,
                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') + 10,
                                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         supp_errors = 'Y')),
         modgaussian = map(data, ~nls_multstart(rate ~ modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                                                data = .x,
                                                iter = 500,
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') - 1,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') +1,
                                                supp_errors = 'Y')),
         boatman = map(data, ~nls_multstart(rate ~ boatman_2017(temp = temp, rmax, tmin, tmax, a, b),
                                            data = .x,
                                            iter = 500,
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') -1,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 1,
                                            supp_errors = 'Y')),
         thomas_2017 = map(data, ~nls_multstart(rate ~ thomas_2017(temp = temp, a, b, c, d, e),
                                                data = .x,
                                                iter = 500,
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') *0.5,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') *1.5,
                                                supp_errors = 'Y')))


#===================
# stack models
d_stack <- gather(d_models, 'model', 'output', 6:ncol(d_models))

# preds
newdata <- tibble(temp = seq(min(d_1$temp), max(d_1$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(output, augment, newdata = newdata)) %>%
  unnest(preds)

# estimate parameters
params <- d_stack %>%
  mutate(., est = map(output, tidy)) %>%
  select(., -c(data, output)) %>%
  unnest(est)

# plot fit
ggplot(d_preds, aes(temp, rate)) +
  geom_point(aes(temp, rate), d_1) +
  geom_line(aes(temp, .fitted, col = model)) +
  facet_wrap(~model, labeller = labeller(model = label_facets_num)) +
  theme_bw(base_size = 16) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  xlab('Temperature (ºC)') +
  ylab('rate') +
  geom_hline(aes(yintercept = 0), linetype = 2)

#use rezende, function (temp, q10, a, b, c)
plot(10:40,rezende_2019(10:40, q10=2.27, a=0.109, b=9.02, c=0.00116), type="b")

nls_multstart(rate ~ rezende_2019(temp = temp, a, q10, b, c),
              data = dat.tpc,
              iter = 500,
              start_lower = get_start_vals(dat.tpc$temp, dat.tpc$rate, model_name = 'rezende_2019') * 0.8,
              start_upper = get_start_vals(dat.tpc$temp, dat.tpc$rate, model_name = 'rezende_2019') * 1.2,
              upper = get_upper_lims(dat.tpc$temp, dat.tpc$rate, model_name = 'rezende_2019'),
              lower = get_lower_lims(dat.tpc$temp, dat.tpc$rate, model_name = 'rezende_2019'),
              supp_errors = 'Y')

