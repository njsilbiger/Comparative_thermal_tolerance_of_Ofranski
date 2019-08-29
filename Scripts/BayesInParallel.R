## R script to construct Bayesian models to compare photosynthesis and respiration rates between Panama and Bermuda
## Created by N. Silbiger


rm(list=ls())

# load libraries###############
library(tidyverse)
library(rjags)
library(lattice)
library(stringi)
library(reshape2)
library(broom)
library(dclone)
library(data.table)
library(gridExtra)
library(tidybayes)
library(gtable)
library(grid)

### Read in Bermuda Data ####
Data<-read.csv('Data/BermudaRates.csv')

#####
# replace negative photosythesis values with a 0. A negative value means there is no P anymore
Data$Pgross_umol.cm2.hr[Data$Pgross_umol.cm2.hr<0]<-0

# take the aboslute value of respiration so that it is positive
Data$Rdark_umol.cm2.hr<-abs(Data$Rdark_umol.cm2.hr)

#make it long format with a colum for photosynthesis and resipiration
Data<-rbind(Data,Data,Data)
Data$PR<-c(rep('Photosynthesis',32),rep('Respiration',32), rep('Net Photosynthesis',32))
Data$PR<-as.factor(Data$PR)
Data$Rate<-c(Data$Pgross_umol.cm2.hr[1:32], Data$Rdark_umol.cm2.hr[1:32], Data$Pnet_umol.cm2.hr[1:32])
Data[,c("Pgross_umol.cm2.hr","Rdark_umol.cm2.hr", "Pnet_umol.cm2.hr")]<-list(NULL)

# add a column with log transformed data for the schoolfield equation
Data$Rate.ln<-log(Data$Rate+1)

Data$Light_Dark<-NA
Data$Light_Dark[Data$PR=='Photosynthesis']<-'Light'
Data$Light_Dark[Data$PR=='Respiration']<-'Dark'
Data$Light_Dark[Data$PR=='Net Photosynthesis']<-'Net'
Data$Light_Dark<-as.factor(Data$Light_Dark)

# join with NEC data
NECData<-read.csv('Data/Bermuda/TAData/NECData3.csv')
NECData$Rate<-NECData$NEC
NECData$Rate.ln<-log(NECData$NEC+1)

#remove stuff I dont need
NECData<-NECData[,c(2:4,6,17,18,31,32,33)]


# BIND THE NEC and PR data
Data<-bind_rows(Data,NECData)

# add a column with just the individual ID
Data<-Data%>%
  separate(col =Fragment.ID, into ='Individual.ID', sep = '-', extra='drop')

############# Read in Panama Data

PData<-read.csv('Data/PanamaRates.csv')
PData$Light_Dark<-as.character(PData$Light_Dark)
PData$Light_Dark[PData$PR=='Net Photosynthesis']<-as.character('Net')
PData$Light_Dark[is.na(PData$Light_Dark)]<-'Net'
PData$Light_Dark<-as.factor(PData$Light_Dark)
# Make a temp category for Panama
PData$Temp.Cat<-round(PData$Temp.C,0)

## combine the Panama and Bermuda data so that I can compare them in the same model

# add a row of locations to the Bermuda Dataset
Data<-Data %>% mutate(Location = 'Bermuda')%>%
  select(Temp.C, Date, Species, Fragment.Number, Genotype, Temp.Cat, PR, Light_Dark, Rate, Location, Individual.ID) # filter out the uncessesary columns and rename them so that they are the same
# put in the actual temp values


AllData<- PData %>%
  select(Temp.C, Date, Species, Fragment.Number, Genotype, Temp.Cat, PR, Light_Dark, umol.cm2.hr)%>% # filter out the uncessesary columns and rename them so that they are the same
  rename(Rate = umol.cm2.hr) %>% ### next make an individual ID
  mutate(Individual.ID=paste0(Species,Genotype))%>% # make an individual ID similar to Bermuda
  mutate(Location = 'Panama')%>% # make a column that says Panama
  mutate(PR = recode(PR,  "Gross Photosynthesis" = 'Photosynthesis')) %>%
  bind_rows(.,Data) %>% # bind the Panama and Bermuda Data
  mutate(unique.ID = paste0(Individual.ID,'.', Location)) %>% # add a unique identifier for each individual
  filter(Species =='OF')%>% # filter out all other species except for OF
  mutate_if(is.character, as.factor) %>%  # convert all the characters to factors
  select(Temp.C, Date, Species, Genotype, Temp.Cat, PR, Light_Dark, Rate, Individual.ID, Location, unique.ID)


#### paramaterize a model that is nor hierarchical #####
# JAGS model name for the broad model. All the JAGS code is in this file
modelBroad = "Scripts/JAGSnestedLocation.R"

# ##### Photosynthesis model for OF ########## 
BBData<-AllData
BBData<-BBData[BBData$Species=='OF' & BBData$PR=='Photosynthesis',]

unq_ind_group <- BBData[,10:11]
unq_ind_group <- unq_ind_group[!duplicated(unq_ind_group),] # this is needed for the heirarchical strucutre it creates a vector of the unique IDs

#make the rates log x+1
BBData$Rate.ln<-log(BBData$Rate+1)

# data list for the model
data.list<-list(
  y=BBData$Rate.ln, # log O2 rate
  x.hat = seq(14+273.15, 40+273.15, 0.1), # for predictions, 14 to 40 degrees C by 0.1 deg C
  W=diag(4), # for scaled wishart distribution
  N = length(BBData$Rate.ln), # first level of hierarchy (all observations)
  temp=BBData$Temp.C+273.15, # temperature data
  L=length(unique(BBData$Location)), # second level of hierarchy (locations)
  loc=as.factor(BBData$Location),
  G = length(unique(BBData$unique.ID)), # third level of hierarchy (genotypes)
  geno  = as.factor(BBData$unique.ID),
  ind_per_group = as.factor(rep(c('Bermuda','Panama'),4))) # vector of individual IDs


y=BBData$Rate.ln # log O2 rate for shiny stan

#varibles to save for inference
variables.zp<-c("lnc.loc","E.loc" ,"Eh.loc",'Topt', 'sigma.obs','sigma', "lnc.geno", 'sigma.geno', 'lnc.pred', 'Tmax')
keep<-c("pvalue.mean","pvalue.sd","pvalue.cv","E.loc" ,"Eh.loc",'Topt', 'y.hat',"y.new", "lnc.loc", 'Tmax','lnc.pred', 'sigma.obs','sigma','sigma.geno') # params to keep for estimates

## MCMC settings
setts <- list('n.iter' = 250000, 'n.thin' = 200, 'n.burn' = 200000)
setts.m <- 10  # this is 1.5M iterations with a 1.5K thinning and 1.2M burn-in

mSetts <- 1
if(mSetts) setts <- lapply(setts, function(v) v * setts.m)
setts$n.chains <- 3 # number of chains
setts$n.adapt <- 500000 # adaptation phase
n.cores <- 3 # number of cores to use in parallel computing
set.seed(129)

# run the bayesian model across multiple cores using the dclone libirary
cl <- makePSOCKcluster(n.cores)
tmp <- clusterEvalQ(cl, library(dclone))
parLoadModule(cl, 'glm')
parLoadModule(cl, 'lecuyer')
parLoadModule(cl, 'dic')
zm.p <- jags.parfit(cl = cl, data = data.list, params = variables.zp, model = modelBroad, 
                    n.chains = setts$n.chains, 
                    n.adapt = setts$n.adapt, 
                    n.update = setts$n.burn,
                    n.iter = setts$n.iter, 
                    thin = setts$n.thin)

stopCluster(cl)

# plot the trace plots
png('Output/traceplots/tracephoto%03d.png', 1500,1500, res = 150)
plot(zm.p)
dev.off()

# extract all the mcmc chains
out.p<-rbindlist(lapply(zm.p, as.data.frame))

# summarize all the data
stat.p<-summary(zm.p)$quantiles

## this is a little easier to extract data for plotting the predictions
# extract information for plotting
jm<-jags.model(modelBroad, data = data.list,  n.chains = setts$n.chains, n.adapt = 200000)
update(jm, n.iter=100000) #remove the burn-in
zp.j<-jags.samples(jm, variable.names = c('y.hat', 'y.new',"pvalue.mean","pvalue.sd","pvalue.cv"), n.iter =100000, n.thin = 1)


## Extract parameters for plotting
#calculate quartiles
q<-c(0.025,0.25,0.5,0.75, 0.975)

# predictions
y.hat.p<-summary(zp.j$y.hat, quantile, q)$stat 
y.new.p<-summary(zp.j$y.new, quantile, q)$stat 

# calculate the quantiles and means f p-values
params<-c('pvalue.mean', 'pvalue.sd', 'pvalue.cv')
p.values<- zp.j %>% 
  map_at(.at = params,function(x){summary(x, mean)$stat})
 
p.values<-cbind(p.values[[1]],p.values[[2]],p.values[[3]])
colnames(p.values)<-params

## predicitive checks
##autocorrelation plot
png('Output/traceplots/AutocorrelationPhoto.png', 3000,1500, res = 150)
acfplot(zm.p)
dev.off()

# save the Rhat statistic
Rhat.p<-gelman.diag(zm.p, multivariate = FALSE)

#DIC.p<-dic.samples(jm, n.iter, thin = thin)
#Neff
Neff.p<-effectiveSize(zm.p)

###### plot the photosynthesis curves #####
pdf('Output/MSFigures/Figure3.pdf', 5,6)

# Photosynthesis
BBData<-AllData
BBData<-BBData[BBData$Species=='OF' & BBData$PR=='Photosynthesis',]
#make the rates log x+1
BBData$Rate.ln<-log(BBData$Rate+1)

mycol1 <- rgb(t(col2rgb('skyblue')), max = 255, alpha = 125, names = "blue50")# transparent blue
mycol2 <- rgb(255,99,71, max = 255, alpha = 125, names = "tomato")# transparent orange
my.cols<-list(mycol1,mycol2)
colors<-c('skyblue','tomato')


par(mfrow=c(1,1))
par(mar=c(4.1,5.1,4.1,2.1)) # add more space on the left

plot(BBData$Temp.C[BBData$Location=='Bermuda'],BBData$Rate.ln[BBData$Location=='Bermuda'],
     ylim = c(-0.5,2),cex.lab = 1.5, cex.axis = 1.5, xlab ='Temperature'~degree~C,
     ylab = expression(paste("Gross Photosynthesis (log ", mu,"mol cm"^{-2}, "hr"^{-1},")")), pch=19, col = 'skyblue', xlim = c(20,40))
points(BBData$Temp.C[BBData$Location=='Panama'],BBData$Rate.ln[BBData$Location=='Panama'], xlab="Temperature"~degree~C, ylab = "Log Rate", ylim = c(-0.3,1.5),pch=19, col = 'tomato')

# fill in lines where there is data
predlow<-which(data.list$x.hat<=273.15+24)
predhigh<-which(data.list$x.hat>=273.15+36)
pred<-which(data.list$x.hat>=273.15+24 & data.list$x.hat<=273.15+36)

# panama went up to 37 deg
predhigh.p<-which(data.list$x.hat>=273.15+36)
pred.p<-which(data.list$x.hat>=273.15+24 & data.list$x.hat<=273.15+37)

pred<-list(pred,pred.p)
# make a polygon

for (i in 1:2){
  polygon(c(rev(data.list$x.hat[pred[[i]]]-273.15), data.list$x.hat[pred[[i]]]-273.15), 
          c(rev(y.hat.p[5,pred[[i]],i]), y.hat.p[1,pred[[i]],i]), col = my.cols[[i]], border = NA)
  lines(data.list$x.hat[pred[[i]]]-273.15,y.hat.p[3,pred[[i]],i], col=colors[i], lwd = 2)
  # add predictions for past where we have data
  for (j in c(1,3,5)){
    lines(data.list$x.hat[-predlow]-273.15,y.hat.p[j,-predlow,i], col=colors[i], lwd = 2, lty=2)
    lines(data.list$x.hat[-predhigh]-273.15,y.hat.p[j,-predhigh,i], col=colors[i], lwd = 2, lty=2)
  }
}
legend('bottomleft', legend = c('Bermuda','Panama'), col = c('skyblue','tomato'), lty=1, bty='n', cex = 1.5)
dev.off()

###rate comparison ######
### make plot of each parameter for both photosynthesis and respiration ###
# extract the values for P and r
E<-rbind(stat.p[1:2,])
Eh<-rbind(stat.p[3:4,])
Fd<-rbind(stat.p[5:6,])-273.15 # put on celcius
Topt<-rbind(stat.p[7:8,])-273.15 # put on celcius
lnc<-rbind(stat.p[9:10,])

# make a list of all the quantiles for easier looping
estimates<-list(E,Eh,Topt,lnc,Fd)
name<-c('E','Eh','Topt',expression(paste('b(T'[c],')')),'CTmax')
## plot the estimates

pdf('Output/MSFigures/FigureS7.pdf', 5,8)
#x labels for the graph
xlabs<-c(c('Energy', 'Energy'),'Temperature'~degree~C,expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})),'Temperature'~degree~C)
# make the plot
par(mfrow=c(3,2))
for (i in 1:length(estimates)){
  plot(estimates[[i]][,3], c(1,2), xlim = c(min(estimates[[i]]), max(estimates[[i]])),
       ylim = c(0.8,2.4),  xlab = xlabs[i], ylab = "", yaxt = 'n', 
       main = name[i], 
       #col = c('blue','lightblue','red','pink'), 
       col = c('skyblue','tomato','skyblue','tomato'),
       pch = c(19),
       cex = 2, cex.lab = 1.5, cex.main = 2
  )
  #error bars
  segments(estimates[[i]][,1], c(1,2),estimates[[i]][,5], c(1,2), col = 'grey') # 95% CI
  segments(estimates[[i]][,2], c(1,2),estimates[[i]][,4], c(1,2), lwd=3) # 50% CI
  if((i %% 2) > 0){ # only add the y labels on the left side
    axis(side=2,at = c(1,2), labels = c('Bermuda','Panama'), cex.axis = 1.5)
  }
}
plot(NA, ylim = c(0,2), xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", bty = "n")
dev.off()


##########Respiration####################
BBData<-AllData
BBData<-BBData[BBData$Species=='OF' & BBData$PR=='Respiration',]

unq_ind_group <- BBData[,10:11]
unq_ind_group <- unq_ind_group[!duplicated(unq_ind_group),] # this is needed for the heirarchical strucutre it creates a vector of the unique IDs

#make the rates log x+1
BBData$Rate.ln<-log(BBData$Rate+1) # might need to change back to 1

# data list for the model
data.list<-list(
  y=BBData$Rate.ln, # log O2 rate
  x.hat = seq(14+273.15, 40+273.15, 0.1), # for predictions, 14 to 40 degrees C by 0.1 deg C
  W=diag(4), # for scaled wishart distribution
  N = length(BBData$Rate.ln), # first level of hierarchy (all observations)
  temp=BBData$Temp.C+273.15, # temperature data
  L=length(unique(BBData$Location)), # second level of hierarchy (locations)
  loc=as.factor(BBData$Location),
  G = length(unique(BBData$unique.ID)), # third level of hierarchy (genotypes)
  geno  = as.factor(BBData$unique.ID),
  ind_per_group = as.factor(rep(c('Bermuda','Panama'),4))) # vector of individual IDs


y=BBData$Rate.ln # log O2 rate for shiny stan

## MCMC settings
setts <- list('n.iter' = 250000, 'n.thin' = 200, 'n.burn' = 200000)
setts.m <- 10  # this is 1.5M iterations with a 1.5K thinning and 1.2M burn-in

mSetts <- 1
if(mSetts) setts <- lapply(setts, function(v) v * setts.m)
setts$n.chains <- 3 # number of chains
setts$n.adapt <- 500000 # adaptation phase
n.cores <- 3 # number of cores to use in parallel computing
#set.seed(129)

cl <- makePSOCKcluster(n.cores)
tmp <- clusterEvalQ(cl, library(dclone))
parLoadModule(cl, 'glm')
parLoadModule(cl, 'lecuyer')
parLoadModule(cl, 'dic')
zm.r <- jags.parfit(cl = cl, data = data.list, params = variables.zp, model = modelBroad, 
                    n.chains = setts$n.chains, 
                    n.adapt = setts$n.adapt, 
                    n.update = setts$n.burn,
                    n.iter = setts$n.iter, 
                    thin = setts$n.thin)

stopCluster(cl)
# plot the trace plots
png('Output/traceplots/traceresp%03d.png', 1500,1500, res = 150)
plot(zm.r)
dev.off()

#extract the quantiles
stat.r<-summary(zm.r)$quantiles

## this is a little easier to extract data for plotting the predictions
# extract information for plotting
jm<-jags.model(modelBroad, data = data.list,  n.chains = setts$n.chains, n.adapt = 200000)
update(jm, n.iter=150000) #remove the burn-in
zp.jr<-jags.samples(jm, variable.names = c('y.hat', 'y.new',"pvalue.mean","pvalue.sd","pvalue.cv"), n.iter =150000, n.thin = 10)


## Extract parameters for plotting
#calculate quartiles
q<-c(0.025,0.25,0.5,0.75, 0.975)

# predictions
y.hat.r<-summary(zp.jr$y.hat, quantile, q)$stat 
y.new.r<-summary(zp.jr$y.new, quantile, q)$stat 

# calculate the quantiles and means f p-values
p.values.r<- zp.jr %>% 
  map_at(.at = params,function(x){summary(x, mean)$stat})

p.values.r<-cbind(p.values.r[[1]],p.values.r[[2]],p.values.r[[3]])
colnames(p.values.r)<-params

########### Plotting Results #######################
### respiration curves
pdf('Output/MSFigures/FigureS2.pdf',6,6)

par(mfrow=c(1,1))
par(mar=c(5.1,5.1,4.1,2.1)) # add more space on the left

plot(BBData$Temp.C[BBData$Location=='Bermuda'],BBData$Rate.ln[BBData$Location=='Bermuda'],
     ylim = c(-0.5,2),cex.lab = 1.5, cex.axis = 1.5, xlab = 'Temperature'~degree~C, 
     ylab = expression(paste("Respiration (log ", mu,"mol cm"^{-2}, "hr"^{-1},")")), pch=15, col = 'skyblue', xlim = c(20,40))
points(BBData$Temp.C[BBData$Location=='Panama'],BBData$Rate.ln[BBData$Location=='Panama'],  ylim = c(-0.3,1.5),pch=15, col = 'tomato')

# make a polygon for the predictions
for (i in 1:2){
  polygon(c(rev(data.list$x.hat[pred[[i]]]-273.15), data.list$x.hat[pred[[i]]]-273.15), 
          c(rev(y.hat.r[5,pred[[i]],i]), y.hat.r[1,pred[[i]],i]), col = my.cols[[i]], border = NA)
  lines(data.list$x.hat[pred[[i]]]-273.15,y.hat.r[3,pred[[i]],i], col=colors[i], lwd = 2)
  # add predictions for past where we have data
  for (j in c(1,3,5)){
    lines(data.list$x.hat[-predlow]-273.15,y.hat.r[j,-predlow,i], col=colors[i], lwd = 2, lty=2)
    lines(data.list$x.hat[-predhigh]-273.15,y.hat.r[j,-predhigh,i], col=colors[i], lwd = 2, lty=2)
  }
}
# add a legend
legend('topleft', c('Bermuda','Panama'), col = c('skyblue','tomato'), lty=1, bty ='n', cex = 1.5)
dev.off()

##autocorrelation plot
png('Output/traceplots/AutocorrelationResp.png', 3000,1500, res = 150)
acfplot(zm.r)
dev.off()

# save the Rhat statistic
Rhat.r<-gelman.diag(zm.r, multivariate = FALSE)
#Neff
Neff.r<-effectiveSize(zm.r)


### make plot of each parameter for both photosynthesis and respiration ###
# extract the values for P and r
E<-rbind(stat.r[1:2,])
Eh<-rbind(stat.r[3:4,])
Fd<-rbind(stat.r[5:6,])-273.15 # put on celcius
Topt<-rbind(stat.r[7:8,])-273.15 # put on celcius
lnc<-rbind(stat.r[9:10,])

# make a list of all the quantiles for easier looping
estimates<-list(E,Eh,Topt,lnc,Fd)

name<-c('E','Eh','Topt',expression(paste('b(T'[c],')')),'CTmax')

#sigmas
sigma.obs<-rbind(stat.p[which(row.names(stat.p)=='sigma.obs'),],stat.r[which(row.names(stat.r)=='sigma.obs'),])
sigma.geno<-rbind(stat.p[which(row.names(stat.p)=='sigma.geno'),],stat.r[which(row.names(stat.r)=='sigma.geno'),])
sigma<-rbind(stat.p[which(row.names(stat.p)=='sigma[1]'):which(row.names(stat.p)=='sigma[4]'),],stat.r[which(row.names(stat.r)=='sigma[1]'):which(row.names(stat.r)=='sigma[4]'),])
# only pull the lnc Sigma
sigma.lnc<-sigma[c(1,5),]
#make a list to make it easier to plot
sigmas<-list(sigma.obs, sigma.geno, sigma.lnc)

pdf('Output/MSFigures/FigureS4.pdf', 5,8)
#x labels for the graph
xlabs<-c(c('Energy', 'Energy'),'Temperature'~degree~C,expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})),'Temperature'~degree~C)
# make the plot
par(mfrow=c(3,2))
for (i in 1:length(estimates)){
  plot(estimates[[i]][,3], c(1,2), xlim = c(min(estimates[[i]]), max(estimates[[i]])),
       ylim = c(0.8,2.4),  xlab = xlabs[i], ylab = "", yaxt = 'n', 
       main = name[i], 
       #col = c('blue','lightblue','red','pink'), 
       col = c('skyblue','tomato','skyblue','tomato'),
       pch = c(19),
       cex = 2, cex.lab = 1.5, cex.main = 2
  )
  #error bars
  segments(estimates[[i]][,1], c(1,2),estimates[[i]][,5], c(1,2), col = 'grey') # 95% CI
  segments(estimates[[i]][,2], c(1,2),estimates[[i]][,4], c(1,2), lwd=3) # 50% CI
  if((i %% 2) > 0){ # only add the y labels on the left side
    axis(side=2,at = c(1,2), labels = c('Bermuda','Panama'), cex.axis = 1.5)
  }
}
plot(NA, ylim = c(0,2), xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", bty = "n")
dev.off()


# 
## Sigmas
png('Output/MSFigures/Figures9a.png', 1500,1500, res = 250)
par(mfrow=c(1,1))
par(mar=c(5.1,12.1,4.1,2.1))
plot(c(sigma.obs[,3], sigma.geno[,3], sigma.lnc[,3]), c(1,1.2,2,2.2,3,3.3), 
     #col = c('blue','red'), 
     cex.lab = 1.5,
     cex.axis = 1.5,
     col = 'black',
     cex = 2, pch = c(19,15), xlim = c(min(data.frame(sigmas)), max(data.frame(sigmas))),
     ylim = c(0.8,3.4), yaxt='n', ylab = "", xlab = "Variance")
segments(c(sigma.obs[,1], sigma.geno[,1], sigma.lnc[,1]), c(1,1.2,2,2.2,3,3.3),c(sigma.obs[,5], sigma.geno[,5], sigma.lnc[,5]), c(1,1.2,2,2.2,3,3.3), col = 'grey') # 95% CI
segments(c(sigma.obs[,2], sigma.geno[,2], sigma.lnc[,2]), c(1,1.2,2,2.2,3,3.3),c(sigma.obs[,4], sigma.geno[,4], sigma.lnc[,4]), c(1,1.2,2,2.2,3,3.3), lwd=3) # 50% CI
#axis(side=2,at = c(1,2,3), labels = c('Observation error','Between genotypes w/i locations','Between locations' ), cex.axis = 1.5, srt = 90)
text(x = par("usr")[1], y = c(1,2,3), labels = c('Observation error','Between genotypes\n
                                                     within locations','Between locations' ), srt = 0, pos = 2, xpd = TRUE, cex = 1.5)
legend('bottomright',c('Gross Photosynthesis','Respiration'), pch = c(19,15), col = 'black', bty='n', cex=1)
dev.off()

# observed vs predicted plots
png('Output/MSFigures/FigureS5.png', 2000,1500, res = 250)

par(mfrow=c(1,2))
y2<-AllData
y2<-y2[y2$Species=='OF' & y2$PR=='Photosynthesis',]
y2$Rate.ln<-log(y2$Rate+1)

# Photosynthesis
plot(y2$Rate.ln,y.new.p[3,], xlab = 'Observed', ylab = 'Predicted', main = 'Photosynthesis')
abline(0,1)

# respiration
plot(y,y.new.r[3,], xlab = 'Observed', ylab = 'Predicted', main = 'Respiration')
abline(0,1)
dev.off()

### make tables of the results
# make a data table of all the checks
diagnostic<-cbind(as.data.frame(Rhat.p[[1]]),as.data.frame(Neff.p),as.data.frame(Rhat.r[[1]]), as.data.frame(Neff.r) )
write.csv(diagnostic, file = "Output/diagnostic.csv")

#confidence intervals of the params
write.csv(stat.p, file = "Output/PhotoResults.csv")
write.csv(stat.r, file = "Output/RespResults.csv")

# p-values
ps<-rbind(p.values,p.values.r)
rownames(ps)<-c('Gross Photosynthesis', 'Respiration')
write.csv(ps,'Output/PValues_locations.csv')


### make some comparison plots using tidybayes ####
# compare groups in tidy-bayes for Photosynthesis
test<-zm.p  %>%
  spread_draws(E.loc[loc], Eh.loc[loc], lnc.loc[loc], Topt[loc], Tmax[loc]) 

test$loc<-ifelse(test$loc==1, 'Ber', 'Pan')

p1<-test %>%  compare_levels(E.loc, by = loc) %>%
  ggplot(aes(x = E.loc, y = loc)) + 
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('E')+
  xlab('Energy')+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_blank())  

# get the means and quantiles of the different comparisons
E.quant<-test %>%  compare_levels(E.loc, by = loc) %>%
  group_by(loc)%>%
  summarise(E.loc = list(enframe(quantile(E.loc, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(E.loc=value)


p2<-test%>%  compare_levels(Eh.loc, by = loc) %>%
  ggplot(aes(x = Eh.loc, y = loc)) + 
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('Eh')+
  xlab('Energy')+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_blank())  
# get the means and quantiles of the different comparisons
Eh.quant<-test %>%  compare_levels(Eh.loc, by = loc) %>%
  group_by(loc)%>%
  summarise(Eh.loc = list(enframe(quantile(Eh.loc, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(Eh.loc=value)


p3<-test%>%  compare_levels(lnc.loc, by = loc) %>%
  ggplot(aes(x = lnc.loc, y = loc)) + 
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('b(Tc)')+
  xlab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_blank()) 

lnc.quant<-test %>%  compare_levels(lnc.loc, by = loc) %>%
  group_by(loc)%>%
  summarise(lnc.loc = list(enframe(quantile(lnc.loc, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(lnc.loc=value)

p4<-test%>%  compare_levels(Topt, by = loc) %>%
  ggplot(aes(x = Topt, y = loc)) +
  #scale_y_continuous(limits=c("Panama - Bermuda"), breaks = 1)+
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('Topt')+
  xlab("Temperature"~degree~C)+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_blank()) 

Topt.quant<-test %>%  compare_levels(Topt, by = loc) %>%
  group_by(loc)%>%
  summarise(Topt = list(enframe(quantile(Topt, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(Topt=value)

p5<-test%>%  compare_levels(Tmax, by = loc) %>%
  ggplot(aes(x = Tmax, y = loc)) +
  #scale_y_continuous(limits=c("Panama - Bermuda"), breaks = 1)+
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('CTmax')+
  xlab("Temperature"~degree~C)+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_blank()) 
 
Tmax.quant<-test %>%  compare_levels(Tmax, by = loc) %>%
  group_by(loc)%>%
  summarise(Tmax = list(enframe(quantile(Tmax, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(Tmax=value)

# put them all in one table
Photoquants<-left_join(E.quant,Eh.quant)%>%
  left_join(.,lnc.quant) %>%
    left_join(.,Topt.quant) %>%
  left_join(.,Tmax.quant)

#put all the photosynthesis in one plots
allP<-grid.arrange(p1,p2,p3,p4,p5, ncol = 2,  top = textGrob("Gross Photosynthesis", gp=gpar(fontsize=15,font=8) ),
                  left = textGrob("Panama - Bermuda", gp=gpar(fontsize=15,font=8, face = 'bold'), rot=90 ))
ggsave(filename = 'Output/MSFigures/Figure4.png', plot = allP, device = 'png', width = 6, height = 7, dpi = 300)
ggsave(filename = 'Output/MSFigures/Figure4.pdf', plot = allP, device = 'pdf', width = 6, height = 7, dpi = 300)


# now respiration
test1<-zm.r  %>%
  spread_draws(E.loc[loc], Eh.loc[loc], lnc.loc[loc], Topt[loc], Tmax[loc]) 
test1$loc<-ifelse(test1$loc==1, 'Ber', 'Pan')


p6<-test1 %>%  compare_levels(E.loc, by = loc) %>%
  ggplot(aes(x = E.loc, y = loc)) + 
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('E')+
  xlab('Energy')+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_blank())  

p7<-test1%>%  compare_levels(Eh.loc, by = loc) %>%
  ggplot(aes(x = Eh.loc, y = loc)) + 
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('Eh')+
  xlab('Energy')+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_blank()) 

p8<-test1%>%  compare_levels(lnc.loc, by = loc) %>%
  ggplot(aes(x = lnc.loc, y = loc)) + 
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('b(Tc)')+
  xlab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_blank()) 

p9<-test1%>%  compare_levels(Topt, by = loc) %>%
  ggplot(aes(x = Topt, y = loc)) +
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('Topt')+
  xlab("Temperature"~degree~C)+
  theme_light()+
  xlim(-3,3)+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_blank())  

p10<-test1%>%  compare_levels(Tmax, by = loc) %>%
  ggplot(aes(x = Tmax, y = loc)) +
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('CTmax')+
  xlab("Temperature"~degree~C)+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_blank()) 
  


# get the means and quantiles of the different comparisons
E.quant.r<-test1 %>%  compare_levels(E.loc, by = loc) %>%
  group_by(loc)%>%
  summarise(E.loc = list(enframe(quantile(E.loc, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(E.loc=value)

Eh.quant.r<-test1 %>%  compare_levels(Eh.loc, by = loc) %>%
  group_by(loc)%>%
  summarise(Eh.loc = list(enframe(quantile(Eh.loc, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(Eh.loc=value)

lnc.quant.r<-test1 %>%  compare_levels(lnc.loc, by = loc) %>%
  group_by(loc)%>%
  summarise(lnc.loc = list(enframe(quantile(lnc.loc, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(lnc.loc=value)

Topt.quant.r<-test1 %>%  compare_levels(Topt, by = loc) %>%
  group_by(loc)%>%
  summarise(Topt = list(enframe(quantile(Topt, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(Topt=value)

Tmax.quant.r<-test1 %>%  compare_levels(Tmax, by = loc) %>%
  group_by(loc)%>%
  summarise(Tmax = list(enframe(quantile(Tmax, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(Tmax=value)

# put them all in one table
Respquants<-left_join(E.quant.r,Eh.quant.r)%>%
  left_join(.,lnc.quant.r) %>%
  left_join(.,Topt.quant.r) %>%
  left_join(.,Tmax.quant.r)

# put all P and R together in one plot
AllP<-grid.arrange(arrangeGrob(p1, p2, p3,p4,p5, top=textGrob("Gross Photosynthesis", gp=gpar(fontsize=15,font=8, face = 'bold') ), ncol = 1), 
                   arrangeGrob(p6, p7,p8,p9,p10, top=textGrob("Respiration", gp=gpar(fontsize=15,font=8, face = 'bold') ), ncol = 1), 
                   ncol = 2,left = textGrob("Panama - Bermuda", gp=gpar(fontsize=15,font=8, face = 'bold'), rot=90 ))
# Just respiration in one plot
allR<-grid.arrange(p6, p7,p8,p9,p10, ncol = 2,  top = textGrob("Dark Respiration", gp=gpar(fontsize=15,font=8) ),
                   left = textGrob("Panama - Bermuda", gp=gpar(fontsize=15,font=8, face = 'bold'), rot=90 ))

#ggsave(filename = 'Output/LocationComparisonAll.png', plot = AllP, device = 'png', width = 10, height = 10, dpi = 300)
ggsave(filename = 'Output/MSFigures/FigureS3.png', plot = allR, device = 'png', width = 10, height = 10, dpi = 300)
ggsave(filename = 'Output/MSFigures/FigureS3.pdf', plot = allR, device = 'pdf', width = 6, height = 7, dpi = 300)

# export the results
write.csv(Photoquants, 'Output/PhotoQuantiles_locations.csv')
write.csv(Respquants, 'Output/RespQuantiles_locations.csv')

