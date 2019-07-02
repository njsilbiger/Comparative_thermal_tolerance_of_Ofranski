## R script to construct Bayesian models to compare photosynthesis, respiration, and calcification rates  to each other in Bermuda
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
# the temp categories in Panama are not all correct
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


#### Run a model comparing all three rates to each other in the same model

BBData<-AllData

BBData<-BBData %>% filter(Species=='OF' & Location=='Bermuda') %>%
  filter(PR=='Photosynthesis' | PR=='Respiration' | PR =='Net Calcification')  %>% 
  droplevels()

unq_ind_group <- BBData[,c(6,11)]
unq_ind_group <- unq_ind_group[!duplicated(unq_ind_group),] # this is needed for the heirarchical strucutre it creates a vector of the unique IDs

#make the rates log x+1
BBData$Rate.ln<-log(BBData$Rate+1)

#There is one coral that is net dissolving and driving the entire pattern.  I think it is a titration error
bad<-which(BBData$Rate.ln< -0.5)
BBData[bad,"Rate.ln"]<-NA


# order the PR levels
# make the ids unique to rate and location
BBData$unique.ID<-factor(paste0(BBData$unique.ID,BBData$PR))


# data list for the model
data.list<-list(
  y=BBData$Rate.ln, # log O2 rate
  x.hat = seq(14+273.15, 40+273.15, 0.1), # for predictions, 14 to 40 degrees C by 0.1 deg C
  W=diag(4), # for scaled wishart distribution
  N = length(BBData$Rate.ln), # first level of hierarchy (all observations)
  temp=BBData$Temp.Cat+273.15, # temperature data
  R=length(unique(BBData$PR)), # second level of hierarchy (rate type)
  rate=as.factor(BBData$PR),
  G = length(unique(BBData$unique.ID)), # third level of hierarchy (genotypes)
  geno  = as.factor(BBData$unique.ID),
  ind_per_group = as.factor(rep(c('Net Calcification','Photosynthesis','Respiration'),4))
  
  ) # vector of individual IDs


y=BBData$Rate.ln # log O2 rate for shiny stan


#### paramaterize a model that is  hierarchical #####
# JAGS model name for the broad model. All the JAGS code is in this file
modelBroad = "Scripts/JAGSnestedRates.R"

#varibles to save for inference
variables.zp<-c("lnc.rate","E.rate" ,"Eh.rate",'Topt', 'sigma.obs','sigma', "lnc.geno", 'sigma.geno', 'lnc.pred', 'Tmax')
keep<-c("pvalue.mean","pvalue.sd","pvalue.cv","E.rate" ,"Eh.rate",'Topt', 'y.hat',"y.new", "lnc.rate", 'Tmax','lnc.pred', 'sigma.obs','sigma','sigma.geno') # params to keep for estimates

## MCMC settings
setts <- list('n.iter' = 150000, 'n.thin' = 100, 'n.burn' = 180000)
setts.m <- 10  # this is 1.5M iterations with a 1.5K thinning and 1.2M burn-in

mSetts <- 1
if(mSetts) setts <- lapply(setts, function(v) v * setts.m)
setts$n.chains <- 3 # number of chains
setts$n.adapt <- 500000 # adaptation phase

n.cores <- 3 # number of cores to use in parallel computing
set.seed(129)

## other settings

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
png('../Output/traceAllRates%03d.png', 1500,1500, res = 150)
plot(zm.p)
dev.off()

# extract all the mcmc chains
out.p<-rbindlist(lapply(zm.p, as.data.frame))

# summarize all the data
stat.p<-summary(zm.p)$quantiles

## this is a little easier to extract data for plotting the predictions
# extract information for plotting
jm<-jags.model(modelBroad, data = data.list,  n.chains = setts$n.chains, n.adapt = 100000)
update(jm, n.iter=100000) #remove the burn-in
zp.j<-jags.samples(jm, variable.names = c('y.hat', 'y.new',"pvalue.mean","pvalue.sd","pvalue.cv"), n.iter =100000, n.thin = 5)


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

#### Plotting ####
pdf('Output/MSFigures/Figure5.pdf', 6,10, useDingbats = FALSE)


mycol1 <- rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50")# transparent blue
mycol2 <- rgb(144, 238, 144, max = 255, alpha = 125, names = "green50")# transparent green
mycol3<-rgb(255,192,203,  max = 255, alpha = 125, names = "pink50")
my.cols<-list(mycol1,mycol2,mycol3)

mycol1 <- rgb(t(col2rgb('skyblue')), max = 255, alpha = 125, names = "blue50")# transparent blue
pt<-c(24,21,22) # points for P, R, and C

par(mfrow=c(2,1))
par(mar=c(2.1,5.1,4.1,5.1))

plot(BBData$Temp.Cat,BBData$Rate.ln, xlab="", cex.lab = 1.5, cex.axis = 1.5,
     ylim = c(-0.5,2),ylab = expression(paste("Rate (log ", mu,"mol cm"^{-2}, "hr"^{-1},")")), pch=19, col = 'white', xlim = c(20,40))

PR<-unique(BBData$PR)[c(3,1,2)]

# fill in lines where there is data
predlow<-which(data.list$x.hat<=273.15+24)
predhigh<-which(data.list$x.hat>=273.15+36)
pred<-which(data.list$x.hat>=273.15+24 & data.list$x.hat<=273.15+36)
# make a polygon

for (i in 1:3){
  points(BBData$Temp.Cat[BBData$PR==PR[i]],BBData$Rate.ln[BBData$PR==PR[i]], bg = 'skyblue', pch = pt[i])
  polygon(c(rev(data.list$x.hat[pred]-273.15), data.list$x.hat[pred]-273.15), 
          c(rev(y.hat.p[3,pred,i]),y.hat.p[1,pred,i]), col = mycol1, border = NA)
  lines(data.list$x.hat[pred]-273.15,y.hat.p[2,pred,i], col='skyblue', lwd = 2)
  # add predictions for past where we have data
  for (j in 1:3){
    lines(data.list$x.hat[-predlow]-273.15,y.hat.p[j,-predlow,i], col='skyblue', lwd = 2, lty=2)
    lines(data.list$x.hat[-predhigh]-273.15,y.hat.p[j,-predhigh,i], col='skyblue', lwd = 2, lty=2)
  }
}
legend('topleft', legend = c('Net Calcification','Gross Photosynthesis','Respiration'), pt.bg = 'skyblue', pch = c(24,21,22), bty='n', cex = 1.5)


## P to R plots
## P/R for each location
photo<-exp(y.hat.p[,,2])-1 # back transporm the data
resp<-exp(y.hat.p[,,3])-1
calc<- exp(y.hat.p[,,1])-1 
calc[which(calc<0)] <-NA#(only take the positive values for calcification)
  
PR<-photo[3,]/resp[3,] #take the median value
CP<-calc[3,]/photo[3,]

par(mar=c(5.1,5.1,4.1,5.1))

plot(data.list$x.hat[pred]-273.15,PR[pred], cex.lab = 1.5,yaxt = 'n',cex.axis = 1.5, type = 'l', xlim = c(20,40),ylim = c(0,3.5), col = 'skyblue', ylab = "", xlab ='Temperature'~degree~C, lwd =2)
abline(h=1, lty=2)
axis(side = 2, cex.axis = 1.5, col = 'skyblue', col.axis = 'skyblue')
mtext(side = 2, line = 3, "Gross Photosynthesis:Respiration", cex = 1.5)

par(new = T)
plot(data.list$x.hat[pred]-273.15,CP[pred], axes=F, xlab=NA, ylab=NA, xlim = c(20,40), type = 'l',ylim = c(0,0.3), col = 'darkblue', lwd =2)
axis(side = 4, cex.axis = 1.5, col = 'darkblue', col.axis = 'darkblue')
mtext(side = 4, line = 3, 'Calcification:Gross Photosynthesis', cex = 1.5)
dev.off()

# pull out the names of the parameters that I want to plot

n<-rownames(stat.p)%in%c('E.rate','Eh.rate','Th.rate','lnc.rate','Topt')

estimates<-stat.p[n]
# create x lables
xlabs<-c(c('Energy', 'Energy'),'Temperature'~degree~C,'log rate','Temperature'~degree~C)
png('Output/MSFigures/FigureS8.png', 1500,2000, res = 250)
E<-stat.p[1:3,]
Eh<-stat.p[4:6,]
Fd<-stat.p[7:9,]-273.15 # put on celcius
Topt<-stat.p[10:12,]-273.15 # put on celcius
lnc<-stat.p[37:39,]

estimates<-list(E,Eh,Topt,lnc,Fd)

name<-c('E','Eh','Topt',expression(paste('b(T'[c],')')),'CTmax')


par(mfrow=c(3,2))
for (j in 1:length(estimates)){
  
  plot(estimates[[j]][,3], c(1:3), xlim = c(min(estimates[[j]]), max(estimates[[j]])),
       ylim = c(0.8,3.2), pch=pt, xlab = xlabs[j], ylab = "", yaxt = 'n', 
       cex.lab = 1.5,
       cex.axis = 1.5,
       col = 'skyblue',
       bg = 'skyblue',
       cex = 2,
       cex.main = 2,
       main = name[j])
  
  segments(estimates[[j]][,1], c(1:3),estimates[[j]][,5],c(1:3), col = 'grey')
  segments(estimates[[j]][,2], c(1:3),estimates[[j]][,4],c(1:3), lwd=3)
  if((j %% 2) > 0){ # only add the y labels on the left side
    axis(side=2,at = c(1:3), labels = c('NC','GP','R'))
  
  }
}
dev.off()

## plot the sigmas
png('Output/MSFigures/FigureS9b.png', 1500,1000, res = 250)

par(mfrow=c(1,1))
par(mar=c(5.1,12.1,4.1,2.1))
sigmas<-cbind(stat.p[45,], stat.p[44,], stat.p[40,])

plot(sigmas[3,1:3], c(1:3), xlim = c(min(sigmas[,1:3]), max(sigmas[,1:3])), xlab = 'Variance', 
     ylim = c(0.8,3.2), pch=19,  ylab = "", yaxt = 'n', cex.lab = 1.5,
     cex.axis = 1.5, cex = 2)
segments(sigmas[1,1:3], c(1:3), sigmas[5,1:3], c(1:3), col = 'grey')
segments(sigmas[2,1:3], c(1:3), sigmas[4,1:3], c(1:3), lwd = 3)
text(x = par("usr")[1], y = c(1,2,3), labels = c('Observation error','Among genotypes\n within functions','Among functions' ),
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

dev.off()


### predictive checks
##autocorrelation plot
png('Output/traceplots/AutocorrelationRates.png', 3000,1500, res = 150)
acfplot(zm.p)
dev.off()

# save the Rhat statistic
Rhat.rate<-gelman.diag(zm.p, multivariate = FALSE)

#Neff
Neff.rate<-effectiveSize(zm.p)

# observed vs predicted plots
png('Output/MSFigures/FiguresS6.png', 1500,1500, res = 250)
par(mfrow=c(1,1))
# All rates
plot(y,y.new.p[3,], xlab = 'Observed', ylab = 'Predicted')
abline(0,1)
dev.off()

### make tables of the results
# make a data table of all the checks
diagnosticR<-cbind(as.data.frame(Rhat.rate[[1]]),as.data.frame(Neff.rate) )
write.csv(diagnosticR, file = "../Output/diagnosticR.csv")

#confidence intervals of the params
write.csv(stat.p, file = "Output/AllRateComparisonResults.csv")

# p-values
write.csv(p.values,'Output/PValuesRateComparison.csv')

### make some comparison plots using tidybayes
# compare groups in tidy-bayes for bermuda
test<-zm.p  %>%
  spread_draws(E.rate[rate], Eh.rate[rate], lnc.rate[rate], Topt[rate], Tmax[rate]) 

#Give the numbers names for comparison
test$rate[which(test$rate==1)]<-'NC'
test$rate[which(test$rate==2)]<-'GP'
test$rate[which(test$rate==3)]<-'R'


p1<-test %>%  compare_levels(E.rate, by = rate) %>%
  ggplot(aes(x = E.rate, y = rate)) + 
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('E')+
  xlab('Energy')+
   theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14)) 

p2<-test %>%  compare_levels(Eh.rate, by = rate) %>%
  ggplot(aes(x = Eh.rate, y = rate)) + 
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('Eh')+
  xlab('Energy')+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14)) 

p3<-test %>%  compare_levels(lnc.rate, by = rate) %>%
  ggplot(aes(x = lnc.rate, y = rate)) + 
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('b(Tc)')+
  xlab(expression(paste("log ", mu,"mol cm"^{-2}, "hr"^{-1})))+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14)) 


p4<-test %>%  compare_levels(Topt, by = rate) %>%
  ggplot(aes(x = Topt, y = rate)) + 
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  xlab("Temperature"~degree~C)+
  ggtitle('Topt')+
 theme_light()+
  xlim(-4,4)+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14)) 


p5<-test %>%  compare_levels(Tmax, by = rate) %>%
  ggplot(aes(x = Tmax, y = rate)) + 
  geom_halfeyeh() +
  geom_vline(xintercept =0, linetype = 2)+
  ylab('')+
  ggtitle('CTmax')+
  xlab("Temperature"~degree~C)+
  theme_light()+
  theme(axis.text= element_text(size = 14),
        plot.title = element_text(size=14, face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 14)) 

allP<-grid.arrange(p1,p2,p3,p4,p5, ncol = 2)
ggsave(filename = 'Output/MSFigures/Figure6.png', plot = allP, device = 'png', width = 8, height = 10, dpi = 300)
ggsave(filename = 'Output/MSFigures/Figure6.pdf', plot = allP, device = 'pdf', width = 8, height = 10, dpi = 300)


# get the means and quantiles of the different comparisons
E.quant.rate<-test %>%  compare_levels(E.rate, by = rate) %>%
  group_by(rate)%>%
  summarise(E.rate = list(enframe(quantile(E.rate, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(E.rate=value)

Eh.quant.rate<-test %>%  compare_levels(Eh.rate, by = rate) %>%
  group_by(rate)%>%
  summarise(Eh.rate = list(enframe(quantile(Eh.rate, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(Eh.rate=value)

lnc.quant.rate<-test %>%  compare_levels(lnc.rate, by = rate) %>%
  group_by(rate)%>%
  summarise(lnc.rate = list(enframe(quantile(lnc.rate, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(lnc.rate=value)

Topt.quant.rate<-test %>%  compare_levels(Topt, by = rate) %>%
  group_by(rate)%>%
  summarise(Topt = list(enframe(quantile(Topt, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(Topt=value)

Tmax.quant.rate<-test %>%  compare_levels(Tmax, by = rate) %>%
  group_by(rate)%>%
  summarise(Tmax = list(enframe(quantile(Tmax, probs=c(0.025,0.5,0.975))))) %>% 
  unnest %>%
  rename(Tmax=value)

# put them all in one table
Ratequants<-left_join(E.quant.rate,Eh.quant.rate)%>%
  left_join(.,lnc.quant.rate) %>%
  left_join(.,Topt.quant.rate) %>%
  left_join(.,Tmax.quant.rate)

write.csv(Ratequants, 'Output/RateComparisonQuantiles.csv')
