#Title: Photosynthesis Irradiance Curves for Panama and Bermuda
#Author: HM Putnam NJ Silbiger
#Edited by: HM Putnam
#Date Last Modified: 20181223


rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("phytotools" %in% rownames(installed.packages()) == 'FALSE') install.packages('phytotools') 

#Read in required libraries
##### Include Versions of libraries
#install_github('colin-olito/LoLinR')
library("devtools")
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library('phytotools')


##### BERMUDA #####
path.p<-"Data/Bermuda/PI_Curve/Run_PI_curve/" #the location of all your respirometry files 

#bring in the files
file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=4)) #generate a 3 column dataframe with specific column names
colnames(Photo.R) <- c("Fragment.ID","Intercept", "µmol.L.sec", "Temp")

#Load Sample meta info Info
Sample.Info <- read.csv(file="Data/Bermuda/PI_Curve/Nubbin_Sample_Info_PI_Curve_Bermuda_QC.csv", header=T) #read sample.info data

#subset the data by light step using time breaks in the data
for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data1 <-read.table(file.path(path.p,file.names[i]),  header=T, sep=",", skip=1, na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data1  <- as.data.frame(cbind(Photo.Data1$Time, Photo.Data1$Value, Photo.Data1$Temp)) #subset columns of interest
  colnames(Photo.Data1) <- c("Time", "Value", "Temp")
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S") #convert time from character to time
  brk <- unique(as.POSIXct(Sample.Info$Start.time,format="%H:%M"))
  Step1 <- subset(Photo.Data1, Photo.Data1$Time > brk[1] & Photo.Data1$Time < brk[2])
  Step2 <- subset(Photo.Data1, Photo.Data1$Time > brk[2] & Photo.Data1$Time < brk[3])
  Step3 <- subset(Photo.Data1, Photo.Data1$Time > brk[3] & Photo.Data1$Time < brk[4])
  Step4 <- subset(Photo.Data1, Photo.Data1$Time > brk[4] & Photo.Data1$Time < brk[5])
  Step5 <- subset(Photo.Data1, Photo.Data1$Time > brk[5] & Photo.Data1$Time < brk[6])
  Step6 <- subset(Photo.Data1, Photo.Data1$Time > brk[6] & Photo.Data1$Time < brk[7])
  Step7 <- subset(Photo.Data1, Photo.Data1$Time > brk[7] & Photo.Data1$Time < brk[8])
  Step8 <- subset(Photo.Data1, Photo.Data1$Time > brk[8] & Photo.Data1$Time < brk[9])
  Step9 <- subset(Photo.Data1, Photo.Data1$Time > brk[9])
  lt.levs <- list(Step1,Step2,Step3,Step4,Step5,Step6,Step7,Step8,Step9) #list levels of segmentation
  
  for(j in 1:length(lt.levs)){    
    Photo.Data <- as.data.frame(lt.levs[j])
    n<-dim(Photo.Data)[1] #identify length of data
    Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
    n<-dim(Photo.Data)[1] #list length of trimmed data
    Photo.Data$sec <- as.numeric(1:n) #set seconds by one from start to finish of run
    Photo.Data$Value <- as.numeric(as.character(Photo.Data$Value)) #save O2 data as numeric
    Photo.Data$Temp <- as.numeric(as.character(Photo.Data$Temp)) #save O2 data as numeric
    
    #Save plot prior to and after data thinning to make sure thinning is not too extreme
    rename <- sub("_.*", "", file.names[i])
    pdf(paste0("Output/Bermuda/PIC_Output/",rename,"_",j,"thinning.pdf"))
    par(omi=rep(0.3, 4)) #set size of the outer margins in inches
    par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'),  axes=FALSE) #plot data as a function of time
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    points(Photo.Data$Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    
    #set data reduction thinning parameter and data usage parameter
    thin <- 20
    alpha <- 0.2
    meth <- "pc"
    #save original unthinned data
    Photo.Data.orig<-Photo.Data
    
    Photo.Data   <-  thinData(Photo.Data.orig , by=thin)$newData1 #thin data 
    Photo.Data $sec <- as.numeric(rownames(Photo.Data)) #maintain numeric values for time
    Photo.Data$Temp <- NA # add a new column to fill with the thinned data
    Photo.Data$Temp <-  thinData(Photo.Data.orig, xy=c(2,3),by=thin)$newData1[,2] #thin data  for the temp values
    plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
    axis(1)
    axis(2, las=1)
    dev.off()
    
    Regs  <-  rankLocReg(xall=Photo.Data$sec, yall=Photo.Data$Value, alpha=alpha, 
                         method=meth, verbose=TRUE) 
    pdf(paste0("Output/Bermuda/PIC_Output/",rename,"_",j,"regression.pdf"))
    plot(Regs)
    dev.off()
    
    s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
    Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
    Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
    Photo.R[j+s[i],4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
  }
}

#view rates
Photo.R
#write.csv(Photo.R,"Data/PI_Curve/Bermuda_PI_Curve_rates.csv")

#Merge rates with sample info
Data <- merge(Photo.R, Sample.Info, by="Fragment.ID")

#correct for the size of the chamber
Data$micromol.s <- Data$µmol.L.sec * Data$Chamber.Vol.L

#calculate the average of the blanks from each time step
blnks <- subset(Data, Sample.Type=="Blank")
blnks <-mean.blnks <- aggregate(micromol.s ~ Light_Level, data=blnks, FUN=mean)
Data <- merge(Data, blnks, by="Light_Level")
colnames(Data)[colnames(Data) == 'micromol.s.x'] <- 'sample.micromol.s'
colnames(Data)[colnames(Data) == 'micromol.s.y'] <- 'blank.micromol.s'

#subtract the average of the blanks from each time step
Data$corr.micromol.s <- Data$sample.micromol.s - Data$blank.micromol.s
Data$micromol.cm2.s <- Data$corr.micromol.s/Data$Surf.Area.cm2
Data$micromol.cm2.h <- Data$micromol.cm2.s*3600

OF1 <- subset(Data, Fragment.Number==6)
OF2 <- subset(Data, Fragment.Number==5)
OF <- subset(Data, Species=="OF")
write.csv(OF,"Data/Bermuda/PI_Curve/Bermuda_PI_Curve_rates_OF.csv")


##### PANAMA #####
path.p<-"Data/Panama/PI_Curve/QC" #the location of all your respirometry files 

#bring in the files
file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
Photo.R <- data.frame(matrix(NA, nrow=length(file.names), ncol=4)) #generate a 3 column dataframe with specific column names
colnames(Photo.R) <- c("Fragment.ID","Intercept", "µmol.L.sec", "Temp")

#Load Sample meta info Info
Sample.Info <- read.csv(file="Data/Panama/PI_Curve/Nubbin_Sample_Info_PI_Curve_Panama_QC.csv", header=T) #read sample.info data

#levs <- unique(Sample.Info$Fragment.ID)

#subset the data by light step using time breaks in the data
for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data <-read.table(file.path(path.p,file.names[i]),  header=T, sep=",", skip=1, na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data  <- as.data.frame(cbind(Photo.Data$Time, Photo.Data$Value,Photo.Data$Temp)) #subset columns of interest
  colnames(Photo.Data) <- c("Time", "Value", "Temp")
  Photo.Data$Time <- as.POSIXct(Photo.Data$Time,format="%H:%M:%S") #convert time from character to time
  
  n<-dim(Photo.Data)[1] #identify length of data
  Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
  n<-dim(Photo.Data)[1] #list length of trimmed data
  Photo.Data$sec <- as.numeric(1:n) #set seconds by one from start to finish of run
  Photo.Data$Value <- as.numeric(as.character(Photo.Data$Value)) #save O2 data as numeric
  Photo.Data$Temp <- as.numeric(as.character(Photo.Data$Temp)) #save O2 data as numeric
  
  #Save plot prior to and after data thinning to make sure thinning is not too extreme
  rename <- sub(".csv", "", file.names[i])
  pdf(paste0("Output/Panama/PIC_Output/",rename,"_thinning.pdf"))
  par(omi=rep(0.3, 4)) #set size of the outer margins in inches
  par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
  plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'),  axes=FALSE) #plot data as a function of time
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data$Value ~ Photo.Data$sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  
  #set data reduction thinning parameter and data usage parameter
  thin <- 20
  alpha <- 0.2
  meth <- "pc"
  #save original unthinned data
  Photo.Data.orig<-Photo.Data
  
  Photo.Data   <-  thinData(Photo.Data.orig , by=thin)$newData1 #thin data 
  Photo.Data $sec <- as.numeric(rownames(Photo.Data)) #maintain numeric values for time
  Photo.Data$Temp <- NA # add a new column to fill with the thinned data
  Photo.Data$Temp <-  thinData(Photo.Data.orig, xy=c(2,3),by=thin)$newData1[,2] #thin data by every 20 points for the temp values
  plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  dev.off()
  
  Regs  <-  rankLocReg(xall=Photo.Data$sec, yall=Photo.Data$Value, alpha=alpha, 
                       method=meth, verbose=TRUE) 
  pdf(paste0("Output/Panama/PIC_Output/",rename,"_regression.pdf"))
  plot(Regs)
  dev.off()
  
  Photo.R[i,2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
  Photo.R[i,1] <- rename #stores the file name in the Date column
  Photo.R[i,4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
}


#view rates
#View(Photo.R)
#write.csv(Photo.R,"Data/Panama/PI_Curve/Panama_PI_Curve_rates.csv")

#Merge rates with sample info
Data <- merge(Photo.R, Sample.Info, by="Fragment.ID")

#correct for the size of the chamber
Data$micromol.s <- Data$µmol.L.sec * Data$Chamber.Vol.L

#calculate the average of the blanks from each time step
blnks <- subset(Data, Sample.Type=="Blank")
blnks <-mean.blnks <- aggregate(micromol.s ~ Light_Level, data=blnks, FUN=mean)
Data <- merge(Data, blnks, by="Light_Level")
colnames(Data)[colnames(Data) == 'micromol.s.x'] <- 'sample.micromol.s'
colnames(Data)[colnames(Data) == 'micromol.s.y'] <- 'blank.micromol.s'

#subtract the average of the blanks from each time step
Data$corr.micromol.s <- Data$sample.micromol.s - Data$blank.micromol.s
Data <- subset(Data, Species!="Blank")
Data$micromol.cm2.s <- Data$corr.micromol.s/Data$Surf.Area.cm2
Data$micromol.cm2.h <- Data$micromol.cm2.s*3600

Data <- Data[with(Data, order(Fragment.Number, Light_Level)), ]

OF1 <- subset(Data, Fragment.Number=="OF1")
OF3 <- subset(Data, Fragment.Number=="OF3")
OF <- subset(Data, Species=="OF")
write.csv(OF,"Data/Panama/PI_Curve/Panama_PI_Curve_rates_OF.csv")


##### PLOTTING CURVES #####
##### Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)

#Plot curves
Ber.data <- read.table("Data/Bermuda/PI_Curve/Bermuda_PI_Curve_rates_OF.csv", header=TRUE, sep=",")
Pan.data <- read.table("Data/Panama/PI_Curve/Panama_PI_Curve_rates_OF.csv", header=TRUE, sep=",")

Ber <- aggregate(micromol.cm2.h ~ Light_Level, data = Ber.data, FUN=mean)
Pan <- aggregate(micromol.cm2.h ~ Light_Level, data = Pan.data, FUN=mean)
Ber.se <- aggregate(micromol.cm2.h ~ Light_Level, data = Ber.data, FUN=std.error)
Pan.se <- aggregate(micromol.cm2.h ~ Light_Level, data = Pan.data, FUN=std.error)

# mean light conditions for the experiment
Pan.light <- mean(c(606,760,600,624,741,622,580,596,560,637,522))
Pan.light.se <- std.error(c(606,760,600,624,741,622,580,596,560,637,522))
Ber.light <- mean(c(530,575,630,525,460,595,630,580,450))
Ber.light.se <- std.error(c(530,575,630,525,460,595,630,580,450))

#Bermuda Data
PAR <- as.numeric(Ber$Light_Level)
Pc <- as.numeric(Ber$micromol.cm2.h)
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(-1, 1.2), cex.lab=0.8,cex.axis=0.8,cex=1, main="A) Bermuda", adj = 0.05) #set plot info
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,start=list(Am=(max(Pc)-min(Pc)),AQY=0.15,Rd=-min(Pc),theta=0.2)) 
my.fit <- summary(curve.nlslrc) #summary of model fit

#draw the curve using the model fit
Ber.curve.fitting <- curve((1/(2*summary(curve.nlslrc)$coef[4,1]))*(summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1]-sqrt((summary(curve.nlslrc)$coef[2,1]*x+summary(curve.nlslrc)$coef[1,1])^2-4*summary(curve.nlslrc)$coef[2,1]*summary(curve.nlslrc)$coef[4,1]*summary(curve.nlslrc)$coef[1,1]*x))-summary(curve.nlslrc)$coef[3,1],lwd=2,col="blue",add=T)
dev.off()

#Amax (max gross photosytnthetic rate) 
Pmax.gross <- my.fit$parameters[1]

#AQY (apparent quantum yield) alpha  
AQY <- my.fit$parameters[2]

#Rd (dark respiration)
Rd <- my.fit$parameters[3]

# Ik light saturation point
Ik <- Pmax.gross/AQY

# Ic light compensation point
Ic <- Rd/AQY

# Net photosynthetic rates
Pmax.net <- Pmax.gross - Rd

#output parameters into a table
Ber.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(Ber.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")


pdf("Output/MSFigures/FigureS1.pdf")
par(mfrow=c(2,1))
#Plot input data and model fit
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(-1,1.2),cex.lab=0.8,cex.axis=0.8,cex=1, main="A) Bermuda", adj = 0.05) #set plot info
segments(PAR,Pc-Ber.se$micromol.cm2.h,PAR,Pc+Ber.se$micromol.cm2.h)
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels
lines(Ber.curve.fitting ,lwd=2,col="blue") #add fit line
abline(v=Ber.PI.Output[5], col="gray", lty = 2)
abline(v=Ber.light, col="black", lty = 2) #add line for TPC light level

# Panama
PAR <- as.numeric(Pan$Light_Level)
Pc <- as.numeric(Pan$micromol.cm2.h)

#fit a model using a Nonlinear Least Squares regression of a non-rectangular hyperbola (Marshall & Biscoe, 1980)
Pan.curve.nlslrc = nls(Pc ~ (1/(2*theta))*(AQY*PAR+Am-sqrt((AQY*PAR+Am)^2-4*AQY*theta*Am*PAR))-Rd,start=list(Am=(max(Pc)-min(Pc)),AQY=0.05,Rd=-min(Pc),theta=0.01)) 
my.fit <- summary(Pan.curve.nlslrc) #summary of model fit

#Amax (max gross photosytnthetic rate) 
Pmax.gross <- my.fit$parameters[1]

#AQY (apparent quantum yield) alpha  
AQY <- my.fit$parameters[2]

#Rd (dark respiration)
Rd <- my.fit$parameters[3]

# Ik light saturation point
Ik <- Pmax.gross/AQY

# Ic light compensation point
Ic <- Rd/AQY

# Net photosynthetic rates
Pmax.net <- Pmax.gross - Rd

#output parameters into a table
Pan.PI.Output <- rbind(Pmax.gross, Pmax.net, -Rd, AQY,Ik,Ic)
row.names(Pan.PI.Output) <- c("Pg.max","Pn.max","Rdark","alpha", "Ik", "Ic")
Pan.PI.Output

#Plot input data and model fit
plot(PAR,Pc,xlab="", ylab="", xlim=c(0,max(PAR)), ylim=c(-2,4),cex.lab=0.8,cex.axis=0.8,cex=1,main="B) Panama", adj = 0.05) #set plot info
segments(PAR,Pc-Pan.se$micromol.cm2.h,PAR,Pc+Pan.se$micromol.cm2.h)
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=1) #add labels
mtext(expression(Rate*" ("*mu*"mol "*O[2]*" "*cm^-2*h^-1*")"),side=2,line=2,cex=1) #add labels
#draw the curve using the model fit
curve((1/(2*summary(Pan.curve.nlslrc)$coef[4,1]))*(summary(Pan.curve.nlslrc)$coef[2,1]*x+summary(Pan.curve.nlslrc)$coef[1,1]-sqrt((summary(Pan.curve.nlslrc)$coef[2,1]*x+summary(Pan.curve.nlslrc)$coef[1,1])^2-4*summary(Pan.curve.nlslrc)$coef[2,1]*summary(Pan.curve.nlslrc)$coef[4,1]*summary(Pan.curve.nlslrc)$coef[1,1]*x))-summary(Pan.curve.nlslrc)$coef[3,1],lwd=2,col="blue",add=T)
abline(v=Pan.PI.Output[5], col="gray", lty = 2) #add Ik line
abline(v=Pan.light, col="black", lty = 2) #add line for TPC light level
dev.off()

# Create Table of PI curve metrics
Ber.PI.Output <- as.data.frame(Ber.PI.Output)
colnames(Ber.PI.Output) <- c("Bermuda")
Pan.PI.Output<- as.data.frame(Pan.PI.Output)
colnames(Pan.PI.Output) <- c("Panama")
PI.Output.Values <- cbind(Ber.PI.Output, Pan.PI.Output)
PI.Output.Values






