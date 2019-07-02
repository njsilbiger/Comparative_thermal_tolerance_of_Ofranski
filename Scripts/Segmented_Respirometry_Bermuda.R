#Title: Photosynthesis and Respiration Calculations
#Author: HM Putnam
#Edited by: NJ Silbiger
#Date Last Modified: 20170924
# This script uses local linear regressions to calculate photosynthesis and respiration rates for data collected in Bermuda

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
library(devtools)
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 


#Read in required libraries
##### Include Versions of libraries
#install_github('colin-olito/LoLinR')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')

#Required Data files

##### PHOTOSYNTHESIS Time 0 #####
path.p<-"Data/Bermuda/Respirometry/" #the location of all your respirometry files 

# bring in the respiration files
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
#basename above removes the subdirectory name from the file
#add file names that include the subdirectory name
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) #list all csv file names in the folder and subfolders

#generate a 3 column dataframe with specific column names
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=5))
colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec","Temp.C","PR")

#Load Sample meta info Info
Sample.Info <- read.csv(file="Data/Bermuda/MetaData/Nubbin_Sample_Info_T0_Bermuda_QC.csv", header=T) #read sample.info data

#load surface area data
#Add in the SA measurements fot the corals measured in ImageJ
SA<-read.csv('Data/Bermuda/MetaData/SA_Master_QC.csv')


#merge the sample info and SA data
Sample.Info<-merge(Sample.Info, SA, by.x = 'Fragment.ID', all.x = TRUE)

# load water chemistry data (to get the volume of the samples)
WaterChem <- read.csv(file="Data/Bermuda/MetaData/WaterChem_TPC Bermuda_QC.csv", header=T) #read sample.info data

#rename the fragment IDto match all the other ones
WaterChem$Fragment.ID<-paste0(WaterChem$Fragment.ID,"_",WaterChem$Date)

#take the average of the two volumes measures for each sample
Chamber.Vol.mL<-aggregate(WaterChem$Volume~WaterChem$Fragment.ID, FUN=mean)
colnames(Chamber.Vol.mL)<-c("Fragment.ID","Volume")

#merge the chamber volumes with all the sample info
Sample.Info<-merge(Sample.Info, Chamber.Vol.mL, by.x = 'Fragment.ID', all.x = TRUE)
# make start and stop times real times
Sample.Info$Start.time <- as.POSIXct(Sample.Info$Start.time,format="%H:%M", tz = "") #convert time from character to time
Sample.Info$Stop.Time <- as.POSIXct(Sample.Info$Stop.Time,format="%H:%M", tz = "") #convert time from character to time

# replace NA in spot test with 0
Sample.Info$Spot_test[is.na(Sample.Info$Spot_test)] <- 0
# sort the data by light and dark
Sample.Info<-Sample.Info[with(Sample.Info, order(Light_Dark, Run, Chamber.Channel)),]

#Add names for photosynthesis or respiration for for loop
PR<-c('Photo','Resp')

for(i in 1:length(file.names.full)) { # for every file in list start at the first and run this following function

    # Kate's computer exported the files differently so I need to import them differently
  #find the lines in sample info that have the same file name that is being brought it
  FRow<-which(Sample.Info$Fragment.ID==strsplit(file.names[i],'.csv'))
  # there are two rows since it is light and dark. Pull the first one
  # Kate's computer was only used on chambers 11-14 so run an if statement to import based on channel or if it is a spot test
  if(Sample.Info[FRow[1],"Chamber.Channel"]>=11 & Sample.Info[FRow[1],"Run"]<=8 || Sample.Info[FRow[1],"Spot_test"]==1){
    Photo.Data1 <-read.csv(file.path(path.p,file.names.full[i]),  header=T) #reads in the data files
    Photo.Data1  <- Photo.Data1[,c(2,7,9)] #subset columns of interest
    # the date is in european format
    # pull out the time
    Photo.Data1$Time<- sapply(strsplit(as.character(Photo.Data1$Date), " "), "[[", 2)
    Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M", tz = "") #convert time from character to time
    # extract the date  
    #rename Oxygen to Value
    colnames(Photo.Data1)[c(2,3)]<-c('Value','Temp')
    #reorder so its similar to all the other files
    Photo.Data1<-Photo.Data1[,c(4,2,3)]
  }else{# the rest of the data loads with this code
  Photo.Data1 <-read.csv(file.path(path.p,file.names.full[i]), skip = 1, header=T) 
  Photo.Data1  <- Photo.Data1[,c(2,9,16)] #subset columns of interest
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
  }
  
  #find the data that corresponds with the light data
  
  Photo<-Photo.Data1[Photo.Data1$Time>Sample.Info[FRow[2],"Start.time"] & Photo.Data1$Time<Sample.Info[FRow[2],"Stop.Time"],]
  Resp<-Photo.Data1[Photo.Data1$Time>Sample.Info[FRow[1],"Start.time"] & Photo.Data1$Time<Sample.Info[FRow[1],"Stop.Time"],]
  
  lt.levs <- list(Photo, Resp) #list levels of segmentation

for(j in 1:length(lt.levs)){    # for photo and resp (1 is light and 2 is dark)
  Photo.Data <- as.data.frame(lt.levs[j])
  
  # when the computer restarted the time data got thrown out of wack so if the computer shut down
  # then resort the data and fill in the mising time stamp with zeros
  if(j==2 & any(diff(Photo.Data$Time) != 1, na.rm=T)){ # if it is respiration and there is any data that is not exactly 1 second intervals
  # first sort the data by time
    Photo.Data<-Photo.Data[order(Photo.Data$Time),]
    
  #create an empty vector of times from start to stop time for dark
  ts <- seq.POSIXt(as.POSIXct(Sample.Info[FRow[1],"Start.time"],'%m/%d/%y %H:%M'),
                   as.POSIXct(Sample.Info[FRow[1],"Stop.Time"],'%m/%d/%y %H:%M'), by="sec")
  # turn it into a dataframe
  df <- data.frame(Time=ts)
  # join the time points which will create NAs at all the missing times
  Photo.Data<- full_join(df,Photo.Data)
  }
  
  #Now add in values for missing time so that regs does a better job of selecting data
  n<-dim(Photo.Data )[1] #identify length of data
  Photo.Data <-Photo.Data [120:(n-3),] #start at data point ~2 minute in to avoid excess noise from start of run and remove last 3 lines containing text
  n<-dim(Photo.Data )[1] #list length of trimmed data
  Photo.Data$sec <- 1:n #set seconds by one from start to finish of run
  
  #Save plot prior to and after data thinning to make sure thinning is not too extreme
  rename <- sub("_.*", "", file.names[i])
  pdf(paste0("Output/Bermuda/Photo_Resp_Output/",rename,"_",PR[j],".pdf"))
  
  par(omi=rep(0.3, 4)) #set size of the outer margins in inches
  par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
  plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot data as a function of time
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  
  #save original unthinned data
  Photo.Data.orig<-Photo.Data
  Photo.Data <-  thinData(Photo.Data ,by=20)$newData1 #thin data by every 20 points for all the O2 values
  Photo.Data$sec <- as.numeric(rownames(Photo.Data )) #maintain numeric values for time
  Photo.Data$Temp<-NA
  Photo.Data$Temp <-  thinData(Photo.Data.orig,xy = c(1,3),by=20)$newData1[,2] #thin data by every 20 points for the temp values
  
  plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
 # dev.off()
  ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
  #option to add multiple outputs method= c("z", "eq", "pc")
  Regs  <-  rankLocReg(xall=Photo.Data $sec, yall=Photo.Data $Value, alpha=0.3, 
                       method="pc", verbose=FALSE) 
  # add the regression data
  plot(Regs)
  dev.off()
  
  s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
  Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
  Photo.R[j+s[i],1] <- rename #stores the file name in the Date column
  #Photo.R[j+s[i],1] <- paste0(rename,"_",j) #stores the file name in the Date column
  Photo.R[j+s[i],4] <- mean(Photo.Data$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
  Photo.R[j+s[i],5] <- PR[j] #stores whether it is photosynthesis or respiration
  
}
}

#save the output in case anything happens
#write.csv(Photo.R, 'Photo.R.csv')

#remove spot test data
Sample.Info<-Sample.Info[Sample.Info$Spot_test==0,]
Sample.Info.SPOT<-Sample.Info[Sample.Info$Spot_test==1,]

#Split up the photostynthesis and respiration data into two dataframes
PHO <- Photo.R[Photo.R$PR=='Photo', ]
RES <- Photo.R[Photo.R$PR=='Resp', ]

#Clean up the frag ID name
PHO$Fragment.ID <- sub("_.*", "", PHO$Fragment.ID)
RES$Fragment.ID <- sub("_.*", "", RES$Fragment.ID)
Sample.Info$Fragment.ID<-sub("_.*", "", Sample.Info$Fragment.ID)

#Load Sample Info
#Convert sample volume to mL
Sample.Info$Vol.L <- Sample.Info$Volume/1000 #calculate volume

#split sample info between photo and resp
Sample.Info.P<-Sample.Info[Sample.Info$Light_Dark=='Light',]
Sample.Info.R<-Sample.Info[Sample.Info$Light_Dark=='Dark',]

#Merge data with sample info
Resp <- merge(RES,Sample.Info.R, by="Fragment.ID" )
Photo <- merge(PHO,Sample.Info.P, by="Fragment.ID")

#Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
Resp$umol.sec <- Resp$umol.L.sec*Resp$Vol.L
Photo$umol.sec <- Photo$umol.L.sec*Photo$Vol.L

## BLANKS HAVE TO BE SPECIFIC TO RESPONSE VARIABLE (I.E., PHOTO OR RESP) AND TEMP
#Calculate the photo resp rate
photo.blnk <- aggregate(umol.sec ~ Sample.Type, data=Photo, mean)
photo.Blanks <- subset(photo.blnk, Sample.Type == "Blank")
Photo$Blank <- photo.Blanks[1,2]

#Calculate resp blank rate
resp.blnk <- aggregate(umol.sec ~ Sample.Type, data=Resp, mean)
resp.Blanks <- subset(resp.blnk, Sample.Type == "Blank")
Resp$Blank <- resp.Blanks[1,2]

#Account for blank rate Subtract Blank by the temperature blank
Resp$umol.sec.corr <- Resp$umol.sec-Resp$Blank
Photo$umol.sec.corr <- Photo$umol.sec-Photo$Blank

#normalize to surface area and h-1
Resp$umol.cm2.hr <- (Resp$umol.sec.corr*3600)/Resp$Surf.Area.cm2
Photo$umol.cm2.hr <- (Photo$umol.sec.corr*3600)/Photo$Surf.Area.cm2

#remove blanks from data set
Pnet <- subset(Photo, Species!= "BK")
Rdark <- subset(Resp, Species!= "BK")

#calculate gross photosynthesis Pnet -- Rdark
NameList <- c("Fragment.ID","Temp.C", "Date", "Chamber.Channel", "Position", "Bin.ID", "Run", "Species", "Genotype", "Fragment.Number", "Temp.Cat", "umol.cm2.hr")
resp.data <- merge(Pnet[,NameList],Rdark[,c("Fragment.ID","umol.cm2.hr")],  by="Fragment.ID")
#rename the columns
names(resp.data)[names(resp.data) == "umol.cm2.hr.x"]<- "Pnet_umol.cm2.hr" 
names(resp.data)[names(resp.data) == "umol.cm2.hr.y"] <- "Rdark_umol.cm2.hr"
#Pnet plus resp (if positive) is pGross
resp.data$Pgross_umol.cm2.hr <- resp.data$Pnet_umol.cm2.hr-resp.data$Rdark_umol.cm2.hr

#Calculate means
AllMeans <- ddply(resp.data, c('Species','Temp.Cat'), summarize,
                  #pnet
                  Pnet.mean= mean(Pnet_umol.cm2.hr, na.rm=T), #mean pnet
                  N = sum(!is.na(Pnet_umol.cm2.hr)), # sample size
                  Pnet.se = sd(Pnet_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                  #Rdark
                  Rdark.mean= mean(Rdark_umol.cm2.hr, na.rm=T), #mean rdark
                  Rdark.se = sd(Rdark_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                  #Pgross
                  Pgross.mean  = mean(Pgross_umol.cm2.hr, na.rm=TRUE),
                  Pgross.se = sd(Pgross_umol.cm2.hr, na.rm=TRUE)/sqrt(N),
                  Temp.mean = mean(Temp.C, na.rm=TRUE),
                  Temp.se = sd(Temp.C, na.rm=TRUE)/sqrt(N)
                  )

#write Results
write.csv(resp.data, file="Data/BermudaRates.csv") # raw data
#write.csv(AllMeans, file="Output/BermudaMeans.csv") # Mean data

